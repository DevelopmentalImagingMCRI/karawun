#!/usr/bin/env python
#
# Copyright 2019 Murdoch Children's Research Institute,
# Melbourne, Australia
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
#
# python tools for creating Brainlab compatible dicoms from nifti and
# mrtrix streamline files.
#
# Author: Richard Beare
#
import pydicom as pydi
import pydicom._storage_sopclass_uids as storage_sopclass
import SimpleITK as sitk
import numpy as np
import time
import os
import uuid
import os.path

import re
import struct
import glob
from . import ciedicom
#  Constants that brainlab uses in the streamline files - not sure
#  if they are important or not

MyNames = True
if MyNames:
    ImplementationClassUID = \
        '2.25.249096206721547513498031017062911078505'
    ImplementationVersionName = 'RBPythonDicom'
    SourceApplicationEntityTitle = 'DevImDicom'
    ManufacturerModelName = 'StreamlineImport'
    DeviceSerialNumber = '1.2.3.4.5.6'
    SoftwareVersions = '0.1'
    Manufacturer = 'DevIm'
    AlgorithmName = 'Fibertracking'
    ContentCreatorName = 'DevIm'
    SeriesDescription = 'MRTrix fiber bundles'

else:
    # versions in the demo export
    # gets overwritten with the pydicom one
    ImplementationClassUID = '1.2.276.0.20.1.1.21.3.4.0'
    ImplementationVersionName = 'DICOMProxy3.4'
    SourceApplicationEntityTitle = 'BRAINLAB'
    ManufacturerModelName = 'Fibertracking'
    DeviceSerialNumber = '52.2.134.22.89.21'
    SoftwareVersions = '1.0.0.76'
    Manufacturer = 'Brainlab'
    AlgorithmName = 'Fibertracking'
    ContentCreatorName = ""
    SeriesDescription = 'Fiber bundles'


# Stuff to do with cie colours

nice_colours_cie = [ciedicom.rgb2DicomLab(col) for
                    col in ciedicom.nice_colours_rgb]


class Error(Exception):
    """Base class for other exceptions"""
    pass


class MissingUIDList(Error):
    """Raised when the required UID list is missing"""
    pass


class RawToLabelImMismatch(Error):
    """Raised when one of the label images doesn't match any of the raw"""
    pass


# mrtrix tckfile stuff
# converts the MIF datatypes
# Bit    bitwise data
# Int8    signed 8-bit (char) integer
# UInt8    unsigned 8-bit (char) integer
# Int16    signed 16-bit (short) integer
# UInt16    unsigned 16-bit (short) integer
# Int16LE    signed 16-bit (short) integer (little-endian)
# UInt16LE    unsigned 16-bit (short) integer (little-endian)
# Int16BE    signed 16-bit (short) integer (big-endian)
# UInt16BE    unsigned 16-bit (short) integer (big-endian)
# Int32    signed 32-bit int
# UInt32    unsigned 32-bit int
# Int32LE    signed 32-bit int (little-endian)
# UInt32LE    unsigned 32-bit int (little-endian)
# Int32BE    signed 32-bit int (big-endian)
# UInt32BE    unsigned 32-bit int (big-endian)
# Float32    32-bit floating-point
# Float32LE    32-bit floating-point (little-endian)
# Float32BE    32-bit floating-point (big-endian)
# Float64    64-bit (double) floating-point
# Float64LE    64-bit (double) floating-point (little-endian)
# Float64BE    64-bit (double) floating-point (big-endian)
#
# into fields suitable for struct.pack and struct.unpack
# returns a tuple,
#    the first element is the endianness,
#        '>' for explicit big-endian
#        '<' for explicit big-endian
#        '' for native
#    the second element is the data type
#        c     char     string of length 1     1
#        b     signed char     integer     1     (3)
#        B     unsigned char     integer     1     (3)
#        h     short     integer     2     (3)
#        H     unsigned short     integer     2     (3)
#        i     int     integer     4     (3)
#        I     unsigned int     integer     4     (3)
#        l     long     integer     4     (3)
#        L     unsigned long     integer     4     (3)
#        q     long long     integer     8     (2), (3)
#        Q     unsigned long long     integer     8     (2), (3)
#        f     float     float     4     (4)
#        d     double


def datatype_get_struct_strings(datatype):
    pat = re.compile(r'(U)?((Int|Float)(\d+))(LE|BE)?')
    mat = pat.match(datatype)

    if mat is None:
        return None
    else:
        T = {'Int8': 'b', 'Int16': 'h', 'Int32': 'i', 'Float32': 'f',
             'Float64': 'd'}
        dataletter = 'u'
        if mat.group(2) in T:
            dataletter = T[mat.group(2)]
            if mat.group(1) is not None:
                dataletter = dataletter.upper()
        size = int(int(mat.group(4)) / 8)
        if mat.group(5) is None:
            endian = ''
        else:
            E = {'LE': '<', 'BE': '>'}
            endian = E[mat.group(5)]
        return dataletter, endian, size


def structstrings_get_datatype(letter, endian):
    T = {'b': 'Int8', 'h': 'Int16', 'i': 'Int32', 'f': 'Float32',
         'd': 'Float64'}

    if letter == letter.upper():
        signedPart = 'U'
    else:
        signedPart = ''

    datatypePart = T[letter.lower()]

    if len(endian) == 1:
        E = {'<': 'LE', '>': 'BE'}
        endianPart = E[endian]
    else:
        endianPart = ''
    return signedPart + datatypePart + endianPart


def load_trackfile(fileName, origVectorMode=False):
    """
    Read an mrtrix .tck file
    :param fileName: the .tck file
    :param origVectorMode: if True, load the file in one go.
    otherwise do by block. Saves a bit of ram
    :return: dict with the following keys
         "init_threshold": the lowest FA to initialise tracks- can
                           be "variable" when tckfiles are combined
         "lmax": the spherical harmonic order
         "max_dist": maximum length of any tracts
         "max_num_attempts": maximum number of tracts sent out
         "max_num_tracks": maximum number of accepted tracts
         "max_trials": number of probabilistic trials at each point
         "method": tractography method
         "min_curv": ???
         "min_dist": minimum track distance
         "no_mask_interp": whether we DONT interpolate the masks
         "sh_precomputed": Legendre polynomials precomputed?
         "source": tensor or CSD image used for tracking
         "step_size": pixel step size
         "stop_when_included": whether we stopped when we went into
                               an include zone - can  be "variable"
                               when tckfiles are combined
         "threshold": FA threshold- can
                      be "variable" when tckfiles are combined
         "unidirectional": whether we sent out in one direction - can
                           be "variable" when tckfiles are combined
         "roi":
            "type (list of strings)": seed or mask
            "file (list of strings)": filenames
        "datatype (dict)":
            "letter": string that has the pack/unpack letter for the datatype
            "endian": '<' for little or '>' for big or '' for native
            "size": the size, in bytes, of each element
         "count": number of accepted tracks
         "total_count": number of tracks sent out
         "tracks (list of numpy arrays)": the track data
    """

    if not os.path.isfile(fileName):
        print("warning: could not find track file: " + fileName)
        return None
    else:
        FID = open(fileName, 'rb')
        # first line should be "mrtrix image"
        firstLine = FID.readline().decode().rstrip()

        if firstLine != 'mrtrix tracks':
            print("file " + fileName + " not mrtrix track file")
            FID.close()
            return None
        else:
            del firstLine

            trackStruct = dict()

            trackStruct['init_threshold'] = None
            trackStruct['lmax'] = None
            trackStruct['max_dist'] = None
            trackStruct['max_num_attempts'] = None
            trackStruct['max_num_tracks'] = None
            trackStruct['max_trials'] = None
            trackStruct['method'] = None
            trackStruct['min_curv'] = None
            trackStruct['min_dist'] = None
            trackStruct['no_mask_interp'] = None
            trackStruct['sh_precomputed'] = None
            trackStruct['source'] = None
            trackStruct['step_size'] = None
            trackStruct['stop_when_included'] = None
            trackStruct['threshold'] = None
            trackStruct['unidirectional'] = None
            trackStruct['mrtrix_version'] = None
            trackStruct['roi'] = dict()
            trackStruct['roi']['type'] = list()
            trackStruct['roi']['file'] = list()
            trackStruct['datatype'] = None
            trackStruct['count'] = None
            trackStruct['total_count'] = None
            trackStruct['tracks'] = list()
            pat = re.compile('^([a-z_]+): (.*)$')

            while True:
                curLine = FID.readline().decode().rstrip()
                if curLine == "END":
                    break
                else:
                    mat = pat.match(curLine)
                    if mat is not None:

                        curKeyword = mat.group(1)
                        curValue = mat.group(2)

                        if curKeyword in ['init_threshold', 'max_dist',
                                          'min_curv', 'min_dist',
                                          'step_size', 'threshold']:
                            try:
                                trackStruct[curKeyword] = float(curValue)
                            except ValueError:
                                trackStruct[curKeyword] = curValue
                        elif curKeyword in ['lmax', 'max_num_attempts',
                                            'max_num_tracks',
                                            'max_trials',
                                            'no_mask_interp',
                                            'sh_precomputed',
                                            'threshold',
                                            'count',
                                            'unidirectional',
                                            'stop_when_included',
                                            'total_count']:
                            try:
                                trackStruct[curKeyword] = int(curValue)
                            except ValueError:
                                trackStruct[curKeyword] = curValue
                        elif curKeyword in ['method', 'source',
                                            'mrtrix_version']:
                            trackStruct[curKeyword] = curValue[:]
                        elif curKeyword == 'datatype':
                            trackStruct['datatype'] = dict()
                            (
                                trackStruct['datatype']['letter'],
                                trackStruct['datatype']['endian'],
                                trackStruct['datatype']['size']
                            ) = datatype_get_struct_strings(curValue)
                        elif curKeyword == 'roi':
                            roimat = re.match(r'(\S+)\s+(\S+)', curValue)
                            if roimat is not None:
                                trackStruct['roi']['type']. \
                                    append(roimat.group(1))
                                trackStruct['roi']['file']. \
                                    append(roimat.group(2))
                            del roimat
                        elif curKeyword == 'file':
                            trackStruct['file'] = curValue[:]
                            filemat = re.match(r'\.\s+(\d+)', curValue)
                            if filemat is not None:
                                trackStruct['dataoffset'] = \
                                    int(filemat.group(1))
                            del filemat
        # if firstLine != 'mrtrix tracks':
        if trackStruct['datatype'] is not None:
            FID.seek(0, 2)  # 2 means the end of the file
            fileSize = FID.tell()

            if (
                    trackStruct['dataoffset'] >= 0 and
                    trackStruct['dataoffset'] < fileSize):
                # 0 means absolute position
                FID.seek(trackStruct['dataoffset'], 0)

                numElements = (
                    (fileSize - trackStruct['dataoffset']) /
                    trackStruct['datatype']['size'])
                if numElements % 3 != 0:
                    print('The data section of ' + fileName +
                          ' does not have a multiple of 3 elements')
                else:
                    trackStruct['tracks'] = list()
                    # original fully vector reading of the file,
                    # read ALL data at once
                    if origVectorMode:
                        T = FID.read()
                        FID.close()
                        T = struct.unpack(
                            trackStruct['datatype']['endian'] +
                            trackStruct['datatype']['letter'] *
                            numElements, T)
                        T = np.reshape(np.array(T),
                                       (3, numElements / 3),
                                       order='F')

                        sepidx = \
                            np.where(np.any(
                                np.logical_not(np.isfinite(T)),
                                axis=0))[0]
                        # print(sepidx)
                        leftIDX = 0

                        for z in range(sepidx.size - 1):
                            trackStruct['tracks'].append(
                                np.array(T[:, leftIDX:sepidx[z]]))
                            leftIDX = sepidx[z] + 1
                            # end original fully vector mode
                    else:
                        # block reading mode
                        numTriplesInFile = int(numElements / 3)
                        # how many triples to read in each time
                        tripleBlockSize = int(5000)

                        numTriplesRead = int(0)

                        # curTrack is used to store tracks
                        # that straddle block boundaries
                        curTrack = list()
                        while numTriplesRead < numTriplesInFile:

                            B = FID.read(
                                tripleBlockSize * 3 *
                                trackStruct['datatype']['size'])
                            curNumTriplesRead = int(
                                len(B) / 3 /
                                trackStruct['datatype']['size'])

                            T = struct.unpack(
                                trackStruct['datatype']['endian'] +
                                trackStruct['datatype']['letter'] *
                                (curNumTriplesRead * 3), B)
                            del B
                            T = np.reshape(
                                np.array(T),
                                (3, curNumTriplesRead), order='F')

                            # print T
                            sepidx = np.where(
                                np.any(np.logical_not(np.isfinite(T)),
                                       axis=0))[0]
                            # allsepidx.append(sepidx + numTriplesRead)
                            numTriplesRead += curNumTriplesRead
                            # no separators in this set, just append
                            # the whole thing to the current track
                            if 0 == np.size(sepidx):
                                curTrack.append(np.array(T))
                            else:
                                # the separator is at the start, push
                                # the current track onto the main list
                                if sepidx[0] == 0:
                                    trackStruct['tracks'].append(
                                        np.concatenate(curTrack,
                                                       axis=1))
                                else:
                                    # otherwise append and
                                    # then add the track
                                    curTrack.append(
                                        T[:, :sepidx[0]])
                                    trackStruct['tracks'].append(
                                        np.concatenate(curTrack,
                                                       axis=1))
                                curTrack = list()

                                if np.size(sepidx) > 1:
                                    # add all the intermediate tracks
                                    for z in range(np.size(sepidx) - 1):
                                        if (
                                            sepidx[z + 1]
                                                > sepidx[z] + 1):
                                            (
                                                trackStruct['tracks'].
                                                append(
                                                    T[:,
                                                      (sepidx[z] + 1):
                                                      sepidx[z + 1]]
                                                )
                                            )
                                if sepidx[-1] < T.shape[1] - 1:
                                    curTrack.append(
                                        T[:, (sepidx[-1] + 1):])
                            del T

                        FID.close()
                    numTracks = len(trackStruct['tracks'])
                    if numTracks != trackStruct['count']:
                        print('Warning: the number of tracks in'
                              ' the file ' + fileName +
                              ', which was ' + str(numTracks) +
                              ', does not match the number specified by'
                              ' the count header field, which was ' +
                              str(trackStruct['count']))

        return trackStruct


def save_trackfile(trstruct, fileName):
    """
    Save a mrtrix track file.
    :param trstruct: (dict)
        MANDATORY FIELDS:
       "datatype" (dict):
            "letter": string that has the pack/unpack letter for the datatype
            "endian": '<' for little or '>' for big or '' for native
            "size": the size, in bytes, of each element
        "count": number of accepted tracks
        "tracks": the track data
        OPTIONAL FIELDS:
         "init_threshold": the lowest FA to initialise tracks
         "lmax": the spherical harmonic order
         "max_dist": maximum length of any tracts
         "max_num_attempts": maximum number of tracts sent out
         "max_num_tracks": maximum number of accepted tracts
         "max_trials": number of probabilistic trials at each point
         "method": tractography method
         "min_curv": ???
         "min_dist": minimum track distance
         "no_mask_interp": whether we DONT interpolate the masks
         "sh_precomputed": Legendre polynomials precomputed?
         "source": tensor or CSD image used for tracking
         "step_size": pixel step size
         "stop_when_included": whether we stopped when we went into
          an include zone
         "threshold": FA threshold
         "unidirectional": whether we sent out in one direction
         "total_count": number of tracks sent out

    :param fileName: destination file.
    """

    # check the must have fields in trstruct
    mustHaveFields = ['datatype', 'tracks']

    for z in range(len(mustHaveFields)):
        if not mustHaveFields[z] in list(trstruct.keys()):
            print("The field: " + mustHaveFields[z] +
                  " was not in the structure")
            return

    if not (trstruct['datatype']['letter'] == 'f' or
            trstruct['datatype']['letter'] == 'd'):
        print("The tracks must have a floating point datatype")
        return

    try:
        FID = open(fileName, 'wb')
    except Exception:
        print("Could not open " + fileName + " for writing")
        return

    FID.write("mrtrix tracks\n".encode())

    for curField in sorted(trstruct.keys()):
        if curField == 'roi':
            for z in range(len(trstruct[curField]['type'])):
                FID.write(
                    ('roi: %s %s\n' %
                     (trstruct[curField]['type'][z],
                      trstruct[curField]['file'][z])).encode()
                )
        elif curField == 'datatype':
            FID.write(
                ('datatype: %s\n' %
                 structstrings_get_datatype(trstruct[curField]['letter'],
                                            trstruct[curField][
                                                'endian'])).encode())
        elif curField != 'tracks':
            FID.write(("%s: %s\n" %
                       (curField, str(trstruct[curField]))).encode())

    if 'count' not in trstruct:
        FID.write(("count: %d\n" % (len(trstruct['tracks']))).encode())
    FID.write("file: ".encode())
    offset = FID.tell()
    offset = offset + 14
    FID.write((". %d\nEND\n" % offset).encode())
    curPos = FID.tell()

    bytesToWrite = offset - curPos

    if bytesToWrite > 0:
        FID.write(('\x00' * bytesToWrite).encode())

    if trstruct['datatype']['letter'] == 'f':
        curPadding = np.single(np.array([np.nan, np.nan, np.nan]))
    else:
        curPadding = np.double(np.array([np.nan, np.nan, np.nan]))

    curPadding.tolist()
    for z in range(len(trstruct['tracks'])):

        if trstruct['datatype']['letter'] == 'f':
            curTrack = np.single(
                trstruct['tracks'][z]).flatten(order='F')
        else:
            curTrack = np.double(
                trstruct['tracks'][z]).flatten(order='F')

        # "vectorised" version
        curTrack = curTrack.tolist()

        FID.write(struct.pack(trstruct['datatype']['endian'] +
                              trstruct['datatype']['letter'] *
                              len(curTrack),
                              *curTrack))
        FID.write(struct.pack(trstruct['datatype']['endian'] +
                              trstruct['datatype']['letter'] *
                              len(curPadding),
                              *curPadding))
        del curTrack

    FID.close()


def get_mrtrixtransform_from_nibabel(NIIHeader):
    """
    Convert between nifti header transforms and mrtrix
    transforms
    :param NIIHeader: an niibabel.load header
    :return: a numpy array  containing the transform.
    """
    QForm, QFormCode = NIIHeader.get_qform(coded=True)
    SForm, SFormCode = NIIHeader.get_sform(coded=True)

    pixDims = NIIHeader.get_header()['pixdim'][1:4]

    if QFormCode <= 0 or SFormCode > 0:
        transform = np.array(SForm)
        voxelSizes = np.sqrt(np.sum(transform[0:3, 0:3] *
                                    transform[0:3, 0:3], axis=0))
        transform[0:3, 0:3] = (transform[0:3, 0:3]
                               / np.atleast_2d(voxelSizes))
    else:
        # use the qcode, just divide by the voxel sizes
        transform = np.array(QForm)
        transform[0:3, 0:3] = (transform[0:3, 0:3]
                               / np.atleast_2d(pixDims))
        pass

    permutation = np.argmax(np.abs(transform[0:3, 0:3] *
                                   transform[0:3, 0:3]), axis=1)

    flip = transform[(np.arange(3), permutation)]
    flip = (flip < 0)
    # print permutation

    if not np.array_equal(permutation, np.array([0, 1, 2])) or \
            np.any(flip):
        origTransform = np.array(transform)
        transform[0:3, 0:3] = transform.take(
            permutation, axis=1).take([0, 1, 2], axis=0)
        for curColumn in range(3):
            dimLength = pixDims[curColumn] * \
                        (NIIHeader.shape[curColumn] - 1)

            if flip[curColumn]:
                transform[0:3, curColumn] = -transform[0:3, curColumn]
                for curRow in range(3):
                    transform[curRow, 3] = \
                        transform[curRow, 3] + \
                        origTransform[curRow, permutation[curColumn]] * \
                        dimLength
    return transform


def tracks_world_to_img(tracksWorldSpace, NIIHeader):
    """
    Transform mrtrix coordinates to nii space
    :param tracksWorldSpace: list of coordinates
    :param NIIHeader: Header returned by nibabel.load
    :return: list of transformed coordinates
    """
    # NIIHeader *MUST* be direct output from nibabel.load
    tracksIMGSpace = list()
    transform = get_mrtrixtransform_from_nibabel(NIIHeader)

    invTransform = np.matrix(np.linalg.inv(transform))

    invTransform = invTransform[0:3, :]

    tracksSZ = np.zeros((len(tracksWorldSpace)), dtype=np.int64)
    for z in range(len(tracksWorldSpace)):
        tracksSZ[z] = tracksWorldSpace[z].shape[1]

    tracksT = np.concatenate(tracksWorldSpace, axis=1)

    tracksT = np.array(
        invTransform *
        np.matrix(
            np.concatenate((tracksT,
                            np.ones((1, tracksT.shape[1]))), axis=0)))
    tracksT = (tracksT
               / np.atleast_2d(NIIHeader.get_header()['pixdim'][1:4]).T)
    tracksT[1] = (NIIHeader.shape[1] - 1) - tracksT[1]

    leftIDX = 0
    tracksIMGSpace = list()
    for z in range(len(tracksWorldSpace)):
        tracksIMGSpace.append(
            tracksT[:, leftIDX:(leftIDX + tracksSZ[z])])
        leftIDX += tracksSZ[z]

    return tracksIMGSpace


# end mrtrix tckfile stuff
##############################################################
# Brainlab dicom streamlines to mrtrix
# Mostly for testing/comparison

def mrt_streamlines(surface):
    def getAllPoints(surface):
        # get the points into an array, with 1 row per point
        if (len(surface.SurfacePointsSequence) > 1):
            raise ValueError(
                'Not expecting more than one SurfacePointSequence')
        p = surface.SurfacePointsSequence[0].NumberOfSurfacePoints
        coordlist = \
            surface.SurfacePointsSequence[0].PointCoordinatesData
        coords = np.array(coordlist)
        coords = coords.reshape((p, 3)).transpose()
        return coords

    def getAllLines(surface):
        # Number of streamlines
        k = len(surface.SurfaceMeshPrimitivesSequence[0].LineSequence)
        if (len(surface.SurfaceMeshPrimitivesSequence) > 1):
            raise ValueError(
                'Not expecting more than one '
                'SurfaceMeshPrimitivesSequence')

        # dump out numpy arrays in a list
        # (note that dicom uses base 1 indexes)
        nn = [np.fromstring(
              surface.SurfaceMeshPrimitivesSequence[0].LineSequence[x].
              PrimitivePointIndexList,
              dtype=np.uint16) - 1 for x in range(k)]
        return nn

    j = getAllPoints(surface)
    g = getAllLines(surface)
    ll = [j[:, m] for m in g]
    ll = [m[:, ] for m in ll]
    return ll


def brainlab_dcm2tck(dcmf, outfile):
    fb = pydi.read_file(dcmf)
    streamlines = [mrt_streamlines(s) for s in fb.SurfaceSequence]
    streamlines = [item for sublist in streamlines for item in sublist]
    tf = np.array([-1, 0, 0, 0, -1, 0, 0, 0, 1])
    tf = tf.reshape(3, 3)
    s2 = [np.array(tf * np.matrix(h)) for h in streamlines]
    datatype = {'endian': '<', 'letter': 'f', 'size': 4}
    tstruct = {'datatype': datatype, 'count': len(streamlines),
               'tracks': s2}
    save_trackfile(tstruct, outfile + ".tck")

###########################################################
# General nifti stuff


def delete_tags(dcm):
    """
    Delete unwanted tags from our sample dicom header

    Fields we don't want to copy
     Rescale intercept, Rescale slope
     Modification time
     modification date
     Image type
     Image orientation
     Slice thickness
     window center
     window width
     SOP instance UID
     PixelData
     Image position patient
     Instance number
     Operator name (retired field)
     SeriesInstanceUID
     StudyIntstanceUID
     Acquisition matrix
     Bits Stored
     High Bit
     Frame of reference

    :param dcm: the pydicom structure
    :return: modified structure
    """

    ignore = [("0x0028", "0x1052"), ("0x0028", "0x1053"),
              ("0x0008", "0x0031"), ("0x0008", "0x0021"),
              ("0x0008", "0x0008"), ("0x0020", "0x0037"),
              ("0x0018", "0x0050"), ("0x0028", "0x1050"),
              ("0x0028", "0x1050"), ("0x0008", "0x0018"),
              ("0x0008", "0x0016"), ("0x7fe0", "0x0010"),
              ("0x0020", "0x0032"), ("0x0020", "0x0013"),
              ("0x0008", "0x1070"), ("0x0020", "0x000e"),
              ("0x0020", "0x000d"), ("0x0018", "0x1310"),
              ("0x0028", "0x0101"), ("0x0028", "0x0102"),
              ("0x0020", "0x0052")]
    for i in ignore:
        tg = pydi.tag.Tag(i)
        if tg in dcm:
            del dcm[tg]

    # Can we delete all the private tags
    private_keys = [k for k in dcm.keys() if k.is_private]
    for k in private_keys:
        del dcm[k]
    return dcm


def dcm_uuid():
    thisuuid = uuid.uuid4()
    uid = '2.25.' + str(thisuuid.int)
    return uid


def str2ds(strings):
    """
    Convert numbers to strings conforming to dicom rules,
    maximum 16 characters.
    :param strings:
    :return: formatted string representation.
    """
    # Decimal string encoding, maximum 16 characters
    h = ['{:+16.8e}'.format(x) for x in strings]
    return h


# Fixed point strings
def str2ds_fixed(strings):
    """
    Convert dicom tag text of fixed point numbers to
    :param strings:
    :return:
    """
    h = ['{:15.10f}'.format(x) for x in strings]
    return h


def mk_file_meta():
    """
    Create the basic pydicom structure for file storage
    :return: a pydicom dataset with a skeleton header
    """
    file_meta = pydi.Dataset()
    # MR Image SOP - Enhanced MR Image SOP
    # contains some extra synchronization
    # fields - use storage_sopclass
    file_meta.MediaStorageSOPClassUID = '1.2.840.10008.5.1.4.1.1.4'
    file_meta.MediaStorageSOPInstanceUID = dcm_uuid()
    # This number uniquely identifies my application.
    # I generated it by calling the UUID tool
    file_meta.ImplementationClassUID = \
        '2.25.249096206721547513498031017062911078505'
    file_meta.ImplementationVersionName = 'RBPythonDicom'
    file_meta.TransferSyntaxUID = pydi.uid.ExplicitVRLittleEndian

    return file_meta


def clone_dcm_meta(dcm):
    """
    Copy an existing pydicom Dataset as a basis for saving
    another image
    :param dcm: the pydicom dataset to be copied
    :return:
    """
    newdcm = pydi.Dataset()
    for k, v in dcm.items():
        newdcm[k] = v
    newdcm.file_meta = mk_file_meta()
    newdcm.is_little_endian = True
    newdcm.is_implicit_VR = False
    newdcm.SOPInstanceUID = newdcm.file_meta.MediaStorageSOPInstanceUID
    newdcm.SOPClassUID = newdcm.file_meta.MediaStorageSOPClassUID
    return newdcm


def scale_slice(imslice, outrange=None):
    """
    Scale a slice for storing as UInt16 dicom
    :param imslice: floating point, 2D simpleITK image
    :param outrange: desired output range, default based on pixel type
    :return: dict containing the following fields
       rescaled_image: simpleitk image with scaling
       origmax: maximum of input
       origmin: minimum of input
       rescaleslope: (origmax - origmin)/outrange
       rescaleintercept: origmin (redundant)
    """

    mm = sitk.MinimumMaximumImageFilter()
    mm.Execute(imslice)
    mx = mm.GetMaximum()
    mn = mm.GetMinimum()
    if outrange is None:
        outrange = float(pow(2, 16) - 1)
    imrange = float(mx - mn)
    rescaleslope = imrange / outrange
    rescaleintercept = mn

    # newslice = (imslice - rescaleintercept)/rescaleslope
    newslice = sitk.ShiftScale(imslice, rescaleintercept,
                               1.0 / rescaleslope)
    newslice = sitk.Cast(newslice, sitk.sitkUInt16)

    return {'rescaled_image': newslice, 'origmax': mx, 'origmin': mn,
            'rescaleslope': rescaleslope,
            'rescaleintercept': rescaleintercept}


def check_isotropy(sitkImage):
    """Return the dimension number which will return a slice with
       isotropic voxels. If the image is isotropic in 3 dimensions,
       return the dimension leading to the smallest number of slices.

       If the image has isotropic axial slices, and the data is
       layed out in MNI space, this function will return 2,
       because sitkImage[:,:, 5] would return an isotropic
       slice and the "5" corresponds to image dimension 2.

       Now rounding to 6 decimal places before selecting the
       plane to use. Rounding is only used in plane selection.
       Original dimensions used in writing dicoms.
       """
    spacing = sitkImage.GetSpacing()
    sparray = np.array(spacing)
    # round to 6 decimal
    spu = np.unique(np.around(sparray, 6))
    if np.unique(spu).ravel().shape[0] == 3:
        message = 'No plane with isotropic voxels - stopping - {},{},{}'.format(sparray[0], sparray[1], sparray[2])
        raise ValueError(message)

    sparray = np.around(sparray, 6)
    sb = sparray == sparray[0]
    if np.all(sb):
        dimensions = sitkImage.GetSize()
        oddoneout = dimensions.index(min(dimensions))
    else:
        if sb.sum() == 2:
            sb = np.logical_not(sb)
        # find the index of True
        oddoneout = np.nonzero(sb)[0].item(0)
    return oddoneout


def get_direction(sitkImage, planeidx):
    k = np.array(sitkImage.GetDirection())
    k = k.reshape([3, 3])
    # plane selection will come into account
    idxs = [0, 1, 2]
    idxs.remove(planeidx)
    idxtuple = (slice(None), idxs)
    direction = list(k[idxtuple].transpose().ravel())
    return direction


def mk_indexing_tuple(index, plane):
    if plane == 2:
        return((slice(None), slice(None), index))
    if plane == 1:
        return((slice(None), index, slice(None)))
    if plane == 0:
        return((index, slice(None), slice(None)))

    raise ValueError("Invalid plane number")


def sitk_nifti_to_dicom(niftifile, dicomfile, dcmprefix, outdir,
                        Description=None, StudyUID=None, FrameUID=None,
                        SeriesNum=None):
    """:param niftifile: (string) path to nifti image
    :param dicomfile: (string) path to a dicom file - used to supply some
                   important tags
    :param dcmprefix: (string) prefix for output filenames.
                  _slicenum.dcm will be appended
    :param outdir: (string) target directory, will be created
                if doesn't exist.
    :param Description: (string) to populate the SeriesDescription field.
    :param StudyUID: (string) Study Instance UID - will generate one
                          if none.
    :param FrameUID: (string) a UID that should be common for all
                 dicoms in a common space
    :param SeriesNum: (int) unique within a study - helps identify
                  files in a series
    :return: dictionary of list of strings containing InstanceUIDs for
        constructing other dicoms, the dicom
        structure containing the common parts, and nifti details
        for matching up with label images.

     Notes :
     Scaling is an issue. Lots of dicom conversion tools don't
     apply scaling, because they aren't
     part of the MRImage standard. Need to decide what
     to do with this.

     Brainlab doesn't cope with dicoms with anisotropic voxels
     within slice. Dicoms we acquire are always isotropic within
     slice, so attempt to figure out automatically which planes to
     write.

     Here is the selection, if in MNI orientation to start with.
     Axial planes : nif[:,:,i], direction = list(k[:,0:2].transpose().ravel())
     Sagittal planes: nif[i,:,:]
     coronal planes: nif[:,i,:]

    """

    # Read in as float and rescale to UInt16, with rescale values
    nif = sitk.ReadImage(niftifile, sitk.sitkFloat32)
    # check the qform
    qformname = nif.GetMetaData('qform_code_name')
    if qformname == 'NIFTI_XFORM_UNKNOWN':
        raise ValueError(niftifile + ': Unknown qform code - stopping')
    
    spacing = nif.GetSpacing()
    oMatrix = nif.GetDirection()
    imsize = nif.GetSize()

    if len(spacing) > 3:
        raise ValueError(niftifile + ': Higher than 3D nifti file - stopping')

    # figure out which plane to write. Aiming for isotropic within plane
    # If image is isotropic, should use minimum number of planes
    try:
        isoidx = check_isotropy(nif)
    except ValueError:
        print("Error processing " + niftifile)
        raise

    otheridx = [0, 1, 2]
    otheridx.remove(isoidx)
    # Load the sample dicom
    t1d = pydi.read_file(dicomfile)

    modification_time = time.strftime("%H%M%S", time.localtime())
    modification_date = time.strftime("%Y%m%d", time.localtime())

    # Columns of this matrix contain the direction cosines
    # matrix is stored in rows
    direction = get_direction(nif, isoidx)
    direction = str2ds(direction)
    sp = nif.GetSpacing()
    sp = str2ds(sp)
    spacing = [sp[i] for i in otheridx]
    slthickness = sp[isoidx]

    # fields we don't want to copy
    # Rescale intercept, Rescale slope
    # Modification time
    # modification date
    # Image type
    #
    # Image orientation
    # Slice thickness
    # window center
    # window width
    # SOP instance UID
    # Delete the things we don't want to copy
    t1d = delete_tags(t1d)
    if StudyUID is None:
        t1d.StudyInstanceUID = dcm_uuid()
    else:
        t1d.StudyInstanceUID = StudyUID
    if SeriesNum is not None:
        t1d.SeriesNumber = SeriesNum
    if Description is not None:
        t1d.SeriesDescription = Description
        t1d.ProtocolName = Description
    # Create the series instance UID
    t1d.SeriesInstanceUID = dcm_uuid()
    # Create the output folder
    os.makedirs(outdir, exist_ok=True)

    # scaling stuff - can't be per slice
    scalednif_dict = scale_slice(nif)
    nifscaled = scalednif_dict['rescaled_image']
    gmx = scalednif_dict['origmax']
    gmn = scalednif_dict['origmin']
    RescaleIntercept = scalednif_dict['rescaleintercept']
    RescaleInterceptDS = str2ds([RescaleIntercept])
    RescaleSlope = scalednif_dict['rescaleslope']
    RescaleSlopeDS = str2ds([RescaleSlope])

    SOPlist = list()
    filenames = list()

    for i in range(imsize[isoidx]):
        # Extract a single slice
        # Query the 3D version to sort out orientation etc
        # image_slice = nifscaled[:,:,i]
        selector = mk_indexing_tuple(i, isoidx)
        image_slice = nifscaled[selector]
        # The input is float, but we need to scale it to Int16 and
        # store the rescale and offset values in the dicom
        # so that the original float values can be recovered
        # U = m*SV + b , where U is read from dicom, SV is stored.
        # we need to reverse this
        mx = gmx
        mn = gmn
        windowcentre = (mx + mn) / 2
        windowwidth = (mx - mn)
        thisslice = clone_dcm_meta(t1d)
        thisslice.ImageType = "DERIVED\\SECONDARY\\OTHER"
        thisslice.SeriesTime = modification_time
        thisslice.SeriesDate = modification_date
        thisslice.WindowCenter = windowcentre
        thisslice.WindowWidth = windowwidth
        thisslice.RescaleIntercept = RescaleInterceptDS
        thisslice.RescaleSlope = RescaleSlopeDS
        thisslice.RescaleType = "US"
        thisslice.SmallestImagePixelValue = int(mn)
        thisslice.LargestImagePixelValue = int(mx)
        thisslice.BitsStored = 16
        thisslice.HighBit = 15

        thisslice.InstanceNumber = i + 1
        # Size stuff
        thisslice.Rows = int(imsize[otheridx[1]])
        thisslice.Columns = int(imsize[otheridx[0]])
        # Position/Orientation
        # plane selection
        corner = [0, 0, 0]
        corner[isoidx] = i
        origin = nif.TransformIndexToPhysicalPoint(corner)
        origin = str2ds(origin)
        thisslice.ImagePositionPatient = origin
        thisslice.ImageOrientationPatient = direction
        thisslice.PixelSpacing = spacing
        thisslice.SliceThickness = slthickness
        if FrameUID is not None:
            thisslice.FrameOfReferenceUID = FrameUID

        # Finally - pixel data
        P = sitk.GetArrayFromImage(image_slice)
        thisslice.PixelData = P.tobytes()
        fname = dcmprefix + "_" + format(i, "04") + ".dcm"
        fname = os.path.join(outdir, fname)
        pydi.filewriter.dcmwrite(fname, thisslice,
                                 write_like_original=False)
        SOPlist.append(thisslice.SOPInstanceUID)
        filenames.append(fname)

    return {'SOPlist': SOPlist, 'dcm': t1d,
            'dcmfiles': filenames, 'spacing': nif.GetSpacing(),
            'size': imsize, 'matrix': oMatrix}


##################################################################
# Streamline stuff

def mk_filemeta_streamlines():
    file_meta = pydi.Dataset()
    # Surface segmentation storage - use storage_sopclass
    file_meta.MediaStorageSOPClassUID = '1.2.840.10008.5.1.4.1.1.66.5'
    file_meta.MediaStorageSOPInstanceUID = dcm_uuid()
    # This number uniquely identifies my application.
    # I generated it by calling the UUID tool
    # Use the same one as for regular images
    file_meta.ImplementationClassUID = ImplementationClassUID
    file_meta.ImplementationVersionName = ImplementationVersionName
    file_meta.TransferSyntaxUID = pydi.uid.ExplicitVRLittleEndian
    file_meta.SourceApplicationEntityTitle = \
        SourceApplicationEntityTitle
    return file_meta


# Label stuff
def mk_filemeta_labelobj():
    file_meta = mk_filemeta_streamlines()
    # Segmentation storage - use storage_sopclass
    file_meta.MediaStorageSOPClassUID = '1.2.840.10008.5.1.4.1.1.66.4'
    file_meta.TransferSyntaxUID = pydi.uid.RLECompressedLosslessSyntaxes
    return file_meta

# Put the mrtrix file details into content description
# (64 characters max)
# filename (without .tck), method, lmax, max_dist


def dicom_date_stamps(newdcm, thefile):
    # Put a files date stamp into dicom fields
    stamp = time.localtime(os.path.getmtime(thefile))
    modification_time = time.strftime("%H%M%S", stamp)
    modification_date = time.strftime("%Y%m%d", stamp)
    newdcm.ContentDate = modification_date
    newdcm.ContentTime = modification_time
    newdcm.InstanceCreationDate = modification_date
    newdcm.InstanceCreationTime = modification_time
    newdcm.SeriesDate = modification_date
    newdcm.SeriesTime = modification_time
    return newdcm


def dicom_patient_stuff(newdcm, templatedcm):
    """

    :param newdcm: the pydicom dataset that we are modifying
    :param templatedcm: the template pydicom dataset (probably loaded
                        from a sample image)
    :return: a modified pydicom dataset with patient details
    """
    newdcm.PatientBirthDate = templatedcm.PatientBirthDate
    newdcm.PatientID = templatedcm.PatientID
    newdcm.PatientName = templatedcm.PatientName
    newdcm.PatientSex = templatedcm.PatientSex
    newdcm.ReferringPhysicianName = templatedcm.ReferringPhysicianName
    newdcm.StudyInstanceUID = templatedcm.StudyInstanceUID
    newdcm.StudyDate = templatedcm.StudyDate
    newdcm.StudyTime = templatedcm.StudyTime
    newdcm.StudyID = templatedcm.StudyID
    return newdcm


def dicom_referenced_series_sequence(UIDlist, SeriesInstance):
    """
    Create the referenced UID list. Assume they are all MRI.
    This indicates which dicoms are related to one another.
    :param UIDlist:
    :param SeriesInstance:
    :return: a modified set of datasets
    """
    #
    # A pydicom sequence is a list of datasets
    def OneRef(uid):
        k = pydi.Dataset()
        k.ReferencedSOPClassUID = "1.2.840.10008.5.1.4.1.1.4"
        k.ReferencedSOPInstanceUID = uid
        return k

    RIS = [OneRef(x) for x in UIDlist]
    m = pydi.Dataset()
    m.SeriesInstanceUID = SeriesInstance
    m.ReferencedInstanceSequence = RIS
    return pydi.Sequence([m])


# SeriesNumber - what to do?
def dicom_fibre_skel():
    """
    Create a dataset skeleton, without patient/doctor details,
    frames of reference
    :return: a skeleton pydicom dataset
    """
    newdcm = pydi.Dataset()
    # These fields copied from a Brainlab example,
    # and some of them are empty
    newdcm.AccessionNumber = b''
    newdcm.ContentCreatorName = ContentCreatorName
    newdcm.ContentLabel = "FIBER"
    # Looks like study date, etc should come from the structurals

    # no need for more than one instance?
    newdcm.InstanceNumber = '1'
    # Some names I've made up
    newdcm.Manufacturer = Manufacturer
    newdcm.ManufacturerModelName = ManufacturerModelName

    newdcm.DeviceSerialNumber = DeviceSerialNumber

    newdcm.Modality = 'SEG'
    # empty in my example
    newdcm.PositionReferenceIndicator = b''
    newdcm.SeriesDescription = SeriesDescription
    newdcm.SeriesInstanceUID = dcm_uuid()
    newdcm.SoftwareVersions = SoftwareVersions
    newdcm.SpecificCharacterSet = 'ISO_IR 100'
    return newdcm


def dicom_label_skel():
    """
    Create a dataset skeleton appropriate for label objects,
    without patient/doctor details, frames of reference
    etc.
    :return: a skeleton pydicom dataset
    """
    newdcm = dicom_fibre_skel()
    newdcm.ContentLabel = "EXT SEG"
    newdcm.SeriesDescription = "objects from workup"
    return newdcm


def mrtrix2surf(tck):
    """
    Convert streamlines to what dicom calls surfaces.
       Lists of locations and lists of indexes
       of connected points.
    :param tck: the track object loaded from mrtrix file
    :return: dict containing the surface and point lists
    """

    # coordinate transform between mrtrix/nifti and dicom
    tf = np.array([-1, 0, 0, 0, -1, 0, 0, 0, 1])
    tf = tf.reshape(3, 3)
    streamlines_dc = [np.array(tf @ np.array(h)) for h in
                      tck['tracks']]
    points_dc = [np.r_[0:x.shape[1]] + 1 for x in streamlines_dc]

    # Maximum of 2^16 points per surface
    slcounts = np.array([h.shape[1] for h in streamlines_dc])
    slcumsum = np.cumsum(slcounts)
    surface_list = list()
    point_list = list()
    last = 1
    while True:
        last = np.argmax(slcumsum > (pow(2, 16) - 1))
        if last == 0:
            last = len(slcumsum)
        surface_list.append(streamlines_dc[0:last])
        point_list.append(points_dc[0:last])

        slcounts = slcounts[last:]
        slcumsum = np.cumsum(slcounts)
        streamlines_dc = streamlines_dc[last:]
        points_dc = points_dc[last:]
        if len(slcumsum) == 0:
            break

    # surface list is now in dicom compatible chunks,
    # each with less than 2^16 points
    # Need to create a list of line segments too.
    chunks = len(surface_list)
    for j in range(chunks):
        lines = len(point_list[j])
        total = 0
        for k in range(lines):
            point_list[j][k] = point_list[j][k] + total
            total = total + len(point_list[j][k])

    return {'sl': surface_list, 'pl': point_list}


def mk_surface_sequence(coords, indexes, chk):
    """
    Create the dicom representation of a tract
    :param coords: locations
    :param indexes: point list (how locations connect together)
    :param chk: the chunk number (to do with the number of points
                dicom can represent at one time)
    :return: pydicom dataset representing a streamline.
    """
    res = pydi.Dataset()
    # Coordinates
    res.SurfacePointsSequence = pydi.Sequence([pydi.Dataset()])
    sz = np.array([x.shape[1] for x in coords])

    res.SurfacePointsSequence[0].NumberOfSurfacePoints = sz.sum()

    combinedpoints = np.concatenate(coords, 1)
    collapsedpoints = combinedpoints.transpose().ravel()
    res.SurfacePointsSequence[0].PointCoordinatesData = list(
        collapsedpoints)

    # private stuff
    tg_0067_0010 = pydi.tag.Tag(0x0067, 0x0010)
    tg_0067_1003 = pydi.tag.Tag(0x0067, 0x1003)
    res.SurfacePointsSequence[0].add_new(tg_0067_0010, 'LO',
                                         'Brainlab-S14-SSO')
    # Eventually we should look these values up from an FA/FOD image,
    # but set them constant for now
    fa = np.array(0.5)
    fa = list(fa.repeat(sz.sum()))
    res.SurfacePointsSequence[0].add_new(tg_0067_1003, 'OF', fa)
    # lines (between points, by index)
    SMPS = pydi.Dataset()
    SMPS.TrianglePointIndexList = b''
    SMPS.EdgePointIndexList = b''
    SMPS.VertexPointIndexList = b''
    SMPS.TriangleStripSequence = pydi.Sequence([])
    SMPS.TriangleFanSequence = pydi.Sequence([])
    SMPS.FacetSequence = pydi.Sequence([])

    def onePtIdxList(idxs):
        PL = pydi.Dataset()
        a = np.array(idxs, dtype=np.uint16)
        PL.PrimitivePointIndexList = a.tobytes()
        return PL

    j = [onePtIdxList(x) for x in indexes]
    SMPS.LineSequence = pydi.Sequence(j)
    res.SurfaceMeshPrimitivesSequence = pydi.Sequence([SMPS])
    # other tags from brainlab
    res.FiniteVolume = 'NO'
    res.Manifold = 'NO'
    res.SurfaceComments = "Chunk " + str(chk + 1)
    res.SurfaceNumber = chk + 1
    # empty normals
    res.SurfacePointsNormalsSequence = pydi.Sequence([])
    res.SurfaceProcessing = 'NO'
    res.RecommendedPresentationType = 'WIREFRAME'
    res.RecommendedDisplayCIELabValue = [52512, 35728, 51251]
    res.RecommendedDisplayGrayscaleValue = 41820
    res.RecommendedPresentationOpacity = 1
    return res


def mk_referenced_surface_sequence(SurfNum, UIDlist):
    """
   Referenced surface sequence - need one for each
     "surface" we create
     part of this is a lot like the ReferencedSeriesSequence
     Seems to be a lot of this kind of book keeping in dicom.
     The UIDlist tends to refer to the first image dataset
     going through conversion.
    :param SurfNum:
    :param UIDlist:
    :return: a pydicom dataset containing the reference structure.
    """
    res = pydi.Dataset()
    res.ReferencedSurfaceNumber = SurfNum
    SSGAI = pydi.Dataset()
    AFCS = pydi.Dataset()
    AFCS.CodeValue = '123101'
    AFCS.CodingSchemeDesignator = 'DCM'
    AFCS.CodeMeaning = 'Neighborhood Analysis'

    SSGAI.AlgorithmFamilyCodeSequence = pydi.Sequence([AFCS])
    SSGAI.AlgorithmVersion = 'undefined'
    # Needs to be Fibertracking
    SSGAI.AlgorithmName = 'Fibertracking'
    res.SegmentSurfaceGenerationAlgorithmIdentificationSequence = \
        pydi.Sequence([SSGAI])

    def OneRef(uid):
        k = pydi.Dataset()
        k.ReferencedSOPClassUID = '1.2.840.10008.5.1.4.1.1.4'
        k.ReferencedSOPInstanceUID = uid
        return k

    RIS = [OneRef(x) for x in UIDlist]
    SSSIS = pydi.Sequence(RIS)
    res.SegmentSurfaceSourceInstanceSequence = SSSIS
    return res


def mk_segment_sequence(description, label, UIDlist, surfaces):
    """
    Create the dicom segmentation structure.
    :param description: something like "mrtrix" to show the method
    :param label: some text, like tract name, or filename
    :param UIDlist: cross referencing stuff
    :param surfaces: the pydicom surface structure
    :return: pydicom dataset containing all the bits and pieces
    (various brainlab codes)
    """
    res = pydi.Dataset()
    ARS = pydi.Dataset()
    # This is always saying "head"
    ARS.CodeValue = 'T-D1100'
    ARS.CodingSchemeDesignator = 'SRT'
    ARS.CodeMeaning = 'Head'
    res.AnatomicRegionSequence = pydi.Sequence([ARS])

    res.SegmentAlgorithmType = 'SEMIAUTOMATIC'
    res.SegmentDescription = description
    res.SegmentLabel = label  # from the file name
    res.SegmentNumber = 1

    SPTCS = pydi.Dataset()
    SPTCS.CodeValue = 'SEG-PT-0080'
    SPTCS.CodingSchemeDesignator = '99BL-GEN'
    SPTCS.CodeMeaning = 'Fiber Structure'

    res.SegmentedPropertyTypeCodeSequence = pydi.Sequence([SPTCS])

    SPCCS = pydi.Dataset()
    SPCCS.CodeValue = 'R-42019'
    SPCCS.CodingSchemeDesignator = 'SRT'
    SPCCS.CodeMeaning = 'Function'
    res.SegmentedPropertyCategoryCodeSequence = pydi.Sequence([SPCCS])

    RSS = [mk_referenced_surface_sequence(idx + 1, UIDlist) for idx in
           range(surfaces)]
    res.ReferencedSurfaceSequence = pydi.Sequence(RSS)
    res.SurfaceCount = surfaces

    # Finally, the dodgy private tags that appear essential,
    # but aren't in the documentation
    # More complex way of adding, as they aren't standard tags
    tg_0067_0010 = pydi.tag.Tag(0x0067, 0x0010)
    tg_0067_1001 = pydi.tag.Tag(0x0067, 0x1001)

    res.add_new(tg_0067_0010, 'LO', 'Brainlab-S14-SSO')

    S_0067_1001 = pydi.Dataset()
    S_0067_1001.CodeValue = 'BL-0006'
    S_0067_1001.CodingSchemeDesignator = '99BL-S14-1'
    S_0067_1001.CodeMeaning = 'Fiber Bundle'
    res.add_new(tg_0067_1001, 'SQ', pydi.Sequence([S_0067_1001]))
    return res


def tck_to_dicom(tckfile, dicomfile, outputfile, seriesNum=0,
                 Description=None, StudyUID=None, UIDlist=None,
                 FrameUID=None):
    """
    Create a brainlab compatible dicom containing streamlines
    from a mrtrix tckfile. Writes to outfile as well as
    returning the structure
    :param tckfile: (string) mrtrix streamline file
    :param dicomfile:  dicom from the newly created MR volume
                            (not the original dicom). This is used
                            to provide the patient information and
                            study ID. StudyInstanceUID will be
                            overwritten if StudyUID is provided
    :param outputfile: target filename
    :param seriesNum: series number - used to distinguish data within
    one set
    :param Description: used to fill SeriesDescription
    :param StudyUID: ID for this study
    :param UIDlist: list of referenced image UIDs
    :param FrameUID: frame of reference stuff
    :return: a pydicom Dataset
    """

    # load the mrtrix file
    tck = load_trackfile(tckfile)
    tckname = os.path.basename(tckfile)
    tckname = os.path.splitext(tckname)[0]
    # collapse it to basic surface structure
    tck_as_surf = mrtrix2surf(tck)

    # make some description strings for content description
    ContDesc = 'mrtrix'
    if tck['mrtrix_version'] is not None:
        ContDesc = tck['mrtrix_version'] + ',' + tck[
            'method'] + ",lmax=" + str(tck['lmax'])
    # create basics of dicom
    dicomtemplate = pydi.read_file(dicomfile)
    fibredcm = dicom_fibre_skel()
    fibredcm = dicom_patient_stuff(fibredcm, dicomtemplate)
    fibredcm = dicom_date_stamps(fibredcm, tckfile)
    fibredcm.file_meta = mk_filemeta_streamlines()
    fibredcm.is_little_endian = True
    fibredcm.is_implicit_VR = False
    fibredcm.SOPInstanceUID = \
        fibredcm.file_meta.MediaStorageSOPInstanceUID
    fibredcm.SOPClassUID = fibredcm.file_meta.MediaStorageSOPClassUID
    fibredcm.SeriesNumber = 777 + seriesNum
    fibredcm.ContentDescription = b''
    if StudyUID is not None:
        fibredcm.StudyInstanceUID = StudyUID
    if UIDlist is not None:
        # Add the references to a dicom volume
        fibredcm.ReferencedSeriesSequence = \
            dicom_referenced_series_sequence(
                UIDlist,
                dicomtemplate.SeriesInstanceUID)
    if FrameUID is not None:
        fibredcm.FrameOfReferenceUID = FrameUID
    if Description is not None:
        fibredcm.SeriesDescription = Description
    # Now turn the surface stuff into dicom format
    u = [
        mk_surface_sequence(tck_as_surf['sl'][M], tck_as_surf['pl'][M], M)
        for M in range(len(tck_as_surf['sl']))]
    fibredcm.SurfaceSequence = pydi.Sequence(u)
    fibredcm.SegmentSequence = pydi.Sequence(
        [mk_segment_sequence(ContDesc, tckname, UIDlist, len(u))])
    fibredcm.NumberOfSurfaces = len(u)

    pydi.filewriter.dcmwrite(outputfile, fibredcm,
                             write_like_original=False)
    return fibredcm


def lookup_cie(labnum):
    global nice_colours_cie
    if labnum > (len(nice_colours_cie) - 1):
        labnum = len(nice_colours_cie) - 1
        print("Error - too many labels")
    return list(nice_colours_cie[labnum])


def process_label_im(im):
    """Create a binary separate image for each label
    cropped in a special way. I can't see how to encode
    a separate image space for each label. Brainlab stores
    a single segmentation per file. We'll do the same.

    Note that Brainlab can read multi-segmentations.

    Returns a dictionary of images and corresponding
    labels (for choosing colours), also a scene bounding
    box.
    """
    # stuff to figure out which way we slice, etc
    isoidx = check_isotropy(im)
    otheridx = [0, 1, 2]
    otheridx.remove(isoidx)

    # direction = get_direction(im, isoidx)
    # sp = im.GetSpacing()
    # sp = str2ds(sp)
    # spacing = [sp[i] for i in otheridx]
    # slthickness = sp[isoidx]

    labstats = sitk.LabelShapeStatisticsImageFilter()
    labstats.Execute(im)
    labels = labstats.GetLabels()
    boxes = [labstats.GetBoundingBox(i) for i in labels]
    # Need to compute bounding box for all labels, as
    # this will set the row/colums
    # boxes are corner and size - this code assumes 3D
    corners = [(x[0], x[1], x[2]) for x in boxes]
    sizes = [(x[3], x[4], x[5]) for x in boxes]

    newcorners = [list(x) for x in corners]
    newsizes = [list(x) for x in sizes]

    ims = [sitk.RegionOfInterest(im, newsizes[i],
                                 newcorners[i]) == labels[i]
           for i in range(len(labels))]

    return({'rois': ims, 'labels': labels,
            'corners': corners,
            'sizes': sizes,
            'original': im})


def count_total_frames(perlabelstuff, isoidx):
    totalframes = 0
    for lab in range(len(perlabelstuff['rois'])):
        sz = perlabelstuff['rois'][lab].GetSize()
        totalframes += sz[isoidx]
    return totalframes


def process_label_imA(im):
    """Crop a label image so that the result contains
    all labels, then return separate images, one for
    each label.
    Returns a dictionary of images and corresponding
    labels (for choosing colours), also a scene bounding
    box. Need to run shape statistics to determine
    the number of labels and the IDs
    """
    # stuff to figure out which way we slice, etc
    isoidx = check_isotropy(im)
    otheridx = [0, 1, 2]
    otheridx.remove(isoidx)

    # direction = get_direction(im, isoidx)
    # sp = im.GetSpacing()
    # sp = str2ds(sp)
    # spacing = [sp[i] for i in otheridx]
    # slthickness = sp[isoidx]

    labstats = sitk.LabelShapeStatisticsImageFilter()
    labstats.Execute(im)
    labels = labstats.GetLabels()
    boxes = [labstats.GetBoundingBox(i) for i in labels]
    # Need to compute bounding box for all labels, as
    # this will set the row/colums
    # boxes are corner and size - this code assumes 3D
    corners = [(x[0], x[1], x[2]) for x in boxes]
    othercorner = [(x[0] + x[3] - 1,
                    x[1] + x[4] - 1,
                    x[2] + x[5] - 1) for x in boxes]
    sizes = [(x[3], x[4], x[5]) for x in boxes]

    all_low_x = [C[0] for C in corners]
    all_low_y = [C[1] for C in corners]
    all_low_z = [C[2] for C in corners]

    low_x = min(all_low_x)
    low_y = min(all_low_y)
    low_z = min(all_low_z)
    lowcorner = (low_x, low_y, low_z)

    all_high_x = [C[0] for C in othercorner]
    all_high_y = [C[1] for C in othercorner]
    all_high_z = [C[2] for C in othercorner]

    high_x = max(all_high_x)
    high_y = max(all_high_y)
    high_z = max(all_high_z)
    highcorner = (high_x, high_y, high_z)

    allsize = (highcorner[0] - lowcorner[0] + 1,
               highcorner[1] - lowcorner[1] + 1,
               highcorner[2] - lowcorner[2] + 1)

    # corners [otheridx] and size[otheridx] should be all the same
    newcorners = [list(x) for x in corners]
    newsizes = [list(x) for x in sizes]

    a = otheridx[0]
    b = otheridx[1]
    for f in range(len(newcorners)):
        newcorners[f][a] = lowcorner[a]
        newcorners[f][b] = lowcorner[b]
        newsizes[f][a] = allsize[a]
        newsizes[f][b] = allsize[b]

    ims = [sitk.RegionOfInterest(im, allsize,
                                 lowcorner) == labels[i]
           for i in range(len(labels))]
    imcrop = sitk.RegionOfInterest(im, allsize, lowcorner)
    return({'rois': ims, 'labels': labels,
            'original': im, 'cropped': imcrop})


def mk_shared_functional_group(nif, SOPList):
    """Create the orientation information that
    goes in the SharedFunctionalGroups sequence
    """
    # set up image geometry details
    isoidx = check_isotropy(nif)
    otheridx = [0, 1, 2]
    otheridx.remove(isoidx)

    direction = get_direction(nif, isoidx)
    direction = str2ds(direction)
    sp = nif.GetSpacing()
    sp = str2ds(sp)
    spacing = [sp[i] for i in otheridx]
    slthickness = sp[isoidx]

    sfgs = pydi.Sequence()
    sfg1 = pydi.Dataset()

    di = pydi.Dataset()
    di.DerivationDescription = 'Segmentation'

    sis = pydi.Sequence()

    def mk_sis(SOPInst):
        si = pydi.Dataset()
        si.ReferencedSOPClassUID = storage_sopclass.MRImageStorage
        si.ReferencedSOPInstanceUID = SOPInst
        ps = pydi.Sequence()
        purpose = pydi.Dataset()
        purpose.CodeValue = '121322'
        purpose.CodingSchemeDesignator = 'DCM'
        purpose.CodeMeaning = 'Source Image for Image Processing Operation'
        ps.append(purpose)
        si.PurposeOfReferenceCodeSequence = ps
        return si

    for g in range(len(SOPList)):
        sis.append(mk_sis(SOPList[g]))

    di.SourceImageSequence = sis

    dcs = pydi.Sequence()
    dc = pydi.Dataset()
    dc.CodeValue = '113076'
    dc.CodingSchemeDesignator = 'DCM'
    dc.CodeMeaning = 'Segmentation'
    dcs.append(dc)

    di.DerivationCodeSequence = dcs
    dis = pydi.Sequence()
    dis.append(di)
    sfg1.DerivationImageSequence = dis

    pos = pydi.Sequence()
    po1 = pydi.Dataset()
    po1.ImageOrientationPatient = direction
    pos.append(po1)
    pms = pydi.Sequence()
    pm1 = pydi.Dataset()
    pm1.SliceThickness = slthickness
    pm1.PixelSpacing = spacing
    pms.append(pm1)
    sfg1.PlaneOrientationSequence = pos
    sfg1.PixelMeasuresSequence = pms
    sfgs.append(sfg1)
    return(sfgs)


def mk_perframe_functional_group(perlabelstuff, UIDlist):
    """key section describing the per label, per slice
    image data that is later run length encoded
    :perlabelstuff: is the subsetted images/bounding boxes
    """

    # Sort out the isotropic slice direction.
    # This repeats code from other places, but
    # saves passing lots of data around
    # Use the original data, not subsets.
    nif = perlabelstuff['original']
    isoidx = check_isotropy(nif)
    otheridx = [0, 1, 2]
    otheridx.remove(isoidx)

    pffgs = pydi.Sequence()
    # one frame content sequence per slice, per label
    lablist = perlabelstuff['labels']

    for labelidx in range(len(lablist)):
        # label_id = lablist[labelidx]

        this_roi = perlabelstuff['rois'][labelidx]
        this_size = this_roi.GetSize()
        this_slices = this_size[isoidx]
        thiscropped = perlabelstuff["rois"][labelidx]
        for thisslice in range(this_slices):
            pffg1 = pydi.Dataset()
            corner = [0, 0, 0]
            corner[isoidx] = thisslice
            origin = thiscropped.TransformIndexToPhysicalPoint(corner)
            fcs = pydi.Sequence()
            fc1 = pydi.Dataset()
            fc1.DimensionIndexValues = [labelidx + 1, thisslice + 1]
            fcs.append(fc1)
            pps = pydi.Sequence()
            pp1 = pydi.Dataset()
            origin = str2ds(origin)
            pp1.ImagePositionPatient = origin

            pps.append(pp1)

            sis = pydi.Sequence()
            si1 = pydi.Dataset()
            # only ever one segment per file
            si1.ReferencedSegmentNumber = labelidx + 1
            sis.append(si1)
            pffg1.PlanePositionSequence = pps
            pffg1.FrameContentSequence = fcs
            pffg1.SegmentIdentificationSequence = sis
            pffgs.append(pffg1)

    return(pffgs)


def mk_rle_data(perlabelstuff):
    """
    Do the run length encoding and encapsulation of
    all labels. Note that the first byte array is
    the Basic Offset Table Item, which appears
    to have a length that depends on the total number
    of frames. Can be confusing if the first one has
    a very different length to the rest.
    """
    nif = perlabelstuff['original']
    isoidx = check_isotropy(nif)
    otheridx = [0, 1, 2]
    otheridx.remove(isoidx)
    labelframe = list()

    for labidx in range(len(perlabelstuff["rois"])):
        roi = perlabelstuff["rois"][labidx]
        sz = roi.GetSize()
        for slce in range(sz[isoidx]):
            selector = mk_indexing_tuple(slce, isoidx)
            roislice = roi[selector]
            slicedat = sitk.GetArrayFromImage(roislice)
            rledat = (pydi.pixel_data_handlers.
                      rle_handler.rle_encode_frame(slicedat))

            labelframe.append(rledat)

    return(pydi.encaps.encapsulate(labelframe,
                                   fragments_per_frame=1, has_bot=True))


def find_match_im(labelfile, nidetails):
    """ Figure out which image the label mask
    was derived from.
    """
    nif = sitk.ReadImage(labelfile, sitk.sitkUInt8)
    qformname = nif.GetMetaData('qform_code_name')
    if qformname == 'NIFTI_XFORM_UNKNOWN':
        raise ValueError(labelfile + ': Unknown qform code - stopping')

    spacing = nif.GetSpacing()
    oMatrix = nif.GetDirection()
    imsize = nif.GetSize()

    def within_tol(t1, t2):
        t1a = np.array(t1)
        t2a = np.array(t2)
        d = np.abs(t1a - t2a)
        return d.max() < 0.0001

    for i in range(len(nidetails)):
        sp = nidetails[i]["spacing"]
        sz = nidetails[i]["size"]
        mt = nidetails[i]["matrix"]
        if (within_tol(spacing, sp) and
           within_tol(imsize, sz) and
           within_tol(oMatrix, mt)):
            return i
    return None


def mk_label_segment_sequence(imname, labs):
    """
    Hardcoded stuff indicating that the mask/labels are imported
    """
    segment_sequence = pydi.Sequence()

    for labnum in range(len(labs)):
        # Segment Sequence: Segment 1
        seg1 = pydi.Dataset()

        # Anatomic Region Sequence
        anatomic_region_sequence = pydi.Sequence()
        seg1.AnatomicRegionSequence = anatomic_region_sequence

        # Anatomic Region Sequence: Anatomic Region 1
        anatomic_region1 = pydi.Dataset()
        anatomic_region1.CodeValue = 'T-D0010'
        anatomic_region1.CodingSchemeDesignator = 'SRT'
        anatomic_region1.CodeMeaning = 'Entire body'
        anatomic_region_sequence.append(anatomic_region1)

        # Segmented Property Category Code Sequence
        seg_property_category_code_sequence = pydi.Sequence()
        seg1.SegmentedPropertyCategoryCodeSequence = (
            seg_property_category_code_sequence)

        # Segmented Property Category Code Sequence:
        # Segmented Property Category Code 1
        seg_property_category_code1 = pydi.Dataset()
        seg_property_category_code1.CodeValue = 'T-D000A'
        seg_property_category_code1.CodingSchemeDesignator = 'SRT'
        seg_property_category_code1.CodeMeaning = 'Anatomical Structure'
        seg_property_category_code_sequence.append(seg_property_category_code1)

        seg1.SegmentNumber = labnum + 1
        seg1.SegmentLabel = imname + " label " + str(labs[labnum])
        seg1.SegmentDescription = ("Label " +
                                   str(labs[labnum]) +
                                   " of " + imname)
        seg1.SegmentAlgorithmType = 'AUTOMATIC'
        seg1.SegmentAlgorithmName = 'Unknown'
        seg1.RecommendedDisplayCIELabValue = lookup_cie(labs[labnum])
        # Segmented Property Type Code Sequence
        segmented_property_type_code_sequence = pydi.Sequence()
        seg1.SegmentedPropertyTypeCodeSequence = (
            segmented_property_type_code_sequence)

        # Segmented Property Type Code Sequence: Segmented Property Type Code 1
        segmented_property_type_code1 = pydi.Dataset()
        segmented_property_type_code1.CodeValue = '111176'
        segmented_property_type_code1.CodingSchemeDesignator = 'DCM'
        segmented_property_type_code1.CodeMeaning = 'Unspecified'
        segmented_property_type_code_sequence.append(
            segmented_property_type_code1)
        segment_sequence.append(seg1)

    return segment_sequence


def sitk_labelnifti_to_dicom(niftifile, dicomfile,
                             outputfile, seriesNum=0,
                             Description=None, StudyUID=None, FrameUID=None,
                             UIDlist=None, SeriesNum=None):
    """
    :param niftifile: (string) path to nifti image - this is expected
                         to be a label image. Each label will be rendered
                         in a different colour in brainlab. Binary
                         images (masks) are also fine.
    :param dicomfile: (string) path to a dicom file - used to supply some
                   important tags
    :param outputprefix: (string) output filename.
    :param Description: (string) to populate the SeriesDescription field.
    :param StudyUID: (string) Study Instance UID - will generate one
                          if none.
    :param FrameUID: (string) a UID that should be common for all
                 dicoms in a common space
    :param SeriesNum: (int) unique within a study - helps identify
                  files in a series
    :param UIDlist: the set of dicom UIDs corresponding to the volume
                    used to create the segmentation
    :return: tuple of list of strings containing InstanceUIDs for
        constructing other dicoms and the dicom
        structure containing the common parts.
    Notes :
    This function creates the special dicom format used by brainlab "objects".
    Useful for rendering masks, such as tumour segmentations,
    functional activation, etc.
    """
    # generate a useful name
    imname = os.path.basename(niftifile)
    imname = os.path.splitext(imname)[0]
    # again in case nifti files are compressed
    imname = os.path.splitext(imname)[0]

    # Read in as float and rescale to UInt16, with rescale values
    nif = sitk.ReadImage(niftifile, sitk.sitkUInt8)
    qformname = nif.GetMetaData('qform_code_name')
    if qformname == 'NIFTI_XFORM_UNKNOWN':
        raise ValueError(niftifile + ': Unknown qform code - stopping')

    spacing = nif.GetSpacing()
    # oMatrix = nif.GetDirection()
    # imsize = nif.GetSize()

    if len(spacing) > 3:
        raise ValueError('Higher than 3D nifti file - stopping')

    try:
        perlabelstuff = process_label_imA(nif)

        # figure out which plane to write. Aiming for isotropic within plane
        # If image is isotropic, should use minimum number of planes
        isoidx = check_isotropy(nif)
    except ValueError:
        print("Error processing " + niftifile)
        raise
    otheridx = [0, 1, 2]
    otheridx.remove(isoidx)

    # modification_time = time.strftime("%H%M%S", time.localtime())
    # modification_date = time.strftime("%Y%m%d", time.localtime())

    # Columns of this matrix contain the direction cosines
    # matrix is stored in rows
    direction = get_direction(nif, isoidx)
    direction = str2ds(direction)
    sp = nif.GetSpacing()
    sp = str2ds(sp)
    spacing = [sp[i] for i in otheridx]
    # slthickness = sp[isoidx]

    # Load the sample dicom
    # create basics of dicom
    dicomtemplate = pydi.read_file(dicomfile)

    labeldcm = dicom_label_skel()
    labeldcm = dicom_patient_stuff(labeldcm, dicomtemplate)
    labeldcm = dicom_date_stamps(labeldcm, niftifile)
    labeldcm.file_meta = mk_filemeta_labelobj()
    labeldcm.is_little_endian = True
    labeldcm.is_implicit_VR = False
    labeldcm.SOPInstanceUID = \
        labeldcm.file_meta.MediaStorageSOPInstanceUID
    labeldcm.SOPClassUID = labeldcm.file_meta.MediaStorageSOPClassUID
    labeldcm.SeriesNumber = 877 + seriesNum
    labeldcm.ImageType = "DERIVED\\PRIMARY"
    labeldcm.SegmentationFractionalType = "PROBABILITY"
    labeldcm.MaximumFractionalValue = 255
    labeldcm.ContentDescription = b'Nifti segmentation objects'
    if StudyUID is not None:
        labeldcm.StudyInstanceUID = StudyUID
    if UIDlist is not None:
        # Add the references to a dicom volume
        labeldcm.ReferencedSeriesSequence = \
            dicom_referenced_series_sequence(
                UIDlist,
                dicomtemplate.SeriesInstanceUID)
    else:
        raise MissingUIDList

    if FrameUID is not None:
        labeldcm.FrameOfReferenceUID = FrameUID
    if Description is not None:
        labeldcm.SeriesDescription = Description

    labeldcm.DimensionOrganizationType = '3D'
    labeldcm.SamplesPerPixel = 1
    labeldcm.PhotometricInterpretation = 'MONOCHROME2'

    cropsize = perlabelstuff["cropped"].GetSize()

    totalframes = count_total_frames(perlabelstuff, isoidx)
    labeldcm.NumberOfFrames = totalframes
    labeldcm.Rows = int(cropsize[otheridx[1]])
    labeldcm.Columns = int(cropsize[otheridx[0]])
    labeldcm.BitsAllocated = 8
    labeldcm.BitsStored = 8
    labeldcm.HighBit = 7
    labeldcm.PixelRepresentation = 0
    labeldcm.LossyImageCompression = '00'
    labeldcm.SegmentationType = 'FRACTIONAL'

    # Dimension organisation
    # - uses macros (0062,0004) - segment number
    #               (0062,0002) - segment sequence
    #                   (tag for a sequence that comes next)
    #               (0020,0032) - image position patient
    #               (0020,9113) - plane position sequence
    dos = pydi.Sequence()
    do1 = pydi.Dataset()
    do1.DimensionOrganizationUID = dcm_uuid()
    dos.append(do1)
    labeldcm.DimensionOrganizationSequence = dos

    dis = pydi.Sequence()
    di1 = pydi.Dataset()
    di1.DimensionOrganizationUID = do1.DimensionOrganizationUID
    di1.DimensionIndexPointer = pydi.tag.Tag(0x0062, 0x0004)
    di1.FunctionalGroupPointer = pydi.tag.Tag(0x0062, 0x0002)
    dis.append(di1)

    di2 = pydi.Dataset()
    di2.DimensionOrganizationUID = do1.DimensionOrganizationUID
    di2.DimensionIndexPointer = pydi.tag.Tag(0x0020, 0x0032)
    di2.FunctionalGroupPointer = pydi.tag.Tag(0x0020, 0x9113)
    dis.append(di2)

    labeldcm.DimensionIndexSequence = dis
    labeldcm.DimensionOrganizationType = "3D"

    # Segmentation sequence
    labeldcm.SegmentSequence = (
        mk_label_segment_sequence(imname, perlabelstuff["labels"]))

    # Shared Functional Groups Sequence - do we need this
    # Brainlab and slicer version have an image orientation patient inside it

    labeldcm.SharedFunctionalGroupsSequence = \
        mk_shared_functional_group(nif, UIDlist)

    # Per-frame functional groups sequence - contains segment number,
    # frame number and position. This is where a lot of the action is.
    labeldcm.PerFrameFunctionalGroupsSequence = (
        mk_perframe_functional_group(perlabelstuff, UIDlist))
    # Pixels at the end - encapsulated form

    labeldcm.PixelData = mk_rle_data(perlabelstuff)
    labeldcm["PixelData"].VR = 'OB'
    labeldcm["PixelData"].is_undefined_length = True

    pydi.filewriter.dcmwrite(outputfile, labeldcm,
                             write_like_original=False)


########################################################################
# Driver scripts to import collections of nifti and tract files


def import_tractography_study(origdcm, niftifiles,
                              tckfiles, labelfiles=None,
                              destdir="./",
                              StudyUID=None, FrameUID=None):
    """
    Create a dicom study for brainlab from nifti and tck files.
    The first nifti file will be used as the reference.
    The files are assumed to be coregistered and a common
    frame of reference (a dicom structure) is created for them.
    Brainlab requires you to accept the coregistration.
    Patient details from origdcm.
    :param origdcm: a dicom used as a template, providing patient
    details etc
    :param niftifiles: a list of nifti filenames
    :param tckfiles: a list of track filenames
    :param destdir: a target directory
    :param StudyUID: a study id. New one will be created if none.
    :param FrameUID: a frame of reference id. New one created
    if none.
    :return:
    """

    # Generate a studyID
    if StudyUID is None:
        StudyUID = dcm_uuid()
    # Generate a FrameUID
    if FrameUID is None:
        FrameUID = dcm_uuid()
    # Generate output names for nifti
    n_bn = [os.path.basename(x) for x in niftifiles]
    n_cn = [os.path.splitext(x)[0] for x in n_bn]
    # again in case nifti files are compressed
    n_cn = [os.path.splitext(x)[0] for x in n_cn]

    # directory names
    n_dir = [os.path.join(destdir, x) for x in n_cn]
    # paste series number?

    nidetails = [
        sitk_nifti_to_dicom(niftifile=niftifiles[idx], dicomfile=origdcm,
                            dcmprefix="IM",
                            outdir=n_dir[idx], Description=n_cn[idx],
                            StudyUID=StudyUID, FrameUID=FrameUID,
                            SeriesNum=idx + 1) for idx in
        range(len(niftifiles))]

    # and for tck files
    t_dir = None
    if tckfiles is not None:
        t_bn = [os.path.basename(x) for x in tckfiles]
        t_cn = [os.path.splitext(x)[0] for x in t_bn]
        t_dir = [os.path.join(destdir, x) for x in t_cn]

        [os.makedirs(x, exist_ok=True) for x in t_dir]

        t_dir = [os.path.join(x, "FT_00.dcm") for x in t_dir]

        tckdetails = [tck_to_dicom(tckfile=tckfiles[idx],
                                   dicomfile=nidetails[0]['dcmfiles'][0],
                                   outputfile=t_dir[idx],
                                   seriesNum=idx + len(nidetails) + 10,
                                   Description=t_cn[idx],
                                   StudyUID=StudyUID,
                                   UIDlist=nidetails[0]['SOPlist'],
                                   FrameUID=FrameUID) for idx in
                      range(len(tckfiles))]

    if labelfiles is not None:
        # now for label images
        ln_bn = [os.path.basename(x) for x in labelfiles]
        ln_cn = [os.path.splitext(x)[0] for x in ln_bn]
        ln_cn = [os.path.splitext(x)[0] for x in ln_cn]

        ln_dir = [os.path.join(destdir, x) for x in ln_cn]

        [os.makedirs(x, exist_ok=True) for x in ln_dir]
        ln_dir = [os.path.join(x, "LB.dcm") for x in ln_dir]

        # figure out which image that we've already converted
        # matches the label image
        # We end up loading the niftifile again.
        for idx in range(len(labelfiles)):
            nif_index = find_match_im(labelfiles[idx], nidetails)

            if nif_index is None:
                raise RawToLabelImMismatch

            sitk_labelnifti_to_dicom(
                labelfiles[idx],
                dicomfile=nidetails[nif_index]['dcmfiles'][0],
                outputfile=ln_dir[idx],
                Description=ln_cn[idx],
                StudyUID=StudyUID,
                UIDlist=nidetails[nif_index]['SOPlist'],
                FrameUID=FrameUID)

    return [n_dir, t_dir]


def get_already_converted_info(origdcmfolder):
    """
    Helper function for appending data to an already converted
    dataset. This function pulls out the frame of reference, study,
    SOP instance for
    # each dicom, series instance
    :param origdcmfolder: Location of converted data
    :return: list of dicom details.
    """
    # order of files matters if we want to match brainlab.
    dcms = os.listdir(origdcmfolder)
    dcms.sort()
    dcms = [os.path.join(origdcmfolder, f) for f in dcms]

    def dcdetails(f):
        db = pydi.read_file(f)
        return {'FrameUID': db.FrameOfReferenceUID,
                'StudyUID': db.StudyInstanceUID,
                'SeriesUID': db.SeriesInstanceUID,
                'SOPInstance': db.SOPInstanceUID,
                'InstanceNumber': db.InstanceNumber}

    alldetails = [dcdetails(x) for x in dcms]
    return alldetails, dcms


def append_imaging_study(origdcmfolder, niftifiles, studystart,
                         destdir="./"):
    """
    Convert more niftis, with FrameID etc coming from
    a previous conversion


    :param origdcmfolder: location of previously converted data
    :param niftifiles: new files to convert
    :param studystart: used for setting the studynumber
    :param destdir: target folder
    :return: list of details of conversion.
    """
    n_bn = [os.path.basename(x) for x in niftifiles]
    n_cn = [os.path.splitext(x)[0] for x in n_bn]
    n_cn = [os.path.splitext(x)[0] for x in n_cn]

    n_dir = [os.path.join(destdir, x) for x in n_cn]

    [os.makedirs(x, exist_ok=True) for x in n_dir]
    alldetails, dcms = get_already_converted_info(origdcmfolder)
    StudyUID = alldetails[0]['StudyUID']
    FrameUID = alldetails[0]['FrameUID']
    # SOPlist = [x['SOPInstance'] for x in alldetails]
    nidetails = [
        sitk_nifti_to_dicom(niftifile=niftifiles[idx], dicomfile=dcms[0],
                            dcmprefix="IM",
                            outdir=n_dir[idx], Description=n_cn[idx],
                            StudyUID=StudyUID, FrameUID=FrameUID,
                            SeriesNum=idx + studystart) for idx in
        range(len(niftifiles))]


def append_tractography_study(origdcmfolder, tckfiles, destdir="./"):
    """
     Put a tck file into an existing dicom study -
    various IDs retrieved from the dicom
    :param origdcmfolder: location of previously converted data
    :param tckfiles: new files to convert
    :param destdir: target folder
    :return: list of pydicom datasets
    """

    t_bn = [os.path.basename(x) for x in tckfiles]
    t_cn = [os.path.splitext(x)[0] for x in t_bn]
    t_dir = [os.path.join(destdir, x) for x in t_cn]

    [os.makedirs(x, exist_ok=True) for x in t_dir]

    t_dir = [os.path.join(x, "FT_00.dcm") for x in t_dir]

    alldetails, dcms = get_already_converted_info(origdcmfolder)

    SOPlist = [x['SOPInstance'] for x in alldetails]
    tckdetails = [tck_to_dicom(tckfile=tckfiles[idx], dicomfile=dcms[0],
                               outputfile=t_dir[idx],
                               seriesNum=idx + 50 + 10,
                               Description="Fiber Bundles",
                               StudyUID=alldetails[0]['StudyUID'],
                               UIDlist=SOPlist,
                               FrameUID=alldetails[0]['FrameUID']) for idx
                  in range(len(tckfiles))]
    print(alldetails)
    return tckdetails


def fix_dwi_shell(origdcmfolder, desiredB, destdir="./", b0thresh=150):
    """
    A test function to modify modern diffusion data with slightly
    varying bvalues to fool brainlab into thinking it is single
    shell. Any bval less than b0thresh is set to 0.
       All others are set to desiredB.
    :param origdcmfolder: multishell, problem dicoms
    :param desiredB: the bvalue we're telling brainlab to use
    :param destdir: output dicom folder
    :param b0thresh: threshold for setting to zero
    :return: list of modified pydicom datasets
    """

    def saveDCMlist(dcmoutfolder, dcmlist):
        os.makedirs(dcmoutfolder, exist_ok=True)
        for x in range(len(dcmlist)):
            NM = "IM" + str(x + 1).zfill(4) + ".dcm"
            NM = os.path.join(dcmoutfolder, NM)
            pydi.filewriter.dcmwrite(NM, dcmlist[x],
                                     write_like_original=False)

    def setB(dcm):
        thisB = int(dcm['0019', '100c'].value)
        newB = 0
        if thisB > b0thresh:
            newB = desiredB
        dcm[0x0019, 0x100c].value = newB
        return dcm

    def changeUIDs(dcm, StudyUID, SeriesUID):
        dcm.file_meta = mk_file_meta()
        dcm.SOPInstanceUID = dcm.file_meta.MediaStorageSOPInstanceUID
        dcm.SOPClassUID = dcm.file_meta.MediaStorageSOPClassUID
        dcm.StudyInstanceUID = StudyUID
        dcm.SeriesInstanceUID = SeriesUID
        return dcm

    fls = glob.glob(os.path.join(origdcmfolder, "*.dcm"))
    fls.sort()
    dd = [pydi.read_file(x) for x in fls]
    bv = [setB(x) for x in dd]
    # Need to change the IDs so that there's no clash
    # with existing data
    # StudyID, SeriesUID
    StudyUID = dcm_uuid()
    SeriesUID = dcm_uuid()
    bv = [changeUIDs(x, StudyUID, SeriesUID) for x in bv]
    saveDCMlist(destdir, bv)
    return bv
