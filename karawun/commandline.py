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

import karawun
import argparse
import sys


def exception_handler(exception_type, exception, traceback):
    # All your trace are belong to us!
    # your format
    print("Exception handler: %s - %s" % (exception_type.__name__, exception))


def is_readable_file(parser, arg):
    try:
        f = open(arg, 'r')
        f.close()
    except Exception:
        raise argparse.ArgumentTypeError("{0} does not exist or is not readable".format(arg))

    return(arg)


parser = argparse.ArgumentParser(description="A tool for creating Brainlab compatible dicom data from "
                                             "mrtrix tract files and nifti images. All image and tck"
                                             "files are assumed to be aligned and will be placed into"
                                             "a dicom frame of reference with the first nifti image"
                                             "as the base for the reference. The alignment must be"
                                             "accepted in brainlab before objects can be viewed on"
                                             "all images.")

parser.add_argument("-d", "--dicom-template",
                    type=lambda x: is_readable_file(parser, x),
                    required=True,
                    help="Dicom file to supply tags, ideally from one "
                         "of the dicoms converted to nifti")

parser.add_argument("-o", "--output-dir", type=str, required=True,
                    help="Output directory for generated dicom")

parser.add_argument("-n", "--nifti", nargs="+",
                    type=lambda x: is_readable_file(parser, x),
                    required=True,
                    help="One or more input nifti files")

parser.add_argument("-t", "--tract-files", nargs="+",
                    type=lambda x: is_readable_file(parser, x),
                    required=False,
                    help="One or more mrtrix streamline files")

parser.add_argument("-l", "--label-files", nargs="+",
                    type=lambda x: is_readable_file(parser, x),
                    required=False,
                    help="One or more nifti label image files")

parser.add_argument("-v", "--verbose",
                    action='store_true',
                    required=False,
                    help="verbose errors, python traceback")


def run_cli(args):
    try:
        karawun.import_tractography_study(origdcm=args.dicom_template,
                                          niftifiles=args.nifti,
                                          tckfiles=args.tract_files,
                                          labelfiles=args.label_files,
                                          destdir=args.output_dir)
    except karawun.RawToLabelImMismatch:
        print("One of the label images is not derived from any of the raw images")
    except karawun.MissingUIDList:
        print("A UID list is required iternally somewhere - this error shouldn't happen")


def import_tractography_cl():
    args = parser.parse_args()
    if not args.verbose:
        sys.excepthook = exception_handler

    run_cli(args)
