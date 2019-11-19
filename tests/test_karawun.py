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

import pytest
import uuid

from unittest import mock
import glob
import karawun
import time
import tempfile
import os
import hashlib
import json

import karawun.karawun

# pytest tests/

# to create the baseline sha object
# pytest --runcreatebaseline -k 'not test_conv1 and not test_conv2' tests/

# Patched timestamp and uuid functions
# for testing


def patchLocalTime(seconds=None):
    return time.struct_time((2019, 8, 1, 12, 5, 20, 0, 0, 0))


class UUIDClass:
    _uuidNum = 0

    def patchdcm_uuid(self):
        roottestuuid = uuid.UUID('8e35875e-b40c-11e9-b0fb-a4c3f0b51dd7')
        thisuuid = uuid.uuid5(roottestuuid, str(self._uuidNum))
        self._uuidNum = self._uuidNum + 1

        uid = '2.25.' + str(thisuuid.int)
        return uid


# sha512 generation tools
def get_sha512(file_path, blocksize=512*16):
    sha = hashlib.sha512()
    with open(file_path, "rb") as f:
        while True:
            chunk = f.read(blocksize)
            if not chunk:
                break
            sha.update(chunk)
    return sha.hexdigest()


def get_all_sha512(location):
    shadict = {}
    for root, subdirs, files in os.walk(location):
        if len(files) > 0:
            thesefiles = [os.path.join(root, f) for f in files]
            shalist = [get_sha512(f) for f in thesefiles]
            # Remove the location component
            relfiles = [os.path.relpath(f, location) for f in thesefiles]
            z = zip(relfiles, shalist)
            shadict.update(dict(z))

    return(shadict)


# test routines
def converter(targetd):
    test_data = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), 'Data')
    t1dcm = os.path.join(test_data, "Dicom",
                         ("1.3.12.2.1107.5.2.43."
                          "167031.2019040213095021814319052.dcm"))
    MR = ['T1brain.nii.gz', 'FLAIRbrain.nii.gz',
          'FAbrain.nii.gz', 'FAbrain_reoriented.nii.gz']
    MR = [os.path.join(test_data, "Tractography", nii) for nii in MR]
    TCK = ['Left_PT_final.tck', 'Right_PT_final.tck']
    TCK = [os.path.join(test_data, "Tractography", tck) for tck in TCK]
    karawun.import_tractography_study(origdcm=t1dcm,
                                      niftifiles=MR,
                                      tckfiles=TCK,
                                      labelfiles=None,
                                      destdir=targetd)

    # test routines
def label_converter(targetd):
    test_data = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), 'Data')
    t1dcm = os.path.join(test_data, "Dicom",
                         ("1.3.12.2.1107.5.2.43."
                          "167031.2019040213095021814319052.dcm"))
    MR = ['T1brain.nii.gz', 'FLAIRbrain.nii.gz',
          'FAbrain.nii.gz']
    MR = [os.path.join(test_data, "Tractography", nii) for nii in MR]
    MR_LAB = ['words.nii.gz', 'globes.nii.gz']
    MR_LAB = [os.path.join(test_data, "Tractography", nii) for nii in MR_LAB]
    
    TCK = ['Left_PT_final.tck', 'Right_PT_final.tck']
    TCK = [os.path.join(test_data, "Tractography", tck) for tck in TCK]
    karawun.import_tractography_study(origdcm=t1dcm,
                                      niftifiles=MR,
                                      tckfiles=TCK,
                                      labelfiles=MR_LAB,
                                      destdir=targetd)


@mock.patch('time.localtime', side_effect=patchLocalTime)
@mock.patch('karawun.karawun.dcm_uuid', side_effect=UUIDClass().patchdcm_uuid)
def test_conv1(dcm_uuid, localtime, tmp_path):
    baselineshaf = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), 'Baseline',
        "image_sha512.json")

    assert os.path.exists(baselineshaf),\
        ("Missing json file - run pytest --runcreatebaseline"
         "-k 'not test_conv1'")
    print("Destination folder = " + str(tmp_path))
    converter(tmp_path)
    this_sha512 = get_all_sha512(tmp_path)
    f = open(baselineshaf)
    baselinesha = json.loads(f.read())
    f.close()
    assert(baselinesha == this_sha512)


@mock.patch('time.localtime', side_effect=patchLocalTime)
@mock.patch('karawun.karawun.dcm_uuid', side_effect=UUIDClass().patchdcm_uuid)
def test_conv2(dcm_uuid, localtime, tmp_path):
    baselineshaf = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), 'Baseline',
        "seg_sha512.json")

    assert os.path.exists(baselineshaf),\
        ("Missing json file - run pytest --runcreatebaseline"
         "-k 'not test_conv1'")
    print("Destination folder = " + str(tmp_path))
    label_converter(tmp_path)
    this_sha512 = get_all_sha512(tmp_path)
    f = open(baselineshaf)
    baselinesha = json.loads(f.read())
    f.close()
    assert(baselinesha == this_sha512)
    

def mkSha(pth, jso):
    baselined = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), 'Baseline')
    baseline_sha = get_all_sha512(pth)
    # save it (json to play more nicely with git)
    json_sha = json.dumps(baseline_sha)
    f = open(os.path.join(baselined, jso), "w")
    f.write(json_sha)
    f.close()


@pytest.mark.createbaseline
@mock.patch('time.localtime', side_effect=patchLocalTime)
@mock.patch('karawun.karawun.dcm_uuid', side_effect=UUIDClass().patchdcm_uuid)
def test_create_baseline_shaA(dcm_uuid, localtime, tmp_path):
    print("Creating baseline checksums for images+tractography")
    converter(tmp_path)
    mkSha(tmp_path, "image_sha512.json")
    

@pytest.mark.createbaseline
@mock.patch('time.localtime', side_effect=patchLocalTime)
@mock.patch('karawun.karawun.dcm_uuid', side_effect=UUIDClass().patchdcm_uuid)
def test_create_baseline_shaB(dcm_uuid, localtime, tmp_path):
    label_converter(tmp_path)
    mkSha(tmp_path, "seg_sha512.json")

    
