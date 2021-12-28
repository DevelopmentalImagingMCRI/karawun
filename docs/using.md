# Using Karawun

*Karawun* provides a single command-line interface, as follows:

```bash
importTractography --dicom-template path/to/a/dicom --nifti T1.nii.gz fa.nii.gz --tract-files left_cst.tck right_cst.tck --output-dir path/to/output/folder
```

or, with label images

```bash
importTractography --dicom-template path/to/a/dicom --nifti T1.nii.gz fa.nii.gz --tract-files left_cst.tck right_cst.tck --label-files lesion.nii.gz white_matter.nii.gz  --output-dir path/to/output/folder
```

## Explanation

`--dicom-template` (required) is a dicom file from which dicom tags are copied,
including all patient information. Ideally this dicom file should be
one of the originals from which the nifti files were derived. All
patient, physician and institutional details are copied from this
template. A single dicom template file is required. A single,
anonymised, dicom is provided with the testing data (see below), that
may be appropriate for other testing purposes, if you do not have
original dicoms available.

`--nifti` (required) are the volume images that will be converted to dicom
images. They are __assumed to be co-registered__, as is typically the
the case for multiple nifti files used in a common image processing
pipeline. The conversion process creates a dicom "frame of reference" to
which all volume images belong. Brainlab requires that the user accept
this registration. Doing so allows overlay of uploaded data in arbitary
combinations.

The names of the nifti files (without suffixes) will be used in the DICOM SeriesDescription tag.

`--tract-files` (optional) are the MRtrix3 .tck files that need to be
displayed in Brainlab. The conversion of these files
creates the DICOM fibre object format that can be viewed as a 3D
object in Brainlab.

`--label-files` (optional) are volume images, also in nifti format, that contain
masks or *label images*. Label images use a constant integer value
to delineate a region. Multiple, non overlapping, regions can
be stored in a nifti volume using different integer values for each
region. *Karawun* currently supports 30 distinct colours, so ensure
that label files only contain values of 30 or below.

The label files must match one of the volume files in terms of voxel spacing
resolution, etc - i.e must have been derived from one of the
images in the `--nifti-files` list.

Label images are explained in more detail in the examples.

The output DICOM files are written to the folder specified by `--output-dir`.
A subfolder is created for each input volume, tract file and label file.
The input file name is used to generate the folder name. For example, a nifti
file named `FA.nii.gz` will be converted to DICOM files named IM_????.dcm in
a folder named `FA`.

A tractfile named `CST.tck` will be converted to `CST/FT_00.dcm` while a label
image named `tumor.nii.gz` will be converted to `tumor/LB.dcm`.


