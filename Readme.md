# Karawun

## Introduction

_Karawun_ converts the results of a diffusion imaging tractography study,
as might be used for surgical planning, to a form that is readable by
the Brainlab software suite.

The planning results will consist of a set of co-registered nifti images
(for example T1, FA, contrast enhanced T1 etc) and associated tract files
created using the [mrtrix3](https://www.mrtrix.org/) tools. _Karawun_ converts
the nifti images to dicom images and the tck files to the dicom object format
used by Brainlab to represent tracts produced by the internal tracking tool.

Tracts generated by mrtrix3 can thus be used in Brainlab in the same
way as tracts generated by Brainlab. This allows you to experiment
with your favourite diffusion sequence or post processing method and
visualize the results in the planning suite.

_Karawun_ is a grass used as a source of fibre by Aboriginal
tribes in south-eastern Australia.

## Quick Installation

_Karawun_ is a python package that can be installed using various
standard python options. It is recommended that some form of virtual environment
is used to isolate _karawun_ and its dependencies from other packages. Installation
with miniconda is illustrated below:

### Miniconda

1. [Install miniconda for python 3.7](https://docs.conda.io/en/latest/miniconda.html)

1. Add package channels:
    ```bash
    conda config --append channels conda-forge \
    --append channels anaconda \
    --append channels SimpleITK
    ```
1. Create a conda environment (using a name of your own choice):
    ```bash
    conda create --name KarawunEnv
    conda install --name KarawunEnv pip
    ```
1. Activate the environment:
    ```bash
    conda activate KarawunEnv
    ```
1. Install _karawun_ (internet access required)
    ```bash
    pip install git+https://github.com/DevelopmentalImagingMCRI/karawun.git@master#egg=karawun
    ```

See below for recommended testing procedure.

## Using
The simplest way of using _karawun_ is via the _importTractograpy_ command, as follows:

```bash
importTractography --dicom-template path/to/a/dicom --nifti T1.nii.gz fa.nii.gz --tract-files left_cst.tck right_cst.tck -o path/to/output/folder
```

### Explanation

The dicom template is a dicom file from which dicom tags are copied,
including all patient information. Ideally this dicom file should be
one of the originals from which the nifti files were derived. All
patient, physician and institutional details are copied from this
template. A single dicom template file is required. A single,
anonymised, dicom is provided with the testing data (see below), that
may be appropriate for other testing purposes, if you do not have
original dicoms available.

The nifti files are the volume images that will be converted to dicom
images. They are assumed to be co-registered, as would typically be
the case for multiple nifti files used in a common image processing
pipeline. The conversion process creates a dicom "frame of reference" to
which all volume images belong. Brainlab requires that the user accept
this registration.

The tract files are the mrtrix .tck files that need to be
displayed/manipulated in Brainlab. The conversion of these files
creates the dicom fibre object format that can be viewed as a 3D
object in Brainlab.

The converted dicoms are placed in the output folder. Separate
subfolders are created for each nifti file and each tck file.

The input file names are used in various dicom description tags.

The frame of reference is based on the first nifti file in the
argument list. Thus the Brainlab "quick viewer" will only allow tracts
to be overlaid on that image. Accepting the registration will allow all
display combinations.

## Typical workflow

A typical workflow involves extracting dicom images from the scanner
system, converting them to nifti formats and peforming image
processing steps using a range of approaches. Examples include
co-registration, tissue classification, brain extraction, activation
detection, distortion correction, etc with tools like FSL, FreeSurfer and
SPM as well as diffusion tractography using mrtrix. The resulting dataset is then
converted to dicom ready for import.

We do not attempt to make the new dicoms align with the original
dicoms from which the nifti were derived. Always treat the images
derived from the external workup as a new and independent set of data.

We anonymise dicoms before import to Brainlab (using
(gdcm)[https://sourceforge.net/p/gdcm/gdcm/]) to avoid any chance of
overwriting other patient data. In theory, copying tags from an
original dicom will make the converted dicoms appear to derive from
the correct patient. Please test very thoroughly before relying on
this facility.

## Longer
install, with testing (recommended)

1. Create and activate an environment, as above, then (download the package sources)[https://github.com/DevelopmentalImagingMCRI/karawun/archive/master.zip] or fetch via git:
    ```bash
    git clone https://github.com/DevelopmentalImagingMCRI/karawun.git
    ```
1. Install the dependencies
    ```bash
    cd karawun
    conda install --name KarawunEnv --file requirements.txt
    ```

1. Run the test:
    ```bash
    python -m pytest -s tests/
    ```

This test uses some data distributed with the package. A successful
result is indicated by 1 passed and 1 skipped test. Success means that
the dicom files created on your system are identical to those created
on the development system. The test destination is displayed as the
test is run.

# Licensing

_Karawun_ is distributed under the Apache License 2.0


# Contributors

Richard Beare, Joseph Yang, Andrew Perry.