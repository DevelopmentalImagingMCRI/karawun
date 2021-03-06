import setuptools
from setuptools import setup

LONG_DESCRIPTION = \
'''Converts nifti files (raw and mask) and mrtrix tck files to DICOM format readable by Brainlab surgical planning and navigation systems. Tck files are converted to 3D objects that can be manipulated by Brainlab tools. Label images are converted to a dicom segmentation format and can also be manipulated as 3D objects in Brainlab'''

setup(
    name='karawun',
    version='0.2.0.1',
    packages=['karawun'],
    python_requires='>=3.6',
    package_dir={'karawun': 'karawun'},
    url='',
    license='',
    author='Richard Beare',
    author_email='richard.beare@mcri.edu.au',
    description=('DICOM image, segmgmentation image and fibre object converter'),
    long_description=(LONG_DESCRIPTION),
    install_requires=["numpy==1.13.0",
                      "pydicom==1.3",
                      "SimpleITK==1.2.0"],
    entry_points={'console_scripts': ['importTractography = karawun.commandline:import_tractography_cl']}
)

