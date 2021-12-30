from setuptools import setup
import versioneer

LONG_DESCRIPTION = """Converts nifti files (raw and mask) and mrtrix tck files to DICOM format
readable by Brainlab surgical planning and navigation systems.
Tck files are converted to 3D objects that can be manipulated by Brainlab tools.
Label images are converted to a dicom segmentation format and can also be
manipulated as 3D objects in Brainlab"""

setup(
    name='karawun',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    packages=['karawun'],
    python_requires='>=3.6',
    package_dir={'karawun': 'karawun'},
    url='https://developmentalimagingmcri.github.io/karawun/',
    license='Apache License, Version 2.0',
    author='Richard Beare',
    author_email='richard.beare@mcri.edu.au',
    description=('DICOM image, segmgmentation image and fibre object converter'),
    long_description=(LONG_DESCRIPTION),
    classifiers=['Programming Language :: Python :: 3',
                 'License :: OSI Approved :: Apache Software License',
                 'Operating System :: OS Independent',
                 'Development Status :: 5 - Production/Stable',
                 'Environment :: Console',
                 'Intended Audience :: Healthcare Industry',
                 'Intended Audience :: Science/Research',
                 'Natural Language :: English',
                 'Topic :: Scientific/Engineering :: Medical Science Apps.'
        ],
    install_requires=["numpy>=1.13.0",
                      "pydicom==1.4.2",
                      "SimpleITK>=1.2.0,<=2.0.2"],
    entry_points={'console_scripts': ['importTractography = karawun.commandline:import_tractography_cl']}
)
