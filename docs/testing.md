# Testing

Regression testing is covered in [Installation](verifying). Those checks
verify that the DICOM files created by an installation are identical
to those created at baseline. Checksums are used to compare
the files to reduce the package size. The test procedures use
python "mocking" facilities to force repeatable creation of timestamps
and UUIDs to simplify the comparison procedures.

## Other test data

An additional set of test data is available on
[github](https://github.com/DevelopmentalImagingMCRI/karawun_test_data.git). Screenshots
created using this test dataset are shown in [](Examples).

This set of test data is derived from an annotated nifti file
distributed with the [FSL
suite](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/). The annotations
indicate Left/Right hemispheres. In addition to the original image, we
provide a range of resliced and reformatted versions in both raw and
label forms, where the label forms include only the Left/Right
markers, and some synthetic streamline data, also showing the
Left/Right markers.

All versions of the images are in the same space, and should remain in
the same space after conversion to DICOM and importation to Brainlab,
irrespective of the data stride pattern or slice thickness.