# FAQ

## Why does my tract file only overlay on one of my volumes?

The tract files are associated with the first image specified in the
`--nifti` argument and by default will
only be overlaid on that volume.

However, all volumes in the `--nifti` list are assigned to a common
frame of reference (encoded in the dicom). It is possible to
accept this registration in the **Image Fusion** module in
Brainlab. Accepting the registation will allow tract and
othe objects to be overlaid on any volume.

## What is wrong with the image brightness?

Karawun does not rescale nifti data. Try rescaling by hand before
converting if default brightness is a problem. Nifti scaled to a
range of 0 to 1 will not display properly. Try multiplying by 10000
before conversion and test whether the situation improves.
