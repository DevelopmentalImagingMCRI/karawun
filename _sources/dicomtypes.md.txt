# DICOM types

Images passed via the `--nifti` flag are converted to DICOM using
the MR modality. All images passed via this flag are assigned to
a common Frame of Reference.

Label images and fiber tracts are both are represented in DICOM using
the Segmentation Information Object Definition (IOD) (SEG modality (tag 0008,0060))

## Label images

Per region text labels (0062,0005) are derived from the
nifti filename and the label number. e.g. in the test dataset the
"globes.nii.gz" file produces segment labels of the form "Label 1 of
globes". The "Segmented Property Category Code Sequence" (0062,0003)
and "Segmented Property Type Code Sequence" (0062,000F) are set to
generic values ("Anatomical Structure" and "Unspecified".

Label pixel values are stored using RLE compression described in [Annex G](https://dicom.nema.org/medical/dicom/current/output/chtml/part05/chapter_G.html)

## Tract files

Fiber objects are represented using the Line Sequence surface mesh primitive of the Surface Mesh Module in  the Surface Information Entity (IE).

## Slicing

3D nifti files are converted to 2D DICOM slices, with slice plane
selected to produce isotropic 2D images. In the case of 3D isotropic images
a slice plane will be selected that minimizes the number of 2D files.

An error is raised if no isotropic slices are possible.