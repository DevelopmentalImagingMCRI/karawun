.. Karawun documentation master file, created by
   sphinx-quickstart on Fri Dec 10 14:28:02 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Karawun's documentation!
===================================

*Karawun* is a python package containing a single command-line tool for
converting diffusion imaging workups into a DICOM form that can be imported
into the Brainlab surgical planning and navigation suite.

*Karawun* offers a unique capability in DICOM conversion, namely creating
two DICOM segmentation file types in addition to the standard MR type.
One of these types is able to represent
streamlines used in tractography while the other is able to represent label
or segmentation images (such as those created when delineating a tumour).
Both types are displayed by Brainlab as 3D :ref:`objects <Examples>`.

Using *Karawun* it is possible to import an analysis
created using open source tools, such as `MRtrix3 <https://www.mrtrix.org/>`__, into Brainlab and view
and manipulate the data as if it had been created within
Brainlab. This allows clinicians to interact with data generated using
the latest experimental tools.


.. toctree::
   :maxdepth: 2
   :caption: Contents: 

   installation.md
   using.md
   importing.md
   screenshots.md
   testing.md
   dicomtypes.md
   whatitsnot.md
   name.md
   citing.md
   faq.md

   
