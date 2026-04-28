This tutorial shows how to create an XSF format file which permits forces
on atoms to be displayed in VESTA.  After running Conquest calculation and
the post-processing to generate the XSF file, some further work is needed
(at present).  The forces, which are printed at the end of the Conquest_out
file, must be added as extra columns in the XSF file.  There are two ways
to do this:

1. Using an appropriate editor, copy the forces from the output file and
   manually paste them into the XSF file
2. Using Unix commands (grep, cut, head, paste, cat) perform the same set
   of operations.  These commands are given at the end of the run script
   file pp_forces_xsf_runCQ.sh

In future, we plan to allow the post processing code to add the forces or
other appropriate vectors to an XSF file automatically.
