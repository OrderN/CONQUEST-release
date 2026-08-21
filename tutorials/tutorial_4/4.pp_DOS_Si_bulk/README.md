This example shows you how to produce the density of states (DOS) for 
bulk Si

When calculating DOS for a solid, a fine k-mesh is required, but a 
self-consistent charge density from a relatively coarse k-mesh can
be used to improve efficiency.  Use the Conquest_input_charge input
file (copied to Conquest_input) to generate the charge density, and
then copy Conquest_input_dos to Conquest_input to generate the 
eigenvalues needed for the DOS calculation.  The same input file
is used for the post-processing code to generate the DOS.
This file has three columns: the energy relative to
the Fermi level (in eV); the DOS; and the integrated DOS.

There is a Jupyter notebook which will produce a plot of the DOS which
shows how to analyse the file.
