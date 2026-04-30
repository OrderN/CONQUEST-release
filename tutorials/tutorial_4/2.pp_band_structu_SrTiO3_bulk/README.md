To perform a band structure plot, two Conquest calculations are needed:

1. To find a self-consistent charge density; this does not need a very fine k-mesh
2. To generate the bands, using a non-self-consistent calculation

The input files provided (Conquest_input_scf and Conquest_input_bnd)
should be copied to Conquest_input in turn before each run.  Once the
second Conquest calculation has been performed, the post-processing
calculation reorders the eigenvalues.dat file into BandStructure.dat
so that bands can be plotted.  The reference band structure was plotted
using Grace.
