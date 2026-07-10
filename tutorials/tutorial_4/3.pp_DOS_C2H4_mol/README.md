This example shows you how to produce the density of states (DOS) for 
a simple molecule, C2H4.

You should first run CONQUEST using the supplied input files, which will 
find the ground state, and output the Kohn-Sham eigenvalues (in the file
eigenvalues.dat).  Running the post-processing code will generate the
DOS file, DOS.dat.  This file has three columns: the energy relative to
the Fermi level (in eV); the DOS; and the integrated DOS.

There is a Jupyter notebook which will produce a plot of the DOS which
shows how to analyse the file.
