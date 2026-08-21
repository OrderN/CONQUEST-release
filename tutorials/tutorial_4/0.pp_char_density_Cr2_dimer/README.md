This example shows you how to produce the cube files from the charge
density calculated at the self-consistent ground state.

The example is a Cr2 dimer, with opposing spins on each Cr atom.
Running the example will generate the files chden_up.001 and 
chden_dn.001.  Running the post-processing code will convert
these files to the CUBE format which can be read by VESTA and
VMD among other utilities.  A plot of the spin difference 
(formed by subtracting the spin down from the spin up density)
visualised in VESTA is given the reference directory.
