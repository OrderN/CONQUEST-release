# CONQUEST: Large-scale DFT calculations

CONQUEST is a DFT code designed for large-scale calculations, with
excellent parallelisation.  It gives a consistent, exact
diagonalisation approach for systems from 1 to 10,000+ atoms, and
brings the possibility of linear scaling calculations on over
1,000,000 atoms.  The code has been demonstrated on scaling to nearly
200,000 cores and 2,000,000 atoms.

## Capabilities

CONQUEST can perform static electronic structure calculations
(including DOS and band structure calculations), structural
relaxation, molecular dynamics (with NVE, NVT and NPT ensembles all
implemented).  It can output energies, forces and stresses as well as
density of states, charge density, orbital density and Tersoff-Hamann
STM images.  The facility to perform delta-SCF and cDFT calculations
is available.  The code can use LibXC for access to a wide variety of
exchange-correlation functionals.

CONQUEST reads pseudopotentials produced by Don Hamann's
[ONCVPSP](http://www.mat-simresearch.com) code, and is fully
compatible with the PseudoDojo database.  It can also read the .ion
files produced by Siesta, using both the pseudopotentials and
pseudo-atomic orbitals (PAOs) in those files.  A CONQUEST code to
generate PAO basis sets is included in the distribution.
