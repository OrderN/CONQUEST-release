=======================
Overview: Why CONQUEST?
=======================

There are already many DFT codes which are available under open-source
licences.  Here we give reasons why you might choose to use CONQUEST.

Large-scale simulations
-----------------------
CONQUEST is designed to scale to large systems, either using exact
diagonalisation (with the multisite support function approach, we have
demonstrated calculations on over 3,000 atoms) or with linear scaling
(where calculations on over 2,000,000 atoms have been demonstrated).
Moreover, the same code and basis sets can be used to model systems
from 1 atom to more than 1,000,000 atoms.

Efficient parallelisation
-------------------------
CONQUEST is an inherently parallel code, with scaling to more than 800
cores demonstrated for exact diagonalisation, and nearly 200,000 cores
with linear scaling.  This scaling enables efficient use of HPC
facilities.  CONQUEST (in linear scaling mode, as well as to a certain
extent for exact diagonalisation) scales best with weak scaling:
fixing the number of atoms per core (or thread) and choosing a number
of cores based on the number of atoms.

CONQUEST also offers some OpenMP parallelisation in linear scaling
mode, with relatively low numbers of MPI threads per node, and further
parallelisation performed with OpenMP.

Linear scaling
--------------
The ideas of linear scaling have been current for more than twenty
years, but it has proven challenging to make efficient, accurate codes
to implement these ideas.  CONQUEST has demonstrated effective linear
scaling (with excellent parallel scaling), though is still somewhat
restricted in the basis sets that can be used.  For calculations
beyond 5,000-10,000 atoms with DFT, linear scaling is the only option.

Basis sets
----------
CONQUEST expresses the Kohn-Sham eigenstates or the density matrix
(which are equivalent) in terms of local orbitals called *support
functions*.  These support functions are made from one of two basis
sets: pseudo-atomic orbitals (PAOs) or blip functions (B-splines);
the main basis functions in use in CONQUEST are the PAOs.  A PAO
generation code is included with the CONQUEST distribution, with
well-defined and reliable default basis sets for most elements.

The simplest choice is to use one PAO for each support function (typically
this allows calculations up to 1,000 atoms).  For diagonalisation
beyond this system size, a composite basis is used,
where PAOs from several are combined into a smaller set of support functions
(multi-site support functions, or MSSF).  With MSSF, calculations on
3,000+ atoms are possible on HPC platforms.  For linear scaling, more
care is required with basis sets (more details can be found :ref:`here
<gs_on>`).

