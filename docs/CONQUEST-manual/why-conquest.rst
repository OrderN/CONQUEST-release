=======================
Overview: Why CONQUEST?
=======================

Large-scale simulations
-----------------------

Efficient parallelisation
-------------------------

Linear scaling
--------------

Basis sets
----------
The main basis set in use in CONQUEST is pseudo-atomic orbitals
(PAOs).  Using these as the primitive basis set, calculations up to
several hundred atoms are standard and a thousand atoms is possible.
To go beyond this with diagonalisation, a composite basis is used,
where PAOs are combined into a smaller set of support functions
(multi-site support functions, or MSSF).  With MSSF, calculations on
10,000-20,000 atoms are possible on HPC platforms.  Beyond this size,
linear scaling must be used, which requires a further simplification
of the support functions.

