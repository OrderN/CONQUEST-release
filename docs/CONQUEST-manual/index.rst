.. CONQUEST

CONQUEST: a local orbital, large-scale DFT code
===============================================

CONQUEST is a massively parallel density functional theory (DFT) code, using a local orbital basis
to represent the Kohn-Sham eigenstates or the density matrix.
CONQUEST can be applied to atoms, molecules, liquids and solids, but
is particularly efficient for large systems.  The code can find the ground
state using exact diagonalisation of the Hamiltonian or via a linear
scaling approach.  The code has demonstrated scaling to over 2,000,000
atoms and 200,000 cores when using linear scaling.
CONQUEST can perform structural relaxation (including unit cell
optimisation) and molecular dynamics (in NVE, NVT and NPT ensembles
with a variety of thermostats).

The main basis set in use in CONQUEST is pseudo-atomic orbitals
(PAOs).  Using these as the primitive basis set, calculations up to
several hundred atoms are standard and a thousand atoms is possible.
To go beyond this with diagonalisation, a composite basis is used,
where PAOs are combined into a smaller set of support functions
(multi-site support functions, or MSSF).  With MSSF, calculations on
10,000-20,000 atoms are possible on HPC platforms.  Beyond this size,
linear scaling must be used, which requires a further simplification
of the support functions.

CONQUEST is distributed freely, under an MIT licence.  We ask that you
acknowledge use of the code by citing appropriate papers **more
details**.  If you have questions or suggestions for developing the
code, please use the GitHub interface.  The developers cannot
guarantee to offer support, though we will try to help.

.. toctree::
   :maxdepth: 1
   :caption: Contents:

   starting
   input-output
   important
   errors
   input_tags

Indices and tables
==================

.. only:: builder_html

   * :ref:`genindex`
   * :ref:`search`

.. only:: not builder_html

   * :ref:`genindex`
