.. _convergence:

=====================
Converging Parameters
=====================

There are various important parameters in CONQUEST that affect the
convergence of the total energy, and need to be tested.  Integrals are
calculated on a grid; the density matrix is found approximately; a
self-consistent charge density is calculated; and support functions
are, in some modes of operations, optimised.  These parameters are
described here.

.. _conv_grid:

Integration Grid
----------------

While many integrals are calculated analytically or on fine grids that
move with the atoms, there are still some integrals that must be found
numerically, and CONQUEST uses an orthorhombic, uniform grid to
evaluate these integrals (this grid is also used for the Fourier
transforms involved in finding the Hartree potential).  The
spacing of the grid will affect the accuracy of the calculation, and
it is important to test the convergence of the total energy with the
grid spacing.

The grid spacing can be set intuitively using an energy (which
corresponds to the kinetic energy of the shortest wavelength wave that
can be represented on the grid).  In atomic units, :math:`E = k^2/2`
with :math:`k = \pi/\delta` for grid spacing :math:`\delta`.  The
cutoff is set with the parameter:

 ::
  
  Grid.GridCutoff E

where ``E`` is an energy in Hartrees.  The grid spacing can also be
set manually, by specifying the number of grid points in each
direction:

 ::

    Grid.PointsAlongX N
    Grid.PointsAlongY N
    Grid.PointsAlongZ N

If setting the grid in this manner, it is important to understand a
little more about the internal workings of CONQUEST.  The grid is divided up into
*blocks* (the default size is 4 by 4 by 4), and the number of grid
points in any direction must correspond to an integer multiple of the
block size in that direction.  The block size can be set by the user:

 ::

    Grid.InBlockX N
    Grid.InBlockY N
    Grid.InBlockZ N

Note that the blocks play a role in parallelisation and memory use, so
that large blocks may require larger memory per process; we recommend
block sizes no larger than 8 grid points in each direction.
There is also, at present, a restriction on the total number of grid
points in anuy direction, that it must have prime factors of only 2, 3 and 5.  This will be
removed in a future release.
    
Go to :ref:`top <convergence>`.

.. _conv_dm:

Finding the density matrix
--------------------------

As discussed in the section on :ref:`finding the ground
state<groundstate>`,
the density matrix is found either
with exact diagonaliation, or the linear scaling approach.  These
two methods require different convergence tests, and are described separately.

.. _conv_dm_bz:

Diagonalisation: Brillouin Zone Sampling
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The sampling of the Brillouin zone must be tested for convergence, and
the parameters are described :ref:`here <gs_diag_bz>`.  The
convergence of charge density will be faster than detailed electronic
structure such as density of states (DOS), and it will be more
accurate for these types of calculations to generate a converged charge
density, and then run non self-consistently (see the section on
:ref:`self consistency <gs_scf>`) with appropriate k-point sampling.

Go to :ref:`top <convergence>`.

.. _conv_on:

Linear Scaling
~~~~~~~~~~~~~~

* L range
* L Tolerance
* Inverse S range
* Inverse S tolerance
* Number of iterations (?)

Go to :ref:`top <convergence>`.

.. _conv_scf:

Self-consistency
----------------

Go to :ref:`top <convergence>`.

.. _conv_suppfunc:

Support Functions
-----------------

Go to :ref:`top <convergence>`.
