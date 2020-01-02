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

The range applied to the density matrix (``DM.L_range``) determines
the accuracy of the calculation, as well as the computational time
required (as the number of non-zero elements will increase based on a
sphere with the radius of the range, the time will increase roughly
proportional to the cube of the range).  In almost all circumstances,
it is best to operate with a range which converges energy
*differences* and forces, rather than the absolute energy.  Testing
for this convergence is an essential part of the preparation for
production calculations.

The tolerance applied to the density matrix optimisation
(``minE.Ltolerance``) must be
chosen to give adequate convergence of the energy and forces.  The
tolerance is applied to the residual in the calculation, defined as:

.. math::

   R = \sqrt{\sum_{i\alpha j\beta} \partial E/\partial L_{i\alpha j\beta}
   \cdot \partial E/\partial L_{i\alpha j\beta} }

The dot product uses the inverse of the overlap matrix as the metric.

The approximate, sparse inversion of the overlap matrix is performed
before the optimisation of the density matrix.  The method used,
Hotelling's method (a version of a Newton-Raphson approach) is
iterative and terminates when the characteristic quantity
:math:`\Omega` increases.  On termination, if :math:`\Omega` is below
the tolerance ``DM.InvSTolerance`` then the inverse is accepted;
otherwise it is set to the identity (the density matrix optimisation
will proceed in this case, but is likely to be inefficient).  We
define:

.. math::

   \Omega = (Tr[I - TS])^2

where :math:`T` is the approximate inverse.  The range for the inverse
must be chosen (``Atom.InvSRange`` in the species block); by default
it is same as the support function range 
(which is then doubled to give the matrix range) but can be
increased.  The behaviour of the inversion with range is not simple,
and must be carefully characterised if necessary.

Go to :ref:`top <convergence>`.

.. _conv_scf:

Self-consistency
----------------

Go to :ref:`top <convergence>`.

.. _conv_suppfunc:

Support Functions
-----------------

Go to :ref:`top <convergence>`.
