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

The standard self-consistency approach uses the Pulay RMM method, and
should be robust in most cases.  It can be monitored via the residual,
which is currently defined as the standard RMS difference in charge
density:

.. math::

   R = \sqrt{\int \mathrm{d}\mathbf{r}\mid \rho^{out}(\mathbf{r}) -
   \rho^{in}(\mathbf{r})\mid^2}

where :math:`\rho^{in}` is the input charge density for an iteration,
and  :math:`\rho^{out}`  is the resulting output charge density.  The
SCF cycle is terminated when this residual is less than the parameter
``minE.SCTolerance``.  The maximum number of iterations is set with
``SC.MaxIters`` (defaults to 50).

There are various further approaches and parameters which can be used
if the SCF cycle is proving hard to converge.  As is standard, the
input for a given iteration is made by combining the charge density
from a certain number of previous steps (``SC.MaxPulay``, default 5).
The balance between input and output charge densities from these
previous steps is set with ``SC.LinearMixingFactor`` (default 0.5;
N.B. for spin polarised calculations,
``SC.LinearMixingFactor_SpinDown`` can be set separately).  Reducing
this quantity may well improve stability, but slow down the rate of
convergence.

Kerker-style preconditioning (damping long wavelength charge
variations) can be selected using ``SC.KerkerPreCondition T`` (this is
most useful in metallic and small gap systems).  The preconditioning
is a weighting applied in reciprocal space:

.. math::

   K = \frac{1}{1+q^2_0/q^2}

where :math:`q_0` is set with ``SC.KerkerFactor`` (default 0.1).
This is often very helpful with slow convergence or instability.

Go to :ref:`top <convergence>`.

.. _conv_suppfunc:

Support Functions
-----------------

The parameters relevant to support functions depend on the basis set
that is used.  In the case of pseudo-atomic orbitals (PAOs), when
support functions are primitive PAOs, the only relevant parameter is
the basis set size, which is set when the ion files are generated.  It
is important to test the accuracy of a given basis set carefully for
the problem that is to be modelled.

When using multi-site support functions (MSSF), the key parameter is
the radius of the MSSF (``Atom.MultisiteRange`` in
the :ref:`atomic specification <input_atomic_spec>` block).
As this is increased, the accuracy of the 
calculation will also increase, but with increased computational
effort.  Full details of the MSSF (and related OSSF) approach are
given in the section on :ref:`multi-site support functions
<basis_mssf>`.

For the blip basis functions, the spacing of the grid where the blips
are defined is key (``Atom.SupportGridSpacing`` in
the :ref:`atomic specification <input_atomic_spec>` block),
and is directly related to an equivalent plane 
wave cutoff (via :math:`k_{bg} = \pi/\delta` and :math:`E_{PW} =
k_{bg}^2/2`, where :math:`\delta` is the grid spacing in Bohr radii
and :math:`E_{PW}` is in Hartrees).  For a particular grid spacing,
the energy will converge monotonically with support function radius
(``Atom.SupportFunctionRange`` in
the :ref:`atomic specification <input_atomic_spec>` block).
A small support function radius will introduce some approximation to
the result, but improve computational performance.  It is vital to
characterise both blip grid spacing and support function radius in any
calculation.  A full discussion of the blip function basis is found
:ref:`here <basis_blips>`.

Go to :ref:`top <convergence>`.
