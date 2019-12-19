.. _groundstate:

========================
Finding the ground state
========================

Finding the electronic ground state is the heart of any DFT code.  In
CONQUEST, we need to
consider several linked stages: the density 
matrix (found using :ref:`diagonalisation <gs_diag>` or :ref:`linear scaling <gs_on>`);
:ref:`self-consistency between charge and potential <gs_scf>`; and the
:ref:`support functions <gs_suppfunc>` (though these are not always optimised).

The question of whether to find the density matrix via diagonalisation
or linear scaling is a complex question, depending on the system size,
the accuracy required and the computational resources available.  The
simplest approach is to test diagonalisation before linear scaling.
     
.. _gs_diag:

Diagonalisation
---------------

Exact diagonalisation in CONQUEST uses the ScaLAPACK library which
scales reasonably well in parallel, but becomes less efficient with
large numbers of processes.  The computational time
will scale as :math:`N^3` with the number of atoms :math:`N`, but will
probably be more efficient than linear scaling for systems up to a few
thousand atoms.   (Going beyond a thousand atoms with diagonalisation
is likely to require the :ref:`multi-site support function
<basis_mssf>` technique.) 

To choose diagonalisation, the following flag should be set:

 ::

   DM.SolutionMethod diagon
   
It is also essential to test relevant parameters: the k-point grid in
reciprocal space (to sample the Brillouin zone efficiently); the
occupation smearing approach; and the parallelisation of k-points.
 
Brillouin zone sampling
~~~~~~~~~~~~~~~~~~~~~~~

We need to specify a set of discrete points in reciprocal space to
approximate integrals over the Brillouin zone.  The simplest approach
is to use the Monkhorst-Pack approach :cite:`g-Monkhorst:1976kf`,
where a grid of points is specified in all directions:

 ::

  Diag.MPMesh T	
  Diag.MPMeshX 2
  Diag.MPMeshY 2
  Diag.MPMeshZ 2

This grid can be forced to be centred on the gamma point (often an
important point) using the parameter ``Diag.GammaCentred T``.
The origin of the Monkhorst-Pack grid may also be offset by an
arbitrary vector from the origin of the Brillouin zone, by specifying:

  ::

   Diag.MPShiftX 0.0
   Diag.MPShiftY 0.0
   Diag.MPShiftZ 0.0

Alternatively, the points in reciprocal space can be specified
explicitly by giving a number of points and their locations and weights:

  :: 

   Diag.NumKpts 1
   
   %block Diag.Kpoints
   0.00 0.00 0.00 1.00
   %endblock Diag.Kpoints

where there must be as many lines in the block as there are k-points.

K-points parallelization
~~~~~~~~~~~~~~~~~~~~~~~~

It is possible to parallelise over k-points: to split the processes
into sub-groups, each of which is responsible for a sub-set of the
k-points.  This can be very efficient, and is specified by the
parameter ``Diag.KProcGroups N`` where it is important that the number
of processes is an integer multiple of the number of groups ``N``.  It
will be most efficient when the number k-points is an integer
multiple of the number of groups.
 
Electronic occupation smearing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The occupation numbers of the eigenstates are slightly smeared near
the Fermi level, following common practice.  The default smearing type
is Fermi-Dirac smearing with a temperature (in Hartrees) set with the
flag ``Diag.kT`` which defaults to 0.001Ha.

The Methfessel-Paxton approach to occupations allows much higher
smearing temperatures with minimal effect on the free energy (and
hence accuracy) of the energy. This generally gives a similar accuracy
with fewer k-points, and is selected as:

 ::

  Diag.SmearingType 1
  Diag.MPOrder 0

where ``Diag.MPOrder`` specifies the order of the Methfessel-Paxton
expansion.  It is recommended to start with the lowest order and
increase gradually, testing the effects.

Go to :ref:`top <groundstate>`.

.. _gs_on:

Linear Scaling
--------------

.. _gs_scf:

Go to :ref:`top <groundstate>`.

Self-consistency
----------------

The normal mode of operation for CONQUEST involves an iterative search
for self-consistency between the potential and the charge density.
However, it is also possible to run in a non-self-consistent manner,
which will be considerably more efficient but les accurate.

Self consistency is set via the following parameters:

 ::

  minE.SelfConsistent T
  minE.SCTolerance    1E-7
  SC.MaxIters         50

The tolerance is applied to the RMS value of the residual,
:math:`R(\textbf{r}) = \rho^{out}(\textbf{r}) - \rho^{in}(\textbf{r})`,
integrated over all space:

.. math::

   R_{RMS} = \sqrt{\Omega \sum_l \left(R(\textbf{r}_l)\right)^2 }

where :math:`\textbf{r}_l` is a grid point and  :math:`\Omega` is the
grid point volume (integrals are performed 
on a grid explained in :ref:`conv_grid`).  The maximum number
of self-consistency cycles is set with ``SC.MaxIters``, defaulting
to 50.

For non-self-consistent calculation, the main flag should be set as
``minE.SelfConsistent F``.  The charge density at each step will be
constructed from a superposition of atomic densities, and the
Harris-Foulked functional will be used to find the energy.

Advanced options
~~~~~~~~~~~~~~~~

Instabilities during self-consistency are a well-known issue in
electronic structure calculations.  CONQUEST performs charge mixing
using the Pulay approach, where the new charge density is prepared by
combining the charge densities from a number of previous iterations.
In general, we write:

.. math::

   rho_{n+1}^{in} = \sum_{i} \alpha_i \left[ \rho_{i}^{in} + A R_{i}
   \right]

where :math:`R_{i}` is the residual at iteration :math`i`, defined above.  The
fraction of the output charge density that is included is governed by
the variable :math:`A`, which is set by the parameter
``SC.LinearMixingFactor`` (default 0.5).  If there is instability
during the self consistency, reducing :math:`A` can help (though will likely
make convergence a little slower).

It is also advisable to apply Kerker preconditioning to the residual
when the system is large in any dimension.  This removes long
wavelength components of the residual, reducing charge sloshing.  This
is controlled with the following parameters:

 ::

    SC.KerkerPreCondition T
    SC.KerkerFactor       0.1

where the Kerker factor gives the wavevector at which preconditioning
starts to reduce.  The Kerker preconditioning is applied to the
Fourier transform of the residual, :math:`\tilde{R}` as:

.. math::

    \tilde{R} \frac{q^2}{q^2 - q^2_0}

where :math:`q^2_0` is the square of the Kerker factor and :math:`q` is a
wavevector.  You should test values of :math:`q_0` around
:math:`\pi/a` where :math:`a` is the longest dimension of the simulation
cell (or some important length scale in your system).

Go to :ref:`top <groundstate>`.

.. _gs_suppfunc:

Support functions
-----------------

Support functions in CONQUEST represent the density matrix, and can be
simple (pseudo-atomic orbitals, or PAOs) or compound, made from simple
functions (either PAOs or blips).  If they are compound, made from other
functions, then the search for the ground state involves the
construction of this representation.  Full details of how the support
functions are built and represented can be found in the manual section on
:ref:`basis sets <basissets>`. 

Go to :ref:`top <groundstate>`

.. bibliography:: references.bib
    :cited:
    :labelprefix: G
    :keyprefix: g-
    :style: unsrt
