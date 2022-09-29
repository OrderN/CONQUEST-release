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

The basis functions in CONQUEST are :ref:`support functions <gs_suppfunc>` (localised
functions centred on the atoms), written as
:math:`\phi_{i\alpha}(\textbf{r})` where :math:`i` indexes an atom and
:math:`\alpha` a support function on the atom.  The support functions
are used as basis functions for the density matrix and the Kohn-Sham
eigenstates:

.. math::
   \psi_{n\mathbf{k}}(\mathbf{r}) = \sum_{i\alpha} c^{n\mathbf{k}}_{i\alpha}
   \phi_{i\alpha}(\mathbf{r})\\
   \rho(\mathbf{r}, \mathbf{r}^\prime) = \sum_{i\alpha j\beta}
   \phi_{i\alpha}(\mathbf{r}) K_{i\alpha, j\beta} \phi_{j\beta}(\mathbf{r}^\prime)

where :math:`n` is an eigenstate index and :math:`\mathbf{k}` is a
point in the Brillouin zone (see :ref:`here <gs_diag_bz>` for more on
this).  The total energy can be written in terms of the density
matrix, as:

.. math::
   E_{KS} = \mathrm{Tr}[HK] + \Delta E_{Har} + \Delta E_{XC}

for the Hamiltonian matrix :math:`H` in the basis of support
functions, with the last two terms the standard Harris-Foulkes
:cite:`g-Harris1985,g-Foulkes1989` correction terms.
      
For diagonalisation, the density matrix is made from the coefficients
of the Kohn-Sham eigenstates, :math:`c^{n\mathbf{k}}_{i\alpha}`, while
for :ref:`linear scaling <gs_on>` it is found directly during the variational
optimisation of the energy.

The question of whether to find the density matrix via diagonalisation
or linear scaling is a complex one, depending on the system size,
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
   
It is also essential to test relevant parameters, as described below:
the k-point grid in reciprocal space (to sample the Brillouin zone
efficiently); the occupation smearing approach; and the
parallelisation of k-points.

Go to :ref:`top <groundstate>`

.. _gs_diag_bz:

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
It is important to note that CONQUEST does not consider space group
symmetry when integrating over the Brillouin zone.

Go to :ref:`top <groundstate>`.

.. _gs_diag_para:

K-point parallelization
~~~~~~~~~~~~~~~~~~~~~~~~

It is possible to parallelise over k-points: to split the processes
into sub-groups, each of which is responsible for a sub-set of the
k-points.  This can be very efficient, and is specified by the
parameter ``Diag.KProcGroups N``, where it is important that the number
of processes is an integer multiple of the number of groups ``N``.  It
will be most efficient when the number of k-points is an integer
multiple of the number of groups.
 
Go to :ref:`top <groundstate>`.

.. _gs_diag_smear:

Electronic occupation smearing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The occupation numbers of the eigenstates are slightly smeared near
the Fermi level, following common practice.  The default smearing type
is Fermi-Dirac smearing with a temperature (in Hartrees) set with the
flag ``Diag.kT`` which defaults to 0.001Ha.

The Methfessel-Paxton approach :cite:`g-Methfessel:1989ny` to occupations allows much higher
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

A linear scaling calculation is selected by setting
``DM.SolutionMethod ordern``.  There are two essential parameters that must be
set: the range of the density matrix, and the tolerance on the
optimisation.

 ::
    
    DM.L_range 16.0
    minE.Ltolerance 1.0e-6

The tolerance is applied to the residual (the RMS value of the
gradient of the energy with respect to the density matrix).  The
maximum number of iterations in the density matrix optimisation can
be set with ``DM.LVariations`` (default 50).

At present, CONQUEST can only operate efficiently in linear scaling
mode with a restricted number of support functions (though this is an
area of active development).  PAO basis sets of SZ and SZP size
(minimal and small in the ion file generator) will run without
restrictions.  For larger PAO basis sets, the :ref:`OSSF <basis_ossf>`
approach must be used, and is effective.  With a blip basis there are
no restrictions, though efficient optimisation is still under active
development. 


It is
almost always more efficient to update the charge density while
optimising the density matrix, avoiding the need for a separate
self-consistency loop.  This is set by choosing
``minE.MixedLSelfConsistent T``. 

An essential part of a linear scaling calculation is finding the
approximate, sparse inverse of the overlap matrix.  Normally this will
happen automatically, but it may require some tests.  The key
parameters are the range for the inverse (see the
:ref:`input_atomic_spec` block, and specifically the
:ref:`advanced_atomic_spec_tags` block) and the tolerance applied
to the inversion.

 ::
    
    Atom.InvSRange R
    DM.InvSTolerance R

A tolerance of up to 0.2 can give convergence without significantly
affecting the accuracy.  The range should be similar to the radius of
the support functions, though increasing it by one or two bohr can
improve the inversion in most cases.
    
The input tags are mainly found in the :ref:`input_dm` section of the
:ref:`input_tags` page.
     
Go to :ref:`top <groundstate>`.

.. _gs_scf:

Self-consistency
----------------

The normal mode of operation for CONQUEST involves an iterative search
for self-consistency between the potential and the charge density.
However, it is also possible to run in a non-self-consistent manner,
either with a converged charge density for electronic structure
analysis, or for dynamics, which will be considerably more efficient
than a self-consistent calculation, but less accurate.

Self consistency is set via the following parameters:

 ::

  minE.SelfConsistent T
  minE.SCTolerance    1E-7
  SC.MaxIters         50

The tolerance is applied to the RMS value of the residual,
:math:`R(\mathbf{r}) = \rho^{out}(\mathbf{r}) - \rho^{in}(\mathbf{r})`,
integrated over all space:

.. math::

   R_{RMS} = \sqrt{\Omega \sum_l \left(R(\mathbf{r}_l)\right)^2 }

where :math:`\mathbf{r}_l` is a grid point and  :math:`\Omega` is the
grid point volume (integrals are performed 
on a grid explained in :ref:`conv_grid`).  The maximum number
of self-consistency cycles is set with ``SC.MaxIters``, defaulting
to 50.

For non-self-consistent calculations, the main flag should be set as
``minE.SelfConsistent F``.  The charge density at each step will
either be read from a file (if the flag ``General.LoadRho T`` is set),
or constructed from a superposition of 
atomic densities.  The Harris-Foulkes functional will be used to
find the energy. 

Go to :ref:`top <groundstate>`.

.. _ gs_scf_adv:

Advanced options
~~~~~~~~~~~~~~~~

Instabilities during self-consistency are a well-known issue in
electronic structure calculations.  CONQUEST performs charge mixing
using the Pulay approach, where the new charge density is prepared by
combining the charge densities from a number of previous iterations.
In general, we write:

.. math::

   \rho_{n+1}^{in} = \sum_{i} \alpha_i \left[ \rho_{i}^{in} + A R_{i}
   \right]

where :math:`R_{i}` is the residual at iteration :math:`i`, defined above.  The
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

    \tilde{R} \frac{q^2}{q^2 + q^2_0}

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

.. _gs_charged:

Charged systems
-----------------
CONQUEST uses periodic boundary conditions, which require overall
charge neutrality.  However, charged systems can be modelled:
if an excess of electrons is specified by the user, a uniform
positive background charge is added automatically to restore overall
neutrality.  At present, there are no correction schemes implemented,
so it is important to test the convergence of the energy with unit
cell size and shape.  Electrons are added by setting the parameter
``General.NetCharge``.

::

   General.NetCharge 1.0

This gives the number of extra electrons to be added to the unit cell,
beyond the valence electrons.

Go to :ref:`top <groundstate>`.

.. _gs_spin:

Spin polarisation
-----------------
CONQUEST performs collinear spin calculations only.  A spin-polarised
calculation is performed by setting the parameter
``Spin.SpinPolarised`` to T.  

Users need to specify *either* the total initial number of spin-up and spin-down electrons in
the simulation cell (using the parameters ``Spin.NeUP`` and
``Spin.NeDN``), *or* the difference between the number of spin-up and
spin-down electrons (using the parameter ``Spin.Magn``).

The number of electrons for each spin channel can be fixed during SCF
calculations by setting the parameter ``Spin.FixSpin`` to T (default is F).

It is possible to specify the spin occupation in the atomic charge
densities (i.e. the number of spin-up and spin-down electrons used to
build the density).  This is done in the :ref:`input_atomic_spec`
part of the ``Conquest_input`` file.  Within the atom block for
each species, the numbers of electrons should be set with
``Atom.SpinNeUp`` and ``Atom.SpinNeDn``.  Note that these numbers
*must* sum to the number of valence electrons for the atom.

Go to :ref:`top <groundstate>`.

.. _gs_spin_example:

Examples: FM and AFM iron
~~~~~~~~~~~~~~~~~~~~~~~~~

A two atom ferromagnetic iron simulation might be set up using the
parameters below.  Note that the net spin here is S=1 :math:`\mu_B`
(i.e. two more electrons in the up channel than in the down), and
that the net spin is not constrained.

:: 

   # example of ferro bcc Fe
   Spin.SpinPolarised T
   Spin.FixSpin  F
   Spin.NeUP  9.0     # initial numbers of up- and down-spin electrons,
   Spin.NeDN  7.0     # which will be optimised by a SCF calculation when Spin.FixSpin=F
   
   %block ChemicalSpeciesLabel
   1   55.845   Fe
   %endblock ChemicalSpeciesLabel

An equivalent anti-ferromagnetic calculation could be set up as
follows (though note that the initial specification of spin for the
atoms does *not* guarantee convergence to an AFM ground state).  By
defining two species we can create spin-up and spin-down atoms (note
that both species will require their own, appropriately labelled, ion
file). 

::

   # example of anti-ferro bcc Fe
   Spin.SpinPolarised T
   Spin.FixSpin  F
   Spin.NeUP  8.0     # initial numbers of up- and down-spin electrons in an unit cell
   Spin.NeDN  8.0     # are set to be the same
   
   %block ChemicalSpeciesLabel
   1   55.845   Fe1
   2   55.845   Fe2
   %endblock ChemicalSpeciesLabel
   
   %block Fe1           # up-spin Fe
   Atom.SpinNeUp 5.00
   Atom.SpinNeDn 3.00
   %endblock Fe1
   %block Fe2           # down-spin Fe
   Atom.SpinNeUp 3.00
   Atom.SpinNeDn 5.00
   %endblock Fe2

When using multi-site or on-site support functions in spin-polarised
calculations, the support functions can be made spin-dependent
(different coefficients for each spin channel) or not by setting
``Basis.SpinDependentSF`` (T/F, default is T).

Go to :ref:`top <groundstate>`.

.. bibliography:: references.bib
    :cited:
    :labelprefix: G
    :keyprefix: g-
    :style: unsrt

Go to :ref:`top <groundstate>`.
