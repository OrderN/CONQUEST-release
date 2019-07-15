.. _input_tags:

==========
Input tags
==========

We have broken down the input tags based on the areas of the code
where they apply.  For each tag, a default is given.  Types of value
are specified as: *integer*;
*real*; *boolean*; or *string* (optA/optB are given for string options).

.. _input_tags_general:

General
-------

General.Title (*string*)
    Title for the calculation

    *default*: none

General.NumberOfSpecies (*integer*)
    Number of species in cell

    *default*: none

General.PseudopotentialType (*string*) siesta/hamann
    Type of pseudopotential (in practice, this defines how the local
    part of the pseudopotential is handled)

    *default*: siest
General.NeutralAtom (*boolean*)
    Use neutral atom potential or not (removes need for Ewald sum)

    *default*: T

General.FunctionalType (*integer*)
    Selects the exchange-correlation functional. If the native
    CONQUEST XC implementation is used, there are three
    parameterisations of the LDA available, as well as three variants
    of the PBE GGA functional, with numbers given below.

    *default*: 3

    =========================================  ======= =======================
    Functional                                 Keyword Ref
    =========================================  ======= =======================
    LDA Perdew-Zunger, no SIC                  1       :cite:`e-Perdew1981`
    LDA Goedecker-Teter-Hutter 96              2       :cite:`e-Goedecker1996`
    LSDA Perdew-Wang 92 (default)              3       :cite:`e-Perdew1992`
    GGA Perdew-Burke-Ernzerhof 96 (PBE)        101     :cite:`e-Perdew1996`
    GGA PBE + Zhang-Yang 98 (revPBE)           102     :cite:`e-Zhang:1998oq`
    GGA PBE + Hammer-Hansen-Norskov 99 (RPBE)  103     :cite:`e-Hammer1999`
    GGA WC                                     104     :cite:`e-Wu:2006cu`
    =========================================  ======= =======================

    At the moment, only LSDA Perdew-Wang 92 and the three GGA
    Perdew-Burke-Ernzerhof functional variants can be used in spin polarised calculations.

    Note that, if the code is compiled with LibXC, the full LibXC
    set of functionals is available, selected with a negative six
    digit number (-XXXCCC or -CCCXXX).

General.EnergyUnits (*string*) Ha/Ry/eV
    **Output only** Chooses units for energy

    *default*: Ha

General.DistanceUnits (*string*) a0/bohr/A
    **Output only** Chooses units for distance (Bohr: a0/bohr or Ångströms: A

    *default*: a0

General.MemoryUnits (*string*) kB/MB/GB
    **Output only** Chooses units for memory use

    *default*: MB

General.PartitionMethod (*string*) File/Hilbert
    Chooses method for partitioning (read from file or use dynamic partitioner
    based on Hilbert curve)

    *default*: Hilbert

    Options:

    -  Hilbert (default) — Automatic partitioning using Hilbert curves;
       safe for initial use though optimum load balancing *not*
       guaranteed
    -  File — Reads a file (name given below in
       Sec. [sec:general-io-flags]; details of how to create in
       Sec. [sec:parall-issu])

General.LoadBalance (*string*) partitions/atoms
    Applies to Hilbert above; chooses whether to distribute atoms or partitions
    evenly between processors (you are *strongly* recommended to use atoms)

    *default*: atoms

General.ManyProcessors (*boolean*)
    Applies to Hilbert above; chooses method for parallelising Hilbert curve work;
    “many” processors here probably means more than two

    *default*: T

General.MaxAtomsPartition (*integer*)
    Applies to Hilbert above; specifies maximum number of atoms
    allowed in a partition; triggers extra level of recursion in
    partitioner

    *default*: 34

General.NPartitions[X/Y/Z] (*integer*)
    Allows the user to specify the number of partitions in x, y and z
    directions

    *default*: 0 (i.e. use Hilbert partitioning, above)

General.NewRun (*boolean*)
    Switches between new run and restart (N.B. restart has *not* been implemented yet)

    *default*: T

General.LoadL (*boolean*)
    Specifies whether to load a previous L matrix from files

    *default*: F

General.LoadRho (*boolean*)
    Specifies whether to load a previous charge density from files

    *default*: F

General.NetCharge (*real*)
    Specifies net charge on unit cell; implemented rather crudely with
    a neutralising background charge assumed. Note that a *positive*
    value indicates *excess* electrons

    *default*: 0.0

General.EwaldAccuracy (*real*)
    Accuracy for ewald sum (in Ha/atom)

    *default*: :math:`10^{-10}`

General.TimeThreshold (*real*)
    Minimum time for a timer to be printed (in seconds)

    *default*: :math:`0.001`

General.vdWDFT (*boolean*)
    Selects vdW DF

    *default*: F

General.DFT\_D2 (*boolean*)
    Selects DFT-D2

    *default*: F

General.MaxTime (*real*)
    Maximum wall time for calculation in seconds. Conquest will exit
    gracefully on completion of an ionic relaxation/MD step

    *default*: 0.0

.. _input_tags_atomic_spec:

Atomic Specification
--------------------

ChemicalSpeciesLabel (*block*)
    Lists all atomic species used in the calculation. Format:

    | ``1 atomic_mass1 element_type1``
    | ``2 atomic_mass2 element_type2``
    | ``...``
    | ``n atomic_mass_n_ element_type_n``

    1-–n are integer numbers used in the coordinate file to identify
    atomic species (see :ref:`coordinate-file` for more details).
    (don’t forget to end the block with %endblock ChemicalSpeciesLabel)

Atom.ValenceCharge (*real*)
    Valence charge of species (e.g. 4 for carbon, 6 for oxygen)

    *default*: 0.0

Atom.NumberOfSupports (*integer*)
    Number of support functions per atom for a species. Don’t confuse
    support functions and PAOs ! Support functions can be expanded in
    a basis set of PAOs or blips

    *default*: none

Atom.SupportFunctionRange (*real*)
    Confinement radius for the support functions for a given species

    *default*: none (though will default to PAO radius if PAO basis is chosen)

Atom.MultisiteRange (*real*)
    Range for multi-site support functions (the PAOs on all atoms
    within this range will be included in the support function)

    *default*: 0.0

Atom.LFDRange (*real*)
    Range for local filter diagonalisation (the Hamiltonian and
    overlap matrix elements from all atoms within this range will be
    included in the cluster diagonalisation)

    *default*: 0.0

.. _io_general_tags:

Input-Output General Tags
-------------------------

IO.Title (*string*)
    Title for run

    *default*: none

IO.Partitions (*string*)
    Name for file containing distribution of partitions over processors
    (generated by accompanying utilities)

    *default*: ``make_prt.dat``

IO.WriteOutToFile (*boolean*)
    Specifies whether the main output information is written to standard output
    or to a file

    *default*: T

IO.OutputFile (*string*)
    Name for the main output file

    *default*: ``Conquest_out``

IO.DumpL (*boolean*)
    Whether to write the auxiliary matrices L to file at each self-consistent steps

    *default*: T

IO.TimingOn (*boolean*)
    Whether time information will be measured and written to output

    *default*: F

IO.TimeAllProcessors (*boolean*)
    Specifies whether time information will be written for all processors or just
    for the input/output process (the default)

    *default*: F

IO.WriteTimeFile (*boolean*)
    Whether time files are written or not. This flag will be ignored if
    ``IO.TimeAllProcessors`` is true, in which case time files are always written.

    *default*: T

IO.TimeFileRoot (*string*)
    Root to be used in the time files, with an extension indicating the processor
    number, e.g. ``.001``

    *default*: ``time``

Atomic Coordinates
------------------

IO.Coordinates (*string*)
    Specifies the file with atomic coordinates. See Section [sec:atomic-positions]
    for details on the file format

    *default*: none

IO.FractionalAtomicCoords (*boolean*)
    Specifies whether fractional or absolute (Cartesian) coordinates are used
    in the coordinate file

    *default*: T

IO.PdbIn (*boolean*)
    Switches between the   coordinate file format (F) and PDB format (T)

    *default*: F

IO.PdbAltLoc (*string*)
    In case of PDB files with multiple locations selects an alternate location.
    Values: A, B, etc., as listed in the pdb file. Note that if the keyword is present
    in the input file but no value is given, only the parts of the system without
    any alternate location specification will be taken into account

    *default*: none

IO.PdbOut (*boolean*)
    Format of the output coordinate file. Writes a PDB file if set to T. In that
    case, either the input must be in pdb format or a PDB “template” file needs to
    be specified (keyword General.PdbTemplate)

    *default*: F

IO.PdbTemplate (*string*)
    A file used as a template for writing out coordinate files in the PDB format,
    i.e., the output file will contain the same information as the template, only
    the atomic coordinates will be overwritten. If the input file is in PDB format,
    it will also be used as the template, although this can still be
    overwritten with this keyword

    *default*: coordinate file

Levels of Output
----------------

The overall level of output is controlled by **IO.Iprint** and can be
fine-tuned with the other IO.Iprint keywords. These are by default set
to the value of iprint, but that will be over-ridden if setting them
explicitly. For instance, IO.Iprint could be set to 0, but IO.Iprint\_MD
could be set to 2 giving more extensive information about atomic
movements but little other information.

N.B. At beta release, these levels of output are still being tuned;
level 0 is reliable, and generally fairly minimal.

IO.Iprint (*integer*)
    The amount of information printed out to the output file
    The larger the value the more detailed the output is.

    | 0 Basic information about the system and the run
    | 1 Breakdown of energies, and details of the SCF cycle
    | 2 Matrix range info, matrix multiplication details (covering set), partition details and general parallelisation info.
    | 3 Subroutines called, messages upon entering/quitting subroutines
    | 4 Details including internal variables of subroutines
    | 5 Don’t do this.


    *default*: 0

IO.Iprint_init (*integer*)
    The initialisation process

    *default*: value of ``IO.Iprint``

IO.Iprint\_mat (*integer*)
    Matrix operations

    *default*: value of ``IO.Iprint``

IO.Iprint\_ops (*integer*)
    Creation of operators H and S

    *default*: value of ``IO.Iprint``

IO.Iprint\_DM (*integer*)
    Density matrix

    *default*: value of ``IO.Iprint``

IO.Iprint\_SC (*integer*)
    Self-consistency

    *default*: value of ``IO.Iprint``

IO.Iprint\_minE (*integer*)
    Energy minimisation

    *default*: value of ``IO.Iprint``

IO.Iprint\_MD (*integer*)
    Molecular dynamics

    *default*: value of ``IO.Iprint``

IO.Iprint\_index (*integer*)
    Indexing routines

    *default*: value of ``IO.Iprint``

IO.Iprint\_gen (*integer*)
    General (not covered by other areas)

    *default*: value of ``IO.Iprint``

IO.Iprint\_pseudo (*integer*)
    Pseudopotentials

    *default*: value of ``IO.Iprint``

IO.Iprint\_basis (*integer*)
    Basis set

    *default*: value of ``IO.Iprint``

IO.Iprint\_intgn (*integer*)
    Integration on the grid (not used at present)

    *default*: value of ``IO.Iprint``

IO.Iprint\_time (*integer*)
    Timing information

    *default*: value of ``IO.Iprint``

Basis Set
---------

Basis.BasisSet (*string*)
    values: blips/PAOs

    Selects the basis set in which to expand the support functions (localised orbitals).

    Options:

    -  PAOs — Pseudo-atomic orbitals :cite:`d-Artacho1999`

    -  blips (default) — B-splines :cite:`d-Hernandez1997`

    *default*: PAOs

Basis.LoadBlip (*boolean*)
    Load blip or PAO coefficients from file. If set to T, for blips the code will
    look for a set of files containing blip coefficients, which is taken to be
    ``blip_coeffs.nnn``, where ``nnn`` is processor number (padded with zeroes);
    for PAOs, the code will look for a *single* file which is ``supp_pao.dat``
    by default, but can be set with ``Basis.SupportPaoFile``

    *default*: F

Basis.SupportPaoFile (*string*)
    Specifies filename for PAO coefficients

    *default*: ``supp_pao.dat``

Basis.UsePulayForPAOs (*boolean*)
    Determines whether to use Pulay DIIS for minimisation of PAO basis coefficients

    *default*: F

Basis.PaoKspaceOlGridspace (*real*)
    Determines the reciprocal-space grid spacing for PAO integrals

    *default*: 0.1

Basis.PaoKspaceOlCutoff (*real*)
    Determines the cutoff for reciprocal-space grid spacing for PAO integrals

    *default*: 1000.0

Basis.PAOs\_StoreAllAtomsInCell (*boolean*)
    Determines whether coefficients for all atoms in cell are stored on each
    processor (improves speed but potentially memory expensive, particularly with
    large systems) or only local atom coefficients (increases communication overhead)

    *default*: T

Basis.SymmetryBreaking (*boolean*)
    Determines whether symmetry-breaking assignment of PAOs to support functions
    is allowed. In general, it is *highly* recommended that all atoms have sufficient
    support functions to span the space of angular momenta used in PAOs
    (i.e. :math:`2l+1` support functions for each :math:`l` channel used for PAOs);
    reducing the number potentially results in symmetry breaking and unphysical behaviour

    *default*: F

Basis.PaoNormFlag (*integer*)
    Determines whether PAOs are normalised

    *default*: 0

Basis.ReadSupportSpec (*boolean*)
    Should a file mapping PAOs to support functions *for initialisation of blip
    coefficients* be read ? See Sec. [sec:other:-basis-spec] for details of file

    *default*: F

Basis.SupportSpecFile (*string*)
    Name of file to be read (see ``Basis.ReadSupportSpec``)

    *default*: ``support.dat``

Basis.TestBasisGradients (*boolean*)
    Chooses whether gradients of energy with respect to basis function coefficients
    should be tested (using numerical vs. analytical gradients). **WARNING :** this
    produces large amounts of data

    *default*: F

Basis.TestBasisGradTot (*boolean*)
    Test total gradient ?

    *default*: F

Basis.TestBasisGradBoth (*boolean*)
    Test both S- and H-derived gradients (i.e. gradients arising from change of
    S or H when support functions vary) ?

    *default*: F

Basis.TestBasisGrad\_S (*boolean*)
    Test S-derived gradient ?

    *default*: F

Basis.TestBasisGrad\_H (*boolean*)
    Test H-derived gradient ?

    *default*: F

Basis.PAOs\_OneToOne (*boolean*)
    Assign PAOs to individual support functions (implies no support function optimisation)

    *default*: F

.. _integration-grid-spacing:

Integration Grid Spacing
------------------------

Grid.GridCutoff (*real*)
    Defines the spacing of the *integration* grid (though for a blip calculation
    must be at least twice as fine as blip grid, and will be adjusted). Note that
    the value chosen will automatically be forced to be a factor of 3, 4 and 5 only
    (to fit with default FFT routines)

    Default: 20 Ha.

Grid.PointsAlong[X/Y/Z] (*integer*)
    Grid points along x (y,z). Overwrites the values set by **Grid.GridCutoff**.
    The default FFT code requires that the number of grid points have prime
    factors of 2, 3 or 5

    *default*: 0

Grid.InBlock[X/Y/Z] (*integer*)
    This is the size of a grid point block (i.e., how many grid points are in one
    block in the x (y,z) direction), which must be 3, 4 or 5

    *default*: 4

Grid.ReadBlocks (*boolean*)
    If specified, the code reads information about blocks from the file make\_blk.dat

    *default*: F

Minimising Energy
-----------------

minE.VaryBasis (*boolean*)
    Chooses whether or not basis coefficients should be varied to minimise the
    total energy

    *default*: F

minE.SelfConsistent (*boolean*)
    Determines whether or not self-consistency cycles are imposed between charge
    density and potential

    *default*: T

minE.MixedLSelfConsistent (*boolean*)
    Determines whether or not to perform self-consistent cycle at the same time
    as energy minimisation with respect to L

    *default*: F

minE.EnergyTolerance (*real*)
    Fractional tolerance for energy on minimisation of support function coefficients

    *default*: 1\ :math:`\times`\ 10\ :math:`^{-5}`

minE.LTolerance (*real*)
    Tolerance on *residual* in O(N) minimisation

    *default*: 1\ :math:`\times`\ 10\ :math:`^{-7}`

minE.SCTolerance (*real*)
    Tolerance on *residual* in self-consistency

    *default*: 1\ :math:`\times`\ 10\ :math:`^{-6}`

minE.SupportVariations (*integer*)
    Maximum number of support-function iterations

    *default*: 20

minE.PreconditionBlips(*boolean*)
    Should blip variation be pre-conditioned? Pre-conditioning is (at present)
    more memory-intensive than it should be, but is efficient

    *default*: F

minE.GlobalTolerance (*boolean*)
    Are the convergence criteria applied to minimisation summed over the whole
    system, or per atom?

    *default*: T

Charge Self-Consistency
-----------------------

SC.LinearMixingSC (*boolean*)
    Should Pulay mixing be used? It is recommended that this is always used

    *default*: T

SC.LinearMixingFactor (*real*)
    Amount of output charge density which is mixed into new charge

    *default*: 0.5

SC.LinearMixingFactor\_SpinDown (*real*)
    Amount of output charge density which is mixed into new charge for spin down channel.

    *default*: value of **SC.LinearMixingFactor**

SC.LinearMixingEnd (*real*)
    Tolerance for end of Pulay mixing

    *default*: self-consistency tolerance

SC.LateStageReset (*integer*)
    If using GR-Pulay, how often is residual calculated fully (rather than interpolated) ?

    *default*: 5

SC.MaxIters (*integer*)
    Maximum self-consistency iterations

    *default*: 50

SC.MaxEarly (*integer*)
    Maximum early-stage iterations

    *default*: 3

SC.MaxPulay (*integer*)
    Number of iterations stored and mixed during Pulay mixing

    *default*: 5

SC.ReadAtomicDensityFile (*string*)
    Filename for radial tables of atomic density (*rarely* used: normally generated from PAOs)

    default:

SC.AtomicDensityFlag (*string*)
    values: pao/read

    Flag determining how atomic densities should be found

    *default*: pao

SC.KerkerPreCondition (*boolean*)
    Flag determining if Kerker precondition is to be used.

    *default*: F

SC.KerkerFactor (*real*)
    Wave-vector magnitude used in Kerker preconditioning, it is :math:`q_0` from
    the factor :math:`q^2 / \left(q^2 + q_0^2\right)`

    *default*: 0.1

SC.WaveDependentMetric (*boolean*)
    Flag determining if wave-dependent metric is to be used in Pulay mixing.

    *default*: F

SC.MetricFactor (*real*)
    Wave-vector magnitude used by wave-dependent metric method, it is :math:`q_1`
    from the factor :math:`\left(q^2 + q_1^2\right) / q^2`.

    *default*: 0.1

Density Matrix
--------------

DM.SolutionMethod (*string*)
    values: ordern/diagon

    Selects the method for finding the ground state density matrix. This can currently
    be either diagonalisation (diagon: minimising the energy with respect to the
    density matrix elements) or an O(N) method (ordern a combination of the
    techniques of Li et al. :cite:`d-Li1993` and Palser and Manolopoulos :cite:`d-Palser1998`.)

    *default*: ordern

DM.L\_range (*real*)
    Cutoff applied to L matrix (total energy will converge with increasing range;
    suggested minimum for O(N) calculations is twice largest support function range;
    see Sec. [sec:find-dens-matr] for more details)

    *default*: 1.0

DM.LVariations (*integer*)
    Maximum number of variations performed in search for ground-state density matrix

    *default*: 50

DM.MaxPulay (*integer*)
    Maximum number of iterations stored for Pulay minimisation

    *default*: 5

DM.MinPulayStepSize (*real*)
    Minimum allowed step size for Pulay minimisation in Energy minimisation stage
    of the calculation. Note that the actual step size is calculated by  automatically,
    but will be constrained within the range defined by ``DM.MinPulayStepSize``
    and ``DM.MaxPulayStepSize``. Not to be confused with the Pulay mixing step
    size for charge self-consistency.

    *default*: 0.001

DM.MaxPulayStepSize (*real*)
    Maximum allowed step size for Pulay minimisation in Energy minimisation stage
    of the calculation. Not to be confused with the Pulay mixing step size
    for charge self-consistency.

    *default*: 0.1

DM.LinTol (*real*)
    Tolerance on linearity required before switching to Pulay minimisation

    *default*: 0.1

DM.InvSTolerance (*real*)
    Tolerance on iterative minimisation to find S\ :math:`^{-1}`. If
    :math:`\Omega = \mathrm{Tr}[(I-TS)^2]/N_{\mathrm{orbitals}}` is above this,
    identity will be used

    *default*: 0.01

DM.InvSMaxSteps (*integer*)
    Sets the maximum number of iterations for finding S\ :math:`^{-1}`

    *default*: 100

DM.InvSDeltaOmegaTolerance (*real*)
    Tolerance which determines when the iterative minimisation to find S\ :math:`^{-1}`
    should finish. :math:`\delta\Omega_n = N_{\mathrm{orbitals}} (\Omega_n - \Omega_{n-1})`,
    where :math:`\Omega` is defined in description for ``DM.InvSTolerance``. This parameter
    differs from ``DM.InvSTolerance`` in that the iterative S\ :math:`^{-1}` finder
    will end iteration when :math:`\delta\Omega` is less than or equal to
    ``DM.InvSDeltaOmegaTolerance``, while ``DM.InvSTolerance`` determines whether
    to reset S\ :math:`^{-1}` to identity (i.e. whether a satisfactory S\ :math:`^{-1}`
    has been found) based on the final :math:`\Omega` produced from the iterative loop

    *default*: 0.0001

DM.ConstantMu (*boolean*)
    Switches between fixed Fermi level (T) and fixed number of electrons (F). You
     are *strongly* recommended to leave at default

    *default*: F

DM.mu (*real*)
    Value of Fermi level for fixed Fermi level calculations

    *default*: 0.0

Diagonalisation
---------------

Diag.NumKpts (*integer*)
    Number of all k-points. No symmetry is applied.

    *default*:

Diag.Kpoints (*block*) 
    Lists fractional coordinates and weights of all k-points: ``x_fract y_fract z_fract weight``
    Generates the Monkhorst-Pack mesh, an equally spaced mesh of k-points.

    *default*:

Diag.MPMesh (*boolean*)
    Switches on/off the Monkhorst-Pack mesh. Note: if this keyword is present in
    the input file, the keyword **Diag.NumKpts** and the block **Kpoints** will
    be ignored.

    *default*:

Diag.MPMesh[X/Y/Z] (*integer*)
    Specifies the number n of k-points along the x(y,z) axis.

    *default*: 1

Diag.ProcRows (*integer*)

    *default*:

Diag.ProcCols (*integer*)

    *default*:

Diag.BlockSizeR (*integer*)

    *default*:

Diag.BlockSizeC (*integer*)
    R ... rows, C ... columns
    These are ScaLAPACK parameters, and can be set heuristically by the code. Blocks
    are sub-divisions of matrices, used to divide up the matrices between processors.
    The block sizes need to be factors of the square matrix size
    (i.e. :math:`\sum_{\mathrm{atoms}}\mathrm{NSF(atom)}`). A value of 64 is considered
    optimal by the ScaLAPACK user’s guide. The rows and columns need to multiply
    together to be less than or equal to the number of processors. If ProcRows
    :math:`\times` ProcCols :math:`<` number of processors, some processors will be left idle.

    *default*:

Diag.MPShift[X/Y/Z] (*real*)
    Specifies the shift *s* of k-points along the x(y,z) axis, in fractional
    coordinates.

    *default*: 0.0

Diag.SmearingType (*integer*)
    Specifies the type of smearing used

    +-----+---------------------+
    | 0   | Fermi-Dirac         |
    +-----+---------------------+
    | 1   | Methfessel-Paxton   |
    +-----+---------------------+

    *default*: 0

Diag.kT (*real*)
    Smearing temperature

    *default*: 0.001

Diag.MPOrder (*integer*)
    Order of Bessel function approximation to delta-function used in Methfessel-Paxton smearing

    *default*: 0

Diag.GaussianHeight (*real*)
    The height of Gaussian function used to determine the width of Methfessel-Paxton
     approximation to delta-function (see section [sec:methf-paxt-smear])

    *default*: 0.1

Diag.EfStepFiness (*real*)
    Parameter controlling the finness of the Fermi energy search step used in
    Methfessel-Paxton smearing method (see section [sec:methf-paxt-smear])

    *default*: 1.0

Diag.NElecLess (*Real*)
    The number of electrons to subtract from the total number of electrons in each
    spin channel, which gives the starting point for searching the lower bound for
    Fermi energy. Used in Methfessel-Paxton smearing method
    (see section [sec:methf-paxt-smear])

    *default*: 10.0

Diag.KProcGroups (*integer*)
    Number of k-point processor groups for k-point parallelisation
    (see section [sec:scal-proc-grid])

    *default*: 1

Diag.ProcRows (*integer*)
    Number of rows in the processor grid for SCALAPACK within each k-point processor
    group (see section [sec:scal-proc-grid])

    *default*: Determined automatically

Diag.ProcCols (*integer*)
    Number of columns in the processor grid for SCALAPACK within each k-point
    processor group (see section [sec:scal-proc-grid])

    *default*: Determined automatically

Moving Atoms
------------
AtomMove.TypeOfRun (*string*)
    values: static/cg/md

    Options:

    static — Single point calculation

    cg — Structure optimisation by conjugate gradients

    md — Velocity Verlet algorithm

    *default*: static

AtomMove.NumSteps (*integer*)
    Maximum number of steps for a structure optimisation or molecular dynamics run

    *default*: 100

AtomMove.MaxForceTol (*real*)
    The structure optimisation will stop when the maximum force component is less
    than **MD.MaxForceTol**

    *default*: 0.0005 Ha/bohr

AtomMove.Timestep (*real*)
    Time step for molecular dynamics

    *default*: 0.5

AtomMove.AppendCoords (*boolean*)
    Chooses whether to append coordinates to ``UpdatedAtoms.dat`` during atomic
    movement (T) or to overwrite (F)

    *default*: T

AtomMove.OutputFreq (*integer*)
    Frequency of output of information. *Not properly implemented*

    *default*: 50

AtomMove.WriteXSF *(boolean*)
    Write atomic coordinates to ``trajectory.xsf`` for ``AtomMove.TypeOfRun = md`` or ``cg``,
    every ``AtomMove.OutputFreq`` steps

    *default*: F

AtomMove.TestForces (*boolean*)
    Flag for testing forces with comparison of analytic and numerical calculations.
    Can produce *large* amounts of output

    *default*: F

AtomMove.TestAllForces (*boolean*)
    Switch to test *all* force contributions or not

    *default*: F

AtomMove.OptCell (*boolean*)
    Turns on conjugate gradient relaxation of the simulation box dimensions a, b
    and c. Note that AtomMove.TypeOfRun must also be set to cg.

    *default*: F

AtomMove.OptCell.EnTol (*real*)
    Sets the energy tolerance to terminate the conjugate gradient relaxation.

    *default*: 1\ :math:`\times`\ 10\ :math:`^{-5}`

AtomMove.OptCell.Constraint (*string*)
    Applies a constraint to the relaxation.

    none: Unconstrained relaxation.

    *Fixing a single cell dimension:*

    a: Fix the x-dimension of the simulation box

    b: Fix the y-dimension of the simulation box

    c: Fix the z-dimension of the simulation box

    *Fixing multiple cell dimensions:*

    any combination of the above separated by a space character. e.g: "a b" fixes
    both the x and y dimensions of the simulation box

    *Fixing Ratios:*

    Any combination of a, b or c separated by a "/" character. e.g "c/a" fixes
    the initial ratio of the z-dimension to the x-direction.

    *Global scaling factor:*

    volume: minimize the total energy by scaling each simulation box dimension by
    the same global scaling factor. Search directions are set by the mean stress.

AtomMove.TestSpecificForce (*integer*)
    Label for which force contribution to test. Note that for PAOs non-local Pulay
    and Hellman-Feynman forces are found together as part of the HF calculation;
    :math:`\phi` Pulay refers to changes in :math:`\phi(\mathbf{r})` when atoms move,
    while S Pulay refers to changes in S when atoms move; see Sec. [sec:forces] and
    references therein for more information. Options:

    1 Total
    2 Total Hellman-Feynman
    3 Total Pulay
    4 Non-SC Correction
    5 Non-local :math:`\phi` Pulay
    6 KE :math:`\phi` Pulay
    7 Local :math:`\phi` Pulay
    8 S Pulay

    *default*: 1

AtomMove.TestForceDirection (*integer*)
    Direction in which atom will be moved (1=x; 2=y; 3=z)

    *default*: 1

AtomMove.TestForceAtom (*integer*)
    Atom to move

    *default*: 1

AtomMove.TestForceDelta (*real*)
    Distance atom will be moved for numerical evaluation of force

    *default*: 10\ :math:`^{-5}` bohr

AtomMove.RestartRun (*boolean*)
    Selects MD restart run (in the case of F, MD step is reset)

    *default*: F

AtomMove.ReuseL (*boolean*)
    Selects the use of L-matrix in MD run

    *default*: F

AtomMove.ReuseInvS (*boolean*)
    Selects the use of T-matrix in MD run

    *default*: F

AtomMove.SkipEarlyDM (*boolean*)
    Selects the skip of earlyDM calculation in MD run

    *default*: F

AtomMove.McWeenyFreq (*integer*)
    Number of steps to apply to McWeeny calculation (with “AtomMove.ReuseL T”)

    *default*:

AtomMove.ExtendedLagrangian (*boolean*)
    Selects XL-BOMD (with “AtomMove.ReuseL T”)

    *default*: F

AtomMove.FixCentreOfMass (*boolean*)
    Remove the centre of mass velocity at every time step

    *default*: T


Molecular Dynamics
------------------

MD.Ensemble (*string*)
    values: nve/nvt/npt

    The molecular dynamics ensemble

    *default*: nve

MD.Thermostat (*string*)
    values: none/nhc/berendsen

    Thermostat type

    ``none``
        No thermostat (used for calculating temperature only)
    ``nhc``
        Nosé-Hoover chain
    ``berendsen``
        Berendsen weak coupling thermostat

    *default*: none

MD.tauT (*real*)
    Coupling time constant for thermostat. Required for Berendsen thermostat, or
    if ``MD.CalculateXLMass = T``

    *default*: 1.0

MD.nNHC (*integer*)
    Number of Nosé-Hoover thermostats in chain

    *default*: 5

MD.CellNHC (*boolean*)
    Use a separate Nosé-Hoover chain for thermostating the unit cell (NPT only)

    *default*: F

MD.NHCMass (*blocks*)
    :math:`<n1> <n2> <n3> \ldots`
    Masses of NHC heat baths

    *default*: 1 1 1 1 1

MD.CellNHCMass (*block*)
    :math:`<n1> <n2> <n3> \ldots`
    Masses of NHC heat baths for unit cell

    *default*: 1 1 1 1 1

MD.Barostat (*string*)
    values: none/berendsen/iso-mttk/ortho-mttk/mttk

    Barostat type

    ``none``
        No barostat (used for calculating pressure only)
    ``berendsen``
        Berendsen weak coupling barostat
    ``iso-mttk``
        MTTK barostat with isotropic cell fluctuations
    ``ortho-mttk``
        MTTK barostat with orthorhombic cell fluctuations
    ``mttk``
        MTTK barostat with fully flexible cell (currently discards off-diagonal
        stress tensor components)

    *default*: none

MD.TargetPressure (*real*)
    External pressure for NPT molecular dynamics

    *default*:

MD.BulkModulusEst (*real*)
    Bulk modulus estimate for system. Only necessary for Berendsen weak pressure
    coupling (``MD.Barostat = berendsen`` or ``MD.BerendsenEquil > 0``)

    *default*: 100

MD.tauP (*real*)
    Coupling time constant for barostat. Required for Berendsen barostat, or if
    MD.CalculateXLMass = T

    *default*: 10.0 (Berendsen) or 100.0 (MTTK)

MD.BoxMass (*real*)
    Mass of box for extended system formalism (MTTK barostats)

    *default*: 1

MD.CalculateXLMass (*boolean*)
    Calculate the mass of the extended system components (thermostats, barostat)
    using the MTTK formulae.

    *default*: F

MD.nYoshida (*integer*)
    values: 1/3/5/7/15/25/125/625

    Order of Yoshida-Suzuki integration

    *default*: 1

MD.nMTS (*integer*)
    Number of time steps in inner loop of MTS scheme

    *default*: 1

MD.BerendsenEquil (*integer*)
    Equilibrate the system for :math:`n` steps using Berendsen weak coupling

    *default*: 0

MD.ThermoDebug (*boolean*)
    Print detailed information about thermostat and extended variables in ``thermostat.dat``

    *default*: F

MD.BaroDebug (*boolean*)
    Print detailed information about barostat and extended variables in ``barostat.dat``

    *default*: F

Spin Polarisation
-----------------

Spin.SpinPolarised (*boolean*)
    Determines if the calculation is spin polarised (collinear) or non-spin polarised.

    *default*: F

Spin.FixSpin (*boolean*)
    Determines if spin populations are to be fixed. Only read if **Spin.FixPolarised** is set.

    *default*: F

Spin.NeUP (*real*)
    Total number of electrons in spin up channel at start of calculation.

    *default*: 0.0

Spin.NeDN (*real*)
    Total number of electrons in spin down channel at start of calculation.

    *default*: 0.0

DeltaSCF
--------

flag\_DeltaSCF (*boolean*)
    Selects delta SCF calculation

    *default*:

DeltaSCF.SourceLevel (*integer*)
    Eigenstate number to remove electron from (source)

    *default*:

DeltaSCF.TargetLevel (*integer*)
    Eigenstate number to promote electron to (target)

    *default*:

DeltaSCF.SourceChannel (*integer*)
    Spin channel for electron source

    *default*:

DeltaSCF.TargetChannel (*integer*)
    Spin channel for electron target

    *default*:

DeltaSCF.SourceNFold (*integer*)
    Allows selection of more than one level for excitation source (N-fold)

    *default*:

DeltaSCF.TargetNFold (*integer*)
    Multiplicity of target (N-fold)

    *default*:

DeltaSCF.LocalExcitation (*boolean*)
    Select an excitation localised on a group of atoms

    *default*:

DeltaSCF.HOMOLimit (*integer*)
    How many states down from HOMO to search for localised excitation

    *default*:

DeltaSCF.LUMOLimit (*integer*)
    How many states up from LUMO to search for localised excitation

    *default*:

DeltaSCF.HOMOThresh (*real*)
    (*please fill in*)

    *default*:

DeltaSCF.LUMOThresh (*real*)
    Threshold for identifying localised excitation (sum over square moduli of coefficients)

    *default*:

Constrained DFT (cDFT)
----------------------

cDFT.Perform\_cDFT (*boolean*)
    Selects cDFT operation

    *default*:

cDFT.Type (*integer*)
    values: 1 or 2

    Selects constraint to be for absolute charge on groups (1) or difference between two groups (2)

    *default*:

cDFT.MaxIterations (*integer*)
    Maximum iterations permitted

    *default*:

cDFT.Tolerance (*real*)
    Tolerance on charge

    *default*:

cDFT.NumberAtomGroups (*integer*)
    Number of groups of atoms

    *default*:

cDFT.AtomGroups (*block*)
    Block with each line specifying: Number of atoms, target charge, label for
    block. For each line, there should be a corresponding block with the appropriate
    label; the block consists of a list of atom numbers for the atoms in the group

vdW-DF
------

vdWDFT.LDAFunctionalType (*string*)
    Selects LDA functional to use with vdW-DF

    *default*:

DFT-D2
------

DFT-D2\_range (*real*)
    DFT-D2 cutoff range (bohr)

    *default*:

XL-BOMD
-------

XL.Kappa (*real*)
    Value of kappa

    *default*:

XL.PropagateX (*boolean*)
    Selects the corrected XL-BOMD

    *default*:

XL.PropagateL (*boolean*)
    Selects the original XL-BOMD

    *default*:

XL.Dissipation (*boolean*)
    Selects the addition of dissipative force

    *default*:

XL.MaxDissipation (*integer*)
    (*please fill in*)

    *default*:

XL.Integrator (*string*)
    Selects the Verlet method or velocity Verlet method

    *default*:

Advanced and obscure tags
-------------------------

.. _advanced_general_tags:

General
*******

General.LoadInvS (*boolean*)
    Selects loading of inverse S matrix from previous step (not
    recommended)

    *default*: F

General.NeutralAtomProjector (*boolean*)
    Selects projector expansion of neutral atom potential; still in
    development.  Only for expert use.  (Allows specification of
    maximum l value for projectors and list of number of projectors
    for each l value.)

    *default*: F

General.PAOFromFiles (*boolean*)
    Allows you to give explicit file name for .ion files in atom block

    *default*: F

General.MaxTempMatrices (*integer*)
    Allows user to increase number of temporary matrices; sometimes
    required for wavefunction output.

    *default*: 100

General.EwaldAccuracy (*real*)
    Accuracy required for Ewald sum

    *default*:1\ :math:`\times`\ 10\ :math:`^{-10}`

General.CheckDFT (*boolean*)
    Calculates DFT energy using output density

    *default*: F

General.AverageAtomicDiameter (*real*)
    Related to space-filling

    *default*: 5.0

General.GapThreshold (*real*)
    Related to space-filling

    *default*: 2.0*(largest support radius)

General.only_Dispersion (*boolean*)
    Selects only DFT\_D2 calculation (no electronic structure etc)

.. _advanced_atomic_spec_tags:

Atomic Specification
********************

Atom.SupportGridSpacing (*real*)
    The spacing of the blip grid (if using). Equivalent (under certain
    circumstances) to a maximum g-vector of
    :math:`\pi`/**SupportGridSpacing**
    plane wave cutoff as region radius and L matrix radius go to infinity. *Not used for PAO
    calculations*. See Sec. [sec:blip-basis] for more
    information. N.B. Grid.GridCutoff will be reset to *at least* half
    SupportGridSpacing if too small.

    *default*: none

Atom.NonLocalFactor  (*real*)
    This is an adjustment factor: the Hamiltonian range is (strictly)
    2 :math:`\times` (support function radius + non-local projector
    radius). However, generally without affecting the results, the
    Hamiltonian range can be set to 2  :math:`\times` (support function
    radius + non\_local\_factor\ :math:`\times` non-local projector radius). If you
    have non\_local\_factor = 1.0 then you get the full range, if 0.0
    then the same range as the S matrix.

    *default*: 0.0

Atom.InvSRange  (*real*)
    Range of inverse S matrix (though actual matrix range is twice
    this for consistency with S matrix range).

    *default*: support function range

Atom.SpinNeUp (*real*)
    Specify the population of spin-up electrons for setting initial
    spin state of atomic densities

    *default*: 0.0

Atom.SpinNeDn (*real*)
    Specify the population of spin-down electrons for setting initial
    spin state of atomic densities

    *default*: 0.0

.. bibliography:: references.bib
    :cited:
    :labelprefix: E
    :keyprefix: e-
    :style: unsrt
