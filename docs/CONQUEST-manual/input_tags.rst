.. _input_tags:

==========
Input tags
==========

We have broken down the input tags based on the areas of the code
where they apply.  For each tag, a default is given.  Types of value
are specified as: *integer*;
*real*; *boolean*; or *string* (optA/optB are given for string options).

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
    Chooses method for partitioning (read from file or use dynamic partitioner based on Hilbert curve)

    *default*: Hilbert

    Options:
    
    -  Hilbert (default) — Automatic partitioning using Hilbert curves;
       safe for initial use though optimum load balancing *not*
       guaranteed
    -  File — Reads a file (name given below in
       Sec. [sec:general-io-flags]; details of how to create in
       Sec. [sec:parall-issu])

General.LoadBalance (*string*) partitions/atoms
    Applies to Hilbert above; chooses whether to distribute atoms or partitions evenly between processors (you are *strongly* recommended to use atoms)

    *default*: atoms

General.ManyProcessors (*boolean*)
    Applies to Hilbert above; chooses method for parallelising Hilbert curve work; “many” processors here probably means more than two

    *default*: T

General.MaxAtomsPartition (*integer*)
    Applies to Hilbert above; specifies maximum number of atoms
    allowed in a partition; triggers extra level of recursion in
    partitioner 

    *default*: 34

General.NPartitionsX/Y/Z (*integer*)
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

Atomic Specification
--------------------

ChemicalSpeciesLabel (*block*)
    Lists all atomic species used in the calculation. Format:
    ``1 atomic_mass1 element_type1 2 atomic_mass2 element_type2 ... n atomic_mass_n_ element_type_n``
    1 – n are integer numbers used in the coordinate file to identify
    atomic species (see Section [sec:atomic-positions] for more details).

%endblock ChemicalSpeciesLabel
    | (don’t forget to end the block)

Atom.ValenceCharge
    | :math:`<f>`
    | Valence charge of species (e.g. 4 for carbon, 6 for oxygen)
    | *default*: 0.0

Atom.NumberOfSupports
    | :math:`<n>`
    | Number of support functions per atom for a species. Don’t confuse support functions and PAOs ! Support functions can be expanded in a basis set of PAOs or blips
    | *default*: none

Atom.SupportFunctionRange
    | :math:`<f>`
    | Confinement radius for the support functions for a given species
    | *default*: none (though will default to PAO radius if PAO basis is chosen)

Atom.SupportGridSpacing
    | :math:`<f>`
    | The spacing of the blip grid (if using). Equivalent (under certain circumstances) to a maximum g-vector of :math:`\pi`/**SupportGridSpacing** plane wave cutoff as region radius and L matrix radius go to infinity. *Not used for PAO calculations*. See Sec. [sec:blip-basis] for more information. N.B. Grid.GridCutoff will be reset to *at least* half SupportGridSpacing if too small.
    | *default*: none

Atom.InvSRange
    | :math:`<f>`
    | Range of inverse S matrix (though actual matrix range is twice this for consistency with S matrix range).
    | *default*: support function range

Atom.NonLocalFactor
    | :math:`<f>`
    | This is an adjustment factor: the Hamiltonian range is (strictly) 2 :math:`\times` (support function radius + non-local projector radius). However, generally without affecting the results, the Hamiltonian range can be set to 2 :math:`\times` (support function radius + non\_local\_factor\*non-local projector radius). If you have non\_local\_factor = 1.0 then you get the full range, if 0.0 then the same range as the S matrix.
    | *default*: 0.0

Advanced and obscure tags
-------------------------

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

    *default*: 1e-10
    
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

.. bibliography:: references.bib
    :cited:
    :labelprefix: E
    :keyprefix: e-
    :style: unsrt
