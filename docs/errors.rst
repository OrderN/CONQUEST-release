.. _errors:

==========================
Errors and troubleshooting
==========================

CONQUEST reports most fatal errors as text rather than as numbered error
codes.  This page lists common messages for which the immediate cause and a
safe first check can be identified.  It is not a complete catalogue of every
internal consistency check in the code.

Where a message contains ``<mode>``, ``<value>`` or ``<filename>``, the angle-
bracketed text stands for a value supplied at runtime; the surrounding text
matches the source diagnostic.

How CONQUEST reports problems
-----------------------------

Fatal errors are written to ``stdout`` and to the configured main output unit
before CONQUEST stops the MPI calculation.  The main output is normally
``Conquest_out``, but its destination can be changed; see :ref:`io_output_main`.
For an early startup failure, the output file might not yet have been opened.
Do not keep only the final MPI abort line: the CONQUEST message immediately
before it is normally more useful.

Warnings issued through CONQUEST's warning routine are written to
``Conquest_warnings`` and, depending on the ``IO.Iprint`` level, can also
appear in the main output.  Check ``Conquest_warnings`` after every run,
including runs that appear to finish normally.

A ``CQ.stop`` file requests that a running calculation stop when CONQUEST
next checks for it.  A stop caused by that file is deliberate and should not
be diagnosed as an error.

First checks after a failed run
-------------------------------

#. Preserve the complete error text, ``Conquest_out``,
   ``Conquest_warnings`` and any scheduler output before restarting.
#. Find the first CONQUEST error or warning, rather than relying only on the
   last message printed by MPI or the scheduler.
#. Check that the run directory contains the intended input, coordinate,
   ion/pseudopotential and restart files.
#. Check coordinate conventions, distance units and species numbers against
   :ref:`io_coords` and :ref:`input_atomic_spec`.
#. Identify whether the problem is an input error, a numerical convergence
   failure, or an external runtime failure before changing tolerances,
   iteration limits or the parallel layout.

Startup and input files
-----------------------

``We need Conquest_input to run !``
   CONQUEST could not open ``Conquest_input`` in the run directory.  Check the
   filename, working directory and read permissions.  The required input files
   and the main input structure are described in :ref:`io_files`.

``No coordinate file specified: please set with IO.Coordinates``
   No coordinate filename was supplied.  Set ``IO.Coordinates`` to the
   intended file.  A PDB input that cannot be opened instead gives
   ``Reading pdb file: file error``.  Check the coordinate filename and the
   selected format against :ref:`io_coords` and :ref:`input_coords`.

``Too few species in ChemicalSpeciesLabel:``
   The ``ChemicalSpeciesLabel`` block contains fewer entries than
   ``General.NumberOfSpecies``.  Make the count and block agree and ensure that
   the block is terminated correctly; see :ref:`input_atomic_spec`.

``Species specified greater than number in input file:``
   A species index in the coordinate file is greater than the number of
   species declared in the input.  Check the coordinate indices and the
   numbered entries in ``ChemicalSpeciesLabel``.  See
   :ref:`input_atomic_spec` for the block format.

``read_pseudopotential: file error``
   A file used by the legacy pseudopotential reader could not be opened.
   Check the filename printed immediately before this error, its location and
   its permissions.  Do not substitute an unrelated pseudopotential merely to
   pass this check; see :ref:`io_ion`, :ref:`input_general` and :ref:`basissets`.

``Functionals differ between pseudopotential files:``
   When ``General.FunctionalType`` is left at its default of zero, CONQUEST
   derives the functional from the first pseudopotential and finds a different
   one in another file.  Use a mutually consistent set of files; see
   :ref:`input_general` and :ref:`io_ion`.

``Functional in input file differs to pseudopotential:``
   ``General.FunctionalType`` differs from the functional stored in a
   pseudopotential.  The calculation stops unless
   ``General.DifferentFunctional`` explicitly permits the mismatch.  Prefer
   matching the input and pseudopotentials; an override does not make an
   incompatible pseudopotential valid.  See :ref:`input_general`.

``Functional in input file differs to pseudopotential but proceeding:``
   The same mismatch was found, but ``General.DifferentFunctional`` permits
   the run to continue with a warning.  Verify that the override is intended
   and that the input and pseudopotentials are physically compatible; see
   :ref:`input_general` and :ref:`io_ion`.

Coordinates and parallel layout
-------------------------------

``Expected fractional coordinates but many are greater than one:``
   This warning means that several values look Cartesian although
   ``IO.FractionalAtomicCoords T`` selects fractional coordinates.

``Expected Cartesian coordinates but many are less than one:``
   This warning is the converse heuristic: several values look fractional
   although ``IO.FractionalAtomicCoords F`` selects Cartesian coordinates.  In
   either case, verify ``IO.FractionalAtomicCoords``,
   ``General.DistanceUnits`` and the coordinate file itself.  The warning is a
   heuristic, so inspect the data rather than changing the flag automatically;
   see :ref:`io_coords` and :ref:`input_coords`.

``Atoms are too close to each other:``
   At least two non-ghost atoms are at or below the configured minimum initial
   separation.

``Atoms are too far apart! Minimum distance is``
   Even the smallest detected separation is above the configured maximum.
   For either distance error, first inspect duplicated atoms, cell dimensions,
   coordinate wrapping, ``IO.FractionalAtomicCoords`` and
   ``General.DistanceUnits``.  Do not disable the check until the geometry and
   units have been verified; see :ref:`io_coords` and :ref:`input_coords`.

``We must have at least one atom per process:``
   The run uses more MPI processes than atoms.  Rerun with no more processes
   than atoms; see :ref:`io_coords` for the coordinate file's atom count.

``More processors than partitions!``
   A manual ``General.NPartitionsX/Y/Z`` product is smaller than the MPI
   process count.  The related messages ``sfc_partitions_module: too few
   partitions created for number of processors`` and
   ``sfc_partitions_module: too few occupied partitions for number of
   processors`` mean that the generated or occupied partition count is also
   insufficient.  Reduce the process count or review the partition settings
   in :ref:`input_general`; do not assume that adding more processes will
   accelerate a small or sparse system.

Diagonalisation and convergence
-------------------------------

``Diag: proc grid product is zero:``
   A manually specified diagonalisation processor grid has a zero dimension.
   Check ``Diag.ProcRows`` and ``Diag.ProcCols`` against :ref:`input_diag`.

``Diag: proc grid product is too large:``
   ``Diag.ProcRows`` multiplied by ``Diag.ProcCols`` exceeds the processes in
   the applicable k-point process group.  Use a valid grid or let CONQUEST
   choose it.  See :ref:`gs_diag_para` and :ref:`input_diag`.

**Messages:** ``block_size_r not a factor of matrix size !``,
``block_size_c not a factor of matrix size !`` or
``Can't find a good block size: please set manually``

These errors concern the ScaLAPACK matrix distribution when matrix padding is
disabled.  Check ``Diag.PaddingHmatrix``, ``Diag.BlockSizeR`` and
``Diag.BlockSizeC`` against the restrictions in :ref:`gs_pad` and
:ref:`input_diag`.  Change them as a consistent layout, not as independent
tuning values.

``Code compiled without ELPA! Set Diag.UseELPA F``
   The input requests ELPA but the executable was built without ELPA support.
   Use ScaLAPACK by disabling ``Diag.UseELPA``, or use an executable built and
   linked with ELPA.  See :ref:`gs_diag_elpa` and :ref:`install_compile`.

``FindEvals: pzhegvx failed for mode <mode> with INFO=<value>``
   The selected generalized eigensolver returned a status that CONQUEST treats
   as fatal.  Despite the ``pzhegvx`` wording, this message can also follow an
   ELPA call.  Preserve the mode, ``INFO`` value and ``Diag.UseELPA`` setting,
   together with the parallel layout and full output.  Do not infer a physical
   cause from the number alone.  Check the diagonalisation setup in
   :ref:`gs_diag` before reporting a reproducible failure.

``ScaLAPACK pzhegvx warning, info=``
   The selected eigensolver returned ``INFO=2`` or ``INFO=4``.  CONQUEST
   permits these statuses to continue and records a warning.  This wording is
   also used when ``Diag.UseELPA`` is enabled, so interpret the status using
   the selected solver.  Check the resulting eigenvalues, occupations and
   convergence, and retain the exact status if reporting unexpected results;
   see :ref:`gs_diag`.

``SelfCon/PulayMixSC_spin: too many SCF iterations:``
   The self-consistency cycle reached ``SC.MaxIters`` without satisfying its
   convergence test, unless ``SC.ContinueOnSCFail`` was selected.  Inspect the
   residual history and verify that the basis, grid and physical setup are
   sensible before increasing the limit.  Mixing and convergence controls are
   discussed in :ref:`conv_scf` and :ref:`gs_scf`.

Atomic movement and restart files
---------------------------------

**Messages:** ``Step too small: safemin2 failed!``,
``Step too small: safemin_cell failed!`` or
``Step too small: safemin_full failed!``

An ionic, cell or combined line search reduced its trial step below the usable
threshold.  Check the starting geometry, forces, preceding energy changes and
line-search messages before changing optimisation controls.  See
:ref:`strucrelax` for the supported relaxation methods and their inputs.

**Messages:** ``Error: Too many SHAKE iterations !`` or
``Error: Too many RATTLE iterations !``

A rigid-bond position or velocity correction did not converge within its
iteration limit.  Check the constraint definitions, the initial constrained
bond lengths and the preceding dynamics before changing tolerances or limits.
See :ref:`moldyn` for molecular-dynamics setup.

``Fail in opening InfoGlobal.dat``
   The restart reader could not open the indexed ``InfoGlobal`` metadata file.
   Despite the diagnostic's ``.dat`` wording, check the actual filename printed
   immediately before it (such as ``InfoGlobal.i00``), the restart settings and
   read permissions.  See :ref:`md_restart`.

**Messages:** ``Fail in opening Lmatrix file in Binary.`` or
``Fail in opening ASCII Lmatrix file.``

The matrix restart reader could not open a file in the selected binary or text
format.  The ``Lmatrix`` wording is also used when reading other matrices,
including ``Kmatrix2``.  Check that the complete indexed matrix-file set is
present and that ``IO.MatrixFile.BinaryFormat.Grab`` matches its format.  Do
not assume the MPI process count must match the previous run: this reader uses
the process count stored in ``InfoGlobal``.  See :ref:`gs_scf_restart` and
:ref:`md_restart`.

``grab_blip_coeffs: failed to open input file <filename>``
   With a blip basis and ``Basis.LoadCoeffs T``, a per-process
   ``blip_coeffs`` file could not be opened.  Check the filename printed by
   CONQUEST, the coefficient files from the same calculation and their read
   permissions; see :ref:`basis_readcoeffs` and :ref:`md_restart`.

External runtime failures
-------------------------

MPI launch failures, scheduler cancellation, out-of-memory termination,
filesystem errors and Unix signals are not CONQUEST error codes, and their
meaning depends on the platform.  Check the scheduler and system logs as well
as the CONQUEST output.  For a segmentation fault seen only with multiple
OpenMP threads, also check the thread-stack guidance in :ref:`install_compile`.

Reporting a reproducible problem
--------------------------------

If the checks above do not identify the problem, follow :ref:`faq_bug` and
include the complete diagnostic and enough material to reproduce the run:

* ``Conquest_input``, the coordinate file and all applicable ion or
  pseudopotential files;
* any restart files needed to reproduce a restart failure;
* ``Conquest_out``, ``Conquest_warnings`` and scheduler output;
* the CONQUEST release, commit or other precise version identifier;
* compiler, numerical-library and MPI implementation versions; and
* MPI process and OpenMP thread counts.

Remove confidential data if necessary, but do not omit the lines immediately
before the failure or silently replace the failing inputs with different
ones.
