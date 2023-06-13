.. _post-proc:

===============================
Post-processing CONQUEST output
===============================

.. _pp_intro:

Introduction
++++++++++++

The utility ``PostProcessCQ`` allows users to post-process the output
of a CONQUEST calculation, to produce structure files, densities of
states, and charge density, band
densities and STM images as CUBE
files (which can be read by the freely available `VESTA
<https://jp-minerals.org/vesta/en/>`_ code).

There are a number of different analyses which can be performed:
coordinate conversion (to formats which can be plotted); conversion of
total charge density to CUBE file format; production of band-resolved
(optionally k-point resolved) densities in CUBE file format; simple
Tersoff-Hamann STM simulation; and calculation of densities of states,
including projected DOS.  You should ensure that all the files
produced during the CONQUEST run are available for the post-processing
(including ``eigenvalues.dat``, ``chden.NNN``, ``make_blk.dat`` or
``hilbert_make_blk.dat`` and ``ProcessNNNNNNNWF.dat`` and
``ProcessSijNNNNNNNWF.dat`` as applicable) as well as the input files.

**Note** that the utility reads the ``Conquest_input`` file, taking some
flags from the CONQUEST run that generated the output, and some
utility-specific flags that are detailed below.

**Note also** that projected DOS, band density and STM simulation are
not at present compatible with multi-site support functions (MSSF),
though we hope to implement this soon.

Go to :ref:`top <post-proc>`.

Coordinate conversion
+++++++++++++++++++++

Set ``Process.Job coo`` to output a coordinate file for further
processing or plotting.  The utility will read the file specified by
``Process.Coordinates`` (which defaults to the file specified by
``IO.Coordinates``).  The output format is selected by specifying the
``Process.CoordFormat`` tag.  The default output format is XYZ (which
adds a ``.xyz`` suffix to the file name) using ``xyz``.  The CASTEP
``.cell`` output format can also be selected using ``cell``.  We plan to
expand this conversion to other formats in the future.

**Note** that for a structural relaxation or molecular dynamics
calculation, if you do not specify ``Process.Coordinates`` then the
``IO.Coordinates`` file, which will be converted, will be the *input*
structure, not the output structure.  Parameters that can be set are:

::

   Process.Coordinates string (default: IO.Coordinates value)
   Process.CoordFormat string (default: xyz; options: xyz, cell)

Go to :ref:`top <post-proc>`.

Charge density
++++++++++++++

Setting ``Process.Job`` to ``cha``, ``chg`` or ``den`` will convert
the files ``chden.NNN`` which are written by CONQUEST to a cube file.
The processing will use the files ``chden.NNN``, ``Conquest_input``
and ``hilbert_make_blk.dat`` or ``raster_make_blk.dat``.  Parameters
that can be set include:

::
   
   Process.ChargeStub string (default: chden)

The ChargeStub simply defines the filename which will be read, and
used for output.

**Note** that to output the ``chden.NNN`` files from CONQUEST, you must
set the flag ``IO.DumpChargeDensity T`` in the CONQUEST run.

Go to :ref:`top <post-proc>`.

Band density
++++++++++++

Setting ``Process.Job`` to ``ban`` produces band densities from wave
function coefficients output by CONQUEST.  The CONQUEST run must have
the following tags set:

::

   IO.outputWF T

A set of bands whose coefficients are output are specified either with
an energy range (the default is to produce *all* bands):

::

   IO.WFRangeRelative T/F
   IO.min_wf_E real (Ha)
   IO.max_wf_E real (Ha)

or with a list of bands:

::

   IO.maxnoWF n

   %block WaveFunctionsOut
   n entries, each a band number
   %endblock

The wavefunction range can be relative to the Fermi level
(``IO.WFRangeRelative T``) otherwise it is absolute.  Either of these
will produce a file containing all eigenvalues at all k-points
(``eigenvalues.dat``) and a series of files containing the
wavefunction expansion coefficients for the selected bands
(``ProcessNNNNNNNWF.dat``).  These files are output as binary
(unformatted) by default (this can be changed by setting
``IO.MatrixFile.BinaryFormat F`` before the CONQUEST run) and will be
read using the same format (it is important to check this!).

From these wavefunction coefficient files, band densities can be
produced in post-processing, using similar tags; either a range:

::

   Process.min_wf_E real (Ha)
   Process.max_wf_E real (Ha)
   Process.WFRangeRelative T/F

or an explicit list of bands:

::

   Process.noWF n

   %block WaveFunctionsProcess
   n entries, each a band number
   %endblock

Note that the bands to be processed must be a subset of the bands
output by CONQUEST.  The bands can be output summed over k-points, or
at individual k-points, by setting ``Process.outputWF_by_kpoint`` to
``F`` or ``T`` respectively.

Go to :ref:`top <post-proc>`.

Tersoff-Hamann STM simulation
+++++++++++++++++++++++++++++

Setting ``Process.Job ter`` will use a very simple Tersoff-Hamann
approach to STM simulation, summing over band densities between the
Fermi level and the bias voltage (this is often surprisingly
accurate).  The following parameters can be set:

::

   STM.BiasVoltage    real (eV)
   STM.FermiOffset    real (eV)
   Process.MinZ       real (Bohr)
   Process.MaxZ       real (Bohr)
   Process.RootFile   string (default: STM)

The ``FermiOffset`` tag allows the user to shift the Fermi level (to simulate
charging or an external field).  The height of the simulation cell
in which the STM image is calculated is set by the ``MinZ`` and
``MaxZ`` tags, and the filename by the ``RootFile`` tag.

Go to :ref:`top <post-proc>`.

Density of states (DOS)
+++++++++++++++++++++++

Setting ``Process.Job dos`` will produce a total density of states
(DOS) for the system, using the eigenvalues output by CONQUEST.  The
following parameters can be set:

::

   Process.min_DOS_E real    (Ha, default lowest eigenvalue)
   Process.max_DOS_E real    (Ha, default highest eigenvalue)
   Process.sigma_DOS real    (Ha, default 0.001)
   Process.n_DOS     integer (default 1001)

The limits for the DOS are set by the first two parameters (note that
CONQUEST will output all eigenvalues, so the limits on these are set
by the eigenspectrum).  The broadening applied to each state is set by
``sigma_DOS``, while the number of bins is set by ``n_DOS``.  The
integrated DOS is also calculated; the user can choose whether this
is the total integrated DOS (i.e. from the lowest eigenvalue,
regardless of the lower limit for DOS) or just the local integrated
DOS (i.e. over the interval specified for the DOS) by setting
``Process.TotalIntegratedDOS`` to ``T`` or ``F``, respectively.

We recommend that, for accurate DOS, CONQUEST should be run
non-self-consistently with a very high k-point density, after reading
in a well-converged input charge density: set ``minE.SelfConsistent
F`` and ``General.LoadRho T`` (which will require that the converged
charge density is written out by CONQUEST by setting ``IO.DumpChargeDensity T``).

Go to :ref:`top <post-proc>`.

Atom-projected DOS
++++++++++++++++++

Setting ``Process.Job pdos`` will produce a total density of states as
above, as well as the density of states projected onto the individual
atoms.  Given support functions :math:`\phi_{i\alpha}(\mathbf{r})`
which are the basis functions of the Kohn-Sham eigenstates
:math:`\psi_{n}(\mathbf{r}) = \sum_{i\alpha}
c^{n}_{i\alpha}\phi_{i\alpha}(\mathbf{r})`, then the projection of a
given state, :math:`n`, onto an atom :math:`i` can be written as
:math:`\sum_{\alpha j\beta} c^{n}_{i\alpha}
S_{i\alpha,j\beta}c^{n\mathbf{k}}_{j\beta}`.  The projected DOS is
constructed using these projections.

If using :ref:`pseudo-atomic orbitals (PAOs) <basis_paos>` as the
basis set, then the atom-projected DOS can be further resolved by
angular momentum (either just :math:`l` or both :math:`l` and
:math:`m`).  If using :ref:`pseudo-atomic orbitals (PAOs)
<basis_paos>` with :ref:`multi-site support functions <basis_mssf>` or
:ref:`blip functions <basis_blips>` then it is not possible to
decompose the DOS any further (in future, it may be possible to
resolve the MSSF coefficients into the individual PAOs, and hence
decompose pDOS by angular momentum).  To output the necessary
coefficients to produce atom-projected DOS, a CONQUEST run must be
performed with the following parameters set:

::

   IO.writeDOS T
   IO.write_proj_DOS T

As for the DOS, very high Brillouin zone sampling is required for
accurate projected DOS, which is most efficiently generated using a
converged charge density and a non-self-consistent calculation with
much higher k-point density.  CONQUEST will produce the wavefunction
files (``ProcessNNNNNNNWF.dat`` and ``ProcessSijNNNNNNNWF.dat``) as
binary (unformatted) by default (change using the flag
``IO.MatrixFile.BinaryFormat F``).

Once the files have been generated by CONQUEST, the output can be
processed by setting the output tag:

::
   
   Process.Job pdos

This is all that is needed for the simplest output.  The number of
bins and smearing of the peaks can be set using:

::

   Process.sigma_DOS 0.002
   Process.n_DOS 10001
   
To resolve the DOS by angular momentum as well as by atom, then the
following flags can be set:

::

   Process.pDOS_l_resolved T
   Process.pDOS_lm_resolved T

Note that only one of these is needed, depending on what level of
resolution is required.  At present, angular momentum resolution is
only available for the PAO basis set (not MSSF or blips) though it
is under development for the MSSF basis (by projection onto the
underlying PAO basis).

The energy range for the projected DOS can
also be specified:

::
   
   Process.min_DOS_E -0.35
   Process.max_DOS_E  0.35
   Process.WFRangeRelative T

where the final tag sets the minimum and maximum values relative to
the Fermi level.

Go to :ref:`top <post-proc>`.

Band structure
++++++++++++++

The band structure of a material can be generated by CONQUEST by performing
a non-self-consistent calculation after reading a well-converged charge density:
set ``minE.SelfConsistent F`` and ``General.LoadRho T`` (which will require that the converged
charge density is written out by CONQUEST by setting ``IO.DumpChargeDensity T``).
The k-points required can be specified as lines of points in k-space;
setting ``Diag.KspaceLines T`` enables the approach, while the number of lines
(e.g. Gamma to L; L to X would be two lines) is set with ``Diag.NumKptLines``
and the number of points along a line with ``Diag.NumKpts``.  The k-point lines
themselves are set with a block labelled ``Diag.KpointLines`` which should have
two lines (starting and finishing k-points) for each k-point line.  So to create
a bandstructure from X-Gamma-L-X (3 lines) with 11 points in each line, you would
use the following input:

::
   Diag.KspaceLines T
   Diag.NumKptLines 3
   Diag.NumKpts 11
   %block Diag.KpointLines
   0.5 0.0 0.0
   0.0 0.0 0.0
   0.0 0.0 0.0
   0.5 0.5 0.5
   0.5 0.5 0.5
   0.5 0.0 0.0
   %endblock

Setting ``Process.Job bst`` and running the post-processing will read the resulting
``eigenvalues.dat`` file, and produce a ``BandStructure.dat`` file.  The x-axis will
be the k-point index by default, but specifying ``Process.BandStrucAxis`` (taking
value ``n`` for number, ``x/y/z`` for a single direction in k-space or ``a`` to
give all k-point coordinates) will allow you to control this.

Go to :ref:`top <post-proc>`.

