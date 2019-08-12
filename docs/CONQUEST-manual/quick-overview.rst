==============
Quick Overview
==============

Setting up a calculation
------------------------

  * Create coordinate file
  * Find ion files
  * Set up the input file (``Conquest_input``)

Coordinates
===========
CONQUEST works with orthorhombic unit cells (i.e. with angles between
lattice vectors at ninety degrees).  The coordinate file is laid out
simply: lattice vectors, number of atoms, atom coordinates (along with
species and movement flags).  Either fractional or Cartesian
coordinates can be read.  CONQUEST also reads and writes PDB format
coordinate files for biomolecular simulations.  More information can
be found in :ref:`coordinate-file`.

Ion files
=========
The ion files follow a format similar to the ion files from Siesta
(and CONQUEST can read Siesta ion files).  A set of default ion files
for different basis set sizes is available **FROM SOMEWHERE**.  A
basis set generation code is included with CONQUEST.  Full details are
found **WHERE**.


Conquest\_input
===============
The ``Conquest_input`` file contains all of the input flags to
control a CONQUEST run.  Simple examples are given in :ref:`examples`
and full documentation of all options can be found in :ref:`input_tags`.

Output from a calculation
-------------------------
The main output from CONQUEST is in a single file, named
``Conquest_out`` by default (this can be changed with the input tag
``IO.OutputFile``, and output to ``stdout`` rather than a file can be
toggled with the ``IO.WriteOutToFile`` boolean input tag).  This file
contains details of the calculation, energies, forces and stresses and
the various electronic structure and atomic movement calculations
performed.  If atoms are moved, new coordinates are output in
``coord_next.dat`` [ **CHECK THIS** ].  Other files that are produced
include:

  * ``input.log`` A log of input options (both set by user and defaults)
  * ``chden.nnn`` The charge density on the grid from process ``nnn``
