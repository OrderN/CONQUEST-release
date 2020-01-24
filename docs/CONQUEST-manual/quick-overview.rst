.. _quick_over:

==============
Quick Overview
==============

.. _setting_up:

Setting up a calculation
------------------------

CONQUEST requires three types of file for a calculation:

  * A coordinate file
  * Ion files (pseudopotentials)
  * The input file (``Conquest_input``)

.. _qo_coords:

Coordinates
===========
CONQUEST works with orthorhombic unit cells (i.e. with angles between
lattice vectors at ninety degrees).  The coordinate file is laid out
simply: lattice vectors, number of atoms, atom coordinates (along with
species and movement flags).  Either fractional or Cartesian
coordinates can be read (the default is fractional; Cartesian
coordinates require a flag to be set in the 
input file).  CONQUEST also reads and writes PDB format
coordinate files for biomolecular simulations.  More information can
be found in :ref:`io_coords`.

.. _qo_ions:

Ion files
=========
The ion files contain the pseudopotentials and pseudo-atomic orbitals
for the elements, and follow a format similar to the ion files from Siesta
(CONQUEST can read Siesta ion files).  A set of default inputs to
generate ion files is available in the directory ``pseudo-and-pao``.
These contain pseudopotentials based on the `PseudoDojo`_ library, and
allow ion files to be produced with the basis set generation code that
is included with CONQUEST in the ``tools/BasisGeneration`` directory.
Full details are found :ref:`here <io_ion>`.

.. _PseudoDojo: https://www.pseudo-dojo.org/

.. _qo_cqinput:

Conquest\_input
===============
The ``Conquest_input`` file contains all of the input flags to
control a CONQUEST run.  At a minimum, the file must specify: the run
type (e.g. ``static`` or ``md``); the coordinate file name; and the number
of species and the ion file names.  For a well characterised
calculation, further options must be given (for instance setting
details for the calculation of the density matrix).  Simple examples
are given in :ref:`examples` and full documentation of all options can
be found in :ref:`input_tags`. 

Go to :ref:`top <quick_over>`

.. _qo_output:

Output from a calculation
-------------------------
The main output from CONQUEST is in a single file, named
``Conquest_out`` by default (this can be changed, and output can be
written to ``stdout`` rather than a file).  This file
contains details of the calculation, energies, forces and stresses and
the various electronic structure and atomic movement calculations
performed.  The most important files that are produced during a run are:

  * ``Conquest_out`` The output file
  * ``Conquest_warnings`` A list of any warnings issued by the code
    (also in ``Conquest_out``)
  * ``coord_next.dat`` The updated set of atomic positions
  * ``conquest.bib`` References suggested for the calculation performed
  * ``input.log`` A log of input options (both set by user and
    defaults)

Other files are produced by different run types, and are discussed
elsewhere.
