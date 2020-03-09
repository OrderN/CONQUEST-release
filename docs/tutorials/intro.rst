.. _intro_tut:

Introductory Tutorials
======================

These introductory tutorials will give you an overview of how to run
CONQUEST, the files and parameter settings required, and what output
to expect.  

.. _intro_one:

Bulk silicon: input, output and SCF
-----------------------------------

We start with a very basic introduction to the input
required for CONQUEST, the output generated, and the self-consistency
(SCF) procedure; it uses the same system (bulk silicon) as the first of the examples
in the manual, but provides more detail.  The files are found in
``docs/tutorials/Introductory_1``. 

CONQUEST requires the following files to run:

* The input file: ``Conquest_input``
* A coordinates file (name set in ``Conquest_input`` with the flag ``IO.Coordinates``; no default)
* Ion files (suffix ``.ion``), which provide the pseudopotentials and
  pseudo-atomic orbitals (PAOs) (names set in ``Conquest_input`` in
  the ``ChemicalSpeciesLabel`` block)

The input file requires the user to provide a certain minimal amount
of information, with most parameters having reasonable defaults set.
The file that is provided for this tutorial gives the most
basic parameters:

::

   # Input/Output
   IO.Title Bulk Si 8 atoms static
   IO.Coordinates ionpos.dat
   
   # General Parameters
   General.NumberOfSpecies 1
   
   %block ChemicalSpeciesLabel
   1  28.0850   Si_SZ
   %endblock

   # Moving Atoms
   AtomMove.TypeOfRun static
   
   # Finding the density matrix
   DM.SolutionMethod diagon
   
   # k-points
   Diag.GammaCentred T
   Diag.MPMesh T
   Diag.MPMeshX 2
   Diag.MPMeshY 2
   Diag.MPMeshZ 2

The most important entries are:

* the coordinate file (``IO.Coordinates``);
* the number of species (``General.NumberOfSpecies``);
* the specification for the species (the block
  ``ChemicalSpeciesLabel`` gives the atomic mass and the ion file name
  for all species);
* the type of run (``AtomMove.TypeOfRun`` which defaults to ``static``)

The Brillouin zone sampling must be investigated carefully, as for
all periodic electronic structure calculations.  The Monkhorst-Pack
mesh (``Diag.MPMesh``) offers a convenient way to do this systematically.
The job title is purely for reference.  Further parameters are
discussed in the next tutorial, with full documentation found in the
:ref:`diagonalisation <input_diag>` section of the manual.

Go to :ref:`top <intro_tut>`.
     
Coordinates
~~~~~~~~~~~

As described in :ref:`the manual <io_coords>`, the coordinate file
requires an orthorhombic simulation cell, with coordinates specified
either in fractional form (the default) or Cartesian (if the tag
``IO.FractionalAtomicCoords F`` is set).  The units for distance are
Bohr radii by default, though can be changed to Angstroms (set
``General.DistanceUnits Ang``).

* The coordinate file ``IO.Coordinates``
* The number of species ``General.NumberOfSpecies``
* The ion files for the species

* The basic input file
* The output
* Changing the output level and destination
* Controlling the SCF (tolerance and iterations, options)

Go to :ref:`top <intro_tut>`.

.. _intro_two:

Bulk silicon: parameters to converge
------------------------------------

* The files that are needed

  * Coordinates
  * Ion files
  * Input file: ``Conquest_input``

* Integration grid
* Brillouin zone sampling
* Possibly basis set size

Go to :ref:`top <intro_tut>`.

.. _intro_three:

Bulk silicon: analysis
----------------------

* The files that are needed

  * Coordinates
  * Ion files
  * Input file: ``Conquest_input``

* Total DOS
* Atom-projected DOS
* Band structure output
* Charge density and bands
* Atomic charges

Go to :ref:`top <intro_tut>`.
