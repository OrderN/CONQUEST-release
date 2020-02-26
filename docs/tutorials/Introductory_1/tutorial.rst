.. Bulk silicon input

Bulk silicon: input, output and SCF
===================================

This tutorial is intended to be a very basic introduction to the input
required for CONQUEST, the output generated, and the self-consistency
(SCF) procedure; it uses the same system as the first of the examples
in the manual, but provides more detail.  The files are found in
``docs/tutorials/Introductory_1``. 

CONQUEST requires the following files to run:

* The input file: ``Conquest_input``
* A coordinates file (name set in ``Conquest_input``; no default)
* Ion files (suffix ``.ion``), which provide the pseudopotentials and
  pseudo-atomic orbitals (PAOs)

The input file requires the user to provide a certain amount of
information.  The minimal file that is provided for this tutorial
gives most of these:

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

The key entries are:

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
discussed in the next tutorial

* The coordinate file ``IO.Coordinates``
* The number of species ``General.NumberOfSpecies``
* The ion files for the species

* The basic input file
* The output
* Changing the output level and destination
* Controlling the SCF (tolerance and iterations, options)
