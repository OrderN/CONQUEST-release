========
Examples
========

All example calculations here use diagonalisation and PAO basis sets
(with a simple one-to-one mapping between PAOs and support functions).

Static calculation
------------------
We will perform a self-consistent electronic structure calculation on
bulk silicon.  The coordinate file that is needed is:

::
   
   10.26  0.00  0.00
    0.00 10.26  0.00
    0.00  0.00 10.26
   8
     0.000 0.000 0.000  1 T T T
     0.500 0.500 0.000  1 T T T
     0.500 0.000 0.500  1 T T T
     0.000 0.500 0.500  1 T T T
     0.250 0.250 0.250  1 T T T
     0.750 0.750 0.250  1 T T T
     0.250 0.750 0.750  1 T T T
     0.750 0.250 0.750  1 T T T

You should save this in an appropriate file (e.g. ``coords.dat``).
The PAO file can be found **somewhere - figure out database**.  The
``Conquest_input`` file requires only a few simple lines at base:

::

   AtomMove.TypeOfRun static
   IO.Coordinates coords.dat
   Grid.GridCutoff  50
   Diag.MPMesh   T
   Diag.MPMeshX  2
   Diag.MPMeshY  2
   Diag.MPMeshZ  2
   General.NumberOfSpecies  1
   %block ChemicalSpeciesLabel
    1 28.086 Si
   %endblock

**NOTE** I anticipate that we will update the CONQUEST input so that
valence charge, number of PAOs and radius are read from the ion file
by default, so I have not included these in the default input. DRB
2019/08/09. **NOTE**

The parameters above should be relatively self-explanatory; the grid
cutoff (in Hartrees) sets the integration grid spacing, and can be
compared to the *charge density* grid cutoff in a plane wave code
(typically four times larger than the plane wave cutoff).  The
Monkhorst-Pack k-point mesh (``Diag.MPMeshX/Y/Z``) is a standard
feature of solid state codes.

The most important parameters set the number of species and give
details of what the species are (``ChemicalSpeciesLabel``).  For each
species label (in this case ``Si``) there should be a corresponding
file with the extension ``.ion`` (again, in this case ``Si.ion``).
CONQUEST will read the necessary information from this file for
default operation, so no further parameters are required.

**Some discussion of the final output - particularly the energy - once
we have decided on the default basis set**

You might like to experiment with the grid cutoff (note that the
number of grid points is proportional to the square root of the energy
while the spacing is proportional to one over this, and
that the computational effort will scale with the *cube* of the number
of grid points) to see how the energy converges; as with all DFT
calculations, you should ensure that you test the convergence with
respect to all parameters.

Relaxation
----------

Atomic Positions
~~~~~~~~~~~~~~~~
I wonder about a molecule here - maybe methane.

Cell Parameters
~~~~~~~~~~~~~~~
We will optimise the lattice constant of the bulk silicon cell that we
studied for the static calculation.  Here we need to change the type
of run, and add one more line:

::

   AtomMove.TypeOfRun cg
   AtomMove.OptCell T

If you run this calculation, you should find a final lattice constant
of **GIVE VALUE IN BOHR** after **NN** iterations.
   
Simple Molecular Dynamics
-------------------------
Maybe we can do some simple NVE or NVT of water.

Tutorials
---------

We recommend that you work through, in order, the tutorials included
in the distribution in the ``tutorials/`` directory
to become familiar with the modes of operation of the code.

Where next?
-----------

While the tutorials have covered the basic operations of Conquest,
there are many more subtle questions and issues.  Detailed
descriptions of the files required and produced by Conquest are given
in :ref:`input-output`, while the input tags themselves are documented
in :ref:`input_tags`.  A set of :ref:`important` give more details
about the operations of the code and the background theory.
