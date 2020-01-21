.. _examples:

====================
Example calculations
====================

All example calculations here use diagonalisation and PAO basis sets
(with a simple one-to-one mapping between PAOs and support functions).

.. _ex_static:

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
The inputs for the ion file can be found in ``pseudo-and-pao/PBE/Si``
(for the PBE functional).  Changing to that directory and running the
``MakeIonFiles`` utility (in ``tools``) will generate the file
``SiCQ.ion``, which should be copied to the run directory, and renamed
to ``Si.ion``. The ``Conquest_input`` file requires only a few simple
lines at base: 

::

   AtomMove.TypeOfRun static
   IO.Coordinates coords.dat
   Grid.GridCutoff  50
   Diag.MPMesh   T
   Diag.GammaCentred T
   Diag.MPMeshX  2
   Diag.MPMeshY  2
   Diag.MPMeshZ  2
   General.NumberOfSpecies  1
   %block ChemicalSpeciesLabel
    1 28.086 Si
   %endblock

The parameters above should be relatively self-explanatory; the grid
cutoff (in Hartrees) sets the integration grid spacing, and can be
compared to the *charge density* grid cutoff in a plane wave code
(typically four times larger than the plane wave cutoff).  The
Monkhorst-Pack k-point mesh (``Diag.MPMeshX/Y/Z``) is a standard
feature of solid state codes; note that the grid can be forced to be
centred on the Gamma point.

The most important parameters set the number of species and give
details of what the species are (``ChemicalSpeciesLabel``).  For each
species label (in this case ``Si``) there should be a corresponding
file with the extension ``.ion`` (again, in this case ``Si.ion``).
CONQUEST will read the necessary information from this file for
default operation, so no further parameters are required.

The output file starts with a summary of the calculation requested,
including parameters set, and gives details of papers that are
relevant to the particular calculation.  After brief details of the
self-consistency, the total energy, forces and stresses are printed,
followed by an estimate of the memory and time required.

You might like to experiment with the grid cutoff (note that the
number of grid points is proportional to the square root of the energy
while the spacing is proportional to one over this, and
that the computational effort will scale with the *cube* of the number
of grid points) to see how the energy converges; as with all DFT
calculations, you should ensure that you test the convergence with
respect to all parameters.

.. _ex_relax:

Relaxation
----------

.. _ex_relax_atoms:

Atomic Positions
~~~~~~~~~~~~~~~~
I wonder about a molecule here - maybe methane.

.. _ex_relax_cell:

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

.. _ex_md:

Simple Molecular Dynamics
-------------------------
Maybe we can do some simple NVE or NVT of water.

.. _ex_tut:

Tutorials
---------

We recommend that you work through, in order, the tutorials included
in the distribution in the ``tutorials/`` directory
to become familiar with the modes of operation of the code.

**NOTE** In the initial pre-release of CONQUEST (January 2020) we have
not included the tutorials; they will be added over the coming months.

.. _ex_next:

Where next?
-----------

While the tutorials have covered the basic operations of Conquest,
there are many more subtle questions and issues, which are given in
the User Guide.
