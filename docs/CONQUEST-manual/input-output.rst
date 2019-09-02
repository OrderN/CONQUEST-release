.. _input-output:

================
Input and output
================

Input files
-----------

Conquest_input
++++++++++++++
All input parameters are specified in the ``Conquest_input`` file,
including the names of the coordinate file and the ion files.  This
file controls the run; there are many sensible default values for
input parameters, but you should ensure that you understand what
they mean.

**Basic set of input tags?**

Full documentation can be found in :ref:`input_tags`.

Ion files
+++++++++

The ion files contain data on the different species being modelled:
valence charge, pseudopotentials, pseudo-atomic orbitals etc.  A
utility for generating these is provided with Conquest, but Siesta ion
files can also be read.  The Conquest utility uses the `ONCVPSP code
from Hamann <http://http://www.mat-simresearch.com>`_.

**More detail on PAO utility**.

**Details of default ion file database?**

.. _coordinate-file:
  
Coordinates
+++++++++++

The coordinates are specified in a separate file with relatively
simple format.  The coordinates can be specified in fractional form
(set the input tag ``IO.FractionalAtomicCoords T``) or cartesian.
Distance units can be Angstroms or Bohr radii (set the input tag
``General.DistanceUnits`` to ``Ang`` or ``bohr`` respectively).  At present,
Conquest only handles *orthorhombic* unit cells.

The coordinate file is formatted as follows:

::
   
   a   0.0 0.0
   0.0 b   0.0
   0.0 0.0 c
   NAtoms
   x y z species MoveX MoveY MoveZ
   .
   .
   .

Note that the flags ``MoveX`` etc take values T/F and indicate whether
atoms are free to move in x, y and z, respectively.  The flag
``species`` is an integer, and selects based on species defined in the
``Conquest_input`` file.

Output files
------------

Main output
+++++++++++

By default, CONQUEST writes output to the ``Conquest_out`` file
(though the flag ``IO.OutputFile`` allows the filename to be
specified, and the flag ``IO.WriteOutToFile`` selects output to file
or ``stdout``).  This file contains all details of the calculation,
including energies, forces and information on the different stages of
the calculation.  The output verbosity is controlled by the
``IO.Iprint`` family of parameters, which allows different levels in
different areas of the code.

Electronic structure
++++++++++++++++++++

Different electronic structure outputs are available; in each case,
the key output flag is given.  Further output flags are described in :ref:`input_tags`.

  * Charge density
  * Band-resolved charge density (``IO.outputWF``)
  * Density of states (``IO.writeDOS``)
  * Atom-projected density of states (``IO.write_proj_DOS``)
  * Atomic charges, using the Becke approach (``IO.AtomChargeOutput``)

The Kohn-Sham eigenvalues are output in the ``Conquest_out`` file.
The charge densities need post-processing to convert from the
standard output format to a file compatible with visualisation
(current supported formats include Gaussian CUBE file and OpenDX
files).

Atomic structure
++++++++++++++++

During structural relaxation and molecular dynamics, the atomic
structure is saved in the output file ``coord_next.dat``.  This is in
the same format as the input.

.. _io_md:

Molecular dynamics
++++++++++++++++++

A molecular dynamics run will generate a number of additional plain text output
files:

  * ``md.stats`` --- summarises thermodynamic quantities at each steps
  * ``md.frames`` --- contains the complete physical state of the system (lattice
    parameters, atomic positions, velocities, forces, stress).
  * ``md.checkpoint`` --- data required for MD restart, namely atomic velocities
    and extended system variables.
  * ``md.positions`` --- Atomic coordinates saved at the moment of checkpointing
  * ``trajectory.xsf`` --- atomic coordinates save in .xsf format, which can be
    visualised using (for example) VMD, if ``AtomMove.WriteXSF`` is true..

Full details are available in :ref:`moldyn`.
