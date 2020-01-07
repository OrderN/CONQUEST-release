.. _groundstate:

========================
Finding the ground state
========================

Finding the electronic ground state is the heart of any DFT code.  In
CONQUEST, we need to consider: the density matrix (found using
diagonalisation or linear scaling); self-consistency; and the support
functions (which can be optimised or not).

.. _gs_diag:

Diagonalisation
---------------

.. _gs_on:

Linear Scaling
--------------

.. _gs_scf:

Self-consistency
----------------

.. _gs_suppfunc:

Support functions
-----------------
Full details of how the support functions are found and represented
can be found in **add this link**.

.. _gs_spin:

Spin polarisation
-----------------
CONQUEST performs collinear spin calculations only.  A spin-polarised
calculation is performed by setting the parameter
``Spin.SpinPolarised`` to T.  

Users need to specify *either* the total initial number of spin-up and spin-down electrons in
the simulation cell (using the parameters ``Spin.NeUP`` and
``Spin.NeDN``), *or* the difference between the number of spin-up and
spin-down electrons (using the parameter ``Spin.Magn``).

The number of electrons for each spin channel can be fixed during SCF
calculations by setting the parameter ``Spin.FixSpin`` to T (default is F).

It is possible to specify the spin occupation in the atomic charge
densities (i.e. the number of spin-up and spin-down electrons used to
build the density).  This is done in the :ref:`input_tags_atomic_spec`
part of the ``Conquest_input`` file.  Within the atom block for
each species, the numbers of electrons should be set with
``Atom.SpinNeUp`` and ``Atom.SpinNeDn``.  Note that these numbers
*must* sum to the number of valence electrons for the atom.

Go to :ref:`top <groundstate>`.

.. _gs_spin_example:

Examples: FM and AFM iron
~~~~~~~~~~~~~~~~~~~~~~~~~

A two atom ferromagnetic iron simulation might be set up using the
parameters below.  Note that the net spin here is S=1 :math:`\mu_B`
(i.e. two more electrons in the up channel than in the down), and
that the net spin is not constrained.

:: 

   # example of ferro bcc Fe
   Spin.SpinPolarised T
   Spin.FixSpin  F
   Spin.NeUP  9.0     # initial numbers of up- and down-spin electrons,
   Spin.NeDN  7.0     # which will be optimised by a SCF calculation when Spin.FixSpin=F
   
   %block ChemicalSpeciesLabel
   1   55.845   Fe
   %endblock ChemicalSpeciesLabel

An equivalent anti-ferromagnetic calculation could be set up as
follows (though note that the initial specification of spin for the
atoms does *not* guarantee convergence to an AFM ground state).  By
defining two species we can create spin-up and spin-down atoms (note
that both species will require their own, appropriately labelled, ion
file). 

::

   # example of anti-ferro bcc Fe
   Spin.SpinPolarised T
   Spin.FixSpin  F
   Spin.NeUP  8.0     # initial numbers of up- and down-spin electrons in an unit cell
   Spin.NeDN  8.0     # are set to be the same
   
   %block ChemicalSpeciesLabel
   1   55.845   Fe1
   2   55.845   Fe2
   %endblock ChemicalSpeciesLabel
   
   %block Fe1           # up-spin Fe
   Atom.SpinNeUp 5.00
   Atom.SpinNeDn 3.00
   %endblock Fe1
   %block Fe2           # down-spin Fe
   Atom.SpinNeUp 3.00
   Atom.SpinNeDn 5.00
   %endblock Fe2

When using multi-site or on-site support functions in spin-polarised
calculations, the support functions can be made spin-dependent
(different coefficients for each spin channel) or not by setting
``Basis.SpinDependentSF`` (T/F, default is T).

Go to :ref:`top <groundstate>`.
