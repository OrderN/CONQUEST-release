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

Spin polarisaion
----------------
Spin-polarised calculation is performed by setting ``Spin.SpinPolarised`` to T. 

Users need to specify the number of spin-up/down electrons in unit cells by ``Spin.NeUP`` and ``Spin.NeDN`` or the difference of spin-up/down electrons by ``Spin.Magn``. When users want to specify the number of spin-up/down electrons for each species, use ``Atom.SpinNeUp`` and ``Atom.SpinNeDn``.

The number of electrons for each spin will be fixed during SCF calculations when we set ``Spin.FixSpin`` to T (default is F).

::

   # example of anti-ferro bcc Fe
   Spin.SpinPolarised T
   Spin.FixSpin  F
   Spin.NeUP  8.0
   Spin.NeDN  8.0
   
   %block ChemicalSpeciesLabel
   1   55.845   Fe1
   2   55.845   Fe2
   %endblock ChemicalSpeciesLabel
   
   %block Fe1
   Atom.SpinNeUp 5.00
   Atom.SpinNeDn 3.00
   %endblock Fe1
   %block Fe2
   Atom.SpinNeUp 3.00
   Atom.SpinNeDn 5.00
   %endblock Fe2

When we use multi-site or on-site support functions in spin-polarised calculations, we can choose whether to use different support function coefficients for each spin by ``Basis.SpinDependentSF`` (T/F, default is T).


