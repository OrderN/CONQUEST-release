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
Support functions are used to express density matrices. The functions are constructed from basis functions, pseudo-atomic orbital (PAO) functions which are given in ``.ion`` files (default) or blip functions. 

When PAOs are used as support functions without any modification, we don't need to set any parameters related to support functions in ``Conquest_input``. 

If we contract multiple PAOs to small number of support functions or use blip functions, we need to provide several input parameters.

Full details of how the support functions are found and represented
can be found in :ref:`basissets`.
