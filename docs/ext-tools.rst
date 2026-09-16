.. _ext-tools:

==============
External tools
==============

.. _et_post_process:

Post-processing for charge density, band density, DOS, STM
----------------------------------------------------------

The utility ``PostProcessCQ`` allows users to post-process the output
of a CONQUEST calculation, to produce the charge density, band
densities, DOS and STM images in useful forms.  It is described fully :ref:`here <post-proc>`.

Atomic Simulation Environment (ASE)
-----------------------------------

.. _et_ase_cq:

`ASE <https://wiki.fysik.dtu.dk/ase>`_ is a set of 
Python tools for setting up, manipulating, running, visualizing and analyzing 
atomistic simulations. ASE contains a CONQUEST interface, also 
called *Calculator* so that it can be used to calculate ``energies``, ``forces`` 
and ``stresses`` as inputs to other calculations such as `Phonon <https://wiki.fysik.dtu.dk/ase/ase/phonons.html#module-ase.phonons>`_ 
or `NEB <https://wiki.fysik.dtu.dk/ase/ase/neb.html#module-ase.neb>`_ that 
are not implemented in CONQUEST. ASE is a versatile tool to manage CONQUEST
calculations without pain either: 

* in a **direct** way where pre-processing, calculation and post-processing are managed on-the-fly by ASE, 
* or in an **indirect** way where the calculation step is performed outside the workflow, ie. on a supercomputer.

The ASE repository containing the Conquest calculator can be
found `here <https://gitlab.com/lionelalexandre/ase-Conquest/-/tree/master?ref_type=heads>`_.
Detailed documentation on how to manage Conquest calculations
with ASE is available :ref:`here <ase-conquest>`.

Go to :ref:`top <ext-tools>`.
