.. _generating_paos:

===============
Generating PAOs
===============

.. _pao_gen_intro:

Introduction
------------

CONQUEST includes a utility for generating the PAO basis files (``MakeIonFiles`` with source code in the directory ``tools/BasisGeneration``), and we also provide pseudopotential files (from the `PseudoDojo <http://www.pseudo-dojo.org/>`_ database).  The input files will generate a reasonable default basis set (though we can offer no guarantees: users must test the accuracy and convergence of their basis sets).  Here we will discuss how to generate both default and custom basis sets.  Full details of the basis sets can be found in a recent paper :cite:`gp-Bowler:2019fv`.

.. _pao_gen_default:

Default basis sets
------------------

To generate basis functions, radii for the PAOs must be specified; by default, the utility will set these radii automatically.  The radii can be set so that the different shells either share the same radii, or share an energy shift associated with confinement.  The default behaviour is to generate basis sets where shells share radii (to change this, add the line ``Atom.Cutoffs energy`` to the input file).  There are four default basis set sizes:

* minimal (single zeta, SZ)
* small   (single zeta and polarisation, SZP)
* medium  (double zeta, single polarisation, DZP)
* large   (triple zeta, triple polarisation, TZTP)

Generally, reasonable results will be obtained with a medium (DZP) basis, though this should always be tested.  Minimal and small basis sets are much faster (and are the only basis sets compatible with :ref:`linear scaling <gs_on>`), but less reliable.  The large basis set will be slower (often it is twice the size of the medium basis set, so diagonalisation will be up to eight times slower) but more reliable and accurate.

We note that Group I and II atoms are a little problematic: the standard approach for most other elements may produce a somewhat limited basis set, so we have created a more accurate, customised input file for these elements (with the exception of Na and Mg, where the pseudopotential does not include l=2 components, so the default approach is all that is possible).  These should be tested carefully.
  
Go to :ref:`top <generating_paos>`.

.. _pao_gen_specify:

Specifying basis sets
---------------------

The generation utility gives the user complete control over the basis sets that are produced.  As an example, we reproduce below the input file for strontium (Sr) and discuss the layout.

::

   %block Sr
   Atom.PseudopotentialFile Sr.in
   Atom.VKBFile Sr.pot
   Atom.Perturbative_Polarised F
   Atom.PAO_N_Shells 5
   Atom.BasisBlock SrBlock
   %endblock
   
   %block SrBlock
   # n, l, number of zetas
   4 0 1
   4 1 1
   5 0 2
   5 1 1
   4 2 1
   # Radii for PAOs (bohr)
   4.0
   5.0
   10.1 5.7
   10.1
   10.1
   %endblock

In this case, we specify five shells (combinations of n and l) via the ``Atom.PAO_N_Shells`` tag, with polarisation functions found simply by solving the Schrodinger equation in the usual way (the alternative, perturbative polarisation, is discussed below).  We must then specify a block that defines these shells (``Atom.BasisBlock SrBlock``).  Within that block, we give the number of zeta functions for each (n,l) pair (specified as a line ``n l nzeta``) followed by the radii for the zeta functions.

Setting radii for the different shells is a complex process which requires considerable time and care, with an extensive literature; we cannot provide significant help, but only make suggestions.  In the first instance, the default radii are a good starting point.  Note that the default setting ``Atom.Cutoffs radii`` averages the radii between shells, while ``Atom.Cutoffs energy`` finds different radii for each shell (so that the energy change due to confinement is the same for all shells).  We recommend starting from one of these sets of radii, and then testing and optimising the radii against some key properties of the system.

A common approach to the generation of polarisation functions (i.e. unoccupied states) is to perturb a valence state (typically the highest energy valence state) to generate a function with angular momentum increased by one; this is the default behaviour.  In this case, the radii for the polarisation state should be the same as the shell being polarised (so for Si, we would perturb the 3p (n=3, l=1) state to get the 3d (n=3, l=2) state), and at present the same number of zeta functions must be specified for the polarisation shell as for the unperturbed shell.  For instance, for Si:

::

   %block Si
   Atom.PseudopotentialFile Si.in
   Atom.VKBFile Si.pot
   Atom.PAO_N_Shells 3
   Atom.BasisBlock SiBlock
   %endblock
   
   %block SiBlock
   # n, l, number of zetas
   3 0 2
   3 1 2
   3 2 2
   # Radii for PAOs (bohr)
   8.0 4.0
   8.0 4.0
   8.0 4.0
   %endblock

The perturbative option can be turned off by specifying `Atom.Perturbative_Polarised F` in the input file.  (Note that in the strontium example above we have specified two polarisation shells, so cannot use the perturbative approach.)

By default, the utility calculates radii which are shared between shells; it is possible to specify instead shared energy shifts using ``Atom.Cutoffs energy``, but this can only be done for valence shells, and so *must* use the perturbative polarisation approach for polarisation functions.


Go to :ref:`top <generating_paos>`.

.. _pao_gen_comp:

Compiling
---------

To compile the code, the same ``system.make`` can be used as is specified for the main code.  Once this is done, simply issue the ccommand ``make`` in the ``tools/BasisGeneration`` directory.  The resulting executable will be placed in the ``bin`` directory.

Go to :ref:`top <generating_paos>`.

.. _pao_gen_pseudo:

Generating new pseudopotentials
-------------------------------

CONQUEST is supplied with a complete set of pseudopotentials for the elements in the `PseudoDojo <http://www.pseudo-dojo.org/>`_ database (covering LDA, PBE and PBEsol exchange-correlation functionals).  In order to generate new pseudopotential files, users will need the `Hamann <http://www.mat-simresearch.com/>`_ pseudopotential code ONCVPSP v3.3.1 (the current release) and the patch file ``Conquest_ONCVPSP_output.patch`` which is in the ``tools`` directory.  After patching and compiling the Hamann code (to patch the code, copy the patch to the ONCVPSP ``src`` directory, and issue the command ``patch -p0 < Conquest_ONCVPSP_output.patch``; we cannot provide any support for this) the ``oncvpsp.x`` utility will generate a file ``VPS.dat`` which should be renamed (something like ``element.pot`` as in the CONQUEST pseudopotential files) and specified in the input file using the ``Atom.VKBFile`` tag.

Go to :ref:`top <generating_paos>`.

.. bibliography:: references.bib
    :cited:
    :labelprefix: GP
    :keyprefix: gp-
    :style: unsrt

Go to :ref:`top <generating_paos>`.
