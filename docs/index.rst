.. CONQUEST

CONQUEST: a local orbital, large-scale DFT code
===============================================

CONQUEST is a local orbital density functional theory (DFT) code,
capable of massively parallel operation with excellent scaling.  It
uses a local orbital basis
to represent the Kohn-Sham eigenstates or the density matrix.
CONQUEST can be applied to atoms, molecules, liquids and solids, but
is particularly efficient for large systems.  The code can find the ground
state using exact diagonalisation of the Hamiltonian or via a linear
scaling approach.  The code has demonstrated scaling to over 2,000,000
atoms and 200,000 cores when using linear scaling, and over 3,400
atoms and 850 cores with exact diagonalisation.
CONQUEST can perform structural relaxation (including unit cell
optimisation) and molecular dynamics (in NVE, NVT and NPT ensembles
with a variety of thermostats).

Getting Started
---------------

* :doc:`why-conquest`
* :doc:`faq`
* :doc:`quick-overview`
* :doc:`installing`
* :doc:`examples`
  
.. toctree::
   :maxdepth: 1
   :hidden:      
   :caption: Getting Started

   why-conquest
   faq
   quick-overview
   installing
   examples

User Guide
----------

* :doc:`input-output`
* :doc:`groundstate`
* :doc:`convergence`
* :doc:`basissets`
* :doc:`strucrelax`
* :doc:`moldyn`
* :doc:`ext-tools`
* :doc:`errors`
* :doc:`input_tags`
  
.. toctree::
   :maxdepth: 1
   :hidden:      
   :caption: User Guide

   input-output
   groundstate
   convergence
   basissets
   strucrelax
   moldyn
   ext-tools
   errors
   input_tags
   generating_paos

Tutorials
---------

* :doc:`tutorials/intro`
* :doc:`tutorials/relax`
* :doc:`tutorials/MD`
* :doc:`tutorials/Basis`
* :doc:`tutorials/Advanced`

.. toctree::
   :maxdepth: 1
   :hidden:      
   :caption: Tutorials

   tutorials/intro
   tutorials/relax
   tutorials/MD
   tutorials/Basis
   tutorials/Advanced

Theory
------

* :doc:`theory-strucrelax`
* :doc:`theory-md`

.. toctree::
   :maxdepth: 1
   :hidden:      
   :caption: Theory

   theory-strucrelax
   theory-md
   
Get in touch
------------

- If you have suggestions for developing the code, please
  use `GitHub issues`_.  The developers cannot guarantee to offer
  support, though we will try to help. 
- Report bugs, or suggest features on `GitHub
  issues`_.  View the source code on the main `GitHub page`_.
- You can ask for help and discuss any problems that you may have on
  the Conquest mailing list (to register for this list, please send an
  email to mlsystem@ml.nims.go.jp with the subject "sub
  conquest-user", though please note that, for a little while, some of the
  system emails may be in Japanese; you will receive a confirmation
  email to which you should simply reply without adding any text).

.. _GitHub issues: https://github.com/OrderN/CONQUEST-release/issues
.. _GitHub page: https://github.com/OrderN/CONQUEST-release

Licence
-------

CONQUEST is available freely under the open source `MIT Licence`__.
We ask that you acknowledge use of the code by citing appropriate
papers, which will be given in the output file (a BiBTeX file
containing these references is also output).  The key CONQUEST
references are:

* D. R. Bowler, T. Miyazaki and M. J. Gillan, J. Phys. Condens. Matter
  **14**, 2781--2798 (2002)
* T. Miyazaki, D. R. Bowler, R. Choudhury
  and M. J. Gillan, J. Chem. Phys. **121**, 6186--6194 (2004)
* D. R. Bowler, R. Choudhury, M. J. Gillan and T. Miyazaki,
  phys. stat. sol. b **243**, 989--1000 (2006)

__ https://choosealicense.com/licenses/mit/
