.. CONQUEST

CONQUEST: a local orbital, large-scale DFT code
===============================================

CONQUEST is a massively parallel density functional theory (DFT) code, using a local orbital basis
to represent the Kohn-Sham eigenstates or the density matrix.
CONQUEST can be applied to atoms, molecules, liquids and solids, but
is particularly efficient for large systems.  The code can find the ground
state using exact diagonalisation of the Hamiltonian or via a linear
scaling approach.  The code has demonstrated scaling to over 2,000,000
atoms and 200,000 cores when using linear scaling.
CONQUEST can perform structural relaxation (including unit cell
optimisation) and molecular dynamics (in NVE, NVT and NPT ensembles
with a variety of thermostats).

Documentation
-------------

**Getting Started**

* :doc:`why-conquest`
* :doc:`faq`
* :doc:`quick-overview`
* :doc:`examples`
* :doc:`installing`
  
.. toctree::
   :maxdepth: 1
   :hidden:      
   :caption: Getting Started

   why-conquest
   faq
   quick-overview
   examples
   installing

**User Guide**

* :doc:`input-output`
* :doc:`important`
* :doc:`errors`
* :doc:`input_tags`
  
.. toctree::
   :maxdepth: 1
   :hidden:      
   :caption: User Guide

   input-output
   important
   errors
   input_tags

Get in touch
------------

- If you have questions or suggestions for developing the code, please
  use **what interface**.  The developers cannot guarantee to offer
  support, though we will try to help. 
- Report bugs, suggest features or view the source code `on GitHub issues`_.

.. _on GitHub issues: http://github.com/OrderN/CONQUEST-release/issues

License
-------

CONQUEST is available freely under the open source `MIT Licence`__.
We ask that you acknowledge use of the code by citing appropriate
papers **details**.

__ https://choosealicense.com/licenses/mit/
