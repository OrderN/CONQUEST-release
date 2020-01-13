==========================
Frequently Asked Questions
==========================

When should I use CONQUEST?
---------------------------
You can use CONQUEST for any DFT simulations that you need to
perform.  It is efficient for small problems, though may not be as
efficient as other codes (e.g. plane wave codes) because it has been
designed for massively parallel operation, which brings some
overhead.  If you need to perform DFT calculations on large systems
(several hundred atoms or beyond) or want to perform highly parallel
calculations, you should definitely consider CONQUEST.

CONQUEST uses `Hamann`_ optimised norm-conserving Vanderbilt (ONCV)
pseudopotentials, which can also be used by `PWSCF`_
and `Abinit`_ which allows direct comparisons between the codes.

.. _Hamann: http://www.mat-simresearch.com
.. _PWSCF: https://www.quantum-espresso.org
.. _Abinit: https://www.abinit.org

When should I use linear scaling?
---------------------------------
You should use linear scaling if you need to model systems with more
than about 5,000 atoms, though gains are often found for smaller
systems (from 1,000 atoms upwards).

Linear scaling calculations offer the prospect of scaling to
significantly larger systems than traditional DFT calculations;
however, they make approximations and require some care and
characterisation.  In particular, instead of solving for eigenvalues
and eigenstates, linear scaling methods solve for the density matrix,
so that energy-resolved information (e.g. DOS and band energies) are
not available.  To enable linear scaling, a range is also imposed on the
density matrix and it is important to test the effect of this range.

Will you implement a specific feature for me?
---------------------------------------------
We cannot guarantee to implement specific features, though we are
always happy to take suggestions.  We also welcome new developers: if
there is something that you would like to see in the code, please do
talk to us about joining the development effort.

How do I report a bug?
----------------------
Please use the `GitHub issues`_ page.  Include details of the compiler
and libraries used, the version of CONQUEST, and the input and output
files (if possible).  We will do our best to check the bug and fix it,
but cannot guarantee to help on any timescale.

.. _GitHub issues: http://github.com/OrderN/CONQUEST-release/issues


How do I get help?
------------------
The Conquest mailing list (**details**) is the best place to get
help.  However, the developers cannot guarantee to answer any
questions, though they will try.  Bug reports should be made through
the `GitHub issues`_ page.

.. _GitHub issues: http://github.com/OrderN/CONQUEST-release/issues
