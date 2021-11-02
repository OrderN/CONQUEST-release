.. _strucrelax:

=====================
Structural relaxation
=====================

This section describes how to find the zero-Kelvin equilibrium atomic structure, given
a starting structure with non-zero forces and/or stresses. CONQUEST
can employ a variety of algorithms to minimise energy with respect to
atomic positions, including: stabilised quasi-Newton method (SQNM); L-BFGS; conjugate gradients; and damped
molecular dynamics (both MDMin and FIRE approaches).  The minimisation
of energy or enthalpy with respect to cell vectors is restricted to
conjugate gradients at present, though L-BFGS will be implemented.

Setting ``AtomMove.WriteXSF T`` for all flavours of optimisation will dump the
trajectory to the file ``trajectory.xsf``, which can be visualised using `VMD
<https://www.ks.uiuc.edu/Research/vmd/>`_ and `XCrysDen <http://http://www.xcrysden.org>`_.
Setting ``AtomMove.AppendCoords T``
will append the structure at each step to ``UpdatedAtoms.dat`` in the format of a
CONQUEST structure input.

For the SQNM, L-BFGS and conjugate gradients relaxations, the progress of the calculation can be
monitored by searching for the word ``GeomOpt``; grepping will print the
following:

::

   $ grep GeomOpt Conquest_out
   GeomOpt - Iter:    0 MaxF:   0.00329282 H:  -0.14168571E+03 dH:   0.00000000
   GeomOpt - Iter:    1 MaxF:   0.00331536 H:  -0.14168995E+03 dH:   0.00424155
   GeomOpt - Iter:    2 MaxF:   0.00350781 H:  -0.14168997E+03 dH:   0.00001651
   GeomOpt - Iter:    3 MaxF:   0.00504075 H:  -0.14169161E+03 dH:   0.00164389
   GeomOpt - Iter:    4 MaxF:   0.00725611 H:  -0.14169172E+03 dH:   0.00010500
   GeomOpt - Iter:    5 MaxF:   0.01134145 H:  -0.14169329E+03 dH:   0.00157361
   GeomOpt - Iter:    6 MaxF:   0.01417229 H:  -0.14169385E+03 dH:   0.00056077
   GeomOpt - Iter:    7 MaxF:   0.01434628 H:  -0.14169575E+03 dH:   0.00190304
   GeomOpt - Iter:    8 MaxF:   0.01711197 H:  -0.14170001E+03 dH:   0.00425400
   GeomOpt - Iter:    9 MaxF:   0.02040556 H:  -0.14170382E+03 dH:   0.00381110
   GeomOpt - Iter:   10 MaxF:   0.01095167 H:  -0.14170752E+03 dH:   0.00370442

In this example, MaxF is the maximum single force component, H is the enthalpy and dH is the
change in enthalpy.

Go to :ref:`top <strucrelax>`.

.. _sr_ions:

Ionic relaxation
----------------

To optimise the ionic positions with respect to the DFT total energy, the
following flags are essential:

::

   AtomMove.TypeOfRun sqnm
   AtomMove.MaxForceTol 5e-4
   AtomMove.ReuseDM T

The parameter ``AtomMove.TypeOfRun`` can take the values ``sqnm``, ``lbfgs`` or
``cg`` for iterative optimisation.  All three algorithms are robust and
relatively efficient in most instances; SQNM :cite:`sr-Schaefer2015` is recommended in most cases,
though if the initial forces are large it may be worth performing quenched
MD to reduce them (see below) before applying SQNM. The
parameter ``AtomMove.MaxForceTol`` specifies the force
convergence criterion in Ha/bohr, i.e. the calculation will terminate
when the largest force component on any atom is below this value.
The parameter
``AtomMove.ReuseDM``  specifies that the density matrix (the K-matrix for
diagonalisation or L-matrix for O(N) calculations) from the
previous step will be used as an initial guess for the SCF cycle after
propagating the atoms; this should generally decrease the number of SCF cycles
per ionic step.

If the self-consistency tolerance is too low, the optimisation may fail to
converge with respect to the force tolerance; this may necessitate a tighter
``minE.SCTolerance`` for diagonalisation (also possibly
``minE.LTolerance`` for O(N) calculations).  A grid which is too
coarse can also cause problems with structural relaxation to high tolerances.

For large initial forces or problematic cases where the relaxation algorithms fail to find a
downhill search direction, it may be worth trying quenched molecular dynamics,
which propagates the equations of motion following a simple NVE
approach, but resets the velocities to zero when the dot product of
force and velocity is zero.

::

   AtomMove.TypeOfRun md
   AtomMove.QuenchedMD T
   AtomMove.MaxForceTol 5e-4
   AtomMove.ReuseDM T

The FIRE algorithm :cite:`sr-Bitzek2006` is a variant of quenched MD
that has been shown to outperform conjugate gradients in some
circumstances. 

::

   AtomMove.TypeOfRun md
   AtomMove.FIRE T
   AtomMove.MaxForceTol 5e-4
   AtomMove.ReuseDM T

Go to :ref:`top <strucrelax>`.

.. _sr_cell:

Simulation cell optimisation
----------------------------

The simulation cell can be optimised with respect to enthalpy *with fixed fractional
coordinates* (``AtomMove.OptCellMethod 1``) using the following input:

::

   AtomMove.TypeOfRun cg
   AtomMove.OptCell T
   AtomMove.OptCellMethod 1
   AtomMove.TargetPressure 1.0
   AtomMove.ReuseDM T
   AtomMove.EnthalpyTolerance 1E-5
   AtomMove.StressTolerance 0.1

Go to :ref:`top <strucrelax>`.

.. _sr_both:

Combined optimisation
---------------------

For simple crystals, the fractional ionic coordinates vary trivially with
changes in the simulation cell lengths; however for more complicated systems such as
molecular crystals and amorphous materials, it is necessary simultaneously relax
the ionic positions and simulation cell lengths (recalling that CONQUEST only
allows *orthorhombic* unit cells). This can be done by setting
``AtomMove.OptCellMethod 3``

::

   AtomMove.TypeOfRun cg
   AtomMove.OptCell T
   AtomMove.OptCellMethod 3
   AtomMove.TargetPressure 1.0
   AtomMove.ReuseDM T
   AtomMove.MaxForceTol 5e-4
   AtomMove.EnthalpyTolerance 1E-5
   AtomMove.StressTolerance 0.1

Note that the enthalpy will generally converge much more rapidly than the force
and stress, and that it may be necessary to tighten ``minE.SCTolerance``
(diagonalisation) or ``minE.LTolerance`` (order(N)) to reach the force
tolerance, if it is even possible.

Due to the nature of the complex partitioning system, large and sudden changes in volume
may cause the calculation to crash, particlularly in the case of combined
optimisation. In such cases, it may help to try ``AtomMove.OptCellMethod 2``,
which uses a simple but robust double-loop minimisation: a full ionic conjugate
gradients relaxation followed by a full simulation cell conjugate gradients
relaxation. This is considerably less efficient, but
may help in particularly problematic cases.

Go to :ref:`top <strucrelax>`.

.. bibliography:: references.bib
    :cited:
    :labelprefix: SR
    :keyprefix: sr-
    :style: unsrt

Go to :ref:`top <strucrelax>`.
