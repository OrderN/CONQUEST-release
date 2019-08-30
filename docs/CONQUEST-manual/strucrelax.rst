.. _strucrelax:

=====================
Structural relaxation
=====================

This section describes how to find the zero-Kelvin equilibrium structure, given
a starting structure with non-zero forces and/or stresses. CONQUEST employs a
conjugate gradients algorithm to minimise energy or enthalpy with respect to
atomic positions and in some cases, cell vectors. Although many codes use a
quasi-Newton algorithm (some variant of BFGS), iteratively improving the Hessian
matrix does not scale well to the kinds of large systems that CONQUEST is
designed to solve. For pathological cases, a damped MD optimiser is also
available.

.. _sr_ions:

Ionic relaxation
----------------

To optimise the ionic positions with respect to the DFT total energy, the
following flags are essential:

::

   AtomMove.TypeOfRun cg
   AtomMove.MaxForceTol 5e-4
   AtomMove.ReuseL T

The first specifies a conjugate gradients relaxation, which is robust and
relatively efficient in most instances. The second specifies the force
convergence criterion, i.e. the calculation will terminate when the *maximum
single component of force* is below this threshold (rather than, for example,
the force residual). The third flag specifies that the K-matrix (for
diagonalisation calculations) or L-matrix (for order(N) calculations) from the
previous step will be used as an initial guess for the SCF cycle after
propagating the atoms. Setting ``AtomMove.WriteXSF T`` for all flavours of
optimisation will dump the trajectory to a .xsf file, which can be visualised
using `VMD <https://www.ks.uiuc.edu/Research/vmd/>`_. Setting
``AtomMove.AppendCoords T`` will the structure at each step to
``UpdatedAtoms.dat`` in the format of a CONQUEST structure input.

If the self-consistency tolerance is too low, the optimisation will fail to
converge with respect to the force tolerance; this may necessitate a tighter
``minE.SCTolerance`` for diagonalisation or ``minE.LTolerance`` for order(N).

For problematic cases where the conjugate gradients algorithm fails to find a
downhill search direction, it may be worth trying quenched molecular dyanamics,
which is akin to a steepest descent algorithm.

::

   AtomMove.TypeOfRun md
   AtomMove.QuenchedMD T
   AtomMove.MaxForceTol 5e-4
   AtomMove.ReuseL T

The FIRE algorithm is a variant of quenched MD that has been shown to outperform
conjugate gradients, and can be switched on using,

::

   AtomMove.TypeOfRun md
   AtomMove.FIRE T
   AtomMove.MaxForceTol 5e-4
   AtomMove.ReuseL T

.. _sr_cell:

Cell optimisation
-----------------

The unit cell can be optimised with respect to enthalpy *with fixed fractional
coordinates* (``AtomMove.OptCellMethod 1``) using the following input:

::

   AtomMove.TypeOfRun cg
   AtomMove.OptCell T
   AtomMove.OptCellMethod 1
   AtomMove.TargetPressure 1.0
   AtomMove.ReuseL T
   AtomMove.EnthalpyTolerance 1E-6
   AtomMove.StressTolerance 0.01

Here, we specify the target pressure in GPa and two new tolerances, the enthalpy
tolerance in Ha and the stress tolerance in GPa.

.. _sr_both:

Combined optimisation
---------------------

For simple crystals, the fractional ionic coordinates vary trivially with
changes in the lattice vectors; however for more complicated systems such as
molecular crystals and amorphous materials, it is necessary simultaneously relax
the ionic positions and lattice vectors. This can be done by setting
``AtomMove.OptCellMethod 3``

::

   AtomMove.TypeOfRun cg
   AtomMove.OptCell T
   AtomMove.OptCellMethod 3
   AtomMove.TargetPressure 1.0
   AtomMove.ReuseL T
   AtomMove.MaxForceTol 5e-4
   AtomMove.EnthalpyTolerance 1E-6
   AtomMove.StressTolerance 0.01

Note that the enthalpy will generally converge much more rapidly than the force
and stress, and that it may be necessary to tighten ``minE.SCTolerance``
(diagonalisation) or ``minE.LTolerance`` (order(N)) to reach the force
tolerance, if it is even possible.

Due to the nature of the complex partitioning system, large and sudden changes in volume
may cause the calculation to crash, particlularly in the case of combined
optimisation. In such cases, it may help to try ``AtomMove.OptCellMethod 2``,
which uses a simple but robust double-loop minimisation: a full ionic conjugate
gradients relaxation for the inner loop and a single cell steepest descent
relaxation for the outer loop. This is considerable less efficient, and is not
guaranteed to converge to the same minimum as ``AtomMove.OptCellMethod 3``, but
may help in particularly problematic cases.
