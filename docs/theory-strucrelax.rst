.. _theory-strucrelax:

=============================
Structural relaxation: Theory
=============================

Structural relaxation involves optimisation of the ionic coordinates,
optimisation of the simulation cell, or both, with respect to the DFT total
energy or the enthalpy if the cell is not fixed.

.. _th_sr_ions:

Ionic relaxation
----------------

.. _th_sr_lbfgs:

L-BFGS:
~~~~~~~

To be written...

.. _th_sr_cg:

Conjugate gradients
~~~~~~~~~~~~~~~~~~~

The most naive geometry optimisation algorithm is steepest descent: we calculate
the gradient of the DFT total energy (i.e. the force) and propagate the system
in the direction of the steepest gradient (the direction of the force vector)
until the energy stops decreasing. We choose the direction (largest gradient in
this case) and perform a line search. This will be sufficient if the potential
energy surface is well-behaved, but in most cases convergence will require many
iterations. Conjugate gradients is a well-established method the improves upon
steepest descent in the choice of search direction. Without going into too much
detail, we choose a new search direction that is orthogonal to all previous
search directions using the *conjugacy ratio* :math:`\beta`. At iteration
:math:`n`, it is given by,

.. math::
   
    \beta_n = \beta_{n-1} + \frac{\mathbf{f}_n^T\mathbf{f}_n}{\mathbf{f}_{n-1}^T\mathbf{f}_{n-1}}

This is the Fletcher-Reeves formulation; note that :math:`\beta_0 = 0`. We can
then construct the search direction at step :math:`n`, :math:`D_n`, 

.. math::

   D_n = \beta_n D_{n-1} + \mathbf{f}_n,

and peform the line minimisation in this direction. This process is repeated
until the maximum force component is below some threshold.

Go to :ref:`top <theory-strucrelax>`.

.. _th_sr_qmd:

Quenched MD
~~~~~~~~~~~

The system is propagated in the direction of steepest descent as determined by
the DFT forces, and the velocity is scaled down as the system approaches its
zero-temperature equilibrium configuration.

Go to :ref:`top <theory-strucrelax>`.

.. _th_sr_fire:

FIRE Quenched MD
~~~~~~~~~~~~~~~~

The system is propagated using the modified equation of motion :cite:`t2-Bitzek2006`,

.. math::

   \mathbf{\dot{v}}(t) = \mathbf{F}(t)/m -
   \gamma(t)|\mathbf{v}(t)|[\mathbf{\hat{v}}(t) - \mathbf{\hat{F}}(t)]

which has the effect of introducing an acceleration in a direction that is
steeper than the current direction of motion. If the power :math:`P(t) =
\mathbf{F}(t)\cdot\mathbf{v}(t)` is positive then the system is moving
"downhill" on the potential energy surface, and the stopping criterion is when
it becomes negative (moving "uphill").

Go to :ref:`top <theory-strucrelax>`.

.. _th_sr_cell_opt:

Cell optimisation
-----------------

When optimising the cell with *fixed fractional ionic coordinates*, the same
conjugate gradients method is used as above, but minimising the enthalpy with
respect to the cell vectors.

Go to :ref:`top <theory-strucrelax>`.

.. _th_sr_comb_opt:

Combined optimisation
---------------------

The ionic and cell degrees of freedom can be relaxed simultaneously by combining
all of their coordinates into a single vector and optimising them with respect
to the enthalpy of the system. However, this atomic forces and total stresses
having numerical values of the same order of magnitude, and changes in ionic
coordinates and cell vectors being of the same order of magnitude. Using the
method of Pfrommer *et al*. :cite:`t2-Pfrommer1997`, the latter can be enforced
by using fractional coordinates for the ionic positions, and fractional lattice
vectors of the form :math:`h = (1 + \epsilon)h_0` where h is the matrix of
lattice vectors, :math:`h_0` is the matrix for some reference configuration and
epsilon is the strain matrix. The "fractional" force on the *i* th atom is then
:math:`\mathbf{F}_i = g\mathbf{f}_i` where :math:`\mathbf{f}_i` is the
DFT-calculated force multiplied by the metric tensor :math:`g = h^Th`. The
"fractional" stress is,

.. math::

   f^{(\epsilon)} = -(\sigma + p\Omega)(1 + \epsilon^T)

where :math:`\sigma` is the DFT-calculated stress, :math:`p` is the target
pressure and :math:`\Omega` is the volume. The resulting vector is optimised
using the same conjugate gradients algorithm as before, minimising the enthalpy.


.. bibliography:: references.bib
    :cited:
    :labelprefix: Tb
    :keyprefix: t2-
    :style: unsrt

Go to :ref:`top <theory-strucrelax>`.


