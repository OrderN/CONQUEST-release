.. _strucrelax:

=====================
Structural relaxation
=====================

Structural relaxation includes optimisation of the ionic coordinates,
optimisation of the simulation cell, and both processes combined in
some fashion.

.. _sr_ions:

Ionic relaxation
----------------

.. _sr_cell:

Cell optimisation
-----------------

Simulation cell optimization is the minimization of the total energy of a
system of electrons and ions with respect to the lattice vectors. CONQUEST,
as currently implemented, can only operate with cubic (:math:`a=b=c`),
tetragonal (:math:`a=b, b\neq c, a\neq c`) and orthorhombic (:math:`a\neq b \neq c`)
cells. This reduces the complexity of the optimization problem to simply

.. math::
   
    min[E_{dft}(a, b, c)]

That is, the values of a, b and c which minimize the the total energy of the
system as predicted by DFT. Note that we do not include the set of ionic positions,
:math:`\{\mathbf{R_{I}} \}`, as a variable here and rather work with fixed ionic
positions.

Unconstrained case
~~~~~~~~~~~~~~~~~~

The first implementation is of an algorithm which allows a, b and c to vary freely.
A widely used method of solving this problem is that of non-linear conjugate gradients.
The algorithm works as follows:

1. Calculate the gradients of the total energy with respect to a, b and c.
These derrivatives are in this case proportional to the diagonal elements of the
stress tensor :math:`\sigma_x` :math:`\sigma_y` and :math:`\sigma_z`.

2. Calcuate the conjugacy update parameter for this step, :math:`\beta_n`
where n is the current step of the algorithm. We use the Flethcher-Reeves version:

.. math::
   
    \beta_n = \beta_{n-1} + \frac{\mathbf{\sigma_n^T}\mathbf{\sigma_n}}{\mathbf{\sigma_{n-1}^T}\mathbf{\sigma_{n-1}}}

Note that :math:`\beta_0 = 0`.

3. Build the search directions

.. math::

   D_{j,n} = \beta_n D_{j, n-1} + \sigma_{j,n}

Where D represents the search direction and the indicies j indicate a cartestion
direction. In this case, j = 1, 2 or 3.

4. Perform line minimizations in these directions:

.. math::

   a_{j,n+1} = a_{j,n} + \alpha_{min} D_{j,n}

Where :math:`a_j` represents a simulation box dimension. e.g. :math:`a_0 = a`,
:math:`a_1 = b`, :math:`a_3 = c`. :math:`\alpha_{min}` is the value of :math:`\alpha`
which minimizes the total energy given the set of :math:`a_j`.

Steps 1-4 are repeated until a small energy residual is reached (energy is currently
used as a surrogate for a stress residual due to known issues).

Constrained case
~~~~~~~~~~~~~~~~~

The current implementation can handle fixing both one and two dimensions
of the simulation cell. The opttimiser will then optimise the dimensions which
are free. What follows is that the vector :math:`\sigma_n` from eq (2) reduces in
stride to accomodate the change. Then, the number of search directions and line
minimisations from eqs (3) and (4) is also reduced. This, in turn, reduces the
complexity of the problem which allows the conjugate gradient algorithm to converge
to a minumum in less iterations.

It is possible to fix ratios of cell dimensions at their initial ratio set in the
coordinate file. Say, for example, we wanted to fix the ratio c/a. What this involves
is appying the constraint

.. math::

   \frac{c_n}{a_n} = \frac{c_{n+1}}{a_{n+1}}

Using the above with eq (4), it can be shown that

.. math::

   c_{n+1} = c_n + \alpha_{min}\frac{c_n}{a_n}D_{n, a}

   .. math::

   a_{n+1} = a_n + \alpha_{min}\frac{a_n}{c_n}D_{n, c}

The third dimension, b, is allowed to vary freely as before.
All other combinations of ratios follow the same logic and lead to similar expressions.


It is possible to scale a.b and c by a global scaling factor by choosing to minimize
the mean stress :math:`\bar{\sigma_n}`. We then build the single search direction
:math:`\bar{D_n}` from eq (3) and perfrom the subsequent line minimizations from eq (4)
using only :math:`\bar{D_n}`.

This approach can be useful optimising cubic phases where
:math:`\sigma_x  = \sigma_y = \sigma_z` but also when investigating
volume changes.
      
.. _sr_both:

Combined optimisation
---------------------
