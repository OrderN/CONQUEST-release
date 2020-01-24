.. _moldyn:

==================
Molecular Dynamics
==================

CONQUEST can perform molecular dynamics both when the density matrix is computed
using diagonalisation and O(N), the latter allowing dynamical simulations of
(but not limited to) tens of thousands of atoms. The equations of motion are
integrated using the velocity Verlet method in the case of the microcanonical
ensemble (NVE), and modifications thereof for the canonical (NVT) and
isobaric-isothermal (NPT) ensembles, the details of which can be found in
:ref:`theory-md`. In addition to converging the parameters for the electronic
structure calculations, the following points must also be considered.

Go to :ref:`top <moldyn>`.

.. _md_scf:

Self-consistency tolerance and XL-BOMD
--------------------------------------

The convergence of the electronic structure is important in MD, as
insufficient convergence can be responsible for "drift" in the
conserved quantity of the dynamics. Although the molecular dynamics
integrators used in CONQUEST are time reversible, *the SCF procedure
is not*. Therefore tight convergence (``minE.SCTolerance`` for
diagonalisation, ``minE.LTolerance`` for linear scaling) is
necessary. In the case of diagonalisation, SCF tolerance of ``1E-6`` is
typically enough to negate the drift. However, extended-Lagrangian
Born-Oppenheimer MD (XL-BOMD) :cite:`md-Niklasson2008`, currently only
implemented for O(N), essentially makes the SCF component of the MD
time-reversible by adding the electronic degrees of freedom to the
Lagrangian, relaxing the constraint on ``minE.LTolerance`` ---
although it is still somewhat dependent on the ensemble.  In the NVE
and NVT ensembles, a L-tolerance of ``1E-5`` has been found to be
sufficient to give good energy conservations, decreasing to ``1E-6``
in the NPT ensemble. The following flags are required for XL-BOMD:

::

   DM.SolutionMethod ordern
   AtomMove.ExtendedLagrangian T

Go to :ref:`top <moldyn>`.

.. _md_restart:

Restarting
----------

Assuming the calculation ended gracefully, it can easily be restarted by
setting,

::

   AtomMove.RestartRun T

This will do several things: it will read the atomic coordinates from
``md.position`` and read the ``md.checkpoint`` file, which contains the
velocities and extended system (Nose-Hoover chain and cell) variables. Depending
on the value of ``DM.SolutionMethod``, it will read the K-matrix files
(``diagon``) or the L-matrix files (``ordern``), and if XL-BOMD is being used,
the X-matrix files. Finally, it will *append* new data to the ``md.stats`` and
``md.frames`` files, but it will overwrite all other files, including
``Conquest_out``. Note that this flag is equivalent to setting the following:

::

   General.LoadL T
   SC.MakeInitialChargeFromK T
   XL.LoadL T

In addition to the files mentioned above, CONQUEST will try to read the K-matrix
from ``Kmatrix2.i00.*`` when using diagonalisation or the L-matrix from
``Lmatrix2.i00.*`` when using O(N), and ``Xmatrix2.i0*.*`` if the
extended-Lagrangian formalism is used. Note that metadata for these files is
stored in ``InfoGlobal.i00.dat`` which is also required when restarting. If the
calculation ended by hitting the walltime limit, the writing of these matrix
files may have been interrupted, rendering them unusable. In this case, the
calculation can be restarted by setting the above flags to ``F`` *after* setting
``AtomMove.RestartRun T``. Setting the flag ``General.MaxTime`` to some number
of seconds less (say 30 minutes) than the calculation wall time limit will force
the calculation to stop gracefully, preventing the aforementioned situation.

Go to :ref:`top <moldyn>`.

.. _md_vis:

Visualising the trajectory
--------------------------

Setting the flag ``AtomMove.WriteXSF T`` dumps the coordinates to the file
``trajectory.xsf`` every ``AtomMove.OutputFreq`` steps. The .xsf file can be
read using `VMD <https://www.ks.uiuc.edu/Research/vmd/>`_. A small VMD script,
``view.vmd`` is included with the code, and can be invoked using,

``vmd -e view.vmd``

assuming the vmd executable is in your path.

Go to :ref:`top <moldyn>`.

.. _md_tdep:

TDEP output
-----------

CONQUEST molecular dynamics data can be used to perform lattice dyanmical
calculations using the `Temperature Dependent Effective Potential (TDEP)
<https://ollehellman.github.io/index.html>`_ code. Setting the flag ``MD.TDEP
T`` will make conquest dump configurations, forces and metadata in a format
readable by TDEP.

Go to :ref:`top <moldyn>`.

.. _md_nonh:

Non-Hamiltonian dynamics
------------------------

.. _md_nvt:

Canonical (NVT) ensemble
++++++++++++++++++++++++

The thermostat is set using the ``MD.Thermostat`` flag, and can take the values
``svr`` (stochastic velocity rescaling) and ``nhc`` (Nose-Hoover
chain). These thermostats generate the correct canonical ensemble
phase space distribution, and both give a conserved quantity that
allows the quality of the dynamics to be monitored.

1. Stochastic velocity rescaling

::

   AtomMove.IonTemperature 300.0
   MD.Ensemble nvt
   MD.Thermostat svr
   MD.tauT 10

While the NHC uses chaotic sensitivity to initial conditions to achieve better
ergodicity, the SVR thermostat :cite:`md-Bussi2007` uses a judiciously chosen stochastic force
coupled to a weak scaling thermostat to correctly generate the
canonical phase space distribution. The ``MD.tauT`` parameter gives
the coupling timescale; the velocity scaling factor is modified by a
factor :math:`\Delta t/\tau`, so a larger :math:`\tau` results in a
more slowly varying temperature.  While some characterisation of the
system is recommended, values of :math:`\tau` around 20--200fs are
reasonable.  To reproduce a simulation, the random number
generator seed can be set with the ``General.RNGSeed <integer>`` flag.

2. Nose-Hoover chain

::

   AtomMove.IonTemperature 300.0
   MD.Ensemble nvt
   MD.Thermostat nhc
   MD.nNHC 5
   MD.nYoshida 5
   MD.tauT 30

When thermostatting using a Nose-Hoover chain :cite:`md-Nose1984,md-Hoover1985,md-Martyna1992`, it may be necessary to set a
couple more flags. ``MD.nNHC`` sets the number of thermostats in the chain (the
default of 5 is generally sensible), and ``MD.nYoshida`` determines the order of
Yoshida-Suzuki integration. This is essentially a higher level integration
scheme that *can* improve energy conservation in cases when rapid changes in the
Nose-Hoover thermostat velocity is causing integration errors. **Note that**
``MD.tauT`` **means something different to the SVR case**. A good guess is
the time period of the highest frequency motion of the system in fs; however, in
the NVT ensemble, the energy conservation is not very sensitive to this value.
The NHC masses can also be set manually using the following block.

::

   MD.CalculateXLMass F
   MD.nNHC 5
   %block MD.NHCmass
     5 1 1 1 1
   %endblock

Go to :ref:`top <moldyn>`.

.. _md_npt:

Isobaric-Isothermal (NPT) ensemble
++++++++++++++++++++++++++++++++++

There is one implemented barostat at present, the extended
system, Parrinello-Rahman :cite:`md-Parrinello1981`. At present the
barostat should be treated as a beta-version implementation, which
will be fully characterised and made robust for the full release of
the code. 

1. Parrinello-Rahman

::

   AtomMove.IonTemperature 300.0
   AtomMove.TargetPressure 10.0
   MD.Ensemble npt
   MD.Thermostat nhc
   MD.Barostat pr
   MD.nNHC 5
   MD.nYoshida 5
   MD.tauT 100
   MD.tauP 200
   MD.PDrag 10.0

The Parrinello-Rahman barostat generates the correct ensemble, but can
be subject to low frequency "ringing" fluctuations in the 
temperature and pressure that can destabilise the system or slow equilibration.
Unlike in the NVT ensemble, this combination of barostat and thermostat is
*very* sensitive to the choice of both ``MD.tauT`` and ``MD.tauP``; note that
their values are somewhat higher in this case, since integration errors in the
NHC tend to be more severe due to coupling of the cell and atomic motions. They
are dependent on the system, so it is advised that you find a combination of
these parameters that gives the best energy conservation. The cell is
thermostatted using a separate Nose-Hoover chain to the atoms by default, but
they can be controlled with the same chain by setting ``MD.CellNHC F``. An *ad
hoc* drag factor specified by ``MD.PDrag`` reduces the thermostat and cell
velocities at every timestep to damp out the ringing fluctuations. In this case,
they are reduced by :math:`10/200 \simeq 5\%`, which strictly speaking breaks the NPT
dynamics, but not significantly, and the stability is significantly improved.

Note that the NPT ensemble can also be generated correctly by thermostatting
using the SVR thermostat, although the meaning of the parameter ``MD.tauT`` is
different in this case, as in NVT dynamics.

Postprocessing tools
--------------------

Details of Python post-processing tools for CONQUEST can be found in :ref:`et_md_scripts`.

Go to :ref:`top <moldyn>`.

.. bibliography:: references.bib
    :cited:
    :labelprefix: MD
    :keyprefix: md-
    :style: unsrt

Go to :ref:`top <moldyn>`.
	    
