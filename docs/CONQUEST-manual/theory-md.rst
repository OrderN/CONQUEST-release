.. _theory-md:

==========================
Molecular Dynamics: Theory
==========================

Microcanonical (NVE) ensemble
-----------------------------

The Hamiltonian for the microcanonical ensemble is,

.. math::
  \mathcal{H} = \sum_{i=1}^N \frac{\mathbf{p}_i^2}{2m_i} + U(\mathbf{r}_i)

where :math:`\mathbf{p}_i` and :math:`\mathbf{r}_i` are the position and momentum of particle :math:`i` and :math:`U` is the DFT total (potential) energy. Hamilton's equations can be solved to give the following equations of motion:

.. math::
  \mathbf{\dot{r}}_i &= \frac{\mathbf{p}_i}{m_i} \\
  \mathbf{\dot{p}}_i &= \frac{\partial U(\mathbf{r}_i)}{\partial\mathbf{r}_i} = \mathbf{F_i}

In order to construct a time-reversible algorithm from these equations, the Liouvillian formulation is employed :cite:`t-Frenkel2002` (trivially, in this case). The Liouville operator :math:`L` can be defined in terms of position and momentum components:

.. math::
  iL = \mathbf{\dot{r}}\frac{\partial}{\partial\mathbf{r}} + \mathbf{\dot{p}}\frac{\partial}{\partial\mathbf{p}} = i(L_r + L_p).

The Liouvillian can be used to construct the classical propagator, which relates the state :math:`f` of the system at time 0 to its state at time :math:`t`:

.. math::
  f[\mathbf{p}^N(t),\mathbf{r}^N(t)] = e^{iLt}f[\mathbf{p}^N(0),\mathbf{r}^N(0)]

Taking the individual position and momentum parts of the Liouvillian :math:`L_r` and :math:`L_p`, it can be shown that applying it to the state :math:`f` result in a simple linear shift in coordinates and a simple linear shift in momentum respectively:

.. math::
  iL_rf(t) &= f[\mathbf{p}^N(0),\mathbf{r}^N(0) + \mathbf{\dot{r}}^N(0)t] \\ 
  iL_pf(t) &= f[\mathbf{p}^N(0) + \mathbf{F}^N(0)t,\mathbf{r}^N(0)]

However, we cannot simply replace :math:`e^{iLt}` with :math:`e^{iL_rt}` because :math:`iL_r` and :math:`iL_p` are non-commuting operators, so we must employ the Trotter-Suzuki identity:

.. math::
  e^{A+B} = \lim_{P\rightarrow\infty}\left(e^{A/2P}e^{B/P}e^{A/2P}\right)^P

Thus for a small enough time step :math:`\Delta t = t/P` and to first order, a discrete time step corresponds to the application of the discrete time propagator :math:`G`,

.. math::
  G(\Delta t) = U_1\left(\frac{\Delta t}{2}\right)U_2\left(\Delta t\right)U_1\left(\frac{\Delta t}{2}\right) = e^{iL_1\frac{\Delta t}{2}}e^{iL_2\Delta t}e^{iL_1\frac{\Delta t}{2}},

which can be shown to be unitary and therefore time-reversible. Applying the operators :math:`U` in the sequence determined by the Trotter decomposition generates the velocity Verlet algorithm, which is used to integrate microcanonical molecular dynamics in CONQUEST. For a detailed derivation of the algorithm, refer to Frenkel & Smit :cite:`t-Frenkel2002`.

Extended Lagrangian Born-Oppenheimer MD (XL-BOMD)
-------------------------------------------------

If the electronic density from the previous ionic step is used as an initila guess for the next SCF cycle, a problem arises because this process breaks the time-reversibility of the dynamics. This is manifested as a gradual drift in the total energy in the case of a NVE simulation, or the conserved quantity in the case of non-Hamiltonian dynamics. The solution proposed by Niklasson :cite:`t-Niklasson2008,t-Niklasson2014` is to introduce auxilliary electronic degrees of freedom into the Lagrangian, which can be propagated via time-reversible integrators.

The extended Lagrangian used in CONQUEST is :cite:`t-Arita2014`,

.. math::
  \mathcal{L}^\mathrm{XBO}\left(\mathbf{X}, \mathbf{\dot{X}}, \mathbf{R}, \mathbf{\dot{R}}\right) = \mathcal{L}^\mathrm{BO}\left(\mathbf{R}, \mathbf{\dot{R}}\right) + \frac{1}{2}\mu\mathrm{Tr}\left[\mathbf{\dot{X}}^2\right] - \frac{1}{2}\mu\omega^2\mathrm{Tr}\left[(\mathbf{LS} - \mathbf{X})^2\right],

where :math:`\mathbf{X}` is a sparse matrix with the same range as :math:`\mathbf{LS}`, :math:`\mu` is the fictitious electron mass and :math:`\omega` is the curvature of the auxiliary harmonic potential. The Euler-Lagrange equations of motion are then,

.. math::
  m_i\mathbf{\ddot{r}_i} &= -\frac{\partial U[{{\mathbf{R;LS}}}]}{\partial\mathbf{r}_i} = \mathbf{F_i} \\
  \mathbf{\ddot{X}} &= \omega^2(\mathbf{LS} - \mathbf{X}),

The first of these is simply Newton's second law, and the velocity update equation of motion in the microcanonical ensemble. The second can be integrated using a time-reversible algorithm, the velocity Verlet scheme in the case of CONQUEST :cite:`t-Arita2014`:

.. math::
  \mathbf{X}(t+\delta t) &= 2\mathbf{X}(t) -\mathbf{X}(t-\delta t) + \delta t^2\omega^2\left[\mathbf{L}(t)\mathbf{S}(t)-\mathbf{X}(t)\right] \\
  &+ a\sum_{m=0}^M c_m\mathbf{X}(t-m\delta t)

i.e. the trajectory of :math:`\mathbf{X}(t)` is time-reversible, and evolves in a harmonic potential centred on the ground state density :math:`\mathbf{L}(t)\mathbf{S}(t)`. The matrix :math:`\mathbf{XS}^{-1}` is a good guess for the :math:`\mathbf{L}` matrix in the Order(N) scheme.

Despite the time-reversitibility, the :math:`\mathbf{X}` matrix tends in practice to gradually drift from the harmonic centre over time, increasing the number of SCF iterations required to reach the minimum over the course of the simulation. To remove such numerical errors, the final dissipative term is included, and is found to have a minimal effect on the time-reversibility. We note that since the auxiliary variable :math:`X` is used to generate an intial guess for the SCF process, it does not appear in the conserved (pseudo-Hamiltonian) quantity for the dynamics.

Non-Hamiltonian dynamics
------------------------

Extended system method
~~~~~~~~~~~~~~~~~~~~~~

Hamiltonian dynamics generally describes systems that are isolated from their surroundings, but in the canonical and isobaric-isothermal ensembles, we need to couple the system to an external heat bath and/or stress. It is possible to model such systems by positing a set of equations of *non-Hamiltonian* equations of motion, and proving that they generate the correct statistical ensemble :cite:`t-Tuckerman2010`. This is the extended system approach: we modify the Hamiltonian to include the thermostat and/or barostat degrees of freedom, derive the (pseudo-) Hamiltonian equations of motion, and demostrate that the correct phase space distribution for the ensemble is recovered.

Canonical (NVT) ensemble
~~~~~~~~~~~~~~~~~~~~~~~~

In the Nose-Hoover formulation :cite:`t-Nose1984,t-Hoover1985`, the Hamiltonian for a system in the canonical ensemble can be written,

.. math::
  \mathcal{H} = \sum_i \frac{1}{2}m_i s^2\mathbf{\dot{r}}_i^2 + U(\mathbf{r}_i) + \frac{1}{2}Q\dot{s}^2 - (n_f + 1)k_B T \ln s,

where :math:`\mathbf{r}_i` and :math:`\mathbf{\dot{r}_i}` are respectively the position and velocity of particle :math:`i`, :math:`U` is the potential energy, in this case the DFT total energy, :math:`s` is a dimensionless quantity that can be interpreted post-hoc as a time step scaling factor, :math:`Q` is the fictitious mass of the heat bath and :math:`n_f` is the number of ionic degrees of freedom. Hamilton's equations can be solved to generate the Nose-Hoover equations of motion. However Martyna *et al*. demonstrate that this method does not generate an ergodic trajectory, and proposed an alternative formulation :cite:`t-Martyna1992` in which the temperature is controlled by a chain of :math:`M` coupled thermostats of mass :math:`Q_k`, notional position :math:`\eta_k` and conjugate momentum :math:`p_{\eta_k}`:

.. math::
  \mathbf{\dot{r}_i} &= \frac{\mathbf{p}_i}{m_i} \\
  \mathbf{\dot{p}_i} &= -\frac{\partial U(\mathbf{r})}{\partial \mathbf{r}_i} - \frac{p_{\eta_1}}{Q_1}\mathbf{p}_i \\
  \dot{\eta}_k &= \frac{p_{\eta_k}}{Q_k} \\
  \dot{p}_{\eta_1} &= \left(\sum_{i=1}^N\frac{\mathbf{p}_i}{m_i} - n_fk_BT\right) - \frac{p_{\eta_{2}}}{Q_{\eta_{2}}}p_{\eta_1} \\
  \dot{p}_{\eta_k} &= \left(\frac{p^2_{\eta_{k-1}}}{Q_{k-1}} - k_BT\right) - \frac{p_{\eta_{k+1}}}{Q_{k+1}}p_{\eta_k} \\
  \dot{p}_{\eta_M} &= \left(\frac{p^2_{\eta_{M-1}}}{Q_{M-1}} - k_BT\right)

The Liouvillian for these equations of motion can be non-uniquely decomposed into components of ionic position (:math:`iL_r`) and momentum (:math:`iL_p`) as in the microcanonical case, the extended Lagrangian (:math:`iL_\mathrm{XL}`, and a Nose-Hoover chain component (:math:`iL_\mathrm{NHC}`)

.. math::
  iL = iL_\mathrm{NHC} + iL_p + iL_{\mathrm{XL}} + iL_r,

which is directly translated into an algorithm with the Trotter-Suzuki expansion,

.. math::
  \exp(iL\Delta t) = &\exp\left(iL_\mathrm{NHC}\frac{\Delta t}{2}\right)\exp\left(iL_p\frac{\Delta t}{2}\right) \times \\
  &\exp\left(iL_\mathrm{XL}\frac{\Delta t}{2}\right)\exp\left(iL_r\Delta t\right)\exp\left(iL_\mathrm{XL}\frac{\Delta t}{2}\right) \times \\
  &\exp\left(iL_p\frac{\Delta t}{2}\right)\exp\left(iL_\mathrm{NHC}\frac{\Delta t}{2}\right)

This is recognisable as the velocity Verlet algorithm with extended Lagrangian integration which can be reduced to a single step, as described in :ref:`Extended Lagrangian Born-Oppenheimer MD (XL-BOMD)`, with a half time step integration of the Nose-Hoover chain equations of motion before and after. For full details of the integration scheme, see Hirakawa *et al*. :cite:`t-Hirakawa2017`.

Isobaric-Isothermal (NPT) ensemble
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Parinello-Rahman equations of motion :cite:`t-Parrinello1981` extend the fixed cell equations of motion to include the cell degrees of freedom in the extended system approach. We use the Martyna-Tobias-Tuckerman-Klein modification :cite:`t-Martyna1996`, which couples the variable cell equations of motion to a Nose-Hoover chain the themrostat the system, recovering the isobaric-isothermal (NPT) ensemble. For an unconstrained cell (i.e. the lattice vectors can change freely), the equations of motion are,

.. math::
  \mathbf{\dot{r}}_i &= \frac{\mathbf{p}_i}{m_i} + \frac{\mathbf{p}_g}{W_g}\mathbf{r}_i \\
  \mathbf{\dot{p}}_i &= \mathbf{F}_i - \frac{\mathbf{p}_g}{W_g}\mathbf{p}_i - \left(\frac{1}{N_f}\right)\frac{\mathrm{Tr}[\mathbf{p}_g]}{W_g}\mathbf{p}_i - \frac{p_\xi}{Q}\mathbf{p}_i \\
  \mathbf{\dot{h}} &= \frac{\mathbf{p}_g\mathbf{h}}{W_g} \\
  \mathbf{\dot{p}_g} &= V(\mathbf{P}_\mathrm{int}-\mathbf{I}P_\mathrm{ext}) + \left[\frac{1}{N_f}\sum_{i=1}^N\frac{\mathbf{p}_i^2}{m_i}\right]\mathbf{I} - \frac{p_\xi}{Q}\mathbf{p}_g \\
  \dot{\xi} &= \frac{p_\xi}{Q} \\
  \mathbf{\dot{p}}_g &= \sum_{i=1}^N\frac{\mathbf{p}_i^2}{m_i} + \frac{1}{W_g}\mathrm{Tr}[\mathbf{p}_g^T\mathbf{p}_g] - (N_f + d^2)kT
   
Here, :math:`\mathbf{r}_i`, :math:`\mathbf{p}_i` and :math:`m_i` are the position, momentum and mass of particle :math:`i` respectively, :math:`\xi`, :math:`p_\xi` and :math:`Q` are the position, momentum and mass of the thermostat, and :math:`\mathbf{h}`, :math:`\mathbf{p}_g` and :math:`W_g` are the matrix of lattice vectors, matrix of cell velocities and cell mass respectively. Note that these equations only include one Nose-Hoover thermostat for simplicity. Conquest uses the Shinoda-Shiga-Mikami splitting of the Liouvillian :cite:`t-Shinoda2004` to propagate the system. The Liouvillian is decomposed as,

.. math::
  iL = iL_r + iL_h + iL_v + iL_\mathrm{bath},

which can be further split,

.. math::
  iL_\mathrm{bath} &= iL_\mathrm{box} + iL_\mathrm{particles} \\
  iL_\mathrm{box} &= iL_\mathrm{vbox} + iL_\xi + iL_{v_{\xi_1}} + iL_{v_{\xi_k}} + iL_{v_{\xi_M}} \\
  iL_\mathrm{particles} &= iL_\mathrm{vpart} + iL_\xi + iL_{v_{\xi_1}} + iL_{v_{\xi_k}} + iL_{v_{\xi_M}}

Using Liouville's theorem, we have,

.. math::
  iL_r &= \sum_{i=1}^N[\mathbf{v}_i + \mathbf{v}_g\mathbf{r}_i]\cdot\nabla_{\mathbf{r}_i} \\
  iL_h &= \sum_{\alpha,\beta}\mathbf{v}_{g,\alpha\beta}\mathbf{h}_{\alpha\beta}\frac{\partial}{\partial\mathbf{h}_{\alpha\beta}} \\
  iL_v &= \sum_{i=1}^N\left(\frac{\mathbf{F}_i}{m_i}\right)\cdot\nabla_{\mathbf{v}_i} \\
  iL_\mathrm{bath} &= iL_\mathrm{vpart} + iL_\mathrm{vbox} + iL_\xi + iL_{v_{\xi_1}} + iL_{v_{\xi_k}} + iL_{v_{\xi_M}} \\
  &= \sum_{i=1}^N\left[-\left\{\mathbf{v}_g + \frac{1}{N_f}\mathrm{Tr}(\mathbf{v}_g) + v_{\xi_1}\right\}\mathbf{v}_i\right]\nabla_{\mathbf{v}_i} \\
  &+ \sum_{\alpha,\beta}\left[\frac{F_\mathrm{box}}{W} - v_{\xi_1}\mathbf{v}_{g,\alpha\beta}\right]\frac{\partial}{\partial\mathbf{v}_{g,\alpha\beta}} \\
  &+ \sum_{k=1}^M v_{\xi_k}\frac{\partial}{\partial\xi_k} \\
  &+ \left[\frac{F_{\mathrm{NHC}_1}}{Q_1} - v_{\xi_1}v_{\xi_2}\right]\frac{\partial}{\partial v_{\xi_1}} \\
  &+ \sum_{k=2}^M\left[\frac{1}{Q_k}(Q_{k-1}v_{\xi_{k-1}}^2 - kT_\mathrm{ext}) - v_{\xi_k}v_{\xi_{k+1}}\right]\frac{\partial}{\partial v_{\xi_k}} \\
  &+ \left[\frac{1}{Q_M}(Q_{M-1}v_{\xi_{M-1}}^2 - kT_\mathrm{ext})\right]\frac{\partial}{\partial v_{\xi_M}}

Here we use :math:`M` heat baths in a Nose-Hoover chain. The Trotter-Suzuki expansion is,

.. math::
  e^{iL\Delta t} = e^{iL_\mathrm{bath}\frac{\Delta t}{2}}e^{iL_v\frac{\Delta t}{2}}e^{iL_h\frac{\Delta t}{2}}e^{iL_r\Delta t}e^{iL_h\frac{\Delta t}{2}}e^{iL_v\frac{\Delta t}{2}}e^{iL_\mathrm{bath}\frac{\Delta t}{2}}.

The Liouvillian for the heat baths can be further expanded:

.. math::
  e^{iL_\mathrm{particles}\frac{\Delta t}{2}} = e^{\left(iL_{v_{\xi_1}} + iL_{v_{\xi_k}} + iL_{v_{\xi_M}}\right)\frac{\Delta t}{4}}e^{\left(iL_\xi + iL_\mathrm{vpart}\right)\frac{\Delta t}{2}}e^{\left(iL_\xi + iL_{v_{\xi_1}} + iL_{v_{\xi_k}} + iL_{v_{\xi_M}}\right)\frac{\Delta t}{4}}

Finally, expanding the first propagator in the previous expression, we have,

.. math::
  e^{\left(iL_{v_{\xi_1}} + iL_{v_{\xi_k}} + iL_{v_{\xi_M}}\right)\frac{\Delta t}{4}} &= e^{-i\left(-v_{\xi_1}v_{\xi_2}\frac{\partial}{\partial \xi_1} - \sum_{k=2}^Mv_{\xi_k}v_{\xi_{k+1}}\frac{\partial}{\partial \xi_k} - v_{\xi_{M-1}}v_{\xi_M}\frac{\partial}{\partial \xi_M}\right)\frac{\Delta t}{8}} \\
  &\times e^{i\left(F_{\mathrm{NHC}_1}\frac{\partial}{\partial v_{\xi_1}} + F_{\mathrm{NHC}_k}\frac{\partial}{\partial v_{\xi_k}} + F_{\mathrm{NHC}_M}\frac{\partial}{\partial v_{\xi_M}}\right)\frac{\Delta t}{4}} \\
  &\times e^{-i\left(-v_{\xi_1}v_{\xi_2}\frac{\partial}{\partial \xi_1} - \sum_{k=2}^Mv_{\xi_k}v_{\xi_{k+1}}\frac{\partial}{\partial \xi_k} - v_{\xi_{M-1}}v_{\xi_M}\frac{\partial}{\partial \xi_M}\right)\frac{\Delta t}{8}}

These expressions are directly translated into the integration algorithm.


Weak coupling thermostat/barostat
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Instead of modifying the Hamiltonian, the Berendsen-type weak coupling method :cite:`t-Berendsen1984` involves coupling the ionic degrees of freedom to a an external temperature and/or pressure bath via "the principle of least local perturbation consistent with the required global coupling." Thermostatting is acheived via a Langevin-type equation of motion, in which the system is globally coupled to a heat bath and subjected to random noise:

.. math::
    m_i\ddot{\mathbf{r}}_i = \mathbf{F}_i + m_i \gamma\left(\frac{T_0}{T}-1\right)\dot{\mathbf{r}}_i,

where :math:`\gamma` is a global friction constant chosen to be the same for all particles. This can be acheived in practice by rescaling the velocities :math:`\mathbf{v}_i \rightarrow \lambda\mathbf{v}_i`, where :math:`\lambda` is,

.. math::
    \lambda = \left[ 1 + \frac{\Delta t}{\tau_T}\left(\frac{T_0}{T}-1\right)\right]^{\frac{1}{2}}

A similar argument can be applied for weak coupling to an external pressure bath. In the isobaric-isoenthalpic ensemble, the velocity of the particles can be expressed,

.. math::
    \dot{\mathbf{r}} = \mathbf{v} - \frac{\beta(P_0 - P)}{3\tau_P}\mathbf{r},

i.e. the fractional coordinates are scaled by a factor determined by the difference between the internal and external pressures, the isothermal compressibility :math:`\beta` and a pressure coupling time constant $\tau_P$. In the isotropic case, the cell scaling factor :math:`\mu` can be expressed,

.. math::
    \mu = \left[ 1 - \frac{\Delta t}{\tau_P}(P_0 - P)\right]^{\frac{1}{3}},

where the compressibility is absorbed into the time time constant :math:`\tau_P`. Allowing for fluctuations of all cell degrees of freedom, the scaling factor becomes,

.. math::
    \mathbf{\mu} = \mathbf{I} - \frac{\beta\Delta t}{3\tau_P}(\mathbf{P}_0 - \mathbf{P})

While trivial to implement and in general stable, the weak-coupling method does not recover the correct phase space distribution for the canonical or isobaric-isothermal ensembles, for which the extended system method is required.

Stochastic velocity rescaling
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Stochastic velocity rescaling (SVR) :cite:`t-Bussi2007` is a modification of the weak coupling method, in which a correctly constructed random force is added to enforce the correct NVT (or NPT) phase space distribution. The kinetic energy is rescaled such that the change in kinetic energy between thermostatting steps is,

.. math::
  dK = (\bar{K} - K)\frac{dt}{\tau} + 2\sqrt{\frac{K\bar{K}}{N_f}}\frac{dW}{\sqrt{\tau}}

where :math:`\bar{K}` is the target kinetic energy (external temperature), :math:`dt` is the time step, :math:`\tau` is the time scale of the thermostat, :math:`N_f` is the number of degrees of freedom and :math:`dW` is a Wiener process. Practically, the particle velocities are rescaled by a factor of :math:`\alpha`, defined via,

.. math::
  \alpha^2 = e^{-\Delta t/\tau} + \frac{\bar{K}}{N_fK}\left(1-e^{-\Delta t/\tau}\right)\left(R_1^2 + \sum_{i=2}^{N_f}R_i^2\right) + 2e^{-\Delta t/2\tau}\sqrt{\frac{\bar{K}}{N_fK}\left(1-e^{-\Delta t/\tau}\right)R_1}

Where :math:`R_i` is a set of :math:`N_f` normally distributed random numbers with unitary variance. This method can be applied to thermostat the NPT ensemble by barostatting the system with the Parinello-Rahman method, and using the above expressions, but with additional :math:`R_i`'s for the cell degrees of freedom, and thermostatting the cell velocities as well as the particle velocities :cite:`t-Bussi2009`.

.. bibliography:: references.bib
    :cited:
    :labelprefix: T
    :keyprefix: t-
    :style: unsrt
