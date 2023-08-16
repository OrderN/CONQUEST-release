.. _theory_energy_force_stress:

=======================================
Background on energy, forces and stress
=======================================

A number of different ways of formulating the energy exist in Conquest at the moment, involving both self-consistent and non-self-consistent densities and potentials, both with and without the neutral atom potential.  With self-consistency, all formulations should give the same result, though numerical issues may give small differences; without self-consistency the Harris-Foulkes functional is more accurate.

Self-consistent calculations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We define the energy in Conquest in two ways that are equivalent at the self-consistent ground state.  The Harris-Foulkes energy is given as:

.. math::
   E_{HF} = 2\mathrm{Tr}\left[KH\right] + \Delta E_{Ha} + \Delta E_{XC} + E_{II}

where the first term is the band structure energy, equivalent to the sum over the energies of the occupied states, the second two terms compensate for double counting and the final term gives the ion-ion interaction:

.. math::
   E_{II} = \frac{1}{2}\left( \sum_{ij} \frac{Z_i Z_j}{\mid \mathbf{R}_i - \mathbf{R}_j \mid} \right)

The Hamiltonian is defined as:

.. math::
   \hat{H} = \hat{T} + \hat{V}_{L} + \hat{V}_{NL} + V_{Ha} + V_{XC}

where the operators are the kinetic energy, the local and non-local pseudopotentials, the Hartree potential, defined as :math:`V_{Ha} = \int d\mathbf{r}^\prime n(\mathbf{r}^\prime)/\mid \mathbf{r} - \mathbf{r}^\prime\mid`, and the exchange-correlation potential.  The alternative form, often known as the DFT energy, is:

.. math::
   E_{DFT} = 2\mathrm{Tr}\left[K(T + V_{L} + V_{NL})\right] + E_{Ha} + E_{XC} + E_{II}

with the Hartree energy defined as usual:

.. math::
   E_{Ha} = \frac{1}{2}\int\int d\mathbf{r}d\mathbf{r}^{\prime} \frac{n(\mathbf{r})n(\mathbf{r}^\prime)}{\mid \mathbf{r} - \mathbf{r}^\prime\mid}

along with the exchange-correlation energy:

.. math::
   E_{XC} = \int d\mathbf{r} \epsilon_{XC}\left[n\right] n(\mathbf{r})

For the Harris-Foulkes and DFT energies to be equal, it is easy to see that the double counting correction terms in the Harris-Foulkes formalism must be:

.. math::
   \Delta E_{Ha} = -E_{Ha} = -\frac{1}{2}\int\int d\mathbf{r}d\mathbf{r}^{\prime} \frac{n(\mathbf{r})n(\mathbf{r}^\prime)}{\mid \mathbf{r} - \mathbf{r}^\prime\mid}

and

.. math::
   \Delta E_{XC} = \int d\mathbf{r} \left(\epsilon_{XC}[n] - V_{XC}[n]\right)n(\mathbf{r})

When calculating forces and stress with self-consistency, we generally use the differentials of the DFT energy rather than the Harris-Foulkes energy; this enables us to separate contributions that are calculated in different ways (in particular on those that are calculated on the integration grid from those that are not).

Go to :ref:`top <theory_energy_force_stress>`.

.. _th_efs_nap:

Neutral atom potential
----------------------

In a DFT code using local orbitals as basis functions, the total energy is most conveniently written in terms of the interaction of neutral atoms: this is simply a reformulation of the total energy which, in particular, reduces the ion-ion interaction to a sum over short-range pair-wise interactions.  The charge density of interest is now the *difference* between the total charge density and a superposition of atomic densities, notated as :math:`\delta n(\mathbf{r}) = n(\mathbf{r}) - \sum_i n_i(\mathbf{r})` for atomic densities :math:`n_i(\mathbf{r})`.  We write:

.. math::
   E_{DFT, NA} = 2\mathrm{Tr}\left[K(T + V_{NA} + V_{NL})\right] + E_{\delta Ha} + E_{XC} + E_{SII}

where the second term is defined as:

.. math::
   E_{\delta Ha} = \frac{1}{2}\int\int d\mathbf{r}d\mathbf{r}^{\prime} \frac{\delta n(\mathbf{r})\delta n(\mathbf{r}^\prime)}{\mid \mathbf{r} - \mathbf{r}^\prime\mid}
   
The final term, the screened ion-ion interaction, is short-ranged, and written as:

.. math::
   E_{SII} = \frac{1}{2}\left( \sum_{ij} \frac{Z_i Z_j}{\mid \mathbf{R}_i - \mathbf{R}_j \mid} - \int d\mathbf{r} n_i(\mathbf{r})V_{Ha,j}(\mathbf{r}) \right)

where :math:`V_{Ha,i}(\mathbf{r})` is the Hartree potential from the atomic density :math:`n_i(\mathbf{r})`.  We define the neutral atom potential for an atom as :math:`V_{NA,i}(\mathbf{r}) = V_{L,i}(\mathbf{r}) + V_{Ha,i}(\mathbf{r})`, combining the local potential and the Hartree potential for the atomic density; the overall neutral atom potential is given as the sum over the atomic densities, :math:`V_{NA}(\mathbf{r}) = \sum_i V_{NA,i}(\mathbf{r})`.  If we write the pseudo-atomic density as :math:`n_{PAD}(\mathbf{r}) = \sum_i n_i(\mathbf{r})` then we can also write :math:`V_{NA}(\mathbf{r}) = V_L(\mathbf{r}) + V_{Ha, PAD}(\mathbf{r})`.  

In this case, we can write the Harris-Foulkes energy as:

.. math::
   E_{HF} = 2\mathrm{Tr}\left[KH\right] + \Delta E_{Ha} + \Delta E_{XC} + E_{SII}

with the Hamiltonian defined as:

.. math::
   \hat{H} = \hat{T} + \hat{V}_{NA} + \hat{V}_{NL} + V_{\delta Ha} + V_{XC}

where :math:`V_{\delta Ha}(\mathbf{r}) = \int d\mathbf{r^\prime} \delta n(\mathbf{r^\prime})/\mid \mathbf{r} - \mathbf{r}^\prime\mid`.  Accordingly, the double counting Hartree correction term has to change:

.. math::
   \Delta E_{Ha} = -E_{\delta Ha} - \int d\mathbf{r} \delta n(\mathbf{r})\sum_i V_{Ha,i}(\mathbf{r}).

Go to :ref:`top <theory_energy_force_stress>`.

.. _th_efs_nsc:   

Non-self-consistent calculations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In non-self-consistent calculations, we use the Harris-Foulkes functional, along with a reasonable guess for the input density, which is normally taken as the superposition of atomic densities, :math:`n_{in}(\mathbf{r})` and write:

.. math::
   E_{NSC} = 2\mathrm{Tr}\left[KH\right] + \Delta E_{Ha}\left[n_{in}\right] + \Delta E_{XC}\left[n_{in}\right] + E_{II}

Notice that we effectively have two densities being used here: :math:`n_{in}` (which is normally the superposition of atomic densities used in the neutral atom case) and and *effective* output density, :math:`n_{out} = \sum_{ij} \phi_i K_{ij} \phi_j` which comes from the band energy (first term); this complicates the calculation of forces and stress compared to the self-consistent case, as we have to consider contributions from both densities.

For the :ref:`neutral atom potential <th_efs_nap>`, :math:`\delta n(\mathbf{r}) = 0` by definition, which also means that :math:`E_{\delta Ha} = 0` and :math:`\Delta E_{Ha}=0`.

Go to :ref:`top <theory_energy_force_stress>`.

.. _th_efs_pcc:

Partial core corrections
~~~~~~~~~~~~~~~~~~~~~~~~

Also known as non-linear core corrections, partial core corrections (PCC) :cite:`efs-Louie:1982aa` add a model core charge to the pseudopotential to allow for the non-linear exchange-correlation interation between core and valence charge (which is linearised in standard pseudopotentials); this generally improves the accuracy of the pseudopotential.  The exchange-correlation potential is evaluated in terms of the combined charge density, :math:`n_v(\mathbf{r}) + n_c(\mathbf{r})` where the valence charge is input or output charge density defined above: :math:`V_{XC}\left[ n_v + n_c \right]`.  The exchange-correlation energy becomes:

.. math::
   E_{XC} = \int d\mathbf{r} \left(n_v(\mathbf{r}) + n_c(\mathbf{r})\right) V_{XC}\left[ n_v + n_c \right] .

Once this change to the charge density has been made, there is no change to the DFT energy.  However, the double counting term for Harris-Foulkes needs redefining, since XC contribution to the band energy is :math:`2Tr[KV_{XC}] = \int d\mathbf{r} n_v(\mathbf{r}) V_{XC}[n_v + n_c]`.  We write:

.. math::
   \Delta E_{XC} &=& \int d\mathbf{r} \left(n_v(\mathbf{r}) + n_c(\mathbf{r})\right)\epsilon_{XC}[n_v + n_c] - \int d\mathbf{r} n_v(\mathbf{r})V_{XC}[n_v + n_c]\\
   &=& \int d\mathbf{r} n_c(\mathbf{r})\epsilon_{XC}[n_v + n_c] + \int d\mathbf{r} \left(\epsilon_{XC}[n_v + n_c] - V_{XC}[n_v + n_c]\right)n_v(\mathbf{r})

There is an extra factor of :math:`\int d\mathbf{r} n_c(\mathbf{r})\epsilon_{XC}[n_v + n_c]` over and above the usual term.

Go to :ref:`top <theory_energy_force_stress>`.

.. _th_efs_fs:

Forces and Stresses
~~~~~~~~~~~~~~~~~~~

It is important that the forces and stresses be the exact derivatives of the energy, for consistency.  In particular, this means that as the energy is calculated in different ways for different contributions, the force or stress contribution must be calculated in the same way.

Go to :ref:`top <theory_energy_force_stress>`.

.. _th_efs_for:

Forces
------

Forces are defined as the change in energy with respect to atomic positions; as the basis functions move with the atoms, these changes will also include Pulay terms.  The forces found in Conquest are documented extensively elsewhere :cite:`efs-Miyazaki2004,efs-Torralba:2009nr` though the changes needed to account for :ref:`PCC <th_efs_pcc>`, particularly in the :ref:`non-self-consistent <th_efs_nsc>` case, have not been published and are given here for completeness.  As well as the Hellmann-Feynman forces (which come from the movement of the local and non-local pseudopotentials with the atoms) we define Pulay forces (divided into two parts, known as phi-Pulay which come from changes in the Hamiltonian matrix, and S-Pulay, which come from changes in the overlap matrix; the phi-Pulay forces are calculated in three contributions, which depend on how the respective parts of the Hamiltonian matrix are calculated: the kinetic energy; the non-local pseudopotential; and the remaining terms which are all found on the integration grid).  The ion-ion interactions also contribute forces.

The inclusion of :ref:`PCC <th_efs_pcc>` adds an extra term to the forces in all calculations, which comes from the change of the core density as the atoms move; the force on atom :math:`i` is given as:

.. math::
   \mathbf{F}^{PCC}_i = -\int d\mathbf{r} \nabla_i n^c_i(\mathbf{r}) V_{XC}[n_v + n_c]

If the :ref:`non-self-consistent <th_efs_nsc>` formalism is used, then a further term is added (the non-self-consistent force changes) to include the gradient of the core charge.  The non-self-consistent force is now written as:

.. math::
   \mathbf{F}^{NSC}_i = -\int d\mathbf{r} V_{\delta Ha}(\mathbf{r}) \nabla_i n^v_i(\mathbf{r}) - \int d\mathbf{r} \delta n(\mathbf{r}) V_{XC}^\prime\left[n^{in}_v + n_c\right] \left( \nabla_i n^v_i(\mathbf{r}) + \nabla_i n^c_i(\mathbf{r}) \right)

where :math:`V^\prime_{XC}` is the derivative of the exchange-correlation potential with respect to charge density.

Go to :ref:`top <theory_energy_force_stress>`.

.. _th_efs_str:

Stress
------

The stress includes all contributions to the change of energy with the lattice constants; the calculation of stress in Conquest is documented in a paper being prepared for publication, but we give a brief overview here.  As Conquest uses orthorhombic cells, only the diagonal stress components (:math:`\sigma_{\alpha\alpha}`) are calculated.

In most cases, forces also contribute to the stress; it is easy to show that the stress contribution is given by:

.. math::
   \sigma_{\alpha\alpha} = \sum_i F_{i\alpha}R_{i\alpha}

where :math:`R_{i\alpha}` is the position of the atom.  As well as these contributions, there are more subtle terms.  Any energies calculated on the grid will contribute to the stress as the integration grid changes with cell size (the stress is simply the energy calculated), and the Hartree potential contributes a term related to the change in the reciprocal lattice vectors (as it is calculated by Fourier transforming the charge density).  If the exchange-correlation functional is a GGA functional, then a further term coming from the change of the gradient of the density with the cell size arises.  (For non-self-consistent calculations this leads to some complications, as this term technically requires both input and output densities; at present, we approximate this as a mixture of the term calculated with input density and the term calculated with output density; the proportion can be adjusted using the parameter ``General.MixXCGGAInOut`` documented in the :ref:`Advanced and obscure tags <advanced_general_tags>` section of the manual, though we do not recommend changing it.)

Go to :ref:`top <theory_energy_force_stress>`.

.. bibliography:: references.bib
    :cited:
    :labelprefix: EFS
    :keyprefix: efs-
    :style: unsrt

Go to :ref:`top <theory_energy_force_stress>`.
