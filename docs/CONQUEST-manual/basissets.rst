.. _basissets:

==========
Basis sets
==========

Basis functions are required to construct :ref:`gs_suppfunc`. 

Two kinds of real-space basis functions are supported in CONQUEST, pseudo-atomic orbitals ``PAOs`` (default) and b-spline function ``blips``, which to use is specified by ``Basis.BasisSet``

::

   Basis.BasisSet PAOs

There are four ways to construct support functions:

* primitive PAOs as support functions (simplest)
* Multi-site support functions from PAOs
* on-site support functions from PAOs
* support functions from blips

For the first case, users don't need to specify any input parameters related to support functions in ``Conquest_input``.

The latter three cases, small number of support functions are constructed by taking linear combinations of the basis functions.
In this case, the number of the support functions is required to be specified by ``Atom.NumberOfSupports``.
The linear-combination coefficients are calculated as exlpained below, or read from the ``SFcoeffmatrix2`` files for PAOs or from the ``blip_coeffs`` files for blips by setting ``Basis.LoadCoeffs T``. 

.. _basis_paos:

Pseudo-atomic orbitals
----------------------

PAOs are the atomic-orbital basis functions found from the pseudo-potentials and consists of radial functions (:math:`\zeta`) multiplied by spherical harmonic functions. 
The minimal PAOs are *single-*:math:`\zeta` *(SZ)* PAOs, in which a radial function is prepared for each spherical harmonic function. The computational cost with SZ PAOs is much lower than larger PAOs, but the accuracy is often insufficient.
In general, the calculation accuracy is improved by using *multiple-*:math:`\zeta` PAOs in which several radial functions are used for each spherical harmonic function, although the systematic improvement is not guaranteed (see also :ref:`basis_blips`). Adding polarization functions is also important for accurate descriptions of electron polarization around the atoms in molecules and solids.

When primitive PAOs are used as the support functions without any modifications, the parameters related to the support functions are automatically set based on the information of the PAOs in ``.ion`` files.


.. _basis_mssf:

Multi-site support functions
----------------------------
When ``Basis.MultisiteSF T`` is set, the *multi-site support functions* are constructed.

::

   Basis.BasisSet PAOs
   Basis.MultisiteSF T

The multi-site support functions are constructed by taking linear combinations of the multiple-:math:`\zeta` PAOs on each atom and its neighbouring atoms in the multi-site range ``Atom.MultisiteRange``.
The number of the support functions is required to be specified by ``Atom.NumberOfSupports``, which should be equal or larger than the number of SZ for the atoms. Setting ``Multisite.nonminimal T`` is required when ``Atom.NumberOfSupports`` is larger than SZ.

There are two methods to determine the linear-combination coefficients, 
the *local filter diagonalization (LFD)* method ``Multisite.LFD`` and the *numerical optimisation* ``minE.VaryBasis``.

LFD
+++++

The coefficients :math:`C` are determined by projecting subspace occupied molecular orbitals :math:`C_{sub}` around each atom to the localized trial vectors :math:`t`,

:math:`C = C_{sub} f(\varepsilon_{sub}) C_{sub}^T S_{sub} t`


The LFD subspace region is determined for each atom with the range ``Atom.LFDRange``, which is required to be equal or larger than ``Atom.MultisiteRange``.
(Currently the largest ``Atom.MultisiteRange`` and ``Atom.LFDRange`` among the atoms are used for every atom.)

::

   Basis.BasisSet PAOs
   Basis.MultisiteSF T
   Multisite.LFD T

   # example of Si
   %block Si
   Atom.NumberOfSupports 4
   Atom.MultisiteRange 8.0
   Atom.LFDRange 8.0
   %endblock


The Fermi function :math:`f` with :math:`\varepsilon_{sub}` ``Multisite.LFD.ChemP`` and :math:`kT` ``Multisite.LFD.kT`` in the equation removes the effects of the subspace molecular orbitals in higher energy region.
In defult, :math:`\varepsilon_{sub}` is automatically set to the mean value of the subspace HOMO and LUMO energies for each subspace. If users want to modify this, set ``Multisite.LFD.UseChemPsub F`` and the :math:`\varepsilon_{sub}` value with ``Multisite.LFD.ChemP``.

For the LFD trial functions :math:`t`, when ``Atom.NumberOfSupports`` is equal to the number of SZ or single-zeta plus polarization (SZP), the PAOs which have the widest radial functions for each spherical harmonic function are chosen as the trial vectors automatically in default.
When ``Atom.NumberOfSupports`` is equal to the number of SZP and ``Multisite.nonminimal.offset`` is set, the other PAOs will have the weight in the trial vectors with the value of ``Multisite.nonminimal.offset``.
The users can also provide the trial vectors from the input file using the ``LFDTrialVector`` block

::

   # Trial vectors of Au (element 1) and O (element 2) atoms.
   # Au: 15 PAOs (DZP) -> 6 support functions, O: 13 PAOs (DZP) -> 4 support functions.
   %block LFDTrialVector
   # species sf npao   s   s   x   y   z  d1  d2  d3  d4  d5  d1  d2  d3  d4  d5 for Au
           1  1   15 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
           1  2   15 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0
           1  3   15 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0
           1  4   15 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0
           1  5   15 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0
           1  6   15 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0
           2  1   13 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
           2  2   13 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
           2  3   13 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0
           2  4   13 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0
   # species sf npao   s   s   x   y   z   x   y   z  d1  d2  d3  d4  d5 for O
   %endblock LFDTrialVector

The first, second and third columns correspond to the indices of species, support functions for each species, and the number of PAOs for each species. The other columns provide the initial values of the trial vectors. For example, in the first line in the above example, the second *s* PAO is chosen as the trial vector for the first support function of Au.


We can use the smearing function to avoid the sudden change of the coefficients at ``Atom.MultisiteRange``. The smearing can be turn on by ``Multisite.Smear``. We can set the smearing-function type ``Multisite.Smear.FunctionType`` (default=1:Fermi-Dirac, 2=Error function), the center position of the function ``Multisite.Smear.Center`` (default is equal to the range of the support functions), offset of the center position ``Multisite.Smear.Shift`` and the width of the Fermi-Dirac function ``Multisite.Smear.Width`` (default=0.1).


The coefficients can be updated by providing new electronic density by the SCF calculation. Therefore, two-step procedure, the SCF calculations and the subsequent update of the coefficients, can be repeated by setting ``Multisite.LFD.Minimise`` until the energy and density converge with the threshold of the total DFT energy Multisite.LFD.Min.ThreshE or the density ``Multisite.LFD.Min.ThreshD``. Since the repeating procedure is not variational, the DFT energy might be increased, especially ``Atom.MultisiteRange`` is not large enough. Therefore, the repeating procedure is treated to be converged when the energy increase is smaller than ``Multisite.LFD.Min.ThreshEnergyRise`` (in default, ten times of ``Multisite.LFD.Min.ThreshE``).
``Multisite.LFD.Min.MaxIteration`` is the maximum iteration number of the repeating procedure.

::

   Multisite.LFD T
   Multisite.LFD.Minimise T
   Multisite.LFD.Min.ThreshE 1.0e-6
   Multisite.LFD.Min.ThreshD 1.0e-6
   Multisite.LFD.Min.MaxIteration 150
   Multisite.LFD.Min.ThreshEnergyRise 1.0


Numerical optimisation
++++++++++++++++++++++++

The linear-combination coefficients are optimised by minimizing the DFT energy with respect to the coefficients. The threshold and the maximum iteration number of the numerical optimisation are specified by ``minE.EnergyTolerance`` and ``minE.SupportVariations``. The optimisation is based on the conjugate gradient (CG) method, and the initial CG step size can be specified by ``minE.InitStep_paomin`` (default is 5.0).

::

   minE.VaryBasis T
   minE.EnergyTolerance 1.0e-6
   minE.SupportVariations 30

The numerical optimisation provides more accurate coefficients than the LFD method but usually more time consuming. Therefore, it is recommended to start from good initial values, for example, the coefficients calculated by LFD. When both ``Multisite.LFD`` (with ``Multisite.LFD.Minimise``) and ``minE.VaryBasis`` are turn on, first the coefficients are calculated by the LFD method (with the LFD repeating procedure) and then optimised numerically. 

::

   Basis.MultisiteSF T
   Multisite.LFD T
   Multisite.LFD.Minimise T
   minE.VaryBasis T

If the users already have some good initial coefficient values as the ``SFcoeffmatrix2`` files, reading the files and performing only the numerical optimisation is also a good choice.

::

   Basis.LoadCoeffs T
   Basis.MultisiteSF T
   Multisite.LFD F
   minE.VaryBasis T





.. _basis_ossf:

On-site support functions
-------------------------

On-site support functions are the linear combinations of the PAOs only on the target atom.
In this case, ``Atom.MultisiteRange`` should be small enough not to include any neighboring atoms.

The coefficient can be determined by the LFD method or the numerical optimisation above. Since the range of on-site support function is small, it is strongly recommended to perform the numerical optimisation subsequently to the LFD calculation to guarantee accuracy. ``Atom.LFDRange`` can contain neighbor atoms to improve the accuracy.

The minimum size of the on-site support functions is SZP, so ``Multisite.nonminimal`` is required to be set to T.

Here, the minimum size of on-site support function is larger than that of multi-site support functions (SZ size), but the order-N calculation is more stable with on-site support functions than with multi-site support functions.

::

   Basis.BasisSet PAOs
   Basis.MultisiteSF T
   Multisite.LFD T
   Multisite.nonminimal T

   minE.VaryBasis T

   # example of Si
   %block Si
   Atom.NumberOfSupports 9
   Atom.MultisiteRange 0.1
   Atom.LFDRange 8.0
   %endblock


.. _basis_blips:

Blips
-----

Blips are useful for very accurate calculations, since the basis set can be systematically improved like planewaves.
However, the calculations are sometimes expensive depending on your parameters, and we are now improving the code for blips. 
Thus, please keep it in mind that the explanations or the keywords in the followings might change in the future.

In the case of blips, each atom has a blip grid, 3D regular grid along :math:`x`, :math:`y`, and :math:`z`,
with the atomic position as its origin. 
The blip grid moves rigidly with the atom, and thus we have a pulay force, as in the PAO case.
With blips, we can systematically improve the basis set, by increasing the support function radius 
and/or reducing the spacing of the blip grids. 
For each species of the atom, we need to provide these two parameters, as well as the number of support functions.
The number of support functions can be the size of a minimal basis set, like multi-site support functions.
(At present, minimum value of blip-grid spacing is used for all species.)

::

	%block **
	Atom.NumberOfSupports                        4
	Atom.SupportFunctionRange                  6.0
	Atom.SupportGridSpacing                    0.3
	%endblock **

You need `*`.ion files of SZ basis sets (``minimal`` in the Basis Generation). 

Blip-grid spacing can be determined from the cutoff energy of pseudo wavefunctions in the planewave calculations.
If you need the cutoff energy :math:`E_{\rm cutoff}` in Hartree, the blip-grid spacing should be 
:math:`\frac{2\pi}{\sqrt{2 E_{\rm cutoff}}}` in bohr.
Note that the grid spacing of integration grids (or FFT grids for the charge density) should be smaller than
the half of the blip grids.

It is essential to optimise the support functions (blip coefficients) in the case of blips, and
you have to set the following keyword.

::

	minE.VaryBasis              T  

You may need to reduce the tolerance and/or increase the number of iterations, to optimise the support functions more.
::

	minE.EnergyTolerance             0.10E-07
	minE.SupportVariations             30 

It is not recommended, but if you would encounter a memory problem for very accurate blip calculations, 
you may need to switch off the preconditinoning procedure for length-scale ill conditioning.

::

	minE.PreconditionBlips              F 


.. _basis_bsse:

Basis Set Superposition Error
----------------------------
The basis set superposition error (BSSE) arises when the two monomer units come closer and the basis set localized on one unit can act as 
diffuse functions for the electrons from the other unit, and therefore could be responsible for the overestimation of the binding energy for the interacting systems. 

To correct this BSSE, the Counterpoise (CP) correction method proposed by "Boys and Bernardi" is used where the artificial 
stabilization is controlled by enabling the atoms to improve their basis sets after borrowing functions of an empty basis set from the ghost atoms. 

Since, the typical interaction energy between two monomers A and B is calculated as:

:math:`E_{AB}^{int} = E_{AB}(AB) - E_A(A) - E_B(B).`

Now, the estimate to the amount of the artificial stabilization of A by the extra basis functions from B is:

:math:`E_{A}^{BSSE} = E_A(AB) - E_A(A),`

where energy of A in its monomer basis is subtracted from energy of A in dimer basis. 

Similarly, for monomer B,

:math:`E_{B}^{BSSE} = E_B(AB) - E_B(B),`

Subtracting the BSSE part of A and B units from the typical interaction energy mentioned above, the counterpoise corrected 
interaction energy without BSSE :math:`(E_{AB}^{CP})` will be:

:math:`E_{AB}^{CP} = E_{AB}^{int} - E_{A}^{BSSE} - E_{B}^{BSSE} = E_{AB}(AB) - E_A(AB) - E_B(AB).`


 
 
Practically, to set up such calculation e.g. to evaluate the energy of A in the dimer basis :math:`(E_A(AB))`, the basis 
functions from B is placed on atomic centers of B, however with zero nuclear charge and mass.  This could be performed in CONQUEST by specifying negative sign to the 
corresponding masses for the ghost atoms(B) in the ``block ChemicalSpeciesLabel`` of the input file:

::

 %block ChemicalSpeciesLabel
   1   1.0100000000 A  A.ion
   2  -1.0100000000 B  B.ion
 %endblock

Similarly, separate calculations are performed for monomer B in dimer basis with ghost atoms on A :math:`(E_B(AB))` and 
also for :math:`E_{AB}(AB)` in complete basis, to get the net interaction energy without BSSE. 


