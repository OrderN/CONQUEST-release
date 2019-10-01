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
The linear-combination coefficients are calculated as exlpained below, or read from the ``SFcoeffmatrix2`` files by setting ``Basis.LoadCoeffs T``.

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
the *local filter diagonalization (LFD)* method ``Multisite.LFD`` and the *numerical optimization* ``minE.VaryBasis``.

LFD
+++++

The coefficients :math:`C` are determined by projecting subspace occupied molecular orbitals :math:`C_{sub}` around each atom to the localized trial vectors :math:`t`,

:math:`C = C_{sub} f(\varepsilon_{sub}) C_{sub}^T S_{sub} t`


The LFD subspace region is determined for each atom with the range ``Atom.LFDRange``, which is required to be equal or larger than ``Atom.MultisiteRange``.
(Currently the largest ``Atom.MultisiteRange`` and ``Atom.LFDRange`` among the atoms are used for every atom.)


The Fermi function :math:`f` with :math:`\varepsilon_{sub}` ``Multisite.LFD.ChemP`` and :math:`kT` ``Multisite.LFD.kT`` in the equation removes the effects of the subspace molecular orbitals in higher energy region.
In defult, :math:`\varepsilon_{sub}` is automatically set to the mean value of the subspace HOMO and LUMO energies for each subspace. If users want to modify this, set ``Multisite.LFD.UseChemPsub F`` and the :math:`\varepsilon_{sub}` value with ``Multisite.LFD.ChemP``.

For the LFD trial functions :math:`t`, when ``Atom.NumberOfSupports`` is equal to the number of SZ or single-zeta plus polarization (SZP), the PAOs which have the widest radial functions for each spherical harmonic function are chosen as the trial vectors automatically in default.
When ``Atom.NumberOfSupports`` is equal to the number of SZP and ``Multisite.nonminimal.offset`` is set, the other PAOs will have the weight in the trial vectors with the value of ``Multisite.nonminimal.offset``.
The users can also provide the trial vectors from the input file as follows:
LFDTrialVector (block)

The coefficients can be updated by providing new electronic density by the SCF calculation. Therefore, two-step procedure, the SCF calculations and the subsequent update of the coefficients, can be repeated by setting ``Multisite.LFD.Minimise`` until the energy and density converge with the threshold of the total DFT energy Multisite.LFD.Min.ThreshE or the density ``Multisite.LFD.Min.ThreshD``. Since the repeating procedure is not variational, the DFT energy might be increased, especially is not large enough. Therefore, the repeating procedure is treated to be converged when the energy increase is smaller than ``Multisite.LFD.Min.ThreshEnergyRise`` (in default, ten times of ``Multisite.LFD.Min.ThreshE``).
``Multisite.LFD.Min.MaxIteration`` is the maximum iteration number of the repeating procedure.

::

   # example of Si
   %block Si
   Atom.NumberOfSupports                        4
   Atom.SupportFunctionRange                  6.0 **(no longer needed?)**
   Atom.MultisiteRange 8.0
   Atom.LFDRange 8.0
   %endblock

Numerical optimization
++++++++++++++++++++++++

The linear-combination coefficients are optimized by minimizing the DFT energy with respect to the coefficients. This method provides more accurate coefficients than the LFD method but usually more time consuming. Therefore, it is recommended to start from good initial values. When both ``Multisite.LFD`` (with ``Multisite.LFD.Minimise``) and ``minE.VaryBasis`` are turn on, first the coefficients are calculated by the LFD method (with the LFD repeating procedure) and then optimized numerically. 
If the users already have some good initial coefficient values as the ``SFcoeffmatrix2`` files, reading the files and performing only the numerical optimization is also a good choice.

::

   Multisite.LFD.ReadTVEC T
   Multisite.LFD F
   minE.VaryBasis T

The threshold and the maximum iteration number of the numerical optimization are specified by ``minE.EnergyTolerance`` and ``minE.SupportVariations``. The optimization is based on the conjugate gradient (CG) method, and the initial CG step size can be specified by ``minE.InitStep_paomin`` (default is 5.0).





.. _basis_ossf:

On-site support functions
-------------------------

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

You need *.ion files of SZ basis sets (``minimal`` in the Basis Generation). 

Blip-grid spacing can be determined from the cutoff energy of pseudo wavefunctions in the planewave calculations.
If you need the cutoff energy :math:`E_{\rm cutoff}` in Hartree, the blip-grid spacing should be 
:math:`\frac{2\pi}{\sqrt{2 E_{\rm cutoff}}}` in bohr.
Note that the grid spacing of integration grids (or FFT grids for the charge density) should be smaller than
the half of the blip grids.

It is essential to optimise the support functions (blip coefficients) in the case of blips, and
you have to set the following keyword.

::

	minE.VaryBasis              T  

You may need to reduce the tolerance and/or increase the number of iterations, to optimize the support functions more.
::

	minE.EnergyTolerance             0.10E-07
	minE.SupportVariations             30 

It is not recommended, but if you would encounter a memory problem for very accurate blip calculations, 
you may need to switch off the preconditinoning procedure for length-scale ill conditioning.

::

	minE.PreconditionBlips              F 

