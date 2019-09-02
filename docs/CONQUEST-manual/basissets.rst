.. _basissets:

==========
Basis sets
==========


.. _basis_paos:

Pseudo-atomic orbitals
----------------------

.. _basis_mssf:

Multi-site support functions
----------------------------

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

