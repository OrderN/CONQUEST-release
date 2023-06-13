.. _elec_struc:

====================
Electronic Structure
====================

Introduction to what we can do with electronic structure; mention
post-processing.  Band structure.

Go to :ref:`top <basissets>`.

.. _band_struc:

Band structure
--------------



Go to :ref:`top <basissets>`.

.. _elec_pol:

Electronic Polarisation
-----------------------

The electronic polarisation (the response of a material to an
external electric field) can be calculated using the approach
of Resta :cite:`es-Resta:1992aa` by setting the tag ``General.CalcPol T``.
The direction in which polarisation is found is set using the tag
``General.PolDir`` (choosing 1-3 gives x, y or z, respectively, while
choosing 0 gives all three directions, though this is normally not
recommended).

The Resta approach is a version of the modern theory of polarisation (MTP)
(perhaps better known in the method of King-Smith and Vanderbilt :cite:`es-KingSmith:1993aa`)
where the polarisation is found as:

.. math::
   \mathbf{P} = -\frac{e\mathrm{L}}{\pi V}\mathrm{Im}\mathrm{ln}\mathrm{det}\mathbf{S}\\
   \mathrm{S}_{mn} = \langle \psi_{m} \vert \exp{i2\pi \mathbf{r}}/L\vert\psi_{n} \rangle

where :math:`\mathrm{L}` is a simulation cell length along an appropriate direction
and :math:`V` is the simulation cell volume.  This approach is only valid in the large
simulation cell limit, with Gamma point sampling.

As with all calculations in the MTP
the only valid physical quantity is a *change* of polarisation between two configurations.
A very common quantity to calculate is the Born effective charge (BEC), which is defined
as :math:`Z^{*}_{\alpha\beta}V\partial P_{\alpha}/\partial u_\alpha` and is most easily
calculated by finding the change in polarisation as one atom (or one sublattice) is
moved a small amount.

Go to :ref:`top <basissets>`.

.. bibliography:: references.bib
    :cited:
    :labelprefix: ES
    :keyprefix: es-
    :style: unsrt

Go to :ref:`top <basissets>`.
