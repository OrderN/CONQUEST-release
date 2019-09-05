.. _groundstate:

========================
Finding the ground state
========================

Finding the electronic ground state is the heart of any DFT code.  In
CONQUEST, we need to consider: the density matrix (found using
diagonalisation or linear scaling); self-consistency; and the support
functions (which can be optimised or not).

.. _gs_diag:

Diagonalisation
---------------
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Along with linear scaling code in CONQUEST to find out the density matrix, the exact diagonalization can also be performed. This option diagonalizes the Hamiltonian matrix and minimizes the energy with respect to the density matrix elements. It will scale with :math:`N^3`, but will probably be more efficient for systems up to a few hundred atoms. 

The minimal parameters for the diagonalization can be:

 ::

   DM.SolutionMethod diagon
   Diag.MPMesh T	
   Diag.MPMeshX 2
   Diag.MPMeshY 2
   Diag.MPMeshZ 2
   Grid.GridCutoff 100 

The explanation for the above is explained in the following points which need to look at while doing diagonalization calculations:

1. k-mesh generation

The k-points is a way to discretize the integral of the Hamiltonian over the Brillouin zone for the calculation of the total energy of the system. ``Diag.NumKpts``  provide the number of Bloch wave-vectors (k-points) to sample for diagonalization while the ``block`` specify positions of k-points listing fractional coordinates and weight of all k-points. 

:: 

   Diag.NumKpts 2
   %block Diag.Kpoints
   0.00 0.00 0.00 1.00
   %endblock Diag.Kpoints

Also, as the k-points may be defined either by specifying a list of k-points as described above or by a Monkhorst-Pack (MP) grid in terms of the dimensions of the k-point mesh. An MP grid is thus specified by just three numbers along each axis. 

 ::

  Diag.MPMesh T	
  Diag.MPMeshX 2
  Diag.MPMeshY 2
  Diag.MPMeshZ 2

The origin of the Monkhorst-Pack grid may be offset by a vector from the origin of the Brillouin zone.

  ::

   Diag.MPShiftX 0
   Diag.MPShiftY 0
   Diag.MPShiftZ 0

Please note that if this keyword is present in the input file, the keyword ``Diag.NumKpts`` and the block Kpoints will be ignored.

2. K-points parallelization

:: 

  Diag.ProcRows  1 
  Diag.ProcCols 4 
  Diag.KProcGroups 4
 
If not set, ``Diag.KProcGroups`` defaults to 1 (no k-point parallelisation).
So, if you are running CONQUEST on n number of MPI processes, and you asked for say, G = 4 k-point groups, and :math:`r*c` SCALAPACK processor grid (r for rows and c for columns) for each group, then :math:`n = G*r*c`.
``Diag.ProcRows``  and ``Diag.ProcCols`` are determined automatically by CONQUEST in default.

3. Electronic smearing

For metallic systems, the problem of the discontinuity in the electron occupation function is solved by applying smearing.  The "softness" of the approximating function is controlled by a "smearing temperature" kT (k is the Boltzmann constant). Generally, smaller the smearing temperature, the better the approximation. And, the effect of smearing is to blur the details of the Fermi surface (by imposing an artificial temperature on the electronic system) . If we use a smaller (less dense) k-point grid, we must use more smearing, but the exact amount of smearing to use is not obvious. The default smearing type is Fermi-Dirac occupation ``Diag.SmearingType 0`` in which the electronic smearing corresponds to a physical thermal distribution at temperature T.

 ::

  Diag.SmearingType 0 
  Diag.kT 0.001

Methfessel-Paxton smearing method allows much higher smearing temperatures with minimal effect on the free energy (and hence accuracy) of the system. So, the use of less k-points and the significant advantage over Fermi-Dirac smearing. ``Diag.MPOrder`` is the order of Bessel function approximation to delta-function 

 ::

  Diag.SmearingType 1
  Diag.MPOrder 0


4. Integration Grid:
 
An energy cutoff used in Hartree units is an intuitive parameter which sets the integration grid spacing. 

 ::
 
  Grid.GridCutoff 100 

For the control of tight grid density, the integration grid density is explicitly chosen by specifying the number of grid points in all three dimensions. 

 :: 
  
  Grid.PointsAlongX 32
  Grid.PointsAlongY 32
  Grid.PointsAlongZ 32


It is to take care that the resulting integration grid spacing depends on the cell size so convergence parameters may not be valid for the same material in a different simulation cell. You are advised to test convergence with respect to the integration grid

.. _gs_on:

Linear Scaling
--------------

.. _gs_scf:

Self-consistency
----------------
++++++++++++++++++++++++++++++++++++++++++++++++++++++
The Conquest can run in self-consistent and non self-consistent diagonalization mode. Self-consistency can be set using:

 ::

  minE.SelfConsistent T
  minE.SCTolerance 1E-7

Here, we chose to perform a charge self-consistent calculation, and set the tolerance of the self-consistency (SC) cycle residual to ``1E-7``. This results in an energy convergence of order ``1E-7`` Ha at the end of the cycle.
We make sure that the maximum number of self-consistency cycles is high enough by setting ``SC.MaxIters`` to a large value. 
We can grep "DFT Total Energy" to obtain a summary of total energy convergence during the SC cycle and when charge self-consistency is achieved, "DFT total energy" should equal the "Harris-Foulkes energy" and any discrepancy is a measure of SC cycle energy convergence. 

In non-self-consistent calculation ``minE.SelfConsistent F``, the charge density is constructed as the superposition of atomic densities, self-consistency between density and potential is not sought and the Harris- Foulkes functional is used for the energy. should be fast and run in a few minutes but when high accuracy is needed, self-consistent is recommended. 

Moreover, the Pulay-SCF iterations are effected by mixing parameter in self-consistent calculations.i It is the amount of output charge density which is mixed into new charge density. 

::

  SC.LinearMixingSC     T 
  SC.LinearMixingFactor 0.3

Using large mixing parameter, SCF iterations may be unstable though fast when stable. So, it is recommended to use small mixing parameter when there is instability in the SCF cycle. 

++++++++++++++++++++++++++++++++++++++++++++++++++++++

.. _gs_suppfunc:

Support functions
-----------------
Full details of how the support functions are found and represented
can be found in **add this link**.
