.. _convergence:

=====================
Converging Parameters
=====================


.. _conv_grid:

Integration Grid
----------------
***************************

To obtain meaningful and accurate results, it is important to choose a fine enough integration grid for either the diagonalization or linear scaling calculations. The higher the ``Grid.GridCutoff``, the finer the grid. 

 ::
  
  Grid.GridCutoff E

We can perform a series of computations, say just to calculate the energies of the system ``Atom.MoveRun static`` by varying Grid.GridCutoff from say 100 to 400 Ha to check the convergence. Setting ``Grid.GridCutoff`` which is in Hartree should be checked. 

***************************


.. _conv_bz:

Brillouin Zone
--------------

The K-points mesh also needs the convergence test for every property need to be studied. So, it can again be performed by increasing the mesh-points in each direction until you see no changes in that structural or energetic property, say energy convergence. So, modify: 


 :: 

  Diag.NumKpts n

Or if Monkhorst-Pack mesh then modify nx, ny and nz, say from ``1  1  1`` to ``10  10  10``:

::
 
 Diag.MPMeshX nx
 Diag.MPMeshY ny
 Diag.MPMeshZ nz

to check the convergence. 

***************************
.. _conv_on:

Linear Scaling
--------------
* L range
* L Tolerance
* Inverse S range
* Inverse S tolerance
* Number of iterations (?)

.. _conv_scf:

Self-consistency
----------------

.. _conv_suppfunc:

Support Functions
-----------------
