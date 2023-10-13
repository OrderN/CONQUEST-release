# Testing different matrix multiplication kernels

This system can be used for profiling different matrix multiplication kernels.
Those can be chosen with the `MULT_KERN` variable in `system.make`.

The additional coordinate files `si_XYZ.xtl` can be used to test weak scaling and 
would work well for increasing the number of nodes: `si_222.xtl` is the same as `coords.dat`
and has 64 atoms. This means it would run well on anywhere from 2MPI/4OpenMP to 8MPI/1OpenMP.
With the rest of the `xtl` files, we double the number of atoms each time, and would need
to double the number of processes.