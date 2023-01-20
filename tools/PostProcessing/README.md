Note to developers:

The Conquest files in this directory are essentially all taken from
../../src with a few changes.

1. io_module.f90 has most routines removed (keeping
read_atomic_positions, write_eigenvalues, banner, get_file_name and
get_file_name_2rank)
2. pseudo_tm_info.f90 has lines referencing modules XC and
sfc_partitions and associated variables (gap_threshold and
functional_*) removed
