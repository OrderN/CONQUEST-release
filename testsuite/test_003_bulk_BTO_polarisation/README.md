# Testing polarisation for BaTiO3

This test finds the electronic ground state for cubic BaTiO3 with one Ba ion displaced a small distance along the z axis, and tests the implementation of the Resta approach to polarisation.  The example provided is completely unrealistic physically (a single five atom unit, run with Gamma point only sampling as required for Resta polarisation), but keeps the runtime relatively short (under 1 minute on an Apple M2 Pro with one process); this will generate a warning (in the file `Conquest_warnings`) but this can be disregarded.  A more realistic though somewhat slower test is provided (see below).

The quantities that are relevant to check are:

 * Harris-Foulkes energy
 * Maximum force
 * Total stress
 * Total polarisation

Note that the last of these is only checked in this test.

## More physically realistic test

The files `Conquest_out.ref_AltRealistic` and `Conquest_input_AltRealistic` form a more physically reasonable test (a 40 atom simulation cell).  If you wish to use this, you should replace `Conquest_out.ref` and `Conquest_input` with these files, though the calculation will take at least four times as long on one process.