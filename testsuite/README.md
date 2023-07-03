# CONQUEST test suite.

This directory currently contains three end-to-end tests

  - `test_001_bulk_Si_1proc_Diag`
  - `test_002_bulk_Si_1proc_OrderN`
  - `test_003_bulk_BTO_polarisation`

for Conquest, and a simple python pytest script for checking the correctness of the outputs.

## Usage

To run the tests

  1. Compile Conquest with `make` in `../src`
  2. Run the `Conquest` executable in the subdirectories named `test_00*`
  3. Check the correctness of the outputs with `pytest`

These steps can be run automatically using the script `run_conquest_test.sh` in this directry.

## Contributing

To add new tests

  1. Add input files and a sample output file (run with `IO.Iprint 0` and named `Conquest_out.ref`) in a new subdirectory under [testsuite](./). The naming convention is test directory names start with `test_` followed by a running index with three digits, e.g. `004`.
  2. Add a test to `TestClass` in [`test_check_output.py`](./test_check_output.py). You can use one of the current tests, named `test_XXX` as templates. Update the directory name passed to `path` and the list of parameters in the `@pytest.mark.parametrize` decorator. The parameters are the fields in `Conquest_out` to be checked against the reference. If necessary, update the `read_conquest_out()` function to parse a new field from the output.
  3. Add it as a new `Run test XXX` step to the [CI workflow](../.github/workflows/makefile.yml)
