# CONQUEST test suite.

This directry currently contains two end-to-end tests

  - `test_001_bulk_Si_1proc_Diag`
  - `test_002_bulk_Si_1proc_OrderN`

for Conquest, and a simple python pytest script for checking the correctness of the outputs.

## Usage

To run the tests

  1. Compile Conquest with `make` in `../src`
  2. Run the `Conquest` executable in the subdirectories named `test_00*`
  3. Check the correctness of the outputs with `pytest`

These steps can be run automatically using the script `run_conquest_test.sh` in this directry.

## Contributing

To add new tests

  1. Add input files and a sample output file (run with `IO.Iprint 0` and named `Conquest_out.ref`) in a new subdirectory under [testsuite](./). The naming convention is test directory names start with `test_` followed by a running index with three digits, e.g. `003`.
  2. Add an entry to the `"test_path"` parameter in [test_check_output.py](./test_check_output.py). The test driver will check for all fields in the output listed in the `"key"` parameter. They are read from the [Conquest_out](test_001_bulk_Si_1proc_Diag/Conquest_out.ref) file by the `read_conquest_out()` function.
  3. Add it as a new `Run test XXX` step to the [CI workflow](../.github/workflows/makefile.yml)
  4. *optional* If a new field in the output needs to be checked, these things are required:
     - Add the logic how to read the line of output to the loop `for line in file.readlines()` in `read_conquest_out()`. The result should be stored in the dictionary `Results`.
     - Add the new key from the dictionary `Results` to the `"key"` parameter for `test_check_outputs()`.
     - By default all keys are checked in all tests. If the new key creates errors in the previous tests (e.g. it is not found in their output), add the offending combination(s) to the `xfail_test` logical variable, so that failures in those tests are expected and don't throw errors. See the `"Not-a-real-key"` key as an example.
  5. *optional* By default the fields are checked to 6 decimal precision. If a custom precision is required, add it to the `custom_precision` dictionary with the corresponding key.
