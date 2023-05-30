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
