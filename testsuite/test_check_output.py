import numpy as np
import os
import pytest

def read_conquest_out(path=".", filename="Conquest_out"):
    '''
    Reads a predefined set of values from a Conquest output file (typically
    named Conquest_out) by matching a string on the line. Currently we read
    Harris-Foulkes energy, Maximum force, Force Residual and Total stress.
    Returns a dictionary with float values.
    '''
    path_to_file = os.path.join(path,filename)
    assert os.path.exists(path_to_file)
    file = open(path_to_file, 'r')
    Results = dict()
    for line in file.readlines():
        if line.find("Harris-Foulkes energy") >= 0:
            Results['Harris-Foulkes energy'] = float(line.split("=")[1].split()[0])
        if line.find("Maximum force") >= 0:
            Results['Max force'] = float(line.split(":")[2].split("(")[0])
        if line.find("Force Residual") >= 0:
            Results['Force residual'] = float(line.split(":")[2].split()[0])
        if line.find("Total stress") >= 0:
            Results['Total stress'] = np.array(line.split(":")[2].split()[:-1], dtype=float)

    return Results

@pytest.mark.parametrize("test_path", ["test_001_bulk_Si_1proc_Diag",
                                       "test_002_bulk_Si_1proc_OrderN"])
@pytest.mark.parametrize("field_name",['Harris-Foulkes energy',
                                       'Max force',
                                       'Force residual',
                                       'Total stress'])
def test_check_outputs(test_path, field_name):
    '''
    Reads a predefined set of results written in Conquest_out files of
    tests 001 and 002 and compares them against results in Conquest_out.ref
    with a tolerance of 4 decimals.
    '''
    ref_result = read_conquest_out(test_path, "Conquest_out.ref")
    test_result = read_conquest_out(test_path, "Conquest_out")

    np.testing.assert_almost_equal(ref_result[field_name],
                                   test_result[field_name],
                                   decimal = 4,
                                   err_msg = test_path+": "+field_name,
                                   verbose = True)
