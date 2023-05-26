import numpy as np
import os
import pytest

def read_conquest_out(path=".", filename="Conquest_out"):

    path_to_file = os.path.join(path,filename)
    assert os.path.exists(path_to_file)
    file = open(path_to_file, 'r')
    for line in file.readlines():
        if line.find("Maximum force") >= 0:
            max_force = float(line.split(":")[2].split("(")[0])
        if line.find("Force Residual") >= 0:
            force_residual = float(line.split(":")[2].split()[0])
        if line.find("Total stress") >= 0:
            total_stress = np.array(line.split(":")[2].split()[:-1], dtype=float)

    return (max_force, force_residual, total_stress)        

@pytest.mark.parametrize("test_path", ["test_001_bulk_Si_1proc_Diag",
                                       "test_002_bulk_Si_1proc_OrderN"])
def test_check_outputs(test_path):

    ref_result = read_conquest_out(test_path, "Conquest_out.ref")
    test_result = read_conquest_out(test_path, "Conquest_out")

    np.testing.assert_almost_equal(ref_result[0], test_result[0],
                                   decimal = 6,
                                   err_msg='max force', verbose=True)
    np.testing.assert_almost_equal(ref_result[1], test_result[1],
                                   decimal = 6,
                                   err_msg='force residual', verbose=True)
    np.testing.assert_almost_equal(ref_result[2], test_result[2],
                                   decimal = 6,
                                   err_msg='total stress', verbose=True)
