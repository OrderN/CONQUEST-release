import numpy as np
import os
import pathlib

def read_conquest_out(path=".", filename="Conquest_out"):
    '''
    Reads a predefined set of values from a Conquest output file (typically
    named Conquest_out) by matching a string on the line. Currently we read
    Harris-Foulkes energy, Maximum force, Force Residual and Total stress.
    Returns a dictionary with float values.
    '''
    path_to_file = os.path.join(path,filename)
    assert os.path.exists(path_to_file), path_to_file
    file = open(path_to_file, 'r')
    Results = dict()
    for line in file.readlines():
        if line.find("Harris-Foulkes energy") >= 0:
            Results['Harris-Foulkes energy'] = float(line.split("=")[1].split()[0])
        if line.find("Maximum force") >= 0:
            Results['Maximum force'] = float(line.split(":")[2].split("(")[0])
        if line.find("Force Residual") >= 0:
            Results['Force Residual'] = float(line.split(":")[2].split()[0])
        if line.find("Total stress") >= 0:
            Results['Total stress'] = np.array(line.split(":")[2].split()[:-1], dtype=float)
        if line.find("Total polarisation") >= 0:
            Results['Total polarisation'] = float(line.split(":")[1].split()[0])

    return Results

def results(path):
    '''
    Reads a result and its reference, selects one value with key and returns it.
    '''
    ref_result = read_conquest_out(path, "Conquest_out.ref")
    test_result = read_conquest_out(path, "Conquest_out")

    return (ref_result, test_result)

def precision(key='_'):
    '''
    Return the relative tolerance used by tests. By default returns 1e-4, but
    takes a key as an argument and you can match it to return a different precision

    For example:

    if(key == 'Special case'):
        return 999.9
    else:
        return 1e-4
    '''
    if (key == 'Harris-Foulkes energy'):
        return 1e-8
    else:
        return 1e-4
        

def testsuite_directory():
    '''
    Return path to testsuite
    '''
    return pathlib.Path().resolve()

def test_conquest_out(test_dir, keys, testsuite_directory):
    path = os.path.join(testsuite_directory, test_dir)
    ref_results, test_results = results(path)
    for key in keys:
            np.testing.assert_allclose(ref_results[key], test_results[key], rtol = precision(key), verbose = True)