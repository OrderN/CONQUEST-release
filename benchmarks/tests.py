import pathlib
import pytest
import os
import sys

# Get the current script's directory
current_dir = os.path.dirname(os.path.abspath(__file__))
# Get the parent directory by going one level up
parent_dir = os.path.dirname(current_dir)
# Add the parent directory to sys.path
sys.path.append(parent_dir)

from conquest_tests import test_conquest_out

class TestClass:
    @pytest.mark.parametrize("test_dir,keys", [
        [
            "test_EXX_isol_C2H4_4proc_PBE0CRI_fullDZP_0.2",
            ['Harris-Foulkes energy','Maximum force','Force Residual','Total stress']
        ],
        [
            "test_EXX_isol_C2H4_4proc_PBE0CRI_fullDZP_0.4",
            ['Harris-Foulkes energy','Maximum force','Force Residual','Total stress']
        ],
        [
            "test_EXX_isol_C2H4_4proc_PBE0CRI_fullDZP_0.4_SCF",
            ['Harris-Foulkes energy','Maximum force','Force Residual','Total stress']
        ],
        [
            "test_EXX_isol_C2H4_4proc_PBE0CRI_fullTZTP_0.6",
            ['Harris-Foulkes energy','Maximum force','Force Residual','Total stress']
        ],
        [
            "test_EXX_isol_C2H4_4proc_PBE0ERI_fullSZP_0.4_SCF",
            ['Harris-Foulkes energy','Maximum force','Force Residual','Total stress']
        ],
        [
            "test_EXX_isol_C2H4_4proc_PBE0ERI_fullSZP_GTO_SCF",
            ['Harris-Foulkes energy','Maximum force','Force Residual','Total stress']
        ],
        [
            "test_EXX_isol_C6H6_4proc_PBE0CRI_fullDZP_0.6",
            ['Harris-Foulkes energy','Maximum force','Force Residual','Total stress']
        ],
        [
            "test_EXX_isol_C6H6_4proc_PBE0CRI_fullSZP_0.6",
            ['Harris-Foulkes energy','Maximum force','Force Residual','Total stress']
        ]
        
    ])
    def test_all(self, test_dir, keys):
        testsuite_directory = pathlib.Path().resolve()
        test_conquest_out(test_dir, keys, testsuite_directory)