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
            "test_001_bulk_Si_1proc_Diag",
            ['Harris-Foulkes energy','Maximum force','Force Residual','Total stress']
        ],
        [
            "test_002_bulk_Si_1proc_OrderN",
            ['Harris-Foulkes energy','Maximum force','Force Residual','Total stress']
        ],
        [
            "test_003_bulk_BTO_polarisation",
            ['Harris-Foulkes energy','Maximum force','Force Residual','Total polarisation']
        ],
        [
            "test_004_isol_C2H4_4proc_PBE0CRI",
            ["Harris-Foulkes energy"]
        ],
        [
            "test_005_isol_C2H4_4proc_PBE0GTO",
            ["Harris-Foulkes energy"]
        ]
    ])
    def test_all(self, test_dir, keys):
        testsuite_directory = pathlib.Path().resolve()
        test_conquest_out(test_dir, keys, testsuite_directory)