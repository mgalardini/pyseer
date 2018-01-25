import unittest
import numpy as np
import pandas as pd
from pyseer.cmdscale import cmdscale


class TestCommandScale(unittest.TestCase):
    input_file = 'tests/distances.tsv.gz'
    Y_file = 'tests/cmdscale.Y.txt.gz'
    e_file = 'tests/cmdscale.e.txt.gz'
    input_data = pd.read_table(input_file, index_col=0)
    Y = np.loadtxt(Y_file)
    e = np.loadtxt(e_file)

    def test_cmdscale(self):
        Y, e = cmdscale(self.input_data)
        self.assertTrue(np.array_equal(Y, self.Y))
        self.assertTrue(np.array_equal(e, self.e))


if __name__ == '__main__':
    unittest.main()
