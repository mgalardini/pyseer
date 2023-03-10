import unittest
import numpy as np
import pandas as pd
from pyseer.cmdscale import cmdscale

precision = 1E-10

class TestCommandScale(unittest.TestCase):
    input_file = 'tests/distances_smaller.tsv.gz'
    Y_file = 'tests/cmdscale.Y.txt.gz'
    e_file = 'tests/cmdscale.e.txt.gz'
    input_data = pd.read_csv(input_file, index_col=0, sep='\t')
    Y = np.loadtxt(Y_file)
    e = np.loadtxt(e_file)
    # ugly hack to take into account minor
    # precision problems between systems
    Y = Y[:, :10]
    e = e[:10]

    def test_cmdscale(self):
        Y, e = cmdscale(self.input_data)
        # ugly hack to take into account minor
        # precision problems between systems
        Y = Y[:, :10]
        e = e[:10]
        self.assertTrue(abs((abs(self.Y) - abs(Y)).max()) < precision)
        self.assertTrue(abs((self.e - e).max()) < precision)


if __name__ == '__main__':
    unittest.main()
