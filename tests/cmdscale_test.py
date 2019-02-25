import unittest
import numpy as np
import pandas as pd
from pyseer.cmdscale import cmdscale


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
        self.assertTrue(abs((self.Y - Y).max()) < 1E-15)
        self.assertTrue(abs((self.e - e).max()) < 1E-15)


if __name__ == '__main__':
    unittest.main()
