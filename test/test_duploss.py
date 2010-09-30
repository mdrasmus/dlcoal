# test dlcoal.duploss

import unittest

import dlcoal

from rasmus.common import *
from rasmus import stats
from rasmus.testing import *

from compbio import coal


class DupLoss (unittest.TestCase):

    def test(self):

        stree = treelib.parse_newick("((A:1000, B:1000):500, (C:700, D:700):800);")
        duprate = .000012
        lossrate = .000011

        


#=============================================================================

show_plots = False
def show_plot():
    if show_plots:
        raw_input()


if __name__ == "__main__":

    if "--" in sys.argv:
        args = sys.argv[sys.argv.index("--")+1:]
        if "plot" in args:
            show_plots = True
        sys.argv = sys.argv[:sys.argv.index("--")]
    
    unittest.main()



