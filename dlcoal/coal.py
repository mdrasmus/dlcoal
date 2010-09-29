
from __future__ import division

from math import *

from rasmus import stats
import compbio.coal
from compbio.coal import *

import dlcoal
from dlcoal.ctypes_export import *


#=============================================================================
# export c functions

ex = Exporter(globals())
export = ex.export


if dlcoal.dlcoalc:

    # replace python function with c
    export(dlcoal.dlcoalc, "prob_coal_counts", c_double,
           [c_int, "a", c_int, "b", c_double, "t", c_double, "n"])
    compbio.coal.prob_coal_counts = prob_coal_counts






