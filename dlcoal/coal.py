
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

    export(dlcoal.dlcoalc, "prob_multicoal_recon_topology", c_double,
           [c_int_p, "ptree", c_int, "nnodes", c_int_p, "recon", 
            c_int_p, "pstree", c_int, "nsnodes", c_double_p, "sdists",
            c_double_p, "popsizes"])


def make_ptree(tree):
    """Make parent tree array from tree"""

    nodes = []
    nodelookup = {}
    ptree = []
    
    def walk(node):
        for child in node.children:
            walk(child)
        nodes.append(node)
    walk(tree.root)
    
    def leafsort(a, b):
        if a.is_leaf():
            if b.isLeaf():
                return 0
            else:
                return -1
        else:
            if b.isLeaf():
                return 1
            else:
                return 0
    
    # bring leaves to front
    nodes.sort(cmp=leafsort)
    nodelookup = {}
    for i, n in enumerate(nodes):
        nodelookup[n] = i
    
    for node in nodes:
        if node == tree.root:
            ptree.append(-1)
        else:
            ptree.append(nodelookup[node.parent])
    
    assert nodes[-1] == tree.root
    
    return ptree, nodes, nodelookup


def make_recon_array(tree, recon, nodes, snodelookup):
    """Make a reconciliation array from recon dict"""
    recon2 = []
    for node in nodes:
        recon2.append(snodelookup[recon[node]])
    return recon2


pmrt = compbio.coal.prob_multicoal_recon_topology
def prob_multicoal_recon_topology(tree, recon, stree, n,
                                  lineages=None, top_stats=None):

    ptree, nodes, nodelookup = make_ptree(tree)
    pstree, snodes, snodelookup = make_ptree(stree)
    recon2 = make_recon_array(tree, recon, nodes, snodelookup)

    popsizes = compbio.coal.init_popsizes(stree, n)
    popsizes2 = [popsizes[snode.name] for snode in snodes]
    sdists = [snode.dist for snode in snodes]
    
    p = dlcoal.dlcoalc.prob_multicoal_recon_topology(
        c_list(c_int, ptree), len(nodes), c_list(c_int, recon2),
        c_list(c_int, pstree), len(snodes), c_list(c_double, sdists),
        c_list(c_double, popsizes2))

    #p2 = pmrt(tree, recon, stree, n, lineages, top_stats)
    #assert abs(p - p2) < .01, (p, p2)
    
    return p

compbio.coal.prob_multicoal_recon_topology = prob_multicoal_recon_topology
