
import random

from rasmus import treelib

from compbio import phylo, birthdeath

import dlcoal


import dlcoal
from dlcoal.ctypes_export import *

#=============================================================================
# export c functions

ex = Exporter(globals())
export = ex.export


if dlcoal.dlcoalc:

    export(dlcoal.dlcoalc, "calcDoomTable", c_int,
           [c_void_p, "tree", c_float, "birth", c_float, "death",
            c_double_p, "doomtable"])
    
    export(dlcoal.dlcoalc, "birthDeathTreePriorFull", c_double,
           [c_void_p, "tree", c_void_p, "stree",
            c_int_p, "recon", c_int_p, "events",
            c_float, "birth", c_float, "death",
            c_double_p, "doomtable"])


def prob_dup_loss(tree, stree, recon, events, duprate, lossrate):
    """Returns the topology prior of a gene tree"""

    if dlcoal.dlcoalc:
        if events is None:
            events = phylo.label_events(tree, recon)

        ptree, nodes, nodelookup = dlcoal.make_ptree(tree)
        pstree, snodes, snodelookup = dlcoal.make_ptree(stree)

        ctree = dlcoal.tree2ctree(tree)
        cstree = dlcoal.tree2ctree(stree)
        recon2 = dlcoal.make_recon_array(tree, recon, nodes, snodelookup)
        events2 = dlcoal.make_events_array(nodes, events)

        doomtable = c_list(c_double, [0] * len(stree.nodes))
        dlcoal.dlcoalc.calcDoomTable(cstree, duprate, lossrate, doomtable)

        p = dlcoal.dlcoalc.birthDeathTreePriorFull(ctree, cstree,
                                    c_list(c_int, recon2), 
                                    c_list(c_int, events2),
                                    duprate, lossrate, doomtable)
        dlcoal.dlcoalc.deleteTree(ctree)
        dlcoal.dlcoalc.deleteTree(cstree)

        return p

    else:
        if "dlcoal_python_fallback" not in globals():
            print >>sys.stderr, "warning: using python code instead of native"
            globals()["dlcoal_python_fallback"] = 1
            # spidir libs
            import spidir
            from spidir import topology_prior
            
        return topology_prior.dup_loss_topology_prior(
            tree, stree, recon, duprate, lossrate,
            events=events)


def sample_dup_times(tree, stree, recon, birth, death,
                     pretime=None, premean=None, events=None):
    """
    Sample duplication times for a gene tree in the dup-loss model
    """

    if events is None:
        events = phylo.label_events(tree, recon)

    # get species tree timestamps
    stimes = treelib.get_tree_timestamps(stree)
    #treelib.check_timestamps(stree, stimes)

    # init timestamps for gene tree
    times = {}


    # set pretimes
    if events[tree.root] != "spec":
        if recon[tree.root] != stree.root:
            # tree root is a dup within species tree
            snode = recon[tree.root]
            start_time = stimes[snode.parent]
            time_span = start_time - stimes[snode]
        else:
            # tree root is a pre-spec dup
            if pretime is None:
                if premean is None:
                    raise Exception("must set pre-mean")

                pretime = 0.0
                while pretime == 0.0:
                    pretime = random.expovariate(1/premean)
            start_time = stimes[stree.root] + pretime
            time_span = pretime

        sample_dup_times_subtree(times, start_time, time_span, tree.root, 
                                 recon, events,
                                 stree, birth, death)

    # set times
    for node in tree.preorder():
        if events[node] == "spec":
            # set speciation time
            times[node] = stimes[recon[node]]


        elif (events[node] == "dup" and
              node.parent is not None and
              recon[node] != recon[node.parent]):
            # set duplication times within duplication subtree
            # node is duproot
            snode = recon[node]
            start_time = stimes[snode.parent]
            time_span = start_time - stimes[snode]
            sample_dup_times_subtree(times, start_time, time_span,
                                     node, 
                                     recon, events,
                                     stree, birth, death)
        elif events[node] == "gene":
            times[node] = 0.0

    return times


def sample_dup_times_subtree(times, start_time, time_span, duproot, 
                             recon, events,
                             stree, birth, death):
    """
    Sample duplication times for only a subtree
    """
    
    # walk duplication subtree
    def walk(dup, parent_time):
        remain = time_span - (start_time - parent_time)

        #while True:
        t = birthdeath.sample_birth_wait_time(1, remain, birth, death)
        times[dup] = parent_time - t
        
        snode = recon[dup]
        for child in dup.children:
            if events[child] == "dup" and recon[child] == snode:
                walk(child, times[dup])
    walk(duproot, start_time)

