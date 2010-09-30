
import random

from rasmus import treelib

from compbio import phylo, birthdeath



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

        while True:
            t = birthdeath.sample_birth_wait_time(1, remain,
                                                  birth, death)
            times[dup] = parent_time - t
            assert t >= 0.0
            
            if times[dup] != parent_time:
                break
            else:
                print t, remain, dup.parent

        snode = recon[dup]
        for child in dup.children:
            if events[child] == "dup" and recon[child] == snode:
                walk(child, times[dup])
    walk(duproot, start_time)

