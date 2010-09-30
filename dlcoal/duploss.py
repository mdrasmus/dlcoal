
import random

from rasmus import treelib

from compbio import phylo, birthdeath




def sample_dup_times(tree, stree, recon, birth, death,
                     pretime=None, premean=None, events=None):
    """
    Sample duplication times for a gene tree in the dup-loss model

    NOTE: Implied speciation nodes must be present
    """

    def gene2species(gene):
        return recon[tree.nodes[gene]].name

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
            time_span = snode.dist
        
        if recon[tree.root] == stree.root:
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
            start_time = times[node] = stimes[recon[node]]
            if node.parent:
                if times[node] > times[node.parent]:
                    print "bad", node.name
                    #raise Exception("bad time")

            # set duplication times within duplication subtree
            for duproot in node.children:
                if events[duproot] == "dup":
                    snode = recon[duproot]
                    time_span = snode.dist

                    #assert start_time - time_span >= stimes[snode], \
                    #       (duproot.name, start_time, time_span, stimes[snode])
                    sample_dup_times_subtree(times, start_time, time_span,
                                             duproot, 
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
    def walk(dup):
        if dup.parent:
            parent_time = times[dup.parent]
        else:
            parent_time = start_time
        remain = time_span - (start_time - parent_time)

        while True:
            t = birthdeath.sample_birth_wait_time(1, remain, birth, death)
            times[dup] = parent_time - t
            assert t >= 0.0
            
            if times[dup] != parent_time:
                break
            else:
                print t, remain, dup.parent
        
        for child in dup.children:
            if events[child] == "dup":
                walk(child)
    walk(duproot)

    
    

