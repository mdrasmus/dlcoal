"""

    Python implementation of topology prior related functions

"""

# python imports
from math import *
import random

# rasmus, compbio, spidir imports
from rasmus import util, stats, treelib
from compbio import birthdeath, phylo
import spidir


def calc_doom_table(tree, birth, death, maxdoom=20, nodelookup=None):
    """Compute a doom table for a species tree 'tree'"""
    
    if nodelookup is None:
        ptree, nodes, nodelookup = spidir.make_ptree(tree)

    doomtable = [0] * len(tree.nodes)
    
    def walk(node):
        if node.is_leaf():
            doomtable[nodelookup[node]] = -util.INF
        else:
            for child in node.children:
                walk(child)
                        
            i = nodelookup[node]
            p = 1.0
            for child in node:
                p *= sum(birthdeath.prob_birth_death1(d, child.dist,
                                                      birth, death) *
                         exp(doomtable[nodelookup[child]]) ** d
                         for d in range(0, maxdoom+1))
            doomtable[i] = util.safelog(p, e, -util.INF)
    walk(tree.root)
    
    return doomtable




def num_redundant_topology(node, gene2species, leaves=None, all_leaves=False):
    """Returns the number of 'redundant' topologies"""

    if leaves is None:
        leaves = node.leaves()
    leaves = set(leaves)
    colors = {}
    nmirrors = [0]

    def walk(node):
        if node in leaves:
            colors[node] = phylo.hash_tree(node, gene2species)
        else:
            # recurse
            for child in node.children:
                walk(child)
            
            childHashes = util.mget(colors, node.children)
            if len(childHashes) > 1 and util.equal(* childHashes):
                nmirrors[0] += 1
            
            childHashes.sort()
            colors[node] = phylo.hash_tree_compose(childHashes)
    walk(node)

    colorsizes = util.hist_dict(util.mget(colors, leaves)).values()

    if all_leaves:
        val = stats.factorial(len(leaves))
    else:
        val = 1
        for s in colorsizes:
            if s > 1:
                val *= stats.factorial(s)
    #print "py val=", val, "nmirrors=", nmirrors[0]
    return val / (2**nmirrors[0])



def get_sub_tree(node, snode, recon, events):
    """Returns the leaves of a duplication subtree"""
    leaves = []

    def walk(node):
        if recon[node] != snode:
            return
        if events[node] != "dup":
            leaves.append(node)
        else:
            for child in node.children:
                walk(child)
    walk(node)

    return leaves
    

def dup_loss_topology_prior(tree, stree, recon, birth, death, maxdoom=20,
                            events=None):
    """
    Returns the log prior of a gene tree topology according to dup-loss model
    """

    def gene2species(gene):
        return recon[tree.nodes[gene]].name

    if events is None:
        events = phylo.label_events(tree, recon)
    leaves = set(tree.leaves())
    phylo.add_implied_spec_nodes(tree, stree, recon, events)
    
    pstree, snodes, snodelookup = spidir.make_ptree(stree)

    # get doomtable
    doomtable = calc_doom_table(stree, birth, death, maxdoom)


    prod = 0.0
    for node in tree:
        if events[node] == "spec":
            for schild in recon[node].children:
                nodes2 = [x for x in node.children if recon[x] == schild]
                if len(nodes2) > 0:
                    node2 = nodes2[0]
                    subleaves = get_sub_tree(node2, schild, recon, events)
                    nhist = birthdeath.num_topology_histories(node2, subleaves)
                    s = len(subleaves)
                    thist = stats.factorial(s) * stats.factorial(s-1) / 2**(s-1)
                    
                    if len(set(subleaves) & leaves) == 0:
                        # internal
                        prod += log(num_redundant_topology(node2, gene2species,
                                                          subleaves, True))
                    else:
                        # leaves
                        prod += log(num_redundant_topology(node2, gene2species,
                                                          subleaves, False))
                    
                else:
                    nhist = 1.0
                    thist = 1.0
                    s = 0

                t = sum(stats.choose(s + i, i) *
                        birthdeath.prob_birth_death1(s + i, schild.dist, birth, death) *
                        exp(doomtable[snodelookup[schild]])**i
                        for i in range(maxdoom+1))
                
                prod += log(nhist) - log(thist) + log(t)

    # correct for renumbering
    nt = num_redundant_topology(tree.root, gene2species)
    prod -= log(nt)
    
    #phylo.removeImpliedSpecNodes(tree, recon, events)
    treelib.remove_single_children(tree)

    return prod




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

    
    


    

