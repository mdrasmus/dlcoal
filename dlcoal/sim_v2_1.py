#!/usr/bin/env python

## second attempt at removing the strong fixation assumption from the D/L sim
##  should work on a species tree
##  makes numerous simplifying assumptions, explained in the doc of sim_tree
##
## note that this is a non-mutating version (c.f. sim_v2.0.py)


from math import *
import random
from rasmus import stats, treelib
from compbio import coal

DEBUGPRINT = True

def debugprint(s):
    if DEBUGPRINT:
        print s


def sim_tree(stree, popsize, freq, dr, lr, freqdup, freqloss, forcetime):
    """
    Runs a relaxed fixation assumption simulation on a species tree.
    Some simplifying assumptions are made for this version of the simulator:
      1) All branches of the species tree have the same population size
      2) All branches of the species tree have the same duplication rate
      3) All branches of the species tree have the same loss rate
      4) All branches of the species tree have the same duplication effect
      5) All branches of the species tree have the same loss effect
      6) All branches of the species tree have the same time between forced
           frequency changes
      7) There is a single allele at the root of the species tree.
    A duplication/loss effect is the change in frequency for either event.
      Appropriate default values for these effects may need to be determined.
    Furture iterations should remove these assumptions by incorporating
      dictionaries to allow values for each branch.
    
    stree is the initial species tree; it may be mutated by the simulator
    popsize is the population size (assmpt. 1)
    freq is the allele frequency (assmpt. 7)
    dr is the duplication rate (in events/myr/indiv(?); assmpt. 2)
    lr is the loss rate (in events/myr/indiv(?); assmpt. 3)
    freqdup is the duplication effect (assmpt. 4)
    freqloss is the loss effect (assmpt. 5)
    forcetime is the maximum time between frequency changes (assmpt. 6)
    """
    
    ## sanity checks before running the simulator; may be removed or relaxed
    treelib.assert_tree(stree)
    assert popsize > 0
    assert 0.0 <= freq and freq <= 1.0
    assert dr >= 0.0
    assert lr >= 0.0
    assert 0.0 <= freqdup and freqdup <= 1.0
    assert 0.0 <= freqloss and freqloss <= 1.0
    assert forcetime >= 0.0
    
    if dr + lr <= 0.0:
        return stree.copy() # no duplications or losses => stree is final tree
    
    # note: the use of < instead of <= is intentional
    #  if lr==0, duprate/fullrate==1, and random() returns from [0.0,1.0)
    def event_is_dup(duprate, fullrate):
        return random.random() < duprate / fullrate
    
    
    def sim_walk(gtree, snode, gnode, p, s_walk_time=0.0, g_walk_time=0.0, \
                    time_until_force=forcetime):
#        debugprint(" Sim on branch" + str(node.name) + " with frequency " + str(p) + " and walk time " + str(walk_time))
#        debugprint(" walking on " + str(gnode.name))
        if p <= 0.0:
            # gnode is 'parent' of extinct node
            #  create new_gnode, set data['freq'] = 0.0
            #  prune at the end
            new_gnode = treelib.TreeNode(gtree.new_name())
            new_gnode.dist = g_walk_time
            new_gnode.data['freq'] = 0.0
            gtree.add_child(gnode, new_gnode)
#            debugprint("   extinction on " + str(gnode.name))
        else: # put everything else in this block to avoid using returns
            p = min(p, 1.0) # sanity check
            eff_dr = dr * p # * popsize #??
            eff_lr = lr * p # * popsize #??
            eff_bothr = eff_dr + eff_lr
            event_time = stats.exponentialvariate(eff_bothr)
            remaining_s_dist = snode.dist - s_walk_time
            if event_time >= min(time_until_force, remaining_s_dist):
                # do not process D/L event; determine whether at force or speciation
                if time_until_force < remaining_s_dist:
                    # force new frequency
                    newp = coal.sample_freq_CDF(p, popsize, forcetime * 1e6)
                      # scale forcetime to years (in myr)
    #                debugprint("   Forced new frequency: " + str(newp))
                    ## TODO: may wish to log newp in node.data
                    new_s_walk_time = s_walk_time + time_until_force
                    new_g_walk_time = g_walk_time + time_until_force
                    sim_walk(gtree, snode, gnode, newp, \
                                s_walk_time=new_s_walk_time, \
                                g_walk_time=new_g_walk_time)
                      # continue walk with new frequency
                      # increase walk_times accordingly
                      # reset time_until_force to forcetime
                else:
                    # speciation event
                    newp = coal.sample_freq_CDF(p, popsize, remaining_s_dist * 1e6)
                      # scale remaining time into years (from myr)
                    new_gnode = treelib.TreeNode(gtree.new_name())
                    new_gnode.dist = g_walk_time + remaining_s_dist
                    new_gnode.data['freq'] = newp
                      # stores frequency of allele at the speciation event
                    gtree.add_child(gnode, new_gnode)
    #                debugprint("   Completed branch; new frequency: " + str(newp))
                    for schild in snode.children:
                        sim_walk(gtree, schild, new_gnode, newp)
    #                return # shouldn't be necessary
            else:
                # process D/L event
                # no WF updates for these events (modelling decision)
                new_s_walk_time = s_walk_time + event_time
                new_g_walk_time = g_walk_time + event_time
                new_time_until_force = time_until_force - event_time
                if event_is_dup(eff_dr, eff_bothr):
                    # perform duplication event
                    new_gnode = treelib.TreeNode(gtree.new_name())
                      # create a node new_gnode for the duplication event
                    new_gnode.dist = new_g_walk_time # set dist to dup
                    new_gnode.data['freq'] = p # set frequency at dup event
    #                debugprint("   Duplication occurred at walk time " + str(new_walk_time))
                    gtree.add_child(gnode, new_gnode)
#                    debugprint("  starting on orig of " + str(new_gnode.name))
                    sim_walk(gtree, snode, new_gnode, p, \
                                s_walk_time=new_s_walk_time, \
                                time_until_force = new_time_until_force)
                      # recurse on remainder of original branch
#                    debugprint("  starting on dup of " + str(new_gnode.name))
                    sim_walk(gtree, snode, new_gnode, freqdup, \
                                s_walk_time=new_s_walk_time, \
                                time_until_force = new_time_until_force)
                      # recurse on dup tree with correct starting frequency
    #                return
                else:
                    # perform loss event
                    newp = max(p - freqloss, 0.0) # sanity check
    #                debugprint("   Loss occurred at walk time " + str(new_walk_time) + " yielding new frequency " + str(newp))
                    sim_walk(gtree, snode, gnode, newp, \
                                s_walk_time=new_s_walk_time, \
                                g_walk_time=new_g_walk_time, \
                                time_until_force=new_time_until_force)
    
    
    # main code
    
    # create new gene tree and simulate its evolution
    gtree = treelib.Tree()
    gtree.make_root()
    gtree.root.dist = 0.0
    gtree.root.data['freq'] = freq
    sim_walk(gtree, stree.root, gtree.root, freq) # should mutate gtree
    
    
    # remove dead branches and single children (inside last method)
    # note that the simplifyRoot argument was added to the treelib methods
    #  so that gtree.root.dist is always equal to 0.0 (and this allows the
    #  root to have a single child)
    #  if this behavior is undesired later, we can simply remove the argument
    #  and the root will be collapsed (and have >0 dist)
    extant_leaves = []
    for leaf in gtree.leaves():
        if leaf.data['freq'] > 0.0:
            extant_leaves.append(leaf.name)
    gtree = treelib.subtree_by_leaf_names(gtree, extant_leaves, \
                                            simplifyRoot=False)
    
    return gtree



if __name__ == "__main__":
    stree = treelib.read_tree('simple.stree')
    popsize = 1e4
    freq = 1e0
    dr = 2.1
    lr = 2.0
    freqdup = .07
    freqloss = .05
    forcetime = 1e0
    gtree = sim_tree(stree, popsize, freq, dr, lr, freqdup, freqloss, forcetime)
    treelib.draw_tree(stree, scale=1)
    print 
    for leaf in stree.leaves():
        print leaf.name, treelib.find_dist(stree, stree.root.name, leaf.name)
    print
    treelib.draw_tree(gtree, scale=1)
    print
    
    if len(gtree.nodes) > 1:
        for leaf in gtree.leaves():
            print leaf.name, treelib.find_dist(gtree,gtree.root.name,leaf.name), leaf.data['freq']
    else:
        print 'only the root remains'
    
#    if len(gtree.nodes) > 1:
#        extant_count = 0
#        for leaf in gtree.leaves():
#            if leaf.data['freq'] > 0.0:
#                extant_count += 1
#                print leaf.name, treelib.find_dist(gtree, gtree.root.name, leaf.name), leaf.data['freq']
#        if extant_count == 0:
#            print 'all leaves are dead'
#    else:
#        print 'only the root remains'
#    
#    print gtree.root.dist
