#!/usr/bin/env python

## second attempt at removing the strong fixation assumption from the D/L sim
##  should work on a species tree
##  makes numerous simplifying assumptions, explained in the doc of sim_tree
##
## note that this is a mutating version
## mutation will be removed in sim_v2.1.py


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
        return stree # no duplications or losses => stree is final tree
    
    # note: the use of < instead of <= is intentional
    #  if lr==0, duprate/fullrate==1, and random() returns from [0.0,1.0)
    def event_is_dup(duprate, fullrate):
        return random.random() < duprate / fullrate
    
    def remove_duds(node):
        if len(node.children) == 0:
            parent = node.parent
            stree.remove(node)
            return remove_duds(parent) if parent else None
#            if parent:
#                remove_duds(parent)
    
    def sim_walk(node, p, walk_time=0.0, time_until_force=forcetime):
        debugprint(" Sim on branch" + str(node.name) + " with frequency " + str(p) + " and walk time " + str(walk_time))
        if p <= 0.0:
            debugprint("  Extinction on branch " + str(node.name))
            parent = node.parent
            stree.remove_tree(node) # extinction event
            remove_duds(parent)
            return
        elif p >= 1.0: # sanity check
            p = 1.0
        eff_dr = dr * p # * popsize #??
        eff_lr = lr * p # * popsize #??
        eff_bothr = eff_dr + eff_lr
        event_time = stats.exponentialvariate(eff_bothr)
        if event_time >= min(time_until_force, node.dist - walk_time): # >= ok?
            # do not process D/L event; determine whether at force or new node
            if time_until_force < node.dist - walk_time:
                # force new frequency
                newp = coal.sample_freq_CDF(p, popsize, forcetime * 1e6)
                  # scale forcetime to years (in myr)
                debugprint("   Forced new frequency: " + str(newp))
                ## TODO: may wish to log newp in node.data
                new_walk_time = walk_time + time_until_force
                return sim_walk(node, newp, walk_time=new_walk_time)
                  # continue walk with new frequency
                  # increase walk_time accordingly
                  # reset time_until_force to forcetime
            else:
                # finish node, determine whether to contine walking on children
                newp = coal.sample_freq_CDF(p, popsize, \
                  (node.dist - walk_time) * 1e6)
                  # scale remaining time into years (from myr)
                node.data['freq'] = newp
                  # stores frequency of allele at the speciation event
                debugprint("   Completed branch; new frequency: " + str(newp))
                return node.recurse(sim_walk, newp)
        else:
            # process D/L event
            # no WF updates for these events (modelling decision)
            new_walk_time = walk_time + event_time
            new_time_until_force = time_until_force - event_time
            if event_is_dup(eff_dr, eff_bothr):
                # perform duplication event
                new_node = treelib.TreeNode(stree.new_name())
                  # create a node new_node for the duplication event
                stree.add_child(node.parent, new_node) # add the dup node
                subtree_copy = treelib.subtree(stree, node)
                  # make a copy of node's subtree (dup tree)
                stree.remove(node) # pull node off of parent node
                stree.add_child(new_node, node) # attach node to dup node
                stree.add_tree(new_node, subtree_copy) # attach dup copy
                new_node.dist = new_walk_time # set dist to dup
                node.dist = node.dist - new_walk_time # correct for dup dist
                subtree_copy.root.dist = node.dist # also correct for dup dist
                new_node.data['freq'] = p # set frequency at dup event
                debugprint("   Duplication occurred at walk time " + str(new_walk_time))
                sim_walk(node, p, time_until_force = new_time_until_force)
                  # recurse on remainder of original branch
                sim_walk(subtree_copy.root, freqdup, time_until_force = new_time_until_force)
                  # recurse on dup tree with correct starting frequency
                return
            else:
                # perform loss event
                newp = p - freqloss
                debugprint("   Loss occurred at walk time " + str(new_walk_time) + " yielding new frequency " + str(newp))
                return sim_walk(node, newp, walk_time=new_walk_time, \
                  time_until_force=new_time_until_force)
    
    def remove_single_child_nodes():
        callagain = False
        for node in stree.preorder():
            if node.parent and len(node.children) == 1:
                parent = node.parent
                child = node.children[0]
                newdist = node.dist + child.dist
                stree.remove(node)
                stree.add_child(parent, child)
                child.dist = newdist
                callagain = True
                break
        if callagain:
            remove_single_child_nodes()
    
    # main code
    sim_walk(stree.root, freq)
    remove_single_child_nodes()
    return stree # poor nomenclature; this will be fixed in v2.1



if __name__ == "__main__":
    stree = treelib.read_tree('simple.stree')
    popsize = 1e4
    freq = 1e0
    dr = 2.1
    lr = 2.0
    freqdup = .05
    freqloss = .05
    forcetime = 1e0
    tree = sim_tree(stree, popsize, freq, dr, lr, freqdup, freqloss, forcetime)
    if tree:
        treelib.draw_tree(tree, scale=1)



### VERSION 1 CODE (for reference)
#
#def sim_branch(N=1e6, T=1e4, S=1e1, dr=.0012, lr=.0011, Fd=None, Fl=None, Fi=[1e0]):
#    """ 
#    Run a Poisson-based D/L simulation on a single branch.
#    N is population size, 
#    T is total time to run the sim, 
#    S is the number of steps in the time period T,
#    dr is the duplication rate (in events/individual/myr), 
#    lr is the loss rate (in events/individual/myr),
#    Fd is the frequency for dups,
#    Fl is the frequency for losses, 
#    Fi is the initial frequency. 
#    """
#    if Fd == None:
#        Fd = 1.0/N
#    if Fl == None:
#        Fl = 1.0/N
#    
#    if type(Fi) != list: # ensure Fi is a list or an int
#        Fi = [Fi] if type(Fi) == int or type(Fi) == float else []
#    
#    if len(Fi) == 0:  # should occur only if the branch goes extinct (and all sub-branches)
#        return None # TODO: may want a different return value
#    
#    if T <= 0.0 or S <= 0.0: # simulation has finished its time
#        return flatten(Fi) # TODO: ensure this is the proper return value
#    
#    for j in xrange(len(Fi)):
#        Fi[j] = min(1.0,max(0.0, Fi[j])) # ensure each frequency falls in [0.0,1.0]
#        Fi[j] = round(N * Fi[j]) / N # kill off unrealistic frequencies (esp. near 0.0)
#    
#    ## kill dead branches
#    j = 0
#    while j < len(Fi):
#        if Fi[j] <= 0.0:
#            Fi = Fi[:j] + Fi[j+1:]
#        else:
#            j += 1
#    if len(Fi) == 0:
#        return None # TODO: may want a different return value
#    
#    debugprint('Time left: ' + str(T))
#    debugprint(' Before step: ' + str(Fi))
#    
#    newFi = sim_step(N,T/S,dr,lr,Fd,Fl,Fi) # run simulation for a time step
#    
#    return sim_branch(N,T-T/S,S-1,dr,lr,Fd,Fl,newFi) # iterate to next step, decreasing total time left


#def sim_step(N, time, dr=.0012, lr=.0011, Fd=1e-5, Fl=1e-5, Fi=[1e0]):
#    """
#    Run a Poisson-based D/L simulation for a single time step.
#    Note that all inputs are assumed to be sanitized (no checks, unlike sim_branch).
#    """
#    for j in xrange(len(Fi)):
#        Fi[j] = sim_step_indiv(N, time, dr, lr, Fd, Fl, Fi[j]) # get frequency after D/L events and WF transforms
#    
#    # TODO: determine if the following steps are needed here...
#    Fi = flatten(Fi)
#    
#    debugprint(" After step: " + str(Fi))
#    
#    j = 0
#    while j < len(Fi):
#        if Fi[j] <= 0.0:
#            Fi = Fi[:j] + Fi[j+1:]
#            debugprint("  Removed a dead branch")
#        else:
#            j += 1
#    
#    return flatten(Fi) # flatten to avoid [foo, [bar, baz]] nesting


#def sim_step_indiv(N, time, dr=.0012, lr=.0011, Fd=1e-5, Fl=1e-5, pi=1e0):
#    """
#    Run a Poisson-based D/L simulation on an individual frequency.
#    All inputs are assumed sanitized
#    """
#    def chooseevent(subrate, fullrate):
#        return random.random() <= subrate / fullrate

#    drpyr = dr / 1e6 # dr in events/myr/indiv; drpyr in events/yr/indiv
#    lrpyr = lr / 1e6 # lr in events/myr/indiv; lrpyr in events/yr/indiv
#    edrpyr = drpyr * pi * N # get population dup rate per year   TODO: may need *N?
#    elrpyr = lrpyr * pi * N # get population loss rate per year  TODO: see note above
#    eventrate = edrpyr + elrpyr # get total event rate per year

#    numlosses = 0

#    clock = time
#    newp = [pi]

#    # simulate Poisson process for D/L events
#    while clock > 0.0:
#        eventtime = stats.exponentialvariate(eventrate) # time to next event
#        if eventtime > clock:
#            clock = 0.0
#            break
#        clock -= eventtime
#        if chooseevent(edrpyr, eventrate): # determine whether event is D or L
#            newp.append(sim_step_indiv(N, clock, dr, lr, Fd, Fl, Fd))
#            debugprint("  Duplication")
#        else:
#            numlosses += 1 # all losses in period grouped together (approx)

#    # poor approx?
#    if numlosses > 0:
#        newp[0] = max(newp[0] - numlosses * Fl, 0.0)
#        debugprint("  Losses: " + str(numlosses))
#    newp[0] = (coal.sample_freq_CDF(newp[0], N, time))

#    return flatten(newp) # needs flattening?





### taken from http://kogs-www.informatik.uni-hamburg.de/~meine/python_tricks
### TODO: rewrite without borrowing code
#def flatten(x):
#    """flatten(sequence) -> list

#    Returns a single, flat list which contains all elements retrieved
#    from the sequence and all recursively contained sub-sequences
#    (iterables).

#    Examples:
#    >>> [1, 2, [3,4], (5,6)]
#    [1, 2, [3, 4], (5, 6)]
#    >>> flatten([[[1,2,3], (42,None)], [4,5], [6], 7, MyVector(8,9,10)])
#    [1, 2, 3, 42, None, 4, 5, 6, 7, 8, 9, 10]
#    """

#    result = []
#    for el in x:
#        #if isinstance(el, (list, tuple)):
#        if hasattr(el, "__iter__") and not isinstance(el, basestring):
#            result.extend(flatten(el))
#        else:
#            result.append(el)
#    return result



#if __name__ == "__main__":
#    N = 1e4
#    T = 1e5
#    S = 1e1
#    dr = .0012
#    lr = .0011
#    Fd = .05 #1.0/N
#    Fl = .05 #1.0/N
#    Fi = .5
#    ans = sim_branch(N=N, T=T, S=S, dr=dr, lr=lr, Fd=Fd, Fl=Fl, Fi=Fi)
#    print ans
