#!/usr/bin/env python
#
# simulation with hemiplasy
#

## second attempt at removing the strong fixation assumption from the D/L sim
##  should work on a species tree
##  makes numerous simplifying assumptions, explained in the doc of sim_tree
##
## this is a non-mutating version of the simulator
##
## this version includes logging of events in the simulated tree
## future versions may remove the logging events from the simulated gene tree,
##  instead keeping the information separate (may not be necessary)

# python imports
import os, sys
import copy
from math import *
import random

# rasmus, compbio imports
from rasmus import util, stats, treelib
from compbio import coal, phylo

# dlcoal imports
import dlcoal
import dlcoal.sim

#=============================================================================
# debugging
DEBUGPRINT = True

def debugprint(s):
    if DEBUGPRINT:
        print s


#=============================================================================
# simulation 

def dlcoal_sims(outdir, nsims, stree, n, duprate, lossrate,
                start=0,
                freq=1.0, freqdup=.05, freqloss=.05, steptime=None,
                nsteps=100,
                full_log=False,
                **options):

    if steptime is None:
        stimes = treelib.get_tree_timestamps(stree)
        steptime = stimes[stree.root] / float(nsteps)
    
    for i in xrange(start, nsims):
        outfile = phylo.phylofile(outdir, str(i), "")
        util.makedirs(os.path.dirname(outfile))
        print "simulating", outfile

        # sample a new tree from DLCoal model
        coal_tree, ex = sample_dlcoal_hem(
            stree, n, duprate, lossrate,
            freq, freqdup, freqloss, steptime,
            keep_extinct=full_log,
            **options)

        # write datastructures
        dlcoal.write_dlcoal_recon(outfile, coal_tree, ex)
        if full_log:
            full_logfile = phylo.phylofile(outdir, str(i), ".locus.info")
            
            full_locus_tree = ex["full_locus_tree"]
            ex2 = generate_extras(stree, full_locus_tree)
            daughters = ex2["daughters"]

            out = open(full_logfile, "w")
            out.write("hem\t%d\n" % 
                      int(is_locus_tree_hemiplasy(full_locus_tree, daughters)))
            out.close()



def sample_dlcoal_hem(stree, n, duprate, lossrate,
                      freq, freqdup, freqloss, steptime,
                      namefunc=lambda x: x,
                      keep_extinct=False,
                      remove_single=True,
                      name_internal="n", minsize=0):
    """Sample a gene tree from the DLCoal model with hemiplasy"""

    # generate the locus tree
    while True:
        locus_tree, locus_extras = sample_locus_tree_hem(
            stree, n, duprate, lossrate, 
            freq, freqdup, freqloss, steptime,
            keep_extinct=keep_extinct)
        if len(locus_tree.leaves()) >= minsize:
            break

    if len(locus_tree.nodes) <= 1: # TODO: check 1 value
        # total extinction
        coal_tree = treelib.Tree()
        coal_tree.make_root()
        coal_recon = {coal_tree.root: locus_tree.root}
        daughters = set()
    else:
        # simulate coalescence
        
        # create new (expanded) locus tree
        logged_locus_tree, logged_extras = locus_to_logged_tree(
            locus_tree, popsize=n)
        daughters = logged_extras[0]
        pops = logged_extras[1]
        log_recon = logged_extras[2]
        
        #treelib.assert_tree(logged_locus_tree)
        
        # removed locus_tree_copy from below
        coal_tree, coal_recon = dlcoal.sim.sample_multilocus_tree(
            logged_locus_tree, n=pops, daughters=daughters,
            namefunc=lambda lognamex: log_recon[lognamex]+'_'+str(lognamex))

        #print set(coal_tree) - set(coal_tree.postorder())
        treelib.assert_tree(coal_tree)
    
        # clean up coal tree
        if remove_single:
            treelib.remove_single_children(coal_tree)
            phylo.subset_recon(coal_tree, coal_recon)


    if name_internal:
        dlcoal.rename_nodes(coal_tree, name_internal)
        dlcoal.rename_nodes(locus_tree, name_internal)


    # store extra information
    extra = {"locus_tree": locus_tree,
             "locus_recon": locus_extras['recon'],
             "locus_events": locus_extras['events'],
             "coal_tree": coal_tree,
             "coal_recon": coal_recon,
             "daughters": daughters}

    if keep_extinct:
        extra["full_locus_tree"] = locus_extras["full_locus_tree"]

    return coal_tree, extra



def sample_locus_tree_hem(stree, popsize, duprate, lossrate,
                          freq=1.0, freqdup=.05, freqloss=.05,
                          steptime=1e6, keep_extinct=False):
    
    """
    Sample a locus tree with birth-death and hemiplasy
    
    
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

    parameters:
    stree is the initial species tree; it may be mutated by the simulator
    popsize is the population size (assmpt. 1)
    freq is the allele frequency (assmpt. 7)
    duprate is the duplication rate (in events/myr/indiv(?); assmpt. 2)
    lossrate is the loss rate (in events/myr/indiv(?); assmpt. 3)
    freqdup is the duplication effect (assmpt. 4)
    freqloss is the loss effect (assmpt. 5)
    forcetime is the maximum time between frequency changes (assmpt. 6)
    
    Returns the locus tree, as well as extra information
    including a reconciliation dictionary and an events dictionary.
    """
    
    ## sanity checks before running the simulator; may be removed or relaxed
    treelib.assert_tree(stree)
    assert popsize > 0
    assert 0.0 <= freq and freq <= 1.0
    assert duprate >= 0.0
    assert lossrate >= 0.0
    assert 0.0 <= freqdup and freqdup <= 1.0
    assert 0.0 <= freqloss and freqloss <= 1.0
    assert steptime > 0.0

    
    # special case: no duplications or losses
    if duprate == 0.0 and lossrate == 0.0:
        locus_tree = stree.copy()
        recon = phylo.reconcile(locus_tree, stree, lambda x: x)
        events = phylo.label_events(locus_tree, recon)

        return locus_tree, {"recon": recon,
                            "events": events,
                            "daughters": set()}
                                
    
    def event_is_dup(duprate, fullrate):
        return random.random() <= duprate / fullrate

    
    def sim_walk(gtree, snode, gparent, p,
                 s_walk_time=0.0, remaining_steptime=steptime,
                 daughter=False):
        """
        eventlog is a log of events along the gtree branch.
        Each entry has the form
          (time_on_branch, event_type, frequency, species_node),
          
        where
           0.0 <= time_on_branch <= branch_node.dist

        event_type is one of
           {'extinction', 'frequency', 'speciation', duplication',
            'loss', 'root', 'gene'},
            
        where 'root' is a unique event not added during the sim_walk process

        frequency is the branch frequency at the event time

        species_node is the name of the node of the species tree branch in
        which the event occurs
        """

        # create new node
        gnode = treelib.TreeNode(gtree.new_name())
        gtree.add_child(gparent, gnode)
        gnode.data = {"freq": p,
                      "log": []}
        eventlog = gnode.data["log"]
        g_walk_time = 0.0
        if daughter:
            eventlog.append((0.0, 'daughter', freqdup, snode.name))
            
        
        # grow this branch, determine next event
        event = None
        while True:
            if p <= 0.0:
                event = "extinct"
                break
            
            # determine remaing time
            remaining_s_dist = snode.dist - s_walk_time
            remaining_time = min(remaining_steptime, remaining_s_dist)

            # sample next dup/loss event
            eff_duprate = duprate * p / freqdup
            eff_lossrate = lossrate * p / freqloss
            eff_bothrate = eff_duprate + eff_lossrate            
            event_time = stats.exponentialvariate(eff_bothrate)

            # advance times
            time_delta = min(event_time, remaining_time)
            s_walk_time += time_delta
            g_walk_time += time_delta

            # sample new frequency
            p = coal.sample_freq_CDF(p, popsize, time_delta)

            # determine event
            if event_time < remaining_time:
                # dup/loss occurs
                if event_is_dup(eff_duprate, eff_bothrate):
                    # dup, stop growing
                    event = "dup"
                    break
                else:
                    # loss, continue growing
                    event = "loss"
                    
            else:
                if remaining_s_dist < remaining_steptime:
                    # we are at a speciation, stop growing
                    event = "spec"
                    break

            # process step
            if event == "loss":
                # LOSS EVENT
                p = max(p - freqloss, 0.0)
                remaining_steptime -= time_delta
                eventlog.append((g_walk_time, 'loss', p, snode.name))
            else:
                # NEXT TIME STEP
                remaining_steptime = steptime
                eventlog.append((g_walk_time, 'frequency', p, snode.name))
                

        # process event
        if event == "extinct":
            # EXTINCTION EVENT (p <= 0)
            gnode.dist = g_walk_time
            gnode.data['freq'] = 0.0
            eventlog.append((g_walk_time, 'extinction', 0.0, snode.name))

        
        elif event == "spec":
            # SPECIATION EVENT
            gnode.dist = g_walk_time
            gnode.data['freq'] = p
                        
            # add speciation event to event log and
            if snode.is_leaf():
                eventlog.append((g_walk_time, 'gene', p, snode.name))
            else:
                eventlog.append((g_walk_time, 'speciation', p, snode.name))
                for schild in snode.children:
                    sim_walk(gtree, schild, gnode, p)


        elif event == "dup":
            # DUPLICATION EVENT
            gnode.dist = g_walk_time
            gnode.data['freq'] = p
            eventlog.append((g_walk_time, 'duplication', p, snode.name))

            # recurse on mother
            sim_walk(gtree, snode, gnode, p, 
                     s_walk_time=s_walk_time, 
                     remaining_steptime=remaining_steptime)

            # recurse on daughter
            sim_walk(gtree, snode, gnode, freqdup, 
                     s_walk_time=s_walk_time, 
                     remaining_steptime=remaining_steptime,
                     daughter=True)

        else:
            raise Exception("unknown event '%s'" % event)
    
    
    # create new gene tree and simulate its evolution
    gtree = treelib.Tree()
    gtree.make_root()
    gtree.root.dist = 0.0
    gtree.root.data['freq'] = freq
    gtree.root.data['log'] = [(0.0, 'speciation', freq, stree.root.name)]

    # simulate locus tree
    sim_walk(gtree, stree.root.children[0], gtree.root, freq)
    sim_walk(gtree, stree.root.children[1], gtree.root, freq)
    
    
    # remove dead branches and single children
    extant_leaves = [leaf.name for leaf in gtree.leaves()
                     if leaf.data['freq'] > 0.0]
    extinctions = [leaf for leaf in gtree.leaves()
                   if leaf.data['freq'] == 0.0]

    if keep_extinct:
        full_gtree = gtree.copy()
        # do deep copy of data
        for node in full_gtree:
            node2 = gtree.nodes[node.name]
            for key, val in node2.data.items():
                node.data[key] = copy.copy(val)
        
    treelib.subtree_by_leaf_names(gtree, extant_leaves, keep_single=True)
    remove_single_children(gtree)

    # determine extra information (recon, events, daughters)
    extras = generate_extras(stree, gtree)

    if keep_extinct:
        extras["full_locus_tree"] = full_gtree
    
    return gtree, extras


def remove_single_children(tree):
    """
    Modified from remove_single_children in treelib.py.
    Added log manipulation.
    """
    
    # find single children
    removed = [node for node in tree
               if len(node.children) == 1 and node.parent]
    
    # actually remove children
    for node in removed:
        newnode = node.children[0]
        
        # add distance
        newnode.dist += node.dist
        
        # update logs
        newnode.data['log'] = node.data['log'][:] + \
                                map(lambda x: (x[0]+node.dist,)+x[1:], \
                                                newnode.data['log'])
        
        # change parent and child pointers
        newnode.parent = node.parent
        index = node.parent.children.index(node)
        node.parent.children[index] = newnode
        
        # remove old node
        del tree.nodes[node.name]

    # do not remove singleton from root (may be added later if desired)
    return removed


def generate_extras(stree, gtree, genename=lambda sp, x: sp + "_" + str(x)):
    """
    generate extra information from the logging in a gene tree

    based on birthdeath.sample_birth_death_gene_tree
    """
    
    recon = {}
    events = {}
    daughters = set()
    
    def walk(gnode):
        (time, event, freq, snname) = gnode.data['log'][-1]
        recon[gnode] = stree.nodes[snname]
        
        if event == 'speciation':
            events[gnode] = 'spec'
            
        elif event == 'duplication':
            events[gnode] = 'dup'
            
            # determine 'daughter' of duplication event
            for gchild in gnode.children:
                if gchild.data['log'][0][1] == 'daughter':
                    daughters.add(gchild)
                    break
                
        elif event == 'gene':
            events[gnode] = 'gene'
            gtree.rename(gnode.name, genename(snname, gnode.name))

        elif event == 'extinction':
            events[gnode] = 'extinction'
            gtree.rename(gnode.name, genename(str(snname), gnode.name))

        else: 
            assert False
        
        for gchild in gnode.children:
            walk(gchild)    
    walk(gtree.root)
    
    return {'recon': recon, 'events': events, 'daughters': daughters}


def is_init_loss(node):
    """
    Returns True if branch contains an 'initial loss', which is a loss
    event when frequency = 1.
    """

    for time, event, freq, snname in node.data["log"]:
        if event == "loss" and freq == 1.0:
            return True
    return False
        

def iter_poly_subtree(locus_tree, daughters):
    """
    Iterates over polymorphic subtress within locus_tree

    Yields (root, leaves)
    """

    for node in locus_tree.preorder():
        if node in daughters or is_init_loss(node):
            yield node, get_poly_subtree_leaves(node, daughters)

def get_poly_subtree_leaves(root, daughters):
    """
    Returns the leaves of a polymorphic subtree rooted at 'root'
    """

    leaves = []

    def walk(node):
        if node.is_leaf():
            leaves.append(node)
        else:
            for child in node.children:
                if child not in daughters and not is_init_loss(child):
                    walk(child)
    walk(root)

    return leaves
    
def is_poly_subtree_hemiplasy(root, leaves):
    """
    Returns True if the polymorphic subtree shows hemiplasy
    """

    extincts = [node.data["freq"] == 0.0 for node in leaves]
    return not util.equal(* extincts)


def is_locus_tree_hemiplasy(locus_tree, daughters):
    """
    Returns True if locus_tree shows hemiplasy
    """

    for root, leaves in iter_poly_subtree(locus_tree, daughters):
        if is_poly_subtree_hemiplasy(root, leaves):
            return True
    return False
    


def locus_to_logged_tree(ltree, popsize=1.0):
    
    ### note: ltree is the locus tree generated by the simulator
    ###       lmctree is the new logged-locus tree (locus-multi-coal-tree)
    
    
    # must have a valid locus tree (of the form from this simulator) passed in
    assert ltree.root.data['log'][-1][1] == 'speciation'
    
    # initializes new logged locus tree
    lmctree = treelib.Tree()
    lmctree.make_root()
    lmctree.root.dist = ltree.root.dist # == 0
#    lmctree.root.data['freq'] = ltree.root.data['freq'] # == starting freq
    
    # initialize the daughters set
    daughters = set()
    
    # initialize the population pops dictionary
    pops = {lmctree.root.name: popsize * ltree.root.data['freq']}
    
    # initialize the newltree to ltree reverse recon
    revrecon = {lmctree.root.name: ltree.root.name}
    
    # initialize the ltree to newltree recon
    recon = {ltree.root.name: lmctree.root.name}
    
    
    # main walking function
    def walk(lnode, lmcparent):
        
        llog = lnode.data['log']
        
        # add a daughter placeholder node to the lmctree 
        #  also include it in the daughters set
        dlog = llog.pop(0) if llog[0][1] == 'daughter' else []
        if dlog:
            # add new placeholder node
            new_lmcnode = treelib.TreeNode(lmctree.new_name())
            new_lmcnode.dist = dlog[0] # == 0.0
#            new_lmcnode.data['freq'] = dlog[2]
            lmctree.add_child(lmcparent, new_lmcnode)
            # add pops entry
            pops[new_lmcnode.name] = popsize * dlog[2]
            # add placeholder node to daughters set
            daughters.add(new_lmcnode)
            # set new parent in the chain
            lmcparent = new_lmcnode
        
        prevdist = 0.0
        
        for i in xrange(len(llog)-1): # perform all actions except for the last node
            # create and add new node
            new_lmcnode = treelib.TreeNode(lmctree.new_name())
            new_lmcnode.dist = llog[i][0] - prevdist
#            new_lmcnode.data['freq'] = llog[i][2]
            lmctree.add_child(lmcparent, new_lmcnode)
            # add pops entry
            pops[new_lmcnode.name] = popsize * llog[i][2]
            # update loop parameters
            prevdist = llog[i][0]
            lmcparent = new_lmcnode
        
        # create and add node for last event (spec/dup)
        new_lmcnode = treelib.TreeNode(lmctree.new_name())
        new_lmcnode.dist = llog[-1][0] - prevdist
#        new_lmcnode.data['freq'] = llog[-1][2]
        lmctree.add_child(lmcparent, new_lmcnode)
        # add pops entry
        pops[new_lmcnode.name] = popsize * llog[-1][2]
        # add revrecon entry
        revrecon[lnode.name] = new_lmcnode.name
        
        # add recon entry
        recon[new_lmcnode.name] = lnode.name
        
        # recurse down the ltree
        for lchild in lnode.children:
            walk(lchild, new_lmcnode)
    
    
    for lchild in ltree.root.children:
        walk(lchild, lmctree.root)
    
    extras = [daughters, pops, recon, revrecon]
    return lmctree, extras




#=============================================================================
# TESTS


def debug_test1():
    stree = treelib.read_tree('../examples/flies.stree')
    for node in stree:
        node.dist *= 1e7 # gen per myr
    popsize = 2e7
    freq = 1e0
    dr = .0012/1e7
    lr = .0006/1e7
    freqdup = freqloss = .05
    forcetime = 1e7
    
    ltree, ex = sample_locus_tree_hem(stree, popsize, freq, dr, lr, freqdup, freqloss, forcetime)
    
    return stree, gtree, ex

def debug_test2():
    stree = treelib.read_tree('examples/flies.stree') # run from ../ of this directory
    for node in stree:
        node.dist *= 1e7 # gen per myr
    popsize = 2e7
    freq = 1e0
    dr = .0012/1e7
    lr = .0006/1e7
    freqdup = freqloss = .05
    forcetime = 1e7
    
#    ltree, ex = sample_locus_tree_hem(stree, popsize, freq, dr, lr, freqdup, freqloss, forcetime)
    
    coal_tree, ex = sample_dlcoal_no_ifix(stree=stree, n=popsize, freq=freq, duprate=dr, lossrate=lr, freqdup=freqdup, freqloss=freqloss, forcetime=forcetime)
    
    treelib.draw_tree(coal_tree, scale=.00000005)



def debug_test3():
    stree = treelib.read_tree('examples/nbin.stree') # run from ../ of this directory
    for node in stree:
        node.dist *= 1e7 # gen per myr
    popsize = 2e7
    freq = 1e0
    dr = .0000012 / 1e7 #.0012/1e7
    lr = .0000011 / 1e7 #.0006/1e7
    freqdup = freqloss = .05
    forcetime = 1e7
    
    for node in stree:
        print node.name, node.dist, len(node.children)
    print
    
    locus_tree, locus_extras = sample_locus_tree_hem(stree, popsize, freq, \
                                                        dr, lr, \
                                                        freqdup, freqloss, \
                                                        forcetime)
    
    for node in locus_tree:
        print node.name, node.dist, len(node.children)
    print
    
    logged_locus_tree, logged_extras = locus_to_logged_tree(locus_tree, popsize)
    daughters = logged_extras[0]
    pops = logged_extras[1]
    
    coal_tree, coal_recon = dlcoal.sample_mutlilocus_tree(logged_locus_tree,
                                    n=pops, daughters=daughters,
                                    namefunc=lambda x: logged_extras[2][x] + '_' + str(x))
    
    #begin debug
    print coal_tree.leaf_names()
    try:
#        print set(coal_tree) - set(coal_tree.postorder())
        treelib.assert_tree(coal_tree)
    except AssertionError:
        print 'assertion error thrown on coal_tree being a proper tree'
        from rasmus import util
        hd= util.hist_dict(x.name for x in coal_tree.postorder())
        for key in hd.keys():
            print key if hd[key]>1 else '',
        print
        print len(coal_tree.nodes) - len(list(coal_tree.postorder()))
    



#def debug_test3():
#    testtree = treelib.parse_newick('(A:1,(C:1)B:1);')
#    for node in testtree:
#        node.dist *= 1e7
#    pops = {}
#    for i in testtree:
#        pops[i.name] = 2e5
#    daughters = set()
#    namefunc = lambda x:x
#    
#    treelib.assert_tree(testtree)
#    
#    coal_tree, coal_recon = dlcoal.sample_multilocus_tree(testtree,
#                                        n=pops,
#                                        daughters=daughters,
#                                        namefunc=namefunc)
#    
#    treelib.assert_tree(coal_tree)


#if __name__ == "__main__":
#    stree = treelib.read_tree('simple.stree')
##    print stree.root.name, stree.root.dist
#    popsize = 1e4
#    freq = 1e0
#    dr = 1.1 / 1e6
#    lr = 1.0 / 1e6
#    freqdup = .07
#    freqloss = .05
#    forcetime = 1e6 #1e0  # updated with species tree branch conversion
#    
#    for node in stree:
#        node.dist *= 1e6
#    
#    coal_tree, coal_extras = sample_dlcoal_no_ifix(stree, popsize, freq, \
#                                dr, lr, freqdup, freqloss,forcetime, minsize=3)
#    locus_tree = coal_extras['locus_tree']
#    
#    for daughter in coal_extras['daughters']:
#        print daughter.name
#    
#    treelib.draw_tree(stree,scale=.000002)
#    print
#    treelib.draw_tree(locus_tree,scale=.000002)
#    print
#    treelib.draw_tree(coal_tree,scale=.000002)
#    print
#    
##    gtree, extras = sample_locus_tree_hem(stree, popsize, freq, dr, lr, freqdup, freqloss, forcetime)
##    treelib.draw_tree(stree, scale=1)
##    print 
##    for leaf in stree.leaves():
##        print leaf.name, treelib.find_dist(stree, stree.root.name, leaf.name)
##    print
##    treelib.draw_tree(gtree, scale=1)
##    print
##    
##    if len(gtree.nodes) > 1:
##        for leaf in gtree.leaves():
##            print leaf.name, treelib.find_dist(gtree,gtree.root.name,leaf.name), leaf.data['freq']
##    else:
##        print 'only the root remains'
##    print
#    
##    for gnode in gtree.nodes:
##        print gnode, extras['recon'][gtree.nodes[gnode]].name

##    print gtree.root.dist, gtree.root.data['log']
##    print len(gtree.root.children), gtree.root.children[0].dist, gtree.root.children[0].data['log']
##    print gtree.leaves()[0].data['log']

##    print gtree.leaves()[0].data['log']
##
##    removed = [node for node in gtree if len(node.children) == 1 and node.parent]
##    print len(removed)
##    ofinterest = [node for node in removed if node.children[0] not in removed]
##    print len(ofinterest)
##    print
##    for i in ofinterest[0].parent.data['log']:
##        print i
##    for i in ofinterest[0].data['log']:
##        print i
##    print len(ofinterest[0].parent.data['log'] + ofinterest[0].data['log'])
##    remove_single_children(gtree)
##    print
##    for i in ofinterest[0].data['log']:
##        print i


''''

    
    def sim_walk(gtree, snode, gnode, p, s_walk_time=0.0, g_walk_time=0.0, 
                 time_until_force=forcetime, eventlog=[]): # eventlog is bug
        """
        eventlog is a log of events along the gtree branch.
        Each entry has the form
          (time_on_branch, event_type, frequency, species_node),
          
        where
           0.0 <= time_on_branch <= branch_node.dist

        event_type is one of
           {'extinction', 'frequency', 'speciation', duplication',
            'loss', 'root', 'gene'},
            
        where 'root' is a unique event not added during the sim_walk process

        frequency is the branch frequency at the event time

        species_node is the name of the node of the species tree branch in
        which the event occurs
        """
        
        if p <= 0.0:
            # EXTINCTION EVENT
            # gnode is 'parent' of extinct node
            #  create new_gnode
            new_gnode = treelib.TreeNode(gtree.new_name())
            new_gnode.dist = g_walk_time
            
            # set new_gnode's frequency
            new_gnode.data['freq'] = 0.0
            gtree.add_child(gnode, new_gnode)
            
            # add extinction event to the event log
            ext_event = (g_walk_time, 'extinction', 0.0, snode.name)
            eventlog.append(ext_event)
            
            # set new_gnode's event log
            new_gnode.data['log'] = eventlog
            return
            
        # extinction does not occur
        eff_duprate = duprate * p # * popsize #??
        eff_lossrate = lossrate * p # * popsize #??
        eff_bothrate = eff_duprate + eff_lossrate

        # sample next dup/loss event
        event_time = stats.exponentialvariate(eff_bothrate)
        remaining_s_dist = snode.dist - s_walk_time

        if event_time >= min(time_until_force, remaining_s_dist):
            # D/L event does not occur in this time period
            # determine whether at force or speciation

            if time_until_force < remaining_s_dist:
                # next time step
                # sample a new frequency

                newp = coal.sample_freq_CDF(p, popsize, forcetime)

                # update walk times
                new_s_walk_time = s_walk_time + time_until_force
                new_g_walk_time = g_walk_time + time_until_force

                # add frequency event to event log
                freq_event = (new_g_walk_time, 'frequency', newp, snode.name)
                eventlog.append(freq_event)

                # continue the walk with a reset forcetime
                sim_walk(gtree, snode, gnode, newp, 
                         s_walk_time=new_s_walk_time, 
                         g_walk_time=new_g_walk_time, 
                         eventlog=eventlog)

            else:
                # SPECIATION EVENT
                # separate into separate root, non-root speciations

                if gnode.data['log'][-1][1] != 'root':
                    # sample a new frequency

                    newp = coal.sample_freq_CDF(
                        p, popsize, remaining_s_dist)

                    # create new_gnode for this event
                    new_gnode = treelib.TreeNode(gtree.new_name())
                    new_g_walk_time = g_walk_time + remaining_s_dist
                    new_gnode.dist = new_g_walk_time

                    # set new node's frequency
                    new_gnode.data['freq'] = newp
                    gtree.add_child(gnode, new_gnode)
                    # add speciation event to event log and
                    # set the new node's log
                    if snode.is_leaf():
                        gene_event = (new_g_walk_time, 'gene',
                                      newp, snode.name)
                        eventlog.append(gene_event)
                        new_gnode.data['log'] = eventlog
                        # end of walk on species branch

                    else:
                        spec_event = (new_g_walk_time,
                                      'speciation', newp, snode.name)
                        eventlog.append(spec_event)
                        new_gnode.data['log'] = eventlog
                        for schild in snode.children:
                            sim_walk(gtree, schild, new_gnode,
                                     newp, eventlog=[])

                else: # gnode is the root
                    spec_event = (0.0, 'speciation', p, snode.name)
                    eventlog = gnode.data['log']
                    eventlog.append(spec_event)
                    gnode.data['log'] = eventlog

                    for schild in snode.children:
                        sim_walk(gtree, schild, gnode, p, eventlog=[])

        else:
            # process D/L event
            # no WF updates for these events (modelling decision)

            new_s_walk_time = s_walk_time + event_time
            new_g_walk_time = g_walk_time + event_time
            new_time_until_force = time_until_force - event_time

            if event_is_dup(eff_duprate, eff_bothrate):
                ## DUPLICATION EVENT
                # create a node new_gnode for the duplication event
                new_gnode = treelib.TreeNode(gtree.new_name())
                new_gnode.dist = new_g_walk_time
                # set new node's frequency
                new_gnode.data['freq'] = p
                gtree.add_child(gnode, new_gnode)
                # add duplication event to event log and
                # set the new node's log
                dup_event = (new_g_walk_time, 'duplication', p, snode.name)
                eventlog.append(dup_event)
                new_gnode.data['log'] = eventlog

                # recurse on remainder of original branch
                sim_walk(gtree, snode, new_gnode, p, \
                            s_walk_time=new_s_walk_time, \
                            time_until_force = new_time_until_force, \
                            eventlog=[])

                # recurse on dup tree with correct starting frequency
                sim_walk(gtree, snode, new_gnode, freqdup, \
                            s_walk_time=new_s_walk_time, \
                            time_until_force = new_time_until_force, \
                            eventlog=[(0.0,'daughter',freqdup,snode.name)])

            else:
                ## LOSS EVENT
                newp = max(p - freqloss, 0.0) # sanity check
                # add loss event to event log
                loss_event = (new_g_walk_time, 'loss', newp, snode.name)
                eventlog.append(loss_event)
                sim_walk(gtree, snode, gnode, newp, \
                            s_walk_time=new_s_walk_time, \
                            g_walk_time=new_g_walk_time, \
                            time_until_force=new_time_until_force, \
                            eventlog=eventlog)
    

'''
