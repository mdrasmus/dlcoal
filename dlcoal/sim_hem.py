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
from math import *
import random

# rasmus, compbio imports
from rasmus import stats, treelib
from compbio import coal, phylo

# dlcoal imports
import dlcoal

#=============================================================================
# debugging
DEBUGPRINT = True

def debugprint(s):
    if DEBUGPRINT:
        print s


#=============================================================================


def sim_DLILS_gene_tree(stree, popsize, freq, dr, lr,
                        freqdup, freqloss, forcetime):
    
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
    
    Update: 30 July 2010
     Will return the gene (locus) tree, as well as extra information including a
      reconciliation dictionary and an events dictionary.
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
                    time_until_force=forcetime, eventlog=[]):
###    Most of the variables are obvious from descriptions in sim_tree or similar.
###    eventlog is a log of events along the gtree branch; each entry has the form
###     (time_on_branch, event_type, frequency, species_node),
###     where 0.0 <= time_on_branch <= branch_node.dist
###     event_type is one of {'extinction', 'frequency', 'speciation', 
###       duplication', 'loss', 'root', 'gene'}, where 'root' is a unique event
###       not added during the sim_walk process
###     frequency is the branch frequency at the event time
###     species_node is the name of the node of the species tree branch 
###       in which the event occurs
        if p <= 0.0:
            ## EXTINCTION EVENT
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
            eventlog = [] # should have no effect; added for debugging on 18 Oct 2010
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
                    ## FREQUENCY UPDATE EVENT
                    # sample a new frequency (note scaling to years from myr) # edit: not any more
                    newp = coal.sample_freq_CDF(p, popsize, forcetime) # * 1e6)
                      # TODO: if we decide not to reset time_until_force at 
                      #  speciation events, the newp generation will need to be
                      #  altered in some form (probably using a new variable)
                    # update walk times
                    new_s_walk_time = s_walk_time + time_until_force
                    new_g_walk_time = g_walk_time + time_until_force
                    # add frequency event to event log
                    freq_event = (new_g_walk_time, 'frequency', newp, snode.name)
                    eventlog.append(freq_event)
                    # continue the walk with a reset forcetime
                    sim_walk(gtree, snode, gnode, newp, \
                                s_walk_time=new_s_walk_time, \
                                g_walk_time=new_g_walk_time, \
                                eventlog=eventlog)
                    eventlog = [] # should have no effect; debug add on 18 Oct 2010
                else:
                    ## SPECIATION EVENT
                    # separate into separate root, non-root speciations
#                    if gnode.parent: # gnode not the root
                    if gnode.data['log'][-1][1] != 'root':
                        # sample a new frequency (note scaling to years from myr) # edit: not any more
                        newp = coal.sample_freq_CDF(p, popsize, remaining_s_dist) # * 1e6)
                        # create new_gnode for this event
                        new_gnode = treelib.TreeNode(gtree.new_name())
                        new_g_walk_time = g_walk_time + remaining_s_dist
                        new_gnode.dist = new_g_walk_time
                        # set new node's frequency
                        new_gnode.data['freq'] = newp
                        gtree.add_child(gnode, new_gnode)
                        # add speciation event to event log and set the new node's log
                        if snode.is_leaf():
                            gene_event = (new_g_walk_time, 'gene', newp, snode.name)
                            eventlog.append(gene_event)
                            new_gnode.data['log'] = eventlog
                            # end of walk on species branch
                            eventlog = [] # should have no effect; debug add on 18 Oct 2010
                        else:
                            spec_event = (new_g_walk_time, 'speciation', newp, snode.name)
                            eventlog.append(spec_event)
                            new_gnode.data['log'] = eventlog
                            for schild in snode.children:
                                sim_walk(gtree, schild, new_gnode, newp, eventlog=[])
                              # TODO: if we decide not to reset time_until_force at
                              #  speciation events, this sim_walk call will need updating
                            eventlog = [] # should have no effect; debug add on 18 Oct 2010
                    else: # gnode is the root
                        spec_event = (0.0, 'speciation', p, snode.name)
                        eventlog = gnode.data['log']
                        eventlog.append(spec_event)
                        gnode.data['log'] = eventlog
#                        ### debug print
#                        print
#                        print 'adding: ', eventlog
#                        ### end debug
                        for schild in snode.children:
                            sim_walk(gtree, schild, gnode, p, eventlog=[])
                        eventlog = [] # should have no effect; debug add on 18 Oct 2010
            else:
                # process D/L event
                # no WF updates for these events (modelling decision)
                new_s_walk_time = s_walk_time + event_time
                new_g_walk_time = g_walk_time + event_time
                new_time_until_force = time_until_force - event_time
                if event_is_dup(eff_dr, eff_bothr):
                    ## DUPLICATION EVENT
                    # create a node new_gnode for the duplication event
                    new_gnode = treelib.TreeNode(gtree.new_name())
                    new_gnode.dist = new_g_walk_time
                    # set new node's frequency
                    new_gnode.data['freq'] = p
                    gtree.add_child(gnode, new_gnode)
                    # add duplication event to event log and set the new node's log
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
                                eventlog=[(0.0,'daughter',freqdup,snode.name)]) # added for daughter detection
                    eventlog = [] # should have no effect; debug add on 18 Oct 2010
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
                    eventlog = [] # should have no effect; debug add on 18 Oct 2010
    
    
    # main code
    
    # create new gene tree and simulate its evolution
    gtree = treelib.Tree()
    gtree.make_root()
    gtree.root.dist = 0.0
    gtree.root.data['freq'] = freq
    root_event = (0.0, 'root', freq, stree.root.name)
    gtree.root.data['log'] = [root_event]
    sim_walk(gtree, stree.root, gtree.root, freq) # should mutate gtree
    
    
#    # remove dead branches and single children (inside last method)
#    # note that the simplifyRoot argument was added to the treelib methods
#    #  so that gtree.root.dist is always equal to 0.0 (and this allows the
#    #  root to have a single child)
#    #  if this behavior is undesired later, we can simply remove the argument
#    #  and the root will be collapsed (and have >0 dist)
    extant_leaves = []
    for leaf in gtree.leaves():
        if leaf.data['freq'] > 0.0:
            extant_leaves.append(leaf.name)
    gtree = treelib.subtree_by_leaf_names(gtree, extant_leaves,
                                          keep_single=True)
    remove_single_children(gtree) # allows for correct logging of events
    extras = generate_extras(stree, gtree)
    return gtree, extras


def remove_single_children(tree):
    """
    Modified from remove_single_children in treelib.py.
    Added log manipulation.
    """
    ###TODO: remove dud events? (remove extinct speciations and duplications?)
    
    # find single children
    removed = [node
               for node in tree
               if len(node.children) == 1 and node.parent]
    
    # actually remove children
    for node in removed:
        newnode = node.children[0]
        
        # add distance
        newnode.dist += node.dist
        
        # update logs
#        for i in xrange(len(newnode.data['log'])):
#            newnode.data['log'][i] = (newnode.data['log'][i][0]+node.dist,) \
#                                        + newnode.data['log'][i][1:]
#        newnode.data['log'] = node.data['log'] + newnode.data['log']
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


## generate extra information from the logging in a gene tree
## based on birthdeath.sample_birth_death_gene_tree
def generate_extras(stree, gtree, genename=lambda sp, x: sp + "_" + str(x)):
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
                    break # end after single child identified as daughter
        elif event == 'gene':
            events[gnode] = 'gene'
            # rename gene to comprehensible name
            gtree.rename(gnode.name, genename(snname, gnode.name))
        else: # extinction? shouldn't encounter this case
            assert False
#            events[gnode] = 'ext?'
        
        for gchild in gnode.children:
            walk(gchild)
    
    walk(gtree.root)
#    walk(gtree.root.children[0])
    
    return {'recon': recon, 'events': events, 'daughters': daughters}


def locus_to_logged_tree(ltree, popsize=1.0):
    
    ### note: ltree is the locus tree generated by the simulator
    ###       lmctree is the new logged-locus tree (locus-multi-coal-tree)
    
    
    # must have a valid locus tree (of the form from this simulator) passed in
    assert ltree.root.data['log'][0][1] == 'root'
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




## based on dlcoal.sample_dlcoal
## uses sim_DLILS_gene_tree to generate a locus tree
## samples a multicoal tree, etc., as in the initial function
def sample_dlcoal_no_ifix(stree, n, freq, duprate, lossrate, freqdup, freqloss,\
                            forcetime, namefunc=lambda x: x, \
                            remove_single=True, name_internal="n", minsize=0):
    """Sample a gene tree from the DLCoal model using the new simulator"""

    # generate the locus tree
    while True:
        locus_tree, locus_extras = sim_DLILS_gene_tree(stree, n, freq, \
                                                        duprate, lossrate, \
                                                        freqdup, freqloss, \
                                                        forcetime)
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
        logged_locus_tree, logged_extras = locus_to_logged_tree(locus_tree, popsize = n)
        daughters = logged_extras[0]
        pops = logged_extras[1]
        log_recon = logged_extras[2]
        
#        treelib.assert_tree(logged_locus_tree)
        
        # removed locus_tree_copy from below
        coal_tree, coal_recon = dlcoal.sample_locus_coal_tree(logged_locus_tree,
                                        n=pops, daughters=daughters,
                                        namefunc=lambda lognamex: log_recon[lognamex] + '_' + str(lognamex))

#        print set(coal_tree) - set(coal_tree.postorder())
        treelib.assert_tree(coal_tree)
    
        # clean up coal tree
        if remove_single:
            treelib.remove_single_children(coal_tree)
            phylo.subset_recon(coal_tree, coal_recon)


    if name_internal:
        dlcoal.rename_nodes(coal_tree, name_internal)
        dlcoal.rename_nodes(locus_tree, name_internal)


    # store extra information
    ### TODO: update this now that we're using logged locus tree, new sample function
    extra = {"locus_tree": locus_tree,
             "locus_recon": locus_extras['recon'],
             "locus_events": locus_extras['events'],
             "coal_tree": coal_tree,
             "coal_recon": coal_recon,
             "daughters": daughters}

    return coal_tree, extra



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
    
    ltree, ex = sim_DLILS_gene_tree(stree, popsize, freq, dr, lr, freqdup, freqloss, forcetime)
    
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
    
#    ltree, ex = sim_DLILS_gene_tree(stree, popsize, freq, dr, lr, freqdup, freqloss, forcetime)
    
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
    
    locus_tree, locus_extras = sim_DLILS_gene_tree(stree, popsize, freq, \
                                                        dr, lr, \
                                                        freqdup, freqloss, \
                                                        forcetime)
    
    for node in locus_tree:
        print node.name, node.dist, len(node.children)
    print
    
    logged_locus_tree, logged_extras = locus_to_logged_tree(locus_tree, popsize)
    daughters = logged_extras[0]
    pops = logged_extras[1]
    
    coal_tree, coal_recon = dlcoal.sample_locus_coal_tree(logged_locus_tree,
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
#    coal_tree, coal_recon = dlcoal.sample_locus_coal_tree(testtree,
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
##    gtree, extras = sim_DLILS_gene_tree(stree, popsize, freq, dr, lr, freqdup, freqloss, forcetime)
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
