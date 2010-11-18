"""

   Code for the DLCoal model
   (duplications, losses, and coalescence)

"""

from __future__ import division

# python libs
import copy
import os
import sys
import random
from itertools import chain, izip
import traceback
from math import *

# import dlcoal C lib
from dlcoal.ctypes_export import *
dlcoalc = load_library(["..", "lib"], "libdlcoal.so")


# add pre-bundled dependencies to the python path,
# if they are not available already
try:
    import rasmus, compbio
except ImportError:
    from . import dep
    dep.load_deps()
    import rasmus, compbio


# rasmus libs
from rasmus import stats, util, treelib

# compbio libs
from compbio import birthdeath, phylo

# dlcoal libs
from . import coal, duploss


#=============================================================================
# export c functions

ex = Exporter(globals())
export = ex.export

if dlcoalc:
    export(dlcoalc, "deleteTree", c_int, [c_void_p, "tree"])
    export(dlcoalc, "makeTree", c_void_p, [c_int, "nnodes",
                                          c_int_p, "ptree"])
    export(dlcoalc, "setTreeDists", c_void_p, [c_void_p, "tree",
                                              c_float_p, "dists"])




#=============================================================================
# probability functions for DLCoal model

def prob_dlcoal_recon_topology(coal_tree, coal_recon,
                               locus_tree, locus_recon, locus_events,
                               daughters,
                               stree, n, duprate, lossrate,
                               pretime=None, premean=None,
                               maxdoom=20, nsamples=100,
                               add_spec=True, info=None):
    """
    Probability of a reconcile gene tree in the DLCoal model.

    coal_tree    -- coalescent tree
    coal_recon   -- reconciliation of coalescent tree to locus tree
    locus_tree   -- locus tree (has dup-loss)
    locus_recon  -- reconciliation of locus tree to species tree
    locus_events -- events dict for locus tree
    stree        -- species tree
    n            -- population sizes in species tree
    duprate      -- duplication rate
    lossrate     -- loss rate

    You must also specify one of the following
    pretime      -- starting time before species tree
    premean      -- mean starting time before species tree

    """

    
    # init popsizes for locus tree
    stree_popsizes = coal.init_popsizes(stree, n)
    popsizes = {}
    for node in locus_tree:
        popsizes[node.name] = stree_popsizes[locus_recon[node].name]
    
    
    # duploss probability
    dl_prob = prob_tree_birth_death(
        locus_tree, stree, locus_recon, locus_events,
        duprate, lossrate, maxdoom=maxdoom)
    
    # daughters probability
    dups = phylo.count_dup(locus_tree, locus_events)
    d_prob = dups * log(.5)
    
    # integrate over duplication times using sampling
    stimes = treelib.get_tree_timestamps(stree)
    prob = prob_locus_coal_recon_topology_samples(
        coal_tree, coal_recon,
        locus_tree, locus_recon, locus_events, popsizes,
        stree, stimes,
        daughters, duprate, lossrate, nsamples,
        pretime, premean)

    
    # logging info
    if info is not None:
        info["duploss_prob"] = dl_prob
        info["daughters_prob"] = d_prob
        info["coal_prob"] = prob
        info["prob"] = dl_prob + d_prob + prob - log(nsamples)
    
    return dl_prob + d_prob + prob - log(nsamples)


def prob_locus_coal_recon_topology_samples(
        coal_tree, coal_recon,
        locus_tree, locus_recon, locus_events, popsizes,
        stree, stimes,
        daughters, duprate, lossrate, nsamples,
        pretime=None, premean=None):
    
    if dlcoalc:
        # sample some reason branch lengths just for logging purposes
        locus_times = duploss.sample_dup_times(
                locus_tree, stree, locus_recon, duprate, lossrate, pretime,
                premean,
                events=locus_events)
        treelib.set_dists_from_timestamps(locus_tree, locus_times)
        #treelib.draw_tree(locus_tree, scale=.0000001)

        # use C code
        return coal.prob_locus_coal_recon_topology_samples(
            coal_tree, coal_recon,
            locus_tree, locus_recon, locus_events, popsizes,
            stree, stimes,
            daughters, duprate, lossrate, nsamples, pretime, premean)
    else:
        # python backup    
        prob = 0.0
        for i in xrange(nsamples):
            # sample duplication times
            locus_times = duploss.sample_dup_times(
                locus_tree, stree, locus_recon, duprate, lossrate, pretime,
                premean,
                events=locus_events)
            treelib.set_dists_from_timestamps(locus_tree, locus_times)

            # coal topology probability
            coal_prob = prob_locus_coal_recon_topology(
                coal_tree, coal_recon, locus_tree, popsizes, daughters)
            
            prob += exp(coal_prob)
        prob = util.safelog(prob / nsamples)

        return prob


def prob_locus_coal_recon_topology(tree, recon, locus_tree, n, daughters):
    """
    Returns the log probability of a reconciled gene tree ('tree', 'recon')
    from the coalescent model given a locus tree 'locus_tree',
    population sizes 'n', and daughters set 'daughters'
    """

    # initialize popsizes, lineage counts, and divergence times
    popsizes = coal.init_popsizes(locus_tree, n)
    lineages = coal.count_lineages_per_branch(tree, recon, locus_tree)
    locus_times = treelib.get_tree_timestamps(locus_tree)


    # calc log probability
    lnp = coal.pmrt(
        tree, recon, locus_tree, popsizes, lineages=lineages)

    def walk(node, gene_counts, leaves):
        if node.is_leaf():
            gene_counts[node.name] = lineages[node][0]
            leaves.add(node)
        else:
            for child in node.children:
                if child in daughters:
                    gene_counts[child.name] = 1
                    leaves.add(child)
                else:
                    walk(child, gene_counts, leaves)

    for daughter in daughters:
        # determine leaves of the coal subtree
        gene_counts = {}
        leaves = set()
        walk(daughter, gene_counts, leaves)

        p = coal.cdf_mrca_bounded_multicoal(
            gene_counts, locus_times[daughter.parent], locus_tree, popsizes,
            sroot=daughter, sleaves=leaves, stimes=locus_times)

        if p == -util.INF:
            return -util.INF

        lnp -= p
    
    return lnp



def prob_tree_birth_death(tree, stree, recon, events, duprate, lossrate,
                          maxdoom=20):
    try:
        return duploss.prob_dup_loss(
            tree, stree, recon, duprate, lossrate,
            maxdoom=maxdoom, events=events)
        
    except Exception, e:
        if "dlcoal_python_fallback" not in globals():
            traceback.print_exception(type(e), e, sys.exc_info()[2],
                                      file=sys.stderr)
            print >>sys.stderr, "warning: using python code instead of native"
            globals()["dlcoal_python_fallback"] = 1
            # spidir libs
            import spidir
            from spidir import topology_prior
            
        return topology_prior.dup_loss_topology_prior(
            tree, stree, recon, duprate, lossrate,
            maxdoom=maxdoom, events=events)
    


#=============================================================================
# sampling from the DLCoal model

def sample_dlcoal(stree, n, duprate, lossrate, namefunc=lambda x: x,
                  remove_single=True, name_internal="n",
                  minsize=0, reject=False):
    """Sample a gene tree from the DLCoal model"""

    # generate the locus tree
    while True:
        # TODO: does this take a namefunc?
        locus_tree, locus_recon, locus_events = \
                    birthdeath.sample_birth_death_gene_tree(
            stree, duprate, lossrate)
        if len(locus_tree.leaves()) >= minsize:
            break

    if len(locus_tree.nodes) <= 1:
        # total extinction
        coal_tree = treelib.Tree()
        coal_tree.make_root()
        coal_recon = {coal_tree.root: locus_tree.root}
        daughters = set()
    else:
        # simulate coalescence
        
        # choose daughter duplications
        daughters = set()
        for node in locus_tree:
            if locus_events[node] == "dup":
                daughters.add(node.children[random.randint(0, 1)])

        if reject:
            # use slow rejection sampling (for testing)
            coal_tree, coal_recon = sample_locus_coal_tree_reject(
                locus_tree, n, daughters=daughters, namefunc=namefunc)
        else:
            coal_tree, coal_recon = sample_locus_coal_tree(
                locus_tree, n, daughters=daughters, namefunc=namefunc)

        # clean up coal tree
        if remove_single:
            treelib.remove_single_children(coal_tree)
            phylo.subset_recon(coal_tree, coal_recon)

    if name_internal:
        rename_nodes(coal_tree, name_internal)
        rename_nodes(locus_tree, name_internal)

    # store extra information
    extra = {"locus_tree": locus_tree,
             "locus_recon": locus_recon,
             "locus_events": locus_events,
             "coal_tree": coal_tree,
             "coal_recon": coal_recon,
             "daughters": daughters}

    return coal_tree, extra


def rename_nodes(tree, prefix="n"):
    """Rename nodes that all names are strings"""
    for node in list(tree.postorder()):
        if isinstance(node.name, int):
            name2 = prefix + str(node.name)
            while name2 in tree.nodes:
                name2 = prefix + str(tree.new_name())
            tree.rename(node.name, name2)


def sample_locus_coal_tree(stree, n, leaf_counts=None,
                           daughters=set(),
                           namefunc=None):
    """
    Returns a gene tree from a locus coalescence process
    n -- population size (int or dict)
         If n is a dict it must map from species name to population size
    """
    
    # initialize vector for how many genes per extant species
    if leaf_counts is None:
        leaf_counts = dict((l, 1) for l in stree.leaf_names())

    # initialize function for generating new gene names
    if namefunc is None:
        spcounts = dict((l, 1) for l in stree.leaf_names())
        def namefunc(sp):
            name = sp + "_" + str(spcounts[sp])
            spcounts[sp] += 1
            return name

    stimes = treelib.get_tree_timestamps(stree)

    # initialize population sizes
    popsizes = coal.init_popsizes(stree, n)

    # init gene counts
    counts = dict((n.name, 0) for n in stree)
    counts.update(leaf_counts)

    # init lineage counts
    lineages = {stree.root: [None, None]}
    for node in stree.leaves():
        lineages[node] = [leaf_counts[node.name], None]
    for node in daughters:
        if node not in lineages:
            lineages[node] = [None, 1]
        else:
            lineages[node][1] = 1
        

    def get_subtree(node, leaves, leaf_counts2):
        """collects info of subtree rooted at node"""
        if node.is_leaf():
            leaves.add(node)
            leaf_counts2[node.name] = leaf_counts[node.name]
        else:
            for child in node.children:
                if child in daughters:
                    leaves.add(child)
                    leaf_counts2[child.name] = 1
                else:
                    get_subtree(child, leaves, leaf_counts2)

    # loop through subtrees
    for snode in chain(daughters, [stree.root]):
        # determine leaves of the coal subtree
        leaves = set()
        leaf_counts2 = {}
        get_subtree(snode, leaves, leaf_counts2)
        
        if snode.parent:
            T  = stimes[snode.parent]
        else:
            T = None

        # calc table
        prob_counts = coal.calc_prob_counts_table(
            leaf_counts2, T, stree, popsizes,
            sroot=snode, sleaves=leaves, stimes=stimes)
        
        # sample lineage counts
        try:
            coal.sample_lineage_counts(snode, leaves, popsizes, stimes, T,
                                       lineages, prob_counts)
        except:
            print snode.name
            treelib.draw_tree_names(stree, maxlen=8)
            util.print_dict(lineages, key=lambda x: x.name)
            raise


    # sample coal times
    tree, recon = coal.coal_cond_lineage_counts(
        lineages, stree.root, set(stree.leaves()),
        popsizes, stimes, None, namefunc)
    
    return tree, recon


def sample_locus_coal_tree_reject(locus_tree, n, leaf_counts=None,
                                  daughters=set(),
                                  namefunc=None):
    
    # use rejection sampling
    while True:
        coal_tree, coal_recon = coal.sample_multicoal_tree(
            locus_tree, n, namefunc=lambda x: x)
        lineages = coal.count_lineages_per_branch(
            coal_tree, coal_recon, locus_tree)
        for daughter in daughters:
            if lineages[daughter][1] != 1:
                break
        else:
            break

    return coal_tree, coal_recon



## added 30 July 2010
##  should mimic dlcoal_sims, but using the new simulator
## TODO: fix this up when code is more integrated
def dlcoal_sims2(outdir, nsims, stree, n, duprate, lossrate, freq=1.0,
                freqdup=.05, freqloss=.05, forcetime=1e6, start=0, minsize=0,
                **options):
    
    import sim_v2_2 # TODO: remove this when code is more integrated
    
    for i in xrange(start, nsims+start): # changed nsims to nsims+start
        outfile = phylo.phylofile(outdir, str(i), "")
        util.makedirs(os.path.dirname(outfile))
        print "simulating", outfile
 
        # sample a new tree from DLCoal model using new simulator
        coal_tree, ex = sim_v2_2.sample_dlcoal_no_ifix(stree, n, freq, \
                            duprate, lossrate, freqdup, freqloss, forcetime, \
                            minsize=minsize, **options)
        write_dlcoal_recon(outfile, coal_tree, ex)



#=============================================================================
# Input/Output


def dlcoal_sims(outdir, nsims, stree, n, duprate, lossrate,
                start=0,
                **options):
    
    for i in xrange(start, nsims):
        outfile = phylo.phylofile(outdir, str(i), "")
        util.makedirs(os.path.dirname(outfile))
        print "simulating", outfile

        # sample a new tree from DLCoal model
        coal_tree, ex = sample_dlcoal(stree, n, duprate, lossrate, **options)
        write_dlcoal_recon(outfile, coal_tree, ex)










def write_dlcoal_recon(filename, coal_tree, extra,
                       exts={"coal_tree": ".coal.tree",
                             "coal_recon": ".coal.recon",
                             "locus_tree": ".locus.tree",
                             "locus_recon": ".locus.recon",
                             "daughters": ".daughters"
                             },
                       filenames={}):
    """Writes a reconciled gene tree to files"""

    # coal
    coal_tree.write(filenames.get("coal_tree", filename + exts["coal_tree"]),
                    rootData=True)
    phylo.write_recon_events(
        filenames.get("coal_recon", filename + exts["coal_recon"]),
        extra["coal_recon"], noevent="none")

    # locus
    extra["locus_tree"].write(
        filenames.get("locus_tree", filename + exts["locus_tree"]),
        rootData=True)
    phylo.write_recon_events(
        filenames.get("locus_recon", filename + exts["locus_recon"]),
        extra["locus_recon"], extra["locus_events"])

    util.write_list(
        filenames.get("daughters", filename + exts["daughters"]),
        [x.name for x in extra["daughters"]])



def read_dlcoal_recon(filename, stree,
                      exts={"coal_tree": ".coal.tree",
                            "coal_recon": ".coal.recon",
                            "locus_tree": ".locus.tree",
                            "locus_recon": ".locus.recon",
                            "daughters": ".daughters"
                            },
                      filenames={}):
    """Reads a reconciled gene tree from files"""

    extra = {}

    # trees
    coal_tree = treelib.read_tree(
        filenames.get("coal_tree", filename + exts["coal_tree"]))
    extra["locus_tree"] = treelib.read_tree(
        filenames.get("locus_tree", filename + exts["locus_tree"]))

    # recons
    extra["coal_recon"], junk = phylo.read_recon_events(
        filenames.get("coal_recon", filename + exts["coal_recon"]),
        coal_tree, extra["locus_tree"])
    extra["locus_recon"], extra["locus_events"] = phylo.read_recon(
        filenames.get("locus_recon", filename + exts["locus_recon"]),
        extra["locus_tree"], stree)


    extra["daughters"] = set(
        extra["locus_tree"].nodes[x] for x in util.read_strings(
        filenames.get("daughters", filename + exts["daughters"])))

    return coal_tree, extra


def read_log(filename):
    """Reads a DLCoal log"""
    stream = util.open_stream(filename)
    for line in stream:
        yield eval(line, {"inf": util.INF})
  

def read_log_all(filename):
    """Reads a DLCoal log"""
    stream = util.open_stream(filename)
    return map(eval, stream)
    



#=============================================================================
# C interface functions

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

    # ensure sort is stable
    def leafsort(a, b):
        if a.is_leaf():
            if b.is_leaf():
                return 0
            else:
                return -1
        else:
            if b.is_leaf():
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


def ptree2ctree(ptree):
    """Makes a c++ Tree from a parent array"""
    pint = c_int * len(ptree)
    tree = makeTree(len(ptree), pint(* ptree))
    return tree


def tree2ctree(tree):
    """Make a c++ Tree from a treelib.Tree data structure"""
    ptree, nodes, nodelookup = make_ptree(tree)
    dists = [x.dist for x in nodes]
    ctree = ptree2ctree(ptree)
    setTreeDists(ctree, c_list(c_float, dists))
    return ctree



def make_recon_array(tree, recon, nodes, snodelookup):
    """Make a reconciliation array from recon dict"""
    recon2 = []
    for node in nodes:
        recon2.append(snodelookup[recon[node]])
    return recon2

def make_events_array(nodes, events):
    """Make events array from events dict"""
    mapping = {"gene": 0,
               "spec": 1,
               "dup": 2}
    return [mapping[events[i]] for i in nodes]



