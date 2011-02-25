
# python imports
import os
import random
from itertools import chain

# rasmus, compbio imports
from rasmus import treelib, util
from compbio import birthdeath, phylo

# dlcoal imports
import dlcoal
from . import coal



def dlcoal_sims(outdir, nsims, stree, n, duprate, lossrate,
                start=0,
                **options):
    
    for i in xrange(start, nsims):
        outfile = phylo.phylofile(outdir, str(i), "")
        util.makedirs(os.path.dirname(outfile))
        print "simulating", outfile

        # sample a new tree from DLCoal model
        coal_tree, ex = sample_dlcoal(stree, n, duprate, lossrate,
                                          **options)
        dlcoal.write_dlcoal_recon(outfile, coal_tree, ex)



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
            coal_tree, coal_recon = sample_multilocus_tree_reject(
                locus_tree, n, daughters=daughters, namefunc=namefunc)
        else:
            coal_tree, coal_recon = sample_multilocus_tree(
                locus_tree, n, daughters=daughters, namefunc=namefunc)

        # clean up coal tree
        if remove_single:
            treelib.remove_single_children(coal_tree)
            phylo.subset_recon(coal_tree, coal_recon)

    if name_internal:
        dlcoal.rename_nodes(coal_tree, name_internal)
        dlcoal.rename_nodes(locus_tree, name_internal)

    # store extra information
    extra = {"locus_tree": locus_tree,
             "locus_recon": locus_recon,
             "locus_events": locus_events,
             "coal_tree": coal_tree,
             "coal_recon": coal_recon,
             "daughters": daughters}

    return coal_tree, extra


def sample_multilocus_tree(stree, n, leaf_counts=None,
                           daughters=set(),
                           namefunc=None):
    """
    Returns a gene tree from a multilocus coalescent process
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


def sample_multilocus_tree_reject(locus_tree, n, leaf_counts=None,
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
