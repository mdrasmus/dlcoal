#
# Python module for SPIDIR library
#


#
# Note, this module requires the rasmus and compbio python modules.
#

import os
import sys
from math import *
from ctypes import *
from spidir.ctypes_export import *

# import spidir C lib
spidir = load_library(["..", "..", "lib"], "libspidir.so")


# add pre-bundled dependencies to the python path,
# if they are not available already
try:
    import rasmus, compbio
except ImportError:
    from . import dep
    dep.load_deps()
    import rasmus, compbio
from compbio import fasta

    


#=============================================================================
# wrap functions from c library

ex = Exporter(globals())
export = ex.export

# typedefs
c_floatlk = c_double
c_floatlk_p = c_double_p


if spidir:
    # common functions
    export(spidir, "gamm", c_double, [c_double, "a"])
    export(spidir, "invgammaPdf", c_double,
           [c_double, "x", c_double, "a", c_double, "b"])
    #export(spidir, "invgammaCdf", c_double,
    #       [c_double, "x", c_double, "a", c_double, "b"])
    #export(spidir, "quantInvgamma", c_double,
    #       [c_double, "p", c_double, "a", c_double, "b"])
    export(spidir, "gammalog", c_double,
           [c_double, "x", c_double, "a", c_double, "b"])
    export(spidir, "gammaPdf", c_double,
           [c_double, "x", c_double, "a", c_double, "b"])
    export(spidir, "gammaDerivX", c_double,
           [c_double, "x", c_double, "a", c_double, "b"])
    export(spidir, "gammaDerivA", c_double,
           [c_double, "x", c_double, "a", c_double, "b"])
    export(spidir, "gammaDerivB", c_double,
           [c_double, "x", c_double, "a", c_double, "b"])
    export(spidir, "gammaDerivV", c_double,
           [c_double, "x", c_double, "v"])
    export(spidir, "gammaDerivV2", c_double,
           [c_double, "x", c_double, "v"])
    export(spidir, "gammaSumPdf", c_double,
           [c_double, "y", c_int, "n", c_float_list, "alpha",
            c_float_list, "beta", 
            c_float, "tol"])
    
    #export(spidir, "negbinomPdf", c_double,
    #       [c_int, "k", c_double, "r", c_double, "p"])
    #export(spidir, "negbinomDerivR", c_double,
    #       [c_int, "k", c_double, "r", c_double, "p"])
    #export(spidir, "negbinomDerivP", c_double,
    #       [c_int, "k", c_double, "r", c_double, "p"])

    #export(spidir, "incompleteGammaC", c_double,
    #       [c_double, "s", c_double, "x"])

    # basic tree functions
    export(spidir, "deleteTree", c_int, [c_void_p, "tree"])
    export(spidir, "makeTree", c_void_p, [c_int, "nnodes",
                                          c_int_p, "ptree"])
    export(spidir, "tree2ptree", c_int, [c_void_p, "tree", c_int_list, "ptree"],
           newname="ctree2ptree")
    export(spidir, "setTreeDists", c_void_p, [c_void_p, "tree",
                                              c_float_p, "dists"])


    # search
    export(spidir, "searchClimb", c_void_p,
           [c_int, "niter", c_int, "quickiter",
            c_int, "nseqs", c_char_p_p, "gene_names", c_char_p_p, "seqs",
            c_int, "nsnodes", c_int_list, "pstree", c_float_list, "sdists",
            c_int_list, "gene2species",
            c_float_list, "sp_alpha", c_float_list, "sp_beta",
            c_float, "gene_rate",
            c_float, "pretime_lambda", c_float, "birth", c_float, "death",
            c_float, "gene_alpha", c_float, "gene_beta",
            c_float_list, "bgfreq", c_float, "kappa",
            c_int, "nsamples", c_int, "approx"])


    # topology prior birthdeath functions
    export(spidir, "inumHistories", c_int, [c_int, "nleaves"])
    export(spidir, "numHistories", c_double, [c_int, "nleaves"])
    #export(spidir, "numTopologyHistories", c_double, [c_void_p, "tree"])
    export(spidir, "birthDeathCount", c_double,
           [c_int, "ngenes", c_float, "time",
            c_float, "birth", c_float, "death"])
    export(spidir, "birthDeathCounts", c_double,
           [c_int, "start", c_int, "end", c_float, "time",
            c_float, "birth", c_float, "death"])
    export(spidir, "birthDeathCountsLog", c_double,
           [c_int, "start", c_int, "end", c_float, "time",
            c_float, "birth", c_float, "death"])
    export(spidir, "birthDeathCountsSlow", c_double,
           [c_int, "start", c_int, "end", c_float, "time",
            c_float, "birth", c_float, "death"])

    # birth death tree counts
    export(spidir, "birthDeathTreeCounts", c_double,
           [c_void_p, "tree", c_int, "nspecies", c_int_list, "counts", 
            c_float, "birth", c_float, "death", c_int, "maxgene",
            c_int, "rootgene", c_void_p, "tab"])
    export(spidir, "birthDeathForestCounts", c_double,
           [c_void_p, "tree", c_int, "nspecies", c_int, "nfams",
            c_int_matrix, "counts", c_int_list, "mult",
            c_float, "birth", c_float, "death", c_int, "maxgene",
            c_int, "rootgene", c_void_p, "tab"])

    export(spidir, "birthDeathCountsML_alloc", c_void_p,
           [c_void_p, "tree", c_int, "nspecies", c_int, "nfams",
            c_int_matrix, "counts", c_int_list, "mult",
            c_float, "birth", c_float, "death", c_float, "step",
            c_int, "maxgene", c_int, "rootgene"])
    export(spidir, "birthDeathCountsML_free", c_void_p,
           [c_void_p, "opt"])
    export(spidir, "birthDeathCountsML_iter", c_int,
           [c_void_p, "opt", c_float_list, "birth", c_float_list, "death",
            c_float_list, "size"])

    # tree topology
    export(spidir, "calcDoomTable", c_int,
           [c_void_p, "tree", c_float, "birth", c_float, "death",
            c_double_p, "doomtable"])
    export(spidir, "birthDeathTreePrior", c_double,
           [c_void_p, "tree", c_void_p, "stree",
            c_int_p, "recon", c_int_p, "events",
            c_float, "birth", c_float, "death",
            c_double_p, "doomtable"])
    export(spidir, "birthDeathTreePriorFull", c_double,
           [c_void_p, "tree", c_void_p, "stree",
            c_int_p, "recon", c_int_p, "events",
            c_float, "birth", c_float, "death",
            c_double_p, "doomtable"])
    export(spidir, "sampleBirthWaitTime", c_double,
           [c_int, "n", c_float, "T", c_float, "birth", c_float, "death"])
    export(spidir, "birthWaitTime", c_double,
           [c_float, "t", c_int, "n", c_float, "T",
            c_float, "birth", c_float, "death"])
    export(spidir, "probNoBirth", c_double,
           [c_int, "n", c_float, "T", c_float, "birth", c_float, "death"])
    export(spidir, "sampleBirthWaitTime1", c_double,
           [c_float, "T", c_float, "birth", c_float, "death"])



    # branch prior functions
    export(spidir, "branchPrior", c_double,
           [c_int, "nnodes", c_int_list, "ptree", c_float_list, "dists",
            c_int, "nsnodes", c_int_list, "pstree", c_float_list, "sdists",
            c_int_list, "recon", c_int_list, "events",
            c_float_list, "sp_alpha", c_float_list, "sp_beta",
            c_float, "generate", 
            c_float, "pretime_lambda", c_float, "dupprob", c_float, "lossprob",
            c_float, "gene_alpha", c_float, "gene_beta",
            c_int, "nsamples", c_int, "approx"])


    # parsimony
    export(spidir, "parsimony", c_void_p,
           [c_int, "nnodes", c_int_list, "ptree", c_int, "nseqs",
            c_char_p_p, "seqs", c_float_list, "dists",
            c_int, "buildAncestral", c_char_p_p, "ancetralSeqs"])

    # sequence likelihood functions
    export(spidir, "makeHkyMatrix", c_void_p,
           [c_float_p, "bgfreq", c_float, "kappa", c_float, "time",
            c_float_p, "matrix"])
    export(spidir, "makeHkyDerivMatrix", c_void_p,
           [c_float_p, "bgfreq", c_float, "kappa", c_float, "time",
            c_float_p, "matrix"])
    export(spidir, "makeHkyDeriv2Matrix", c_void_p,
           [c_float_p, "bgfreq", c_float, "kappa", c_float, "time",
            c_float_p, "matrix"])
    export(spidir, "branchLikelihoodHky", c_floatlk,
           [c_floatlk_p, "probs1", c_floatlk_p, "probs2",
            c_int, "seqlen", c_float_p, "bgfreq", c_float, "kappa",
            c_float, "time"])
    export(spidir, "branchLikelihoodHkyDeriv", c_floatlk,
           [c_floatlk_p, "probs1", c_floatlk_p, "probs2",
            c_int, "seqlen", c_float_p, "bgfreq",
            c_float, "kappa", c_float, "time"])       
    export(spidir, "branchLikelihoodHkyDeriv2", c_floatlk,
           [c_floatlk_p, "probs1", c_floatlk_p, "probs2",
            c_int, "seqlen", c_float_p, "bgfreq",
            c_float, "kappa", c_float, "time"])
    export(spidir, "mleDistanceHky", c_float,
           [c_floatlk_p, "probs1", c_floatlk_p, "probs2",
            c_int, "seqlen", c_float_p, "bgfreq",
            c_float, "kappa", c_float, "t0", c_float, "t1"])
    export(spidir, "calcSeqProbHky", c_floatlk,
           [c_void_p, "tree", c_int, "nseqs", c_char_p_p, "seqs",
            c_float_p, "bgfreq", c_float, "kappa"])
    export(spidir, "findMLBranchLengthsHky", c_floatlk,
           [c_int, "nnodes", c_int_p, "ptree", c_int, "nseqs",
            c_char_p_p, "seqs", c_float_p, "dists",
            c_float_p, "bgfreq", c_float, "kappa",
            c_int, "maxiter", c_int, "parsinit"])


    # training functions
    export(spidir, "train", c_void_p,
           [c_int, "ntrees", c_int, "nspecies",
            c_int_list, "gene_sizes",
            c_float_matrix, "lengths",
            c_float_list, "times",
            c_float_list, "sp_alpha", c_float_list, "sp_beta",
            c_float_list, "gene_alpha", c_float_list, "gene_beta",
            c_int, "nrates", c_int, "max_iter"])
    export(spidir, "allocRatesEM", c_void_p,
           [c_int, "ntrees", c_int, "nspecies", c_int, "nrates",
            c_int_p, "gene_sizes",
            c_float_p_p, "lengths", c_float_p, "times",
            c_float_p, "sp_alpha", c_float_p, "sp_beta", 
            c_float, "gene_alpha", c_float, "gene_beta"])
    export(spidir, "freeRatesEM", c_void_p, [c_void_p, "em"])
    export(spidir, "RatesEM_Init", c_void_p, [c_void_p, "em"])
    export(spidir, "RatesEM_EStep", c_void_p, [c_void_p, "em"])
    export(spidir, "RatesEM_MStep", c_void_p, [c_void_p, "em"])
    export(spidir, "RatesEM_likelihood", c_float, [c_void_p, "em"])
    export(spidir, "RatesEM_getParams", c_void_p,
           [c_void_p, "em", c_float_p, "params"])




#=============================================================================
# additional python interface

def read_params(filename):
    """Read SPIDIR model parameters to a file"""
    
    infile = file(filename)
    params = {}
    
    for line in infile:
        tokens = line.split("\t")
        key = tokens[0]
        values = tokens[1:]
        if key[0].isdigit():
            key = int(key)
        params[key] = map(float, values)
        
    return params


def write_params(filename, params):
    """Write SPIDIR model parameters to a file"""
    
    out = open(filename, "w")
    for key, value in params.iteritems():
        out.write("\t".join(map(str, [key] + value)) + "\n")
    out.close()

    

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


def ptree2tree(ptree, genes):
    """Create a Tree object from a ptree array"""

    from rasmus.treelib import Tree, TreeNode

    tree = Tree()
    nodes = [TreeNode(gene) for gene in genes] + \
            [TreeNode(i) for i in xrange(len(genes)-1)]
    
    for i, p in enumerate(ptree):
        if p == -1:
            tree.root = nodes[i]
        else:
            tree.add_child(nodes[p], nodes[i])
    
    return tree

    


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


def ctree2tree(ctree, genes):
    """Makes a treelib.Tree from a c++ Tree"""
    nnodes = 2*len(genes) - 1
    ptree = [0] * nnodes
    ctree2ptree(ctree, ptree)
    return ptree2tree(ptree, genes)
    
        

def make_gene2species_array(genes, stree, snodelookup, gene2species):
    """Make a gene2species array"""
    gene2speciesarray = []

    for g in genes:
        gene2speciesarray.append(snodelookup[stree.nodes[gene2species(g)]])
    for i in xrange(len(genes)-1):
        gene2speciesarray.append(-1)
    
    return gene2speciesarray


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


def search_climb(genes, align, stree, gene2species,
                 params, birth, death, pretime,
                 bgfreq=[.25,.25,.25,.25], kappa=1.0,
                 niter=50, quickiter=100,                 
                 nsamples=100, branch_approx=True):
    """Search for a MAP gene tree"""

    nseqs = len(align)
    calign = c_list(c_char_p, align)
    gene_names = c_list(c_char_p, genes)
    pstree, snodes, snodelookup = make_ptree(stree)    
    nsnodes = len(snodes)
    sdists = [x.dist for x in snodes]
    
    smap = make_gene2species_array(genes, stree, snodelookup, gene2species)

    sp_alpha = [params[x.name][0] for x in snodes]
    sp_beta = [params[x.name][1] for x in snodes]

    gene_alpha, gene_beta = params["baserate"]
    gene_rate = -1
    
    ctree = searchClimb(niter, quickiter,
                        nseqs, gene_names, calign,
                        nsnodes, pstree, sdists,
                        smap,
                        sp_alpha, sp_beta, gene_rate,
                        pretime, birth, death,
                        gene_alpha, gene_beta,
                        bgfreq, kappa,
                        nsamples, branch_approx)

    nnodes = 2 * nseqs - 1
    tree = ctree2tree(ctree, genes)
    
    deleteTree(ctree)

    return tree

    
    
def calc_joint_prob(align, tree, stree, recon, events, params,
                    birth, death, pretime,
                    bgfreq, kappa,
                    nsamples=100, branch_approx=True,
                    terms=False):
    """Calculate the joint probability of a gene tree"""
    
    branchp = branch_prior(tree, stree, recon, events,
                                  params, birth, death, pretime,
                                  nsamples, approx=branch_approx)
    topp = calc_birth_death_prior(tree, stree, recon,
                                  birth, death, 
                                  events=events)
    seqlk = calc_seq_likelihood_hky(tree, align, bgfreq, kappa)
    
    if terms:
        return branchp, topp, seqlk
    else:
        return branchp + topp + seqlk


def find_parsimony(tree, align):

    ptree, nodes, nodelookup = make_ptree(tree)
    nnodes = len(nodes)
    dists = [x.dist for x in nodes]

    nseqs = len(align)
    seqlen = len(align.values()[0])
    leaves = [x for x in nodes if x.is_leaf()]
    calign = (c_char_p * nseqs)(* [align[x.name] for x in leaves])
    cancestral = (c_char_p * (nseqs-1))(* ["-"*seqlen for i in xrange(nseqs-1)])

    #print ">>>", list(cancestral)
    parsimony(nnodes, ptree, nseqs, calign, dists, True, cancestral)
    #print list(cancestral)
    
    ancestral = fasta.FastaDict()
    for i, key in enumerate(node.name for node in nodes if not node.is_leaf()):
        ancestral[key] = cancestral[i]

    return ancestral


#=============================================================================
# topology prior

def calc_birth_death_prior(tree, stree, recon, birth, death, 
                           events=None):
    """Returns the topology prior of a gene tree"""

    from rasmus.bio import phylo
    
    if events is None:
        events = phylo.label_events(tree, recon)

    ptree, nodes, nodelookup = make_ptree(tree)
    pstree, snodes, snodelookup = make_ptree(stree)

    ctree = tree2ctree(tree)
    cstree = tree2ctree(stree)
    recon2 = make_recon_array(tree, recon, nodes, snodelookup)
    events2 = make_events_array(nodes, events)

    doomtable = c_list(c_double, [0] * len(stree.nodes))
    calcDoomTable(cstree, birth, death, doomtable)
    
    p = birthDeathTreePriorFull(ctree, cstree,
                                c_list(c_int, recon2), 
                                c_list(c_int, events2),
                                birth, death, doomtable)
    deleteTree(ctree)
    deleteTree(cstree)

    return p


def birth_death_tree_counts(stree, counts, birth, death,
                            maxgene=50, rootgene=1):
    
    if birth == death:
        birth = 1.01 * death

    ctree = tree2ctree(stree)
    nspecies = len(counts)
    
    prob = birthDeathTreeCounts(ctree, nspecies, counts,
                                birth, death, maxgene, rootgene, 0)
    deleteTree(ctree)

    return prob


def birth_death_forest_counts(stree, counts, birth, death,
                              maxgene=50, rootgene=1, mult=None):
    
    if birth == death:
        birth = 1.01 * death

    if mult is None:
        hist = {}
        for row in counts:
            row = tuple(row)
            hist[row] = hist.get(row, 0) + 1
        counts, mult = zip(*hist.items())
        counts = map(list, counts)
        mult = list(mult)
        

    ctree = tree2ctree(stree)

    nfams = len(counts)
    nspecies = len(counts[0])
    
    logl = birthDeathForestCounts(ctree, nspecies, nfams, counts, mult,
                                  birth, death, maxgene, rootgene, 0)
    deleteTree(ctree)

    return logl

def birth_death_counts_ml_alloc(stree, counts, birth0, death0, step,
                                maxgene=50, rootgene=1, mult=None):
    
    if birth0 == death0:
        birth0 = 1.01 * death0

    if mult is None:
        hist = {}
        for row in counts:
            row = tuple(row)
            hist[row] = hist.get(row, 0) + 1
        counts, mult = zip(*hist.items())
        counts = map(list, counts)
        mult = list(mult)

    ctree = tree2ctree(stree)

    nfams = len(counts)
    nspecies = len(counts[0])
    
    opt = birthDeathCountsML_alloc(ctree, nspecies, nfams, counts, mult,
                                   birth0, death0, step, maxgene, rootgene)
    return (ctree, opt)

def birth_death_counts_ml_free(opt):
    birthDeathCountsML_free(opt[1])
    deleteTree(opt[0])


def birth_death_counts_ml_iter(opt):

    birth = [0.0]
    death = [0.0]
    size = [0.0]
    status = birthDeathCountsML_iter(opt[1], birth, death, size)

    return status, size[0], (birth[0], death[0])


#=============================================================================
# branch prior

def branch_prior(tree, stree, recon, events, params, birth, death,
                 pretime_lambda=1.0,
                 nsamples=1000, approx=True,
                 generate=None):
    """Returns the branch prior of a gene tree"""

    ptree, nodes, nodelookup = make_ptree(tree)
    pstree, snodes, snodelookup = make_ptree(stree)
    recon2 = make_recon_array(tree, recon, nodes, snodelookup)
    events2 = make_events_array(nodes, events)

    dists = [x.dist for x in nodes]
    sdists = [x.dist for x in snodes]

    nnodes = len(nodes)
    nsnodes = len(snodes)

    sp_alpha = [params[x.name][0] for x in snodes]
    sp_beta = [params[x.name][1] for x in snodes]

    if generate is None:
        generate = -1
    
    p = branchPrior(nnodes, ptree, dists, nsnodes, pstree, sdists,
                    recon2, events2,
                    sp_alpha, sp_beta,
                    generate,
                    pretime_lambda, birth, death,
                    params["baserate"][0],
                    params["baserate"][1],
                    nsamples,
                    approx)
    return p


#=============================================================================
# sequence likelihood

def make_hky_matrix(bgfreq, kappa, t):
    """
    Returns a HKY matrix

    bgfreq -- the background frequency A,C,G,T
    kappa  -- is transition/transversion ratio
    """
    
    matrix = [0.0] * 16
    matrix = c_list(c_float, matrix)
    makeHkyMatrix(c_list(c_float, bgfreq), kappa, t, matrix)
    return [matrix[0:4],
            matrix[4:8],
            matrix[8:12],
            matrix[12:16]]

def make_hky_deriv_matrix(bgfreq, kappa, t):
    """
    Returns a HKY Derivative matrix

    bgfreq -- the background frequency A,C,G,T
    kappa  -- is transition/transversion ratio
    """
    
    matrix = [0.0] * 16
    matrix = c_list(c_float, matrix)
    makeHkyDerivMatrix(c_list(c_float, bgfreq), kappa, t, matrix)
    return [matrix[0:4],
            matrix[4:8],
            matrix[8:12],
            matrix[12:16]]

def make_hky_deriv2_matrix(bgfreq, kappa, t):
    """
    Returns a HKY 2nd Derivative matrix

    bgfreq -- the background frequency A,C,G,T
    kappa  -- is transition/transversion ratio
    """
    
    matrix = [0.0] * 16
    matrix = c_list(c_float, matrix)
    makeHkyDeriv2Matrix(c_list(c_float, bgfreq), kappa, t, matrix)
    return [matrix[0:4],
            matrix[4:8],
            matrix[8:12],
            matrix[12:16]]


def branch_likelihood_hky(probs1, probs2, seqlen, bgfreq, kappa, time):
    return branchLikelihoodHky(c_list(c_floatlk, probs1),
                               c_list(c_floatlk, probs2),
                               seqlen, c_list(c_float, bgfreq), kappa, time)

def branch_likelihood_hky_deriv(probs1, probs2, seqlen, bgfreq, kappa, time):
    return branchLikelihoodHkyDeriv(
        c_list(c_floatlk, probs1), c_list(c_floatlk, probs2), seqlen, 
        c_list(c_float, bgfreq), kappa, time)

def branch_likelihood_hky_deriv2(probs1, probs2, seqlen, bgfreq, kappa, time):
    return branchLikelihoodHkyDeriv2(
        c_list(c_floatlk, probs1), c_list(c_floatlk, probs2), seqlen, 
        c_list(c_float, bgfreq), kappa, time)

def mle_distance_hky(probs1, probs2, seqlen, bgfreq, kappa, t0, t1):    
    return mleDistanceHky(c_list(c_floatlk, probs1), c_list(c_floatlk, probs2),
                          seqlen, c_list(c_float, bgfreq), kappa,
                          t0, t1)


def calc_seq_likelihood_hky(tree, align, bgfreq, kappa):

    ptree, nodes, nodelookup = make_ptree(tree)
    leaves = [x for x in nodes if x.is_leaf()]
    calign = (c_char_p * len(align))(* [align[x.name] for x in leaves])    
    ctree = tree2ctree(tree)
    
    l = calcSeqProbHky(ctree, len(align), calign, c_list(c_float, bgfreq),
                       kappa)

    deleteTree(ctree)

    return l


def find_ml_branch_lengths_hky(tree, align, bgfreq, kappa, maxiter=20,
                               parsinit=True):

    ptree, nodes, nodelookup = make_ptree(tree)
    leaves = [x for x in nodes if x.is_leaf()]
    calign = (c_char_p * len(align))(* [align[x.name] for x in leaves])
    dists = c_list(c_float, [n.dist for n in nodes])

    l = findMLBranchLengthsHky(len(ptree), c_list(c_int, ptree),
                               len(align), calign,
                               dists, c_list(c_float, bgfreq), kappa,
                               maxiter, int(parsinit))
    
    for i, node in enumerate(nodes):
        node.dist = dists[i]
    
    return l


#=============================================================================
# training

def mean(vals):
    return sum(vals) / float(len(vals))
    
def variance(vals):
    """Variance"""
    u = mean(vals)
    return sum((x - u)**2 for x in vals) / float(len(vals)-1)


def read_length_matrix(filename, minlen=.0001, maxlen=1.0,
                       nooutliers=True):
    """Read a length matrix made by spidir-prep"""

    from rasmus import util

    dat = [line.rstrip().split("\t") for line in open(filename)]
    species = dat[0][2:]
    lens = util.map2(float, util.submatrix(dat, range(1, len(dat)),
                                           range(2, len(dat[0]))))
    gene_sizes = map(int, util.cget(dat[1:], 1))
    files = util.cget(dat[1:], 0)

    if nooutliers:
        treelens = map(sum, lens)
        m = mean(treelens)
        ind = util.find(lambda x: x<5*m, treelens)
        files, gene_sizes, lens, treelens = [util.mget(x, ind) for x in
                                             files, gene_sizes, lens, treelens]



    for row in lens:
        for i in xrange(len(row)):
            if row[i] < minlen:
                row[i] = minlen

    
    return species, lens, gene_sizes, files



def train_params(gene_sizes, length_matrix, times, species,
                 nrates=10, max_iter=10):

    ntrees = len(length_matrix)
    nspecies = len(length_matrix[0])

    sp_alpha = [1.0] * nspecies
    sp_beta = [1.0] * nspecies
    gene_alpha = [1.0]
    gene_beta = [1.0]

    train(ntrees, nspecies, gene_sizes, length_matrix,
          times,
          sp_alpha, sp_beta, gene_alpha, gene_beta,
          nrates, max_iter)

    params = {}
    params["baserate"] = [gene_alpha[0], gene_beta[0]]
    for i, sp in enumerate(species):
        params[sp] = [sp_alpha[i], sp_beta[i]]

    return params

def alloc_rates_em(gene_sizes, length_matrix, times, species, nrates):

    ntrees = len(length_matrix)
    nspecies = len(length_matrix[0])

    assert len(gene_sizes) == ntrees
    assert len(times) == nspecies
    assert len(species) == nspecies

    sp_alpha = [1.0] * nspecies
    sp_beta = [1.0] * nspecies
    gene_alpha = 1.0
    gene_beta = 1.0
    
    return allocRatesEM(ntrees, nspecies, nrates,
                        c_list(c_int, gene_sizes),
                        c_matrix(c_float, length_matrix),
                        c_list(c_float, times),
                        c_list(c_float, sp_alpha),
                        c_list(c_float, sp_beta), 
                        gene_alpha, gene_beta)


def free_rates_em(em):
    freeRatesEM(em)


def rates_em_get_params(em, species, stree=None, lens=None, times=None):
    c_params = c_list(c_float, [0.0] * (2*len(species) + 2))    
    RatesEM_getParams(em, c_params)

    params = {"baserate": [c_params[0], c_params[1]]}
    for i, sp in enumerate(species):
        if isinstance(sp, basestring) and sp.isdigit():
            sp = int(sp)
        params[sp] = c_params[2+2*i:4+2*i]

    if stree and lens and times:
        treelens = map(sum, lens)
        m = mean(treelens)
        grates = [i/m for i in treelens]
        rates = [i/j/g for row, g in zip(lens, grates)
                 for i, j in zip(row, times)]
        mu = mean(rates)
        sigma2 = variance(rates)
        params[stree.root.name] = [mu*mu/sigma2, mu/sigma2]

    return params
    

