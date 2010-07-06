#!/usr/bin/env python
# make simulated data for the 12 flies

# imports
import os
from rasmus import treelib, util
from compbio import dlcoal


#=============================================================================
# parameters

# species tree and mapping function (gene name -> species name)
stree = treelib.read_tree("config/flies.stree")
def gene2species(gene):
    return gene.split("_")[0]


gen_per_myr = 1e6 * 10.0 # generations per million years

# (convert branch lengths to generations)
for node in stree:
    node.dist *= gen_per_myr

# effective population size
n = int(1e7) * 2

# dup/loss rates expressed as events per generation
duprate = .0012 / gen_per_myr
lossrate = .0011 / gen_per_myr


#=============================================================================
# simulate
outdir = "data/flies"
dlcoal.dlcoal_sims(outdir, 500, stree, n, duprate, lossrate)



