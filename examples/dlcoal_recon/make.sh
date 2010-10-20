#!/usr/bin/bash


# make simulations
#   Creates dir 'data/flies' with subdirs 0-499.  Each subdir is one simulated 
#   gene family.  The following files will be created:
#      X.coal.tree   -- coalescence tree (newick)
#      X.coal.recon  -- reconciliation from coal tree to locus tree
#      X.locus.tree  -- locus tree
#      X.daughters   -- daughter edges in locus tree
#      X.locus.recon -- reconciliation from locus tree to species tree
dlcoal_sim -i 1 -s config/flies.stree -S config/flies.smap \
    -n 1e7 -D .0012 -L .0011 -g .1 -o data/flies


# view a simulation
view_recon -s config/flies.stree -g .1 data/flies/0/0


# get real relations
#   creates the file 'data/flies.rel.txt' which contains every duplication,
#   loss, and ortholog relationship in the 500 trees
phylofiles -e data/flies .locus.tree | tree-relations \
    -d -s config/flies.stree -S config/flies.smap \
    -T .tree -R .recon > data/flies.rel.txt




#=============================================================================
# run mpr
#   this will create the files
#      X.mpr.recon -- reconciliation from coal tree to species tree
phylofiles -e data/flies .coal.tree | (while read x; do
echo $x
mpr -s config/flies.stree -S config/flies.smap \
  -I .coal.tree -O .mpr $x
done)

# test to see if MPR completed (should return 500 files)
phylofiles -e data/flies .mpr.recon | wc -l


# get mpr relations and compare them to the truth
#   results will be written to 'data/flies.mpr.rel-summary.txt'
phylofiles -e data/flies .coal.tree | tree-relations \
    -d -s config/flies.stree -S config/flies.smap \
    -T .coal.tree -R .mpr.recon > data/flies.mpr.rel.txt
tree-relations-cmp data/flies.rel.txt data/flies.mpr.rel.txt > \
    data/flies.mpr.rel-summary.txt
cat data/flies.mpr.rel-summary.txt



#=============================================================================
# run dlcoal_recon
#   this will create the files
#      X.dlcoal.coal.tree   -- a copy of X.coal.tree
#      X.dlcoal.coal.recon  -- inferred recon from coal to locus tree
#      X.dlcoal.locus.tree  -- inferred locus tree
#      X.dlcoal.daughters   -- inferred daughters (not implemented yet)
#      X.dlcoal.locus.recon -- inferred recon from locus to species tree
phylofiles -n data/flies .dlcoal.locus.recon | (while read x; do
    x=${x/.dlcoal.locus.recon/.coal.tree}
    echo $x
    dlcoal_recon -i 100 --nsamples 1 \
        -s config/flies.stree -S config/flies.smap \
        -n 1e7 -D .0012 -L .0011 -g .1 -I .coal.tree -O .dlcoal $x
done)

# test to see if DLCoalRecon completed (should report 500 files)
phylofiles -e data/flies .dlcoal.locus.recon | wc -l

# view a dlcoal reconstruction
view_recon -s config/flies.stree -g .1 data/flies/0/0.dlcoal

# run dlcoal once
dlcoal_recon -i 100 --nsamples 1 \
        -s config/flies.stree -S config/flies.smap \
        -n 1e7 -D .0012 -L .0011 -g .1 -I .coal.tree -O .dlcoal \
    data/flies/0/0.coal.tree



# get dlcoal relations and compare them to the truth
#   results will be written to 'data/flies.dlcoal.rel-summary.txt'
phylofiles -e data/flies .dlcoal.locus.tree | tree-relations \
    -d -s config/flies.stree -S config/flies.smap \
    -T .tree -R .recon > data/flies.dlcoal.rel.txt
tree-relations-cmp data/flies.rel.txt data/flies.dlcoal.rel.txt > \
    data/flies.dlcoal.rel-summary.txt
cat data/flies.dlcoal.rel-summary.txt



#=============================================================================
# Additional info
#
# reconciliation file format:
# The format is tab-delimited with 3 columns (gene node, species node, event)
# 1. 'gene node' specifies a node in the gene tree (i.e. 1st tree)
# 2. 'species node' specifies a node in the species tree (i.e. 2nd tree)
# 3. 'event' specifies whether the gene node is an extant gene ('gene'),
#    duplication ('dup'), speciation ('spec'), or nothing ('none').

