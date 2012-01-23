#!/usr/bin/bash
#
# This is an example of how to use DLCoal to reconcile gene families.
#
# In this example, we reconcile several simulated fly gene families
#
# Don't execute this script all at once.  Instead try copying and pasting
# each command to the command line one at a time in order to learn how the
# commands work.
#
# For documentation on all of DLCoal's file formats, see doc/dlcoal_manual.html
#

#=============================================================================
# setup/install

# Make sure tools are compiled and installed before running the commands in 
# this tutorial.  This can be done with this command:
cd ..
make install
cd examples

# or you can install spimap into the prefix of your choice
cd ..
make install prefix=YOUR_INSTALL_PREFIX
cd examples

# or you can run from the source directory by setting these environment 
# variables:
cd ..
make
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:../lib
cd examples


#=============================================================================
# make simulation data 

mkdir -p data/flies
mkdir -p data/flies-rel


#   Creates dir 'data/flies' with subdirs 0-499.  Each subdir is one simulated 
#   gene family.  The following files will be created:
#      X.coal.tree   -- coalescence tree (newick)
#      X.coal.recon  -- reconciliation from coal tree to locus tree
#      X.locus.tree  -- locus tree
#      X.daughters   -- daughter edges in locus tree
#      X.locus.recon -- reconciliation from locus tree to species tree
../bin/dlcoal_sim -i 100 -s config/flies.stree  \
    -n 1e7 -D .0012 -L .0011 -g .1 -o data/flies

# or you can untar a pre-made simulated dataset
tar -C data -xvf data/flies.tar.gz


# view a simulation
#  From top to bottom: species tree, locus tree, gene tree.
../bin/view_recon -s config/flies.stree -g .1 data/flies/0/0



# get real relations
#   creates the file 'data/flies.rel.txt' which contains every duplication,
#   loss, and ortholog relationship in the 500 trees
../bin/phylofiles -e data/flies .locus.tree | ../bin/tree-relations \
    -d -s config/flies.stree -S config/flies.smap \
    -T .tree -R .recon > data/flies-rel/flies.rel.txt



#=============================================================================
# run MPR (maximum parsimonious reconciliation)


#   this will create the files
#      X.mpr.recon -- reconciliation from coal tree to species tree
../bin/phylofiles -e data/flies .coal.tree | (while read x; do
echo $x
../bin/mpr -s config/flies.stree -S config/flies.smap \
  -I .coal.tree -O .mpr $x
done)

# test to see if MPR completed (should return 100 files)
../bin/phylofiles -e data/flies .mpr.recon | wc -l


# get mpr relations and compare them to the truth
#   results will be written to 'data/flies.mpr.rel-summary.txt'
../bin/phylofiles -e data/flies .coal.tree | ../bin/tree-relations \
    -d -s config/flies.stree -S config/flies.smap \
    -T .coal.tree -R .mpr.recon > data/flies-rel/flies.mpr.rel.txt
../bin/tree-relations-cmp data/flies-rel/flies.rel.txt \
    data/flies-rel/flies.mpr.rel.txt > \
    data/flies-rel/flies.mpr.rel-summary.txt
cat data/flies-rel/flies.mpr.rel-summary.txt


# you should see something like:
#
# dup actual:    46
# dup pred:      147
# dup sn:        0.891304347826
# dup ppv:       0.278911564626
# loss actual:   34
# loss pred:     339
# loss sn:       0.941176470588
# loss ppv:      0.094395280236
# orth actual:    6856
# orth pred:      5755
# orth sn:        0.839410735123
# orth ppv:       1.0


#=============================================================================
# run dlcoal_recon

# example of running dlcoal for one gene family
../bin/dlcoal_recon -i 100 --nsamples 1 \
        -s config/flies.stree -S config/flies.smap \
        -n 1e7 -D .0012 -L .0011 -g .1 -I .coal.tree -O .dlcoal \
    data/flies/0/0.coal.tree

# example of running dlcoal on all 100 gene trees.
# this will create the files:
#      X.dlcoal.coal.tree   -- a copy of X.coal.tree
#      X.dlcoal.coal.recon  -- inferred recon from coal to locus tree
#      X.dlcoal.locus.tree  -- inferred locus tree
#      X.dlcoal.daughters   -- inferred daughters (not implemented yet)
#      X.dlcoal.locus.recon -- inferred recon from locus to species tree
../bin/phylofiles -n data/flies .dlcoal.locus.recon | (while read x; do
    x=${x/.dlcoal.locus.recon/.coal.tree}
    echo $x
    ../bin/dlcoal_recon -i 100 --nsamples 1 \
        -s config/flies.stree -S config/flies.smap \
        -n 1e7 -D .0012 -L .0011 -g .1 -I .coal.tree -O .dlcoal $x
done)

# test to see if DLCoalRecon completed (should report 100 files)
../bin/phylofiles -e data/flies .dlcoal.locus.recon | wc -l



# get dlcoal relations and compare them to the truth
#   results will be written to 'data/flies.dlcoal.rel-summary.txt'
../bin/phylofiles -e data/flies .dlcoal.locus.tree | ../bin/tree-relations \
    -d -s config/flies.stree -S config/flies.smap \
    -T .tree -R .recon > data/flies-rel/flies.dlcoal.rel.txt
../bin/tree-relations-cmp data/flies-rel/flies.rel.txt \
    data/flies-rel/flies.dlcoal.rel.txt > \
    data/flies-rel/flies.dlcoal.rel-summary.txt
cat data/flies-rel/flies.dlcoal.rel-summary.txt

# you should see something like:
#
# dup actual:     46
# dup pred:       49
# dup sn:         0.95652173913
# dup ppv:        0.897959183673
# loss actual:    34
# loss pred:      43
# loss sn:        0.970588235294
# loss ppv:       0.767441860465
# orth actual:    6856
# orth pred:      6758
# orth sn:        0.985414235706
# orth ppv:       0.999704054454




