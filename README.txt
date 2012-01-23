DLCoal (duplications, losses, and coalescence)
http://compbio.mit.edu/dlcoal/
Matthew Rasmussen

=============================================================================
ABOUT

DLCoal is a software package containing the DLCoalRecon program as well as
several other useful utilities for working with gene trees.

DLCoalRecon is a reconciliation program that maps a gene tree to species tree
in order to determine how gene duplications and losses have occurred.  
DLCoalRecon's unique feature is that it can perform reconciliation despite
the presence of incomplete linegae sorting.

DLCoal citation: 
Rasmussen, Kellis.  A unified model of gene duplication, loss, and coalescence using a locus tree. Genome Research. 2012.


=============================================================================
DEPENDENCIES

DLCoal has the following requirements:

- GNU Scientific library (GSL) http://www.gnu.org/software/gsl/
- Python (2.5 or greater) http://python.org/


=============================================================================
INSTALL

NOTE: Makefile installation will work best on UNIX or CYGWIN (Windows).


To compile the DLCoal library use the Makefile.

    make

Once compiled, to install the DLCoal programs (default install in /usr) use:

    make install

To specify your own installation path use:
    
    make install prefix=/usr/local

DLCoal can also run directly from the source directory.  Simply add the
bin/ directory to your PATH or create symlinks to the scripts within bin/
to any directory on your PATH.


=============================================================================
USAGE

Running dlcoal_recon with no arguments will print out its command-line usage:


Usage: dlcoal_recon [options] GENE_TREE1 [GENE_TREE2 ...]

Options:
  -h, --help            show this help message and exit
  -s SPECIES_TREE, --stree=SPECIES_TREE
                        species tree file in newick format (myr)
  -S GENE_TO_SPECIES_MAP, --smap=GENE_TO_SPECIES_MAP
                        gene to species map
  -n POPULATION_SIZE, --popsize=POPULATION_SIZE
                        Effective population size
  -D DUPLICATION_RATE, --duprate=DUPLICATION_RATE
                        rate of a gene duplication (dups/gene/myr)
  -L LOSS_RATE, --lossrate=LOSS_RATE
                        rate of gene loss (losses/gene/myr)
  -g GENRATION_TIME, --gentime=GENRATION_TIME
                        generation time (years)
  -i ITERATIONS, --iter=ITERATIONS
                        number of search iterations

  File extensions:
    -I INPUT_EXT, --inext=INPUT_EXT
                        input file extension (default='')
    -O OUTPUT_EXT, --outext=OUTPUT_EXT
                        output file extension (default='.dlcoal')

  Miscellaneous:
    --nprescreen=NUM_PRESCREENS
                        number of prescreening iterations
    --nsamples=NUM_SAMPLES
                        number of samples for dup-loss integration
                        (default=100)
    --init-locus-tree=TREE_FILE
                        initial locus tree for search
    -x RANDOM_SEED, --seed=RANDOM_SEED
                        random number seed
    -l, --log           if given, output debugging log


#=============================================================================
# Examples

See examples/make.sh for an example of how to use each program
in the DLCoal package.


#=============================================================================
# Documentation

See doc/dlcoal-manual.html for further documentation of the software and
its associated file-formats.

