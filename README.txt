SPIMAP (Spieces informed max a poseriori)
http://compbio.mit.edu/spimap/
Matthew Rasmussen

=============================================================================
ABOUT

SPIMAP is a phylogenetic program that uses a species tree to aide in
reconstructing gene trees.  It uses code from the SPIDIR phylogenetic library.

SPIMAP citation: 
Rasmussen, Kellis.  A Bayesian Approach for Fast and Accurate Gene-tree 
Reconstruction. Mol Biol Evol. 2011 Jan;28(1):273-90. Epub 2010 Jul 25.

SPIDIR citation:
Rasmussen, Kellis. Accurate gene-tree reconstruction by learning
gene- and species-specific substitution rates across multiple complete genomes.
Genome Research. 2007

This package includes the C++ source of the SPIMAP program and SPIDIR library
as well as several library interfaces for C and python.


=============================================================================
DEPENDENCIES

SPIMAP has the following requirements:

- GNU Scientific library (GSL) http://www.gnu.org/software/gsl/
- Python (2.5 or greater) http://python.org/


=============================================================================
INSTALL

NOTE: Makefile installation will work best on UNIX or CYGWIN.


To compile the SPIMAP stand-alone program use the Makefile.

    make

To compile the SPIDIR C-library use:
    
    make lib  

Once compiled, to install the SPIMAP program (installs by default in /usr) use:

    make install

To specify your own installation path use:
    
    make install prefix=/usr/local

To use the training scripts (bin/spimap-prep-duploss, etc), python must be
installed.


=============================================================================
USAGE

Running spimap with no arguments will print out its command-line usage:

Usage: spimap [OPTION]

  -a,--align  <alignment fasta>
    sequence alignment in fasta format

  -S,--smap  <species map>
    gene to species map

  -s,--stree  <species tree>
    species tree file in newick format

  -p,--param  <params file>
    substitution rate parameters file

  -o,--output  <output filename prefix>
    prefix for all output filenames

Sequence evolution model
  -k,--kappa  <transition/transversion ratio>
    used for HKY model (default=estimate)

  -f,--bgfreq  <A freq>,<C ferq>,<G freq>,<T freq>
    background frequencies (default: estimate)

Dup/loss evolution model
  -D,--duprate  <duplication rate>
    rate of a gene duplication (default=0.1)

  -L,--lossrate  <loss rate>
    probability of loss (default=0.1)

  -P,--pretime  <pre-speciation time parameter>
    lambda param of pre-speciation distribution (default=1.0)

Search
  -i,--niter  <# iterations>
    number of iterations

  --quickiter  <quick iterations>
    number of subproposals (default=50)

  -b,--boot  <# bootstraps>
    number of bootstraps to perform (default: 1)

Information
  -V,--verbose  <verbosity level>
    verbosity level 0=quiet, 1=low, 2=medium, 3=high

  --log  <log filename>
    log filename.  Use '-' to display on stdout.

  -v,--version  
    display version information

  -h,--help  
    display help information

  --help-debug  
    display help information about debug options


#=============================================================================
# Examples

See examples/analyze-fungi.sh for an example of how to use each program
in the SPIMAP package.

