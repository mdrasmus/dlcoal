/*=============================================================================

  Matt Rasmussen
  Copyright 2007-2011

  Parsimony phylogeny algorithm

=============================================================================*/



#ifndef SPIDIR_PARSIMONY_H
#define SPIDIR_PARSIMONY_H

#include "Tree.h"

namespace spidir {

extern "C" {

void parsimony(int nnodes, int *ptree, int nseqs, char **seqs, float *dists,
               bool buildAncestral=false, char **ancetralSeqs=NULL);

} // extern "C"

void parsimony(Tree *tree, int nseqs, char **seqs,
               bool buildAncestral=false, char **ancetralSeqs=NULL);

} // namespace spidir

#endif
