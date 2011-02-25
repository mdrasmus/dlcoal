/*=============================================================================

  Matt Rasmussen
  Copyright 2010-2011

  Estimating birth-death rates from gene counts (implementing Hahn et al 2005)

=============================================================================*/


#ifndef SPIDIR_BIRTHDEATH_ML_H
#define SPIDIR_BIRTHDEATH_ML_H

#include "Tree.h"

namespace spidir {

extern "C" {


double birthDeathTreeCounts(Tree *tree, int nspecies, int *counts, 
                            float birth, float death, int maxgene,
                            int rootgene, double **tab=NULL);

double birthDeathForestCounts(Tree *tree, int nspecies, int nfams,
                              int **counts, int *mult,
                              float birth, float death, int maxgene,
                              int rootgene, double **tab=NULL);


}

} // namespace spidir

#endif // SPIDIR_BIRTHDEATH_ML_H
