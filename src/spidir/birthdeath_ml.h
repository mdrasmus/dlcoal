#ifndef SPIDIR_BIRTHDEATH_ML_H
#define SPIDIR_BIRTHDEATH_ML_H

#include "Tree.h"
#include "spidir.h"

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
