/*=============================================================================

  Matt Rasmussen
  Copyright 2007-2011

  Gene tree topology prior of the SPIMAP model

=============================================================================*/


#ifndef SPIDIR_TOP_PRIOR_H
#define SPIDIR_TOP_PRIOR_H

#include "Tree.h"

namespace spidir {

extern "C" {

void calcDoomTable(Tree *tree, float birth, float death, double *doomtable);

void getSpecSubtree(Node *node, Node *snode, int *recon, int *events,
                    ExtendArray<Node*> &nodes);

double birthDeathTreePrior(Tree *tree, Tree *stree, int *recon, 
                          int *events, float birth, float death,
                          double *doomtable);

double birthDeathTreePriorFull(Tree *tree, Tree *stree, int *recon, 
                              int *events, float birth, float death,
                              double *doomtable);



}

} // namespace spidir

#endif // SPIDIR_TOP_PRIOR_H
