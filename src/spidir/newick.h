/*=============================================================================

  Matt Rasmussen
  Copyright 2007-2011

  Newick tree reading/writing

=============================================================================*/

#ifndef SPIDIR_NEWICK_H
#define SPIDIR_NEWICK_H

#include <stdlib.h>
#include <string>

#include "Tree.h"

namespace spidir {

Tree *readNewickTree(FILE *infile, Tree *tree=NULL);
Tree *readNewickTree(const char *filename, Tree *tree=NULL);

void writeNewickNode(FILE *out, Node *node, int depth, bool oneline=false);
void writeNewickTree(FILE *out, Tree *tree, int depth, bool oneline=false);
bool writeNewickTree(const char *filename, Tree *tree, bool oneline=false);


} // namespace spidir


#endif // SPDIR_NEWICK_H

