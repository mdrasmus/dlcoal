/*=============================================================================

  Matt Rasmussen
  Copyright 2007-2011

  Visualize trees

=============================================================================*/

#ifndef SPIDIR_TREEVIS_H
#define SPIDIR_TREEVIS_H

#include <stdlib.h>
#include <string>


using namespace std;

namespace spidir {

void displayTree(Tree *tree, FILE *outfile=stdout, 
                 float xscale=20.0, int yscale=2);
void displayTreeMatrix(Tree *tree, float xscale, int yscale, 
                       char ***matrix, int *nrows, int *ncols);


} // namespace spidir


#endif // SPDIR_TREEVIS_H

