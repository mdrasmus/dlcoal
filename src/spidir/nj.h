/*=============================================================================

  Matt Rasmussen
  Copyright 2007-2011

  Neighbor-joining algorithm

=============================================================================*/


#ifndef SPIDIR_NJ_H
#define SPIDIR_NJ_H


namespace spidir {

void neighborjoin(int ngenes, float **distmat, int *ptree, float *branches);

} // namespace spidir

#endif // SPIDIR_NJ_H
