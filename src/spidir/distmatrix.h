/*=============================================================================

  Matt Rasmussen
  Copyright 2007-2011

  Distance matrix functions

=============================================================================*/


#ifndef DISTMATRIX_H
#define DISTMATRIX_H

using namespace std;

namespace spidir {


void calcDistMatrix(int nseqs, int seqlen, char **seqs, float **distmat);
bool writeDistMatrix(const char *filename, int ngenes, float **dists, 
                     string *names);


} // namespace spidir

#endif // DISTMATRIX_H
