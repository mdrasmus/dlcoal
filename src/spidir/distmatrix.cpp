/*=============================================================================

  Matt Rasmussen
  Copyright 2010-2011

  Distance matrix functions

=============================================================================*/


// c++ headers
#include <stdlib.h>
#include <stdio.h>
#include <string>

// spidir headers
#include "distmatrix.h"

namespace spidir {


// calculate the pairwise distances between sequences
// NOTE: simple version implemented first
void calcDistMatrix(int nseqs, int seqlen, char **seqs, float **distmat)
{
    for (int i=0; i<nseqs; i++) {
        distmat[i][i] = 0.0;
    
        for (int j=i+1; j<nseqs; j++) {
            float changes = 0.0;
            int len = 0;
            
            for (int k=0; k<seqlen; k++) {
                if (seqs[i][k] != '-' && seqs[j][k] != '-') {
                    len++;                
                    if (seqs[i][k] != seqs[j][k]) {
                        changes += 1;
                    }
                }
            }
            
            if (len > 0) {
                distmat[i][j] = changes / len;
                distmat[j][i] = changes / len;
            } else {
                distmat[i][j] = 1.0;
                distmat[j][i] = 1.0;
            }
        }
    }
}


bool writeDistMatrix(const char *filename, int ngenes, float **dists, 
                     string *names)
{
    FILE *stream = NULL;
    
    if ((stream = fopen(filename, "w")) == NULL) {
        fprintf(stderr, "cannot open '%s'\n", filename);
        return false;
    }
    
    // print number of genes
    fprintf(stream, "%d\n", ngenes);
    
    for (int i=0; i<ngenes; i++) {
        fprintf(stream, "%s ", names[i].c_str());
        
        for (int j=0; j<ngenes; j++) {
            fprintf(stream, "%f ", dists[i][j]);
        }
        fprintf(stream, "\n");
    }
    
    fclose(stream);
    return true;
}


} // namespace spidir
