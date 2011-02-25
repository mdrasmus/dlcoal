/*=============================================================================

  Matt Rasmussen
  Copyright 2007-2011

  SPIMAP model parameters

=============================================================================*/


#ifndef SPIDIR_MODEL_PARAMS_H
#define SPIDIR_MODEL_PARAMS_H

#include <string>

using namespace std;

namespace spidir {


// forward declaration
class SpeciesTree;

// spidir parameters
class SpidirParams
{
public:
    SpidirParams(int size, string *_names, 
                 float *_sp_alpha, float *_sp_beta, 
		 float _gene_alpha, float _gene_beta,
		 float _pretime_lambda=1.0) :
        nsnodes(size),
        gene_alpha(_gene_alpha),
        gene_beta(_gene_beta),
	pretime_lambda(_pretime_lambda)
    {
        names = new string [nsnodes];
        sp_alpha = new float [nsnodes];
        sp_beta = new float [nsnodes];
        
        for (int i=0; i<nsnodes; i++) {
            if (_names)
                names[i] = _names[i];
            sp_alpha[i] = _sp_alpha[i];
            sp_beta[i] = _sp_beta[i];
        }
    }
    
    ~SpidirParams()
    {
        delete [] names;
        delete [] sp_alpha;
        delete [] sp_beta;
    }
    
    // sorts the parameters to match newnames
    bool order(SpeciesTree *tree);    

    int nsnodes;
    string *names;
    float *sp_alpha;
    float *sp_beta;
    float gene_alpha;
    float gene_beta;
    float pretime_lambda;
};


SpidirParams *readSpidirParams(const char* filename);


} // namespace spidir

#endif
