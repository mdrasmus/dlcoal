#ifndef SPIDIR_H
#define SPIDIR_H

#include <string>

using namespace std;

namespace spidir {


// convert dna characters into standard numbers
extern const int dna2int[256];

// convert standard numbers to dna characters
extern const char *int2dna;

// base numbers
enum {
    DNA_A = 0,
    DNA_C = 1,
    DNA_G = 2,
    DNA_T = 3,
    DNA_PURINE,
    DNA_PRYMIDINE
};


// get the base type of a nucleotide
extern int dnatype[];



// compute background frequencies
void computeBgfreq(int nseq, char **seqs, float *bgfreq);

//=============================================================================
// Distance Matrices

void calcDistMatrix(int nseqs, int seqlen, char **seqs, float **distmat);
bool writeDistMatrix(const char *filename, int ngenes, float **dists, 
                     string *names);





//=============================================================================
// SPIDIR parameters


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
