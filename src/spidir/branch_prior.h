/*=============================================================================

  Matt Rasmussen
  Copyright 2010-2011

  Branch length prior of the SPIMAP model

=============================================================================*/


#ifndef SPIDIR_BRANCH_PRIOR_H
#define SPIDIR_BRANCH_PRIOR_H

#include "Tree.h"
#include "model_params.h"
#include "birthdeath.h"

namespace spidir {

typedef void (*geneRateCallback) (float generate, Tree *tree, void *userdata);

extern "C" {



// fractional branches
enum {
    FRAC_NONE,
    FRAC_DIFF,
    FRAC_PARENT,
    FRAC_NODE,
    FRAC_ONE
};



// Branch distribution parameters for one branch
class BranchParams
{
public:
    BranchParams(float _alpha=-1.0, float _beta=-1.0) :
        alpha(_alpha),
        beta(_beta)
    {}
    
    bool isNull()
    {
        return alpha == -1.0;
    }
    
    float alpha;
    float beta;
};


class BranchPart
{
public:
    BranchPart(int species=-1, int frac=-1) :
	species(species),
	frac(frac)
    {}

    int species;
    int frac;
};

// Reconciliation parameters
class ReconParams
{
public:
    ReconParams(int nnodes, SpidirParams *params) :
        nnodes(nnodes),
	params(params),
        unfold(-1),
        unfolddist(0)
    {
	parts = new ExtendArray<BranchPart> [nnodes];	
        midpoints = new float [nnodes];
        
        freebranches = new bool [nnodes];
    }
    
    ~ReconParams()
    {
	delete [] parts;
        delete [] midpoints;
        
        delete [] freebranches;
    }
    
    int nnodes;
    SpidirParams *params;

    ExtendArray<BranchPart> *parts;

    float pretime;
    float *midpoints;

    bool *freebranches;
    int unfold;
    float unfolddist;
};




double branchPrior(int nnodes, int *ptree, float *dists,
                   int nsnodes, int *pstree, float *sdists,
                   int *recon, int *events,
                   float *sp_alpha, float *sp_beta, float generate, 
                   float predupprob=1.0, float dupprob=1.0, float lossprob=1.0,
                   float gene_alpha=0, float gene_beta=0,
                   int nsamples=1000, bool approx=true);

} // extern "C"

double branchPrior(Tree *tree,
                   SpeciesTree *stree,
                   int *recon, int *events, SpidirParams *params,
                   float generate, 
                   float predupprob, float dupprob, float lossprob,
                   int nsamples=1000, bool approx=true);

// get nodes in preorder (starting with given node)
//void getSubtree(int **ftree, int node, int *events, ExtendArray<int> *subnodes);

bool getSubtree(Node *node, int *events, ExtendArray<Node*> *subnodes);

BranchParams getBranchParams(int node, Tree *tree, ReconParams *reconparams);

// Reconcile a branch to the species tree
void reconBranch(int node, 
		 Tree *tree, SpeciesTree *stree, 
		 int *recon, int *events, 
                 SpidirParams *params,
                 ReconParams *reconparams);

void determineFreeBranches(Tree *tree, SpeciesTree *stree, 
                           int *recon, int *events, float generate,
                           int *unfold, float *unfolddist, bool *freebranches);

void setRandomMidpoints(int root, Tree *tree, SpeciesTree *stree,
                        Node **subnodes, int nsubnodes, 
                        int *recon, int *events, 
                        ReconParams *reconparams,
			float birth, float death);


//=============================================================================
// old code

float maxPosteriorGeneRate(Tree *tree, SpeciesTree *stree,
                           int *recon, int *events, SpidirParams *params);

float maxPosteriorGeneRate(int nnodes, int *ptree, float *dists,
                           int nsnodes, int *pstree, 
                           int *recon, int *events,
                           float *mu, float *sigma,
                           float alpha, float beta);


void samplePosteriorGeneRate(Tree *tree,
                             int nseqs, char **seqs,
                             const float *bgfreq, float ratio,
                             SpeciesTree *stree,
                             int *gene2species,
                             SpidirParams *params,
                             int nsamples,
                             geneRateCallback callback,
                             void *userdata);


void samplePosteriorGeneRate(Tree *tree,
                             int nseqs, char **seqs, 
                             const float *bgfreq, float ratio,
                             SpeciesTree *stree,
                             int *recon, int *events, SpidirParams *params,
                             int nsamples,
                             geneRateCallback callback,
                             void *userdata);


void generateBranchLengths(int nnodes, int *ptree, 
                           int nsnodes, int *pstree,
                           int *recon, int *events,
                           float *mu, float *sigma,
                           float alpha, float beta,
                           float *dists);

void generateBranchLengths(Tree *tree,
                           SpeciesTree *stree,
                           int *recon, int *events,
                           SpidirParams *params,
                           float generate=-1.0, 
                           int subnode=-1, int subchild=-1);





} // namespace spidir

#endif // SPIDIR_BRANCH_PRIOR_H
