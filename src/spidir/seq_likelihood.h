/*=============================================================================

  Matt Rasmussen
  Copyright 2007-2011

  Computing sequence likelihood
  Maximum Likelihood Branch Length Estimation

=============================================================================*/


#ifndef SPIDIR_SEQ_LIKELIHOOD_H
#define SPIDIR_SEQ_LIKELIHOOD_H


#include "Tree.h"

namespace spidir {

typedef double floatlk;

double findMLBranchLengthsHky(Tree *tree, int nseqs, char **seqs, 
                              const float *bgfreq, float kappa, 
                              int maxiter=100, 
                              double minlen=.0001, double maxlen=10);

template <class Model>
floatlk getTotalLikelihood(ExtendArray<floatlk*> &lktable, Tree *tree, 
                           int nseqs, int seqlen, char **seqs, Model &model,
                           const float *bgfreq);

extern "C" {

void makeHkyMatrix(const float *bgfreq, float ratio, float t, float *matrix);

void makeHkyDerivMatrix(const float *bgfreq, float ratio, float t, float *matrix);

floatlk branchLikelihoodHky(floatlk *probs1, floatlk *probs2, int seqlen, 
                            const float *bgfreq, float kappa, float t);

floatlk branchLikelihoodHkyDeriv(floatlk *probs1, floatlk *probs2, int seqlen, 
                                 const float *bgfreq, float kappa, float t);

floatlk branchLikelihoodHkyDeriv2(floatlk *probs1, floatlk *probs2, 
                                  int seqlen, 
                                  const float *bgfreq, float kappa, float t);

floatlk mleDistanceHky(floatlk *probs1, floatlk *probs2, int seqlen, 
                       const float *bgfreq, float kappa,
                       float t0, float t1);

floatlk calcSeqProbHky(Tree *tree, int nseqs, char **seqs, 
                     const float *bgfreq, float kappa);

floatlk findMLBranchLengthsHky(int nnodes, int *ptree, int nseqs, char **seqs, 
                               float *dists, const float *bgfreq, float kappa, 
                               int maxiter, bool parsinit=false);

double findMLKappaHky(Tree *tree, int nseqs, char **seqs, 
                      const float *bgfreq, float minkappa, float maxkappa,
                      float kappastep);

} // extern "C"

} // namespace spidir

#endif // SPIDIR_SEQ_LIKELIHOOD_H
