/*=============================================================================

  Matt Rasmussen
  Copyright 2007-2011

  Train parameters of branch length prior in SPIMAP model

=============================================================================*/


#ifndef SPIDIR_TRAIN_H
#define SPIDIR_TRAIN_H


namespace spidir {


extern "C" {

void train(int ntrees, int nspecies, float **lengths, float *times,
           float *sp_alpha, float *sp_beta, float *gene_alpha, float *gene_beta,
           int nrates, int max_iter)

} // extern "C"


} // namespace spidir


