#ifndef SPIDIR_SEQ_LIKELIHOOD_H
#define SPIDIR_SEQ_LIKELIHOOD_H


namespace spidir {


extern "C" {

void train(int ntrees, int nspecies, float **lengths, float *times,
           float *sp_alpha, float *sp_beta, float *gene_alpha, float *gene_beta,
           int nrates, int max_iter)

} // extern "C"


} // namespace spidir


