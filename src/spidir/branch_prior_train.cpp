/*=============================================================================

  Matt Rasmussen
  Copyright 2007-2011

  Train parameters of branch length prior in SPIMAP model

=============================================================================*/

// c++ headers
#include <math.h>
#include <time.h>

// 3rd party
#include <gsl/gsl_roots.h>
#include <gsl/gsl_multimin.h>

// spidir headers
#include "common.h"
#include "gamma.h"
#include "Matrix.h"
#include "parsimony.h"
#include "seq_likelihood.h"
#include "Tree.h"


namespace spidir {



class RatesEM
{
public:

    RatesEM(int ntrees, int nspecies, int nrates, 
            int *gene_sizes,
            float **lengths, float *times,
            float *sp_alpha, float *sp_beta, 
            float gene_alpha, float gene_beta) : 
        ntrees(ntrees),
        nspecies(nspecies), 
        gene_sizes(gene_sizes),
        lengths(lengths),
        times(times),
        sp_alpha(sp_alpha),
        sp_beta(sp_beta),
        gene_alpha(gene_alpha),
        gene_beta(gene_beta),
        nrates(nrates),
        gtab(ntrees, nrates),
        pgtab(ntrees, nrates)
    {
        // allocate gene rate optimizer
	sol_gene_rate = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);

        // setup optimizer for gene rates
	opt_gene_rate.function = &gene_rate_f;
        opt_gene_rate.params = this; 

        // setup optimizer for species rates
        const int ndim = 2;
        sol_sp_rate = gsl_multimin_fdfminimizer_alloc(
            gsl_multimin_fdfminimizer_vector_bfgs2, ndim);
        opt_sp_rate.f = &sp_rate_f;
        opt_sp_rate.df = &sp_rate_df;
        opt_sp_rate.fdf = &sp_rate_fdf;
        opt_sp_rate.n = ndim;
        opt_sp_rate.params = this;
    }


    ~RatesEM()
    {
	//gsl_root_fdfsolver_free(sol_gene_rate);
	gsl_root_fsolver_free(sol_gene_rate);
        gsl_multimin_fdfminimizer_free(sol_sp_rate);
    }

    // gene rates function and derivative
  
    static double gene_rate_f(double x, void *params)
    {
        RatesEM *em = (RatesEM*) params;
        double gene_nu = x;

        // clamp gamma params
        if (gene_nu < .001)
            gene_nu = .001;

        double b = gene_nu;
        double a = gene_nu + 1;

        double sum = 0.0;
        for (int j=0; j<em->ntrees; j++) {
            for (int k=0; k<em->nrates; k++) {
                double g = em->gtab[j][k];
                double igamma = invgammaPdf(g, a, b);
                double dgamma = invgammaDerivG(g, b);
                
                if (isinf(dgamma) || isnan(dgamma) ||
                    isinf(igamma) || isnan(igamma)) {
                    fprintf(stderr, "igamma=%f, dgamma=%f, b=%f\n", 
                            igamma, dgamma, b);
                }

                if (!isnan(igamma) && igamma != 0.0) {
                    sum += em->pgtab[j][k] * dgamma / igamma;
                }
            }
        }
        assert(!isnan(sum) && !isinf(sum));

        return sum;
    }

    static double gene_rate_df(double x, void *params)
    {
        RatesEM *em = (RatesEM*) params;
        double gene_nu = x;
        // clamp gamma params
        if (gene_nu < .001)
            gene_nu = .001;
        
        double b = gene_nu;
        double a = gene_nu + 1;

        double sum = 0.0;
        for (int j=0; j<em->ntrees; j++) {
            for (int k=0; k<em->nrates; k++) {
                double g = em->gtab[j][k];
                double igamma = invgammaPdf(g, a, b);
		double dgamma = invgammaDerivG(g, b);
                double dgamma2 = invgammaDerivG2(g, b);

                if (!isnan(igamma) && igamma != 0.0) {
                    sum += em->pgtab[j][k] * 
                        (igamma * dgamma2 - dgamma*dgamma)
			 / (igamma*igamma);
                }
            }
        }

        assert(!isnan(sum) && !isinf(sum));

        return sum;
    }


    static void gene_rate_fdf(double x, void *params, 
                              double *f, double *df)
    {
        RatesEM *em = (RatesEM*) params;
        double gene_nu = x;

        // clamp gamma params
        if (gene_nu < .001)
            gene_nu = .001;

        double b = gene_nu;
        double a = gene_nu + 1;
        
        double sum = 0.0;
        double sum2 = 0.0;
        for (int j=0; j<em->ntrees; j++) {
            for (int k=0; k<em->nrates; k++) {
                double g = em->gtab[j][k];
                double igamma = invgammaPdf(g, a, b);
		double dgamma = invgammaDerivG(g, b);
                
                // TODO: substitute a better test
                if (igamma != 0.0) {
		    sum += em->pgtab[j][k] * 
                        invgammaDerivG(g, b) / igamma;
		    sum2 += em->pgtab[j][k] * 
                        (igamma * invgammaDerivG2(g, b) - dgamma*dgamma)
			 / (igamma*igamma);
                }
            }
        }

        assert(!isnan(sum) && !isinf(sum));
        assert(!isnan(sum2) && !isinf(sum2));

        // set return
        *f = sum;
        *df = sum2;
    }
    
    
    
    //======================================================
    // species rates function and derivative

    // f(a_i, b_i) = - sum_j sum_k pgtab_jk 
    //                       log(NB(c_ij, a_i, b_i  / (g_jk d_j t_i + b_i)))
    static double sp_rate_f(const gsl_vector *x, void *params)
    {
        RatesEM *em = (RatesEM*) params;
        double sp_alpha_i = gsl_vector_get(x, 0);
        double sp_beta_i = gsl_vector_get(x, 1);
        int i = em->cur_species;

        double sum = 0.0;
        for (int j=0; j<em->ntrees; j++) {
            // skip nulls
            if (em->lengths[j][i] < 0)
                continue;

            for (int k=0; k<em->nrates; k++) {
                double lngamma = gammalog(em->lengths[j][i], sp_alpha_i, 
                                sp_beta_i / (em->gtab[j][k] *
                                             em->times[i]));
                if (!isnan(lngamma))
                    sum += em->pgtab[j][k] * lngamma;
            }
        }

        return -sum;
    }
    
    // d/d a_i f(a_i, b_i) = - sum_j sum_k pgtab_jk 
    //                (NB'_r(c_ij, a_i, b_i  / (g_jk d_j t_i + b_i)) /
    //                    NB(c_ij, a_i, b_i  / (g_jk d_j t_i + b_i)))
    static void sp_rate_df(const gsl_vector *x, void *params, gsl_vector *df)
    {
        RatesEM *em = (RatesEM*) params;
        double sp_alpha_i = gsl_vector_get(x, 0);
        double sp_beta_i = gsl_vector_get(x, 1);
        int i = em->cur_species;

        double alpha_sum = 0.0;
        double beta_sum = 0.0;
        for (int j=0; j<em->ntrees; j++) {
            double l = em->lengths[j][i];
            // skip nulls
            if (l < 0)
                continue;


            for (int k=0; k<em->nrates; k++) {
                double g = em->gtab[j][k];
                double bgt = sp_beta_i / (g * em->times[i]);
                double lngamma = gammalog(l, sp_alpha_i, bgt);
                double gamma = exp(lngamma);

                if (gamma != 0.0) {
                    alpha_sum += em->pgtab[j][k] * 
                        gammaDerivA(l, sp_alpha_i, bgt) / gamma;
                    beta_sum += em->pgtab[j][k] *
                        gammaDerivB(l, sp_alpha_i, bgt) /
                        (gamma * g * em->times[i]);
                }
            }
        }

        gsl_vector_set(df, 0, -alpha_sum);
        gsl_vector_set(df, 1, -beta_sum);
    }

    // f(a_i, b_i) = - sum_j sum_k pgtab_jk 
    //                       log(NB(c_ij, a_i, b_i  / (g_jk d_j t_i + b_i)))
    // d/d a_i f(a_i, b_i) = - sum_j sum_k pgtab_jk 
    //                (NB'_r(c_ij, a_i, b_i  / (g_jk d_j t_i + b_i)) /
    //                    NB(c_ij, a_i, b_i  / (g_jk d_j t_i + b_i)))
    static void sp_rate_fdf(const gsl_vector *x, void *params, 
                              double *f, gsl_vector *df)
    {
        RatesEM *em = (RatesEM*) params;
        double sp_alpha_i = gsl_vector_get(x, 0);
        double sp_beta_i = gsl_vector_get(x, 1);  
        int i = em->cur_species;
   
        double sum = 0.0;
        double alpha_sum = 0.0;
        double beta_sum = 0.0;
        for (int j=0; j<em->ntrees; j++) {
            double l = em->lengths[j][i];
            // skip nulls
            if (l < 0)
                continue;

            for (int k=0; k<em->nrates; k++) {
                double g = em->gtab[j][k];
                double bgt = sp_beta_i / (g * em->times[i]);
                double lngamma = gammalog(l, sp_alpha_i, bgt);
                double gamma = exp(lngamma);

                sum += em->pgtab[j][k] * lngamma;
                
                if (gamma != 0.0) {
                    alpha_sum += em->pgtab[j][k] *
                        gammaDerivA(l, sp_alpha_i, bgt) / gamma;
                    beta_sum += em->pgtab[j][k] * 
                        gammaDerivB(l, sp_alpha_i, bgt) /
                        (gamma * g * em->times[i]);
                }
            }
        }

        // set return
        *f = -sum;
        gsl_vector_set(df, 0, -alpha_sum);
        gsl_vector_set(df, 1, -beta_sum);
    }


    float likelihood()
    {
        double logl = 0.0;
	
        for (int j=0; j<ntrees; j++) {
            double sum = 0.0;
            for (int k=0; k<nrates; k++) {
                double prod = 0.0;
                for (int i=0; i<nspecies; i++) {
                    double l = lengths[j][i];
                    // skip nulls
                    if (l < 0)
                        continue;
                    prod += gammalog(l, sp_alpha[i], 
				     sp_beta[i] / (times[i] * gtab[j][k]));
                }
                sum += pgtab[j][k] * exp(prod);
            }
            logl += log(sum);
        }

        return logl;
    }

    //====================

    void init_params()
    {
        ExtendArray<float> vec(ntrees);
        ExtendArray<float> treelens(ntrees);

        // compute tree lengths and their mean
        float meantreelens = 0.0;
        for (int j=0; j<ntrees; j++) {
            treelens[j] = 0.0;
            for (int i=0; i<nspecies; i++)
                if (lengths[j][i] > 0)
                    treelens[j] += lengths[j][i];
            meantreelens += treelens[j];
        }
        meantreelens /= ntrees;
        
        // compute variance of tree lengths
        float vartreelens = variance(treelens.get(), ntrees, meantreelens);

        // initialize gene rate parameter
        gene_nu = (1.0 / vartreelens) + 1.0;

        // initialize species rate parameters
        for (int i=0; i<nspecies; i++) {
            float sum = 0.0;
            for (int j=0; j<ntrees; j++) {
                if (lengths[j][i] > 0) {
                    vec[j] = lengths[j][i] / 
                        (treelens[j] / meantreelens) / times[i];
                    sum += vec[j];
                }
            }
            
            float mu = sum / ntrees;
            float sigma = stdev(vec.get(), ntrees);
            
            sp_alpha[i] = mu*mu / sigma / sigma;
            sp_beta[i] = mu / sigma / sigma;
        }
    }

    inline double gene_post(double g, double A, double B)
    {
        return loginvgammaPdf(g, A, B);
    }
    

    double find_upper_g(double m, double A, double B,
                        double tol1=.05, double tol2=.01)
    {        
        double fm = gene_post(m, A, B);
        double top = 2*m;
        double bot = m;
        tol2 = log(tol2);
        tol1 = log(tol1);

        // extent top
        while (true) {
            double ftop = gene_post(top, A, B);
            if (ftop - fm <= tol2)
                break;
            top *= 2;
        }

        // binary search
        while (true) {
            double u = (top + bot) / 2.0;
            double fu = gene_post(u, A, B);

            if (fu - fm > tol1)
                bot = u;
            else if (fu - fm < tol2)
                top = u;
            else
                return u;
        }
    }


    double find_lower_g(double m, double A, double B,
                       double tol1=.05, double tol2=.01)
    {
        double fm = gene_post(m, A, B);
        double top = m;
        double bot = 0;
        tol2 = log(tol2);
        tol1 = log(tol1);

        // binary search
        while (true) {
            double u = (top + bot) / 2.0;
            double fu = gene_post(u, A, B);

            if (fu - fm > tol1)
                top = u;
            else if (fu - fm < tol2)
                bot = u;
            else
                return u;
        }
    }


    // populate gene rate posteriors
    void EStep()
    {
        // temp variables for PDF of posterior gene rate
        float x[nrates+1];
        double y[nrates+1];

        // determine commonly used coefficients
        double A2 = gene_alpha;
        for (int i=0; i<nspecies; i++)
            A2 += sp_alpha[i];


        for (int j=0; j<ntrees; j++) {

            // determine commonly used coefficients
            double B = gene_beta;
            double A = A2;
            for (int i=0; i<nspecies; i++) {
                const double l = lengths[j][i];
                if (l > 0)                     
                    B += sp_beta[i] * l / times[i];
                else
                    A -= sp_alpha[i];
            }

            // find main range of gene rates
            double mid = B / (A+1);
            double top = find_upper_g(mid, A, B);
            double bot = find_lower_g(mid, A, B);

            int half_nrates = (nrates + 1) / 2;
            double step1 = (mid - bot) / half_nrates;
            double step2 = (top - mid) / (nrates + 1 - half_nrates);

            // compute x, y for posterior gene rate PDF
            for (int k=0; k<half_nrates; k++) {
                x[k] = bot + step1 * k;
                y[k] = invgammaPdf(x[k], A, B);
            }
            for (int k=half_nrates; k<nrates+1; k++) {
                x[k] = mid + step2 * (k-half_nrates);
                y[k] = invgammaPdf(x[k], A, B);
            }

            // compute gtab and pgtab
            for (int k=0; k<nrates; k++) {
                gtab[j][k] = (x[k] + x[k+1]) / 2.0;
                pgtab[j][k] = (y[k] + y[k+1]) * (x[k+1] - x[k]) / 2.0;
            }
        }
    }
    

    // maximize each model parameter given the hidden data estimated from
    // last iteration
    void MStep()
    {
        // optimization config
        double step_size = .01;
        double tol = .1;
        const double epsabs = .01;
        gsl_vector *init_x = gsl_vector_alloc(2);
        int status;
	double low = .0001, high = gene_nu * 20;
	int iter = 0;
	const int maxiter = 10;
        
        // optimize gene rate parameters
	gsl_root_fsolver_set(sol_gene_rate, &opt_gene_rate, low, high);

	status = GSL_CONTINUE;
	for (iter=0; iter<maxiter && status==GSL_CONTINUE; iter++) {
            // do one iteration
	    status = gsl_root_fsolver_iterate(sol_gene_rate);
            if (status)
                break;

            // check convergence
	    low = gsl_root_fsolver_x_lower(sol_gene_rate);
	    high = gsl_root_fsolver_x_upper(sol_gene_rate);
	    status = gsl_root_test_interval(low, high, 0, 0.01);
	}
	gene_nu = gsl_root_fsolver_root(sol_gene_rate);
	gene_alpha = gene_nu + 1;
        gene_beta = gene_nu;
        fprintf(stderr, "nu = %f (iter=%d)\n", gene_nu, iter);
        
        
        // optimize each species rate parmater set
        for (int i=0; i<nspecies; i++) {
            cur_species = i;
            gsl_vector_set(init_x, 0, sp_alpha[i]);
            gsl_vector_set(init_x, 1, sp_beta[i]);
            gsl_multimin_fdfminimizer_set(sol_sp_rate, &opt_sp_rate, init_x, 
                                          step_size, tol);            
            status = GSL_CONTINUE;

            
	    for (iter=0; iter<maxiter && status==GSL_CONTINUE; iter++) {
                // do one iteration
                status = gsl_multimin_fdfminimizer_iterate(sol_sp_rate);
                if (status)
                    break;        
                // get gradient
                status = gsl_multimin_test_gradient(sol_sp_rate->gradient, epsabs);
            }

            double lk = likelihood();
            fprintf(stderr, "species %d %d %f\n", i, iter, lk);

            sp_alpha[i] = gsl_vector_get(sol_sp_rate->x, 0);
            sp_beta[i] = gsl_vector_get(sol_sp_rate->x, 1);

            //printf("sp[%d] = (%f, %f)\n", i, sp_alpha[i], sp_beta[i]);
        }

        gsl_vector_free(init_x);
    }

    // data
    int ntrees;
    int nspecies;
    int *gene_sizes;
    float **lengths;

    // given fixed parameters
    float *times;

    // model parameters
    float *sp_alpha;
    float *sp_beta;
    float gene_alpha;
    float gene_beta;
    float gene_nu;
    
    //protected:

    // hidden data
    int nrates;
    Matrix<float> gtab;
    Matrix<double> pgtab;

    // optimizers
    gsl_multimin_fdfminimizer *sol_sp_rate;
    gsl_root_fsolver *sol_gene_rate;
    gsl_function opt_gene_rate;
    gsl_multimin_function_fdf opt_sp_rate;
    int cur_species;
};


extern "C" {

void train(int ntrees, int nspecies, int *gene_sizes,
           float **lengths, float *times,
           float *sp_alpha, float *sp_beta, 
           float *gene_alpha, float *gene_beta,
           int nrates, int max_iter)
{
    RatesEM em(ntrees, nspecies, nrates, gene_sizes, lengths, times,
               sp_alpha, sp_beta, *gene_alpha, *gene_beta);
    
    
    // make initial guess for model parameters
    em.init_params();
    
    // iterate until convergence
    for (int iter=0; iter<max_iter; iter++) {
        em.EStep();
        em.MStep();
    }
      
}


RatesEM *allocRatesEM(int ntrees, int nspecies, int nrates,
                      int *gene_sizes,
                      float **lengths, float *times,
                      float *sp_alpha, float *sp_beta, 
                      float gene_alpha, float gene_beta)
{
    int *gene_sizes2 = new int [ntrees];
    for (int j=0; j<ntrees; j++)
        gene_sizes2[j] = gene_sizes[j];

    // copy lengths
    float **lengths2 = new float* [ntrees];
    for (int j=0; j<ntrees; j++) {
        lengths2[j] = new float [nspecies];
        for (int i=0; i<nspecies; i++)
            lengths2[j][i] = lengths[j][i];
    }

    float *times2 = new float [nspecies];
    float *sp_alpha2 = new float [nspecies];
    float *sp_beta2 = new float [nspecies];

    for (int i=0; i<nspecies; i++) {
        times2[i] = times[i];
        sp_alpha2[i] = sp_alpha[i];
        sp_beta2[i] = sp_beta[i];
    }

    return new RatesEM(ntrees, nspecies, nrates, gene_sizes2,
                       lengths2, times2,
		       sp_alpha2, sp_beta2, gene_alpha, gene_beta);
}


void freeRatesEM(RatesEM *em)
{
    delete [] em->gene_sizes;

    for (int j=0; j<em->ntrees; j++)
        delete [] em->lengths[j];
    delete [] em->lengths;

    delete [] em->times;
    delete [] em->sp_alpha;
    delete [] em->sp_beta;
    
    delete em;
}


void RatesEM_Init(RatesEM *em)
{
    em->init_params();
}


void RatesEM_EStep(RatesEM* em)
{
    em->EStep();
}

void RatesEM_MStep(RatesEM* em)
{
    em->MStep();
}

float RatesEM_likelihood(RatesEM *em)
{
    return em->likelihood();
}


void RatesEM_getParams(RatesEM *em, float *params)
{
    params[0] = em->gene_nu + 1; //em->gene_alpha;
    params[1] = em->gene_nu; //em->gene_beta;

    for (int i=0; i<em->nspecies; i++) {
        params[2+2*i] = em->sp_alpha[i];
        params[2+2*i+1] = em->sp_beta[i];
    }
}



} // extern "C"



} // spidir
