/*=============================================================================

  Matt Rasmussen
  Copyright 2010-2011

  Estimating birth-death rates from gene counts (implementing Hahn et al 2005)

=============================================================================*/

// c/c++ includes
#include <math.h>

// 3rd party
#include <gsl/gsl_roots.h>
#include <gsl/gsl_multimin.h>

// use old simplex method as backup if it new one does not exist
#ifndef gsl_multimin_fminimizer_nmsimplex2
#define gsl_multimin_fminimizer_nmsimplex2 gsl_multimin_fminimizer_nmsimplex
#endif


// spidir includes
#include "common.h"
#include "Matrix.h"
#include "Tree.h"
#include "phylogeny.h"
#include "birthdeath.h"


namespace spidir
{

extern "C" {


//=============================================================================
// gene counts on species tree


double birthDeathTreeCounts(Tree *tree, int nspecies, int *counts, 
                            float birth, float death, int maxgene,
                            int rootgene, double **tab)
{
    // set up dynamic table
    bool cleanup = false;
    if (!tab) {
        // allocate dynamic table
        tab = allocMatrix<double>(tree->nnodes, maxgene);
        cleanup = true;
    }

    // initialize leaves
    for (int i=0; i<nspecies; i++) {
        for (int j=0; j<maxgene; j++)
            tab[i][j] = -INFINITY;
        tab[i][counts[i]] = 0.0;
    }

    // perform post order traversal of tree
    ExtendArray<Node*> postnodes(0, tree->nnodes);
    getTreePostOrder(tree, &postnodes);
    for (int a=0; a<tree->nnodes; a++) {
        const Node *node = postnodes[a];
        //const int i = node.name;

        // skip leaves
        if (node->isLeaf())
            continue;

        for (int j=0; j<maxgene; j++) {

            // compute product over children
            double prod = 0.0;        
            for (int ci=0; ci<node->nchildren; ci++) {
                const int c = node->children[ci]->name;
                const double t = node->children[ci]->dist;
                double sum = -INFINITY;
                
                for (int j2=0; j2<maxgene; j2++) {
                    if (j < 20 && j2 < 20)
                        sum = logadd(sum, log(birthDeathCounts(j, j2, t, 
                                                               birth, death)) +
                                     tab[c][j2]);
                    else
                        sum = logadd(sum, birthDeathCountsLog(j, j2, t, 
                                                              birth, death) +
                                     tab[c][j2]);
                }
                prod += sum;
            }

            tab[node->name][j] = prod;
        }
    }

    double prob = tab[tree->root->name][rootgene];

    // cleanup
    if (cleanup)
        freeMatrix(tab, tree->nnodes);

    return prob;
}



double birthDeathForestCounts(Tree *tree, int nspecies, int nfams,
                              int **counts, int *mult,
                              float birth, float death, int maxgene,
                              int rootgene, double **tab)
{
    // set up dynamic table
    bool cleanup = false;
    if (!tab) {
        tab = allocMatrix<double>(tree->nnodes, maxgene);
        cleanup = true;
    }

    double logl = 0.0;

    // loop through families
    for (int i=0; i<nfams; i++) {
        int top = 0;
        for (int j=0; j<nspecies; j++) {
            if (counts[i][j] > top) top = counts[i][j];
        }

        int maxgene2 = top * 2;
        if (maxgene2 < 10) maxgene2 = 10;
        if (maxgene2 > maxgene) maxgene2 = maxgene;
        
        logl += mult[i] * birthDeathTreeCounts(tree, nspecies, counts[i], 
                                               birth, death, maxgene2,
                                               rootgene, tab);
    }

    // cleanup
    if (cleanup)
        freeMatrix(tab, tree->nnodes);

    return logl;
}


class BirthDeathCountsML
{
public:
    BirthDeathCountsML(Tree *tree, int nspecies, int nfams,
                       int **counts, int *mult, 
                       double birth, double death, 
                       double step,
                       int maxgene,
                       int rootgene=1) :
        iter(0),
        tree(tree),
        nspecies(nspecies),
        nfams(nfams),
        counts(counts),
        mult(mult),
        birth(birth),
        death(death),
        maxgene(maxgene),
        rootgene(rootgene)
    {

        // allocate dynamic table
        tab = allocMatrix<double>(tree->nnodes, maxgene);

        // allocate optimizer
        opt = gsl_multimin_fminimizer_alloc(
             gsl_multimin_fminimizer_nmsimplex2, 2);
        func.f = &function;
        func.n = 2;
        func.params = this;
                
        // init optimizer
        gsl_vector *init_x = gsl_vector_alloc(2);
        gsl_vector *step_size = gsl_vector_alloc(2);
        
        gsl_vector_set(init_x, 0, birth);
        gsl_vector_set(init_x, 1, death);
        gsl_vector_set(step_size, 0, step);
        gsl_vector_set(step_size, 1, step);
        gsl_multimin_fminimizer_set(opt, &func, init_x, step_size);
        gsl_vector_free(init_x);

        /*
        gsl_root_fsolver *sol_gene_rate = 
            gsl_root_fsolver_alloc(gsl_root_fsolver_brent);

        // setup optimizer for gene rates
        gsl_function opt_gene_rate;
        opt_gene_rate.function = &gene_rate_f;
        opt_gene_rate.params = NULL;
        */
    }

    ~BirthDeathCountsML()
    {
        gsl_multimin_fminimizer_free(opt);
        freeMatrix(tab, tree->nnodes);

        /*
        gsl_root_fsolver_free(sol_gene_rate);
        gsl_multimin_fdfminimizer_free(sol_sp_rate);
        */
    }

    static double function(const gsl_vector *x, void *params)
    {
        double birth = gsl_vector_get(x, 0);
        double death = gsl_vector_get(x, 1);

        double bpenalty = 0.0;
        double dpenalty = 0.0;

        if (birth < 0) {
            bpenalty = -birth;
            birth = 0.0000001;
        }

        if (death < 0) {
            dpenalty = -death;
            death = 0.0000002;
        }
        
        BirthDeathCountsML *p = (BirthDeathCountsML*) params;

        double prob = -birthDeathForestCounts(p->tree, p->nspecies, p->nfams,
                                              p->counts, p->mult,
                                              birth, death, p->maxgene,
                                              p->rootgene, p->tab);

        // constraint penalties
        prob += exp(bpenalty) - 1.0 + exp(dpenalty) - 1.0;

        return prob;
    }


    int iterate()
    {
        int status;
        
        // do one iteration
        status = gsl_multimin_fminimizer_iterate(opt);
        birth = gsl_vector_get(opt->x, 0);
        death = gsl_vector_get(opt->x, 1);
        if (status)
            return status;

        double epsabs = min(fabs(birth * .001), fabs(death * .001));
        //double epsabs = .01;
        double size = gsl_multimin_fminimizer_size(opt);
        
        // get gradient
        status = gsl_multimin_test_size(size, epsabs);
        
        birth = gsl_vector_get(opt->x, 0);
        death = gsl_vector_get(opt->x, 1);
        
        return status;
    }    

    void getBirthDeath(float *_birth, float *_death)
    {
        *_birth = birth;
        *_death = death;
    }

    double getSize()
    {
        return gsl_multimin_fminimizer_size(opt);
    }


    int iter;
    Tree *tree;
    int nspecies;
    int nfams;
    int **counts;
    int *mult;
    double birth;
    double death;
    int maxgene;
    int rootgene;
    double **tab;

    gsl_multimin_fminimizer *opt;
    gsl_multimin_function func;
};

void *birthDeathCountsML_alloc(Tree *tree, int nspecies, int nfams,
                               int **counts, int *mult,
                               float birth, float death, float step,
                               int maxgene,
                               int rootgene)
{
    // handle errors myself
    gsl_set_error_handler_off();

    // copy arrays to separate memory
    int **counts2 = allocMatrix<int>(nfams, nspecies);
    int *mult2 = new int [nfams];

    for (int i=0; i<nfams; i++) {
        for (int j=0; j<nspecies; j++) {
            counts2[i][j] = counts[i][j];
        }
        mult2[i] = mult[i];
    }

    return new BirthDeathCountsML(tree, nspecies, nfams,
                                  counts2, mult2,
                                  birth, death, step,
                                  maxgene, rootgene);
}

void birthDeathCountsML_free(void *opt)
{
    BirthDeathCountsML* opt2 = (BirthDeathCountsML*) opt;
    freeMatrix(opt2->counts, opt2->nfams);
    delete [] opt2->mult;
    delete opt2;
}


int birthDeathCountsML_iter(void *opt, float *birth, float *death, float *size)
{
    BirthDeathCountsML* opt2 = (BirthDeathCountsML*) opt;
    int status = opt2->iterate();
    opt2->getBirthDeath(birth, death);
    *size = opt2->getSize();
    return status;
}

/*
double birthDeathCountsML(Tree *tree, int nspecies, int nfams,
                          int **counts, int *mult,
                          float *birth, float *death, int maxgene,
                          int rootgene)
{ 
}
*/



}

} // namespace spidir
