/*=============================================================================

  Matt Rasmussen
  Copyright 2007-2011

  Computing sequence likelihood
  Maximum Likelihood Branch Length Estimation

=============================================================================*/

// c++ headers
#include <math.h>
#include <time.h>

// 3rd party
#include <gsl/gsl_roots.h>
#include <gsl/gsl_multimin.h>

// spidir headers
#include "common.h"
#include "hky.h"
#include "logging.h"
#include "Matrix.h"
#include "parsimony.h"
#include "roots.h"
#include "seq.h"
#include "seq_likelihood.h"
#include "Tree.h"




namespace spidir {



// prototype
template <class Model>
void calcLkTableRow(int seqlen, Model &model,
		    floatlk *lktablea, floatlk *lktableb, floatlk *lktablec, 
		    float adist, float bdist);

template <class Model>
void calcLogLkTableRow(int seqlen, Model &model,
                       floatlk *lktablea, floatlk *lktableb, floatlk *lktablec,
                       float adist, float bdist);


//=============================================================================


template <class Model, class DModel>
class DistLikelihoodDeriv
{
public:
    DistLikelihoodDeriv(int seqlen, 
			Model *model, DModel *dmodel) :
        probs1(NULL),
        probs2(NULL),
        seqlen(seqlen),
        bgfreq(NULL),
        model(model),
	dmodel(dmodel)
    {
	probs3 = new floatlk [seqlen*4];
	probs4 = new floatlk [seqlen*4];
    }
    
    ~DistLikelihoodDeriv()
    {
	delete [] probs3;
	delete [] probs4;
    }

    void set_params(floatlk *_probs1, floatlk *_probs2, const float *_bgfreq)
    {
        probs1 = _probs1;
        probs2 = _probs2;
        bgfreq = _bgfreq;
    }

    double operator()(float t)
    {
        // trivial case
        if (t < 0)
            return INFINITY;
	

	// g(t, j)
	calcLkTableRow(seqlen, *model, probs1, probs2, probs3, t, 0);

	// g'(t, j)
	calcDerivLkTableRow(seqlen, *dmodel, probs1, probs2, probs4, t);


	// interate over sequence
	double dlogl = 0.0;
	for (int j=0; j<seqlen; j++) {
	    double sum1 = 0.0, sum2 = 0.0;
	    for (int k=0; k<4; k++) {
		sum1 += bgfreq[k] * probs3[matind(4,j,k)];
		sum2 += bgfreq[k] * probs4[matind(4,j,k)];
	    }
	    dlogl += sum2 / sum1;
	}
        
	return dlogl;
    }
    
    floatlk *probs1;
    floatlk *probs2;
    floatlk *probs3;
    floatlk *probs4;
    int seqlen;
    const float *bgfreq;
    Model *model;
    DModel *dmodel;
};



template <class Model, class DModel, class D2Model>
class DistLikelihoodDeriv2
{
public:
    DistLikelihoodDeriv2(int seqlen, 			
                         Model *model, DModel *dmodel, D2Model *d2model) :
        seqlen(seqlen),
        model(model),
	dmodel(dmodel),
        d2model(d2model)
    {
	probs3 = new floatlk [seqlen*4];
	probs4 = new floatlk [seqlen*4];
        probs5 = new floatlk [seqlen*4];
    }
    
    ~DistLikelihoodDeriv2()
    {
	delete [] probs3;
	delete [] probs4;
        delete [] probs5;
    }

    void set_params(floatlk *_probs1, floatlk *_probs2, const float *_bgfreq)
    {
        probs1 = _probs1;
        probs2 = _probs2;
        bgfreq = _bgfreq;
    }


    double operator()(float t)
    {
        // trivial case
        if (t < 0)
            return -INFINITY;
	

	// g(t, j)
	calcLkTableRow(seqlen, *model, probs1, probs2, probs3, t, 0);

	// g'(t, j)
	calcDerivLkTableRow(seqlen, *dmodel, probs1, probs2, probs4, t);

	// g''(t, j)
	calcDerivLkTableRow(seqlen, *d2model, probs1, probs2, probs5, t);


	// interate over sequence
	double d2logl = 0.0;
	for (int j=0; j<seqlen; j++) {
	    double g = 0.0, dg = 0.0, d2g = 0.0;
	    for (int k=0; k<4; k++) {
		g += bgfreq[k] * probs3[matind(4,j,k)];
		dg += bgfreq[k] * probs4[matind(4,j,k)];
                d2g += bgfreq[k] * probs5[matind(4,j,k)];
	    }
	    d2logl += - dg*dg/(g*g) + d2g/g;
	}
        
	return d2logl;
    }
    
    floatlk *probs1;
    floatlk *probs2;
    floatlk *probs3;
    floatlk *probs4;
    floatlk *probs5;
    int seqlen;
    const float *bgfreq;
    Model *model;
    DModel *dmodel;
    D2Model *d2model;
};


// Find the maximum likelihood estimate (MLE) of the distance (time) between 
// two sequences (represented probabilistically as probs1 and probs2)
//
//  bgfreq = background frequency of bases
//  model = sequence evolution model
//
template <class Model>
float mleDistance(floatlk *probs1, floatlk *probs2, int seqlen, 
                  const float *bgfreq, Model &model, 
                  float t0=.001, float t1=1, float step=.0001)
{
    typename Model::Deriv *dmodel = model.deriv();

    DistLikelihoodDeriv<Model, typename Model::Deriv> df(
        probs1, probs2, seqlen, bgfreq, &model, dmodel);
    float mle = bisectRoot(df, t0, t1, .0001);

    delete dmodel;
    return mle;
}

extern "C" {

/*

  lk = prod_j sum_k bgfreq[k] * lktable[root][j,k] 
     = prod_j sum_k bgfreq[k] * (sum_x P(x|k, t_a) lktable[a][j,x]) *
                                (sum_y P(y|k, t_b) lktable[b][j,y])

 */
floatlk branchLikelihoodHky(floatlk *probs1, floatlk *probs2, int seqlen, 
                            const float *bgfreq, float kappa, float t)
{

    HkyModel hky(bgfreq, kappa);

    // allocate precomputed probability terms
    floatlk *probs3 = new floatlk [seqlen*4];
    
    calcLkTableRow(seqlen, hky, probs1, probs2, probs3, 0, t);

    floatlk logl = 0.0;

    // interate over sequence
    for (int j=0; j<seqlen; j++) {
	floatlk sum = 0.0;
	for (int k=0; k<4; k++)
	    sum += bgfreq[k] * probs3[matind(4,j,k)];
	logl += log(sum);
    }

    // free probability table
    delete [] probs3;

    return logl;
}


/*

  lk = prod_j sum_k bgfreq[k] * lktable[root][j,k] 
     = prod_j sum_k bgfreq[k] * (sum_x P(x|k, t_a) lktable[a][j,x]) *
                                (sum_y P(y|k, t_b) lktable[b][j,y])

 */
floatlk branchLogLikelihoodHky(floatlk *probs1, floatlk *probs2, int seqlen, 
                               const float *bgfreq, float kappa, float t)
{

    HkyModel hky(bgfreq, kappa);

    // allocate precomputed probability terms
    floatlk *probs3 = new floatlk [seqlen*4];
    
    calcLogLkTableRow(seqlen, hky, probs1, probs2, probs3, 0, t);

    floatlk logl = 0.0;

    // interate over sequence
    for (int j=0; j<seqlen; j++) {
	floatlk sum = -INFINITY;
	for (int k=0; k<4; k++)
	    sum = logadd(sum, log(bgfreq[k]) + probs3[matind(4,j,k)]);
	logl += sum;
    }

    // free probability table
    delete [] probs3;
    
    return logl;
}



floatlk branchLikelihoodHkyDeriv(floatlk *probs1, floatlk *probs2, int seqlen, 
                                 const float *bgfreq, float kappa, float t)
{
    HkyModel hky(bgfreq, kappa);
    HkyModelDeriv dhky(bgfreq, kappa);
    DistLikelihoodDeriv<HkyModel, HkyModelDeriv> df(seqlen, &hky, &dhky);
    df.set_params(probs1, probs2, bgfreq);
    return df(t);
}

floatlk branchLikelihoodHkyDeriv2(floatlk *probs1, floatlk *probs2, 
                                  int seqlen, 
                                  const float *bgfreq, float kappa, float t)
{
    HkyModel hky(bgfreq, kappa);
    HkyModelDeriv dhky(bgfreq, kappa);
    HkyModelDeriv2 d2hky(bgfreq, kappa);
    DistLikelihoodDeriv2<HkyModel, HkyModelDeriv, HkyModelDeriv2> d2f
	(seqlen, &hky, &dhky, &d2hky);
    d2f.set_params(probs1, probs2, bgfreq);
    return d2f(t);
}

floatlk mleDistanceHky(floatlk *probs1, floatlk *probs2, int seqlen, 
                       const float *bgfreq, float kappa,
                       float t0, float t1)
{
    HkyModel hky(bgfreq, kappa);
    HkyModelDeriv dhky(bgfreq, kappa);
    DistLikelihoodDeriv<HkyModel, HkyModelDeriv> df(seqlen, &hky, &dhky);
    df.set_params(probs1, probs2, bgfreq);
    
    return bisectRoot(df, t0, t1, .0001);
}

} // extern "C"


//=============================================================================
// Dynamic programming of conditional likelihood


class LikelihoodTable 
{
public:

    LikelihoodTable(int nnodes, int seqlen) :
        nnodes(nnodes),
        seqlen(seqlen)
    {
        // allocate conditional likelihood dynamic programming table
        lktable = new floatlk* [nnodes];
        for (int i=0; i<nnodes; i++)
            lktable[i] = new floatlk [4 * seqlen];
    }

    ~LikelihoodTable()
    {
        // cleanup
        for (int i=0; i<nnodes; i++)
            delete [] lktable[i];
        delete [] lktable;
    }

    floatlk **lktable;

    int nnodes;
    int seqlen;
    
};



/*

    From: Inferring Phylogenies. Felsenstein. p 254.
    
    Baisc recursion of conditional likelihoods
        lktable[c][j,k] = (sum_x P(x|k, t_a) lktable[a][j,x]) *
                          (sum_y P(y|k, t_b) lktable[b][j,y])
    
    where c is a node, with child nodes a and b
          t_a, t_b are the branch lengths of branches a and b
          j indexes sites
          k,x,y are DNA bases (e.g. A=1, C=2, G=3, T=4)

*/


// conditional likelihood recurrence
template <class Model>
void calcLkTableRow(int seqlen, Model &model,
		    floatlk *lktablea, floatlk *lktableb, floatlk *lktablec, 
		    float adist, float bdist)
{
    float atransmat[16];
    float btransmat[16];
    
    // build transition matrices
    model.getMatrix(adist, atransmat);
    model.getMatrix(bdist, btransmat);    
    
    // iterate over sites
    for (int j=0; j<seqlen; j++) {
        const floatlk *terma = &lktablea[matind(4, j, 0)];
        const floatlk *termb = &lktableb[matind(4, j, 0)];
        
        for (int k=0; k<4; k++) {
            const float *aptr = &atransmat[4*k];
            const float *bptr = &btransmat[4*k];
            
            // sum_x P(x|k, t_a) lktable[a][j,x]
            const floatlk prob1 = aptr[0] * terma[0] +
                                  aptr[1] * terma[1] +
                                  aptr[2] * terma[2] +
                                  aptr[3] * terma[3];

            // sum_y P(y|k, t_b) lktable[b][j,y]
            const floatlk prob2 = bptr[0] * termb[0] +
                                  bptr[1] * termb[1] +
                                  bptr[2] * termb[2] +
                                  bptr[3] * termb[3];
            

            lktablec[matind(4, j, k)] = prob1 * prob2;
        }
    }
}


// conditional likelihood recurrence
template <class Model>
void calcLogLkTableRow(int seqlen, Model &model,
                       floatlk *lktablea, floatlk *lktableb, floatlk *lktablec,
                       float adist, float bdist)
{
    float atransmat[16];
    float btransmat[16];
    
    // build transition matrices
    model.getMatrix(adist, atransmat);
    model.getMatrix(bdist, btransmat);

    for (int i=0; i<16; i++) {
        atransmat[i] = log(atransmat[i]);
        btransmat[i] = log(btransmat[i]);
    }

    // iterate over sites
    for (int j=0; j<seqlen; j++) {
        const floatlk *terma = &lktablea[matind(4, j, 0)];
        const floatlk *termb = &lktableb[matind(4, j, 0)];
        
        for (int k=0; k<4; k++) {
            const float *aptr = &atransmat[k*4];
            const float *bptr = &btransmat[k*4];
            
            // sum_x P(x|k, t_a) lktable[a][j,x]
            floatlk prob1 = aptr[0] + terma[0];
            prob1 = logadd(prob1, aptr[1] + terma[1]);
            prob1 = logadd(prob1, aptr[2] + terma[2]);
            prob1 = logadd(prob1, aptr[3] + terma[3]);

            // sum_y P(y|k, t_b) lktable[b][j,y]
            floatlk prob2 = bptr[0] + termb[0];
            prob2 = logadd(prob2, bptr[1] + termb[1]);
            prob2 = logadd(prob2, bptr[2] + termb[2]);
            prob2 = logadd(prob2, bptr[3] + termb[3]);

            lktablec[matind(4, j, k)] = prob1 + prob2;
        }
    }
}



// conditional likelihood recurrence
template <class DModel>
void calcDerivLkTableRow(int seqlen, DModel &dmodel,
			 floatlk *lktablea, floatlk *lktableb, 
                         floatlk *lktablec, 
			 float adist)
{
    float btransmat[16];
    
    // build transition matrix
    dmodel.getMatrix(adist, btransmat);
    
    // iterate over sites
    for (int j=0; j<seqlen; j++) {
        const floatlk *terma = &lktablea[matind(4, j, 0)];
        const floatlk *termb = &lktableb[matind(4, j, 0)];
        
        for (int k=0; k<4; k++) {
            //float *aptr = &atransmat[4*k];
            const float *bptr = &btransmat[4*k];
            
            // sum_y P(y|k, t_b) lktable[b][j,y]
            const floatlk prob2 = bptr[0] * termb[0] +
                                  bptr[1] * termb[1] +
                                  bptr[2] * termb[2] +
                                  bptr[3] * termb[3];
            
            lktablec[matind(4, j, k)] = terma[k] * prob2;
        }
    }
}


// initialize the condition likelihood table
template <class Model>
void calcLkTable(floatlk** lktable, Tree *tree, 
                 int nseqs, int seqlen, char **seqs, Model &model)
{
    // recursively calculate cond. lk. of internal nodes
    ExtendArray<Node*> nodes(0, tree->nnodes);
    getTreePostOrder(tree, &nodes);
    
    for (int l=0; l<nodes.size(); l++) {
        Node *node = nodes[l];
        int i = node->name;
        
        if (node->isLeaf()) {
            // initialize leaves from sequence
        
            // iterate over sites
            for (int j=0; j<seqlen; j++) {
                int base = dna2int[int(seqs[i][j])];

                if (base == -1) {
                    // handle gaps
                    lktable[i][matind(4, j, 0)] = 1.0;
                    lktable[i][matind(4, j, 1)] = 1.0;
                    lktable[i][matind(4, j, 2)] = 1.0;
                    lktable[i][matind(4, j, 3)] = 1.0;
                } else {
                    // initialize base
                    lktable[i][matind(4, j, 0)] = 0.0;
                    lktable[i][matind(4, j, 1)] = 0.0;
                    lktable[i][matind(4, j, 2)] = 0.0;
                    lktable[i][matind(4, j, 3)] = 0.0;

                    lktable[i][matind(4, j, base)] = 1.0;
                }
            }
        } else {
            // compute internal nodes from children
            Node *node1 = node->children[0];
            Node *node2 = node->children[1];
            
            calcLkTableRow(seqlen, model, 
                           lktable[node1->name], 
			   lktable[node2->name], 
			   lktable[node->name],
                           node1->dist, node2->dist);
        }
    }
}


// initialize the condition likelihood table
template <class Model>
void calcLogLkTable(floatlk** lktable, Tree *tree, 
                    int nseqs, int seqlen, char **seqs, Model &model)
{
    // recursively calculate cond. lk. of internal nodes
    ExtendArray<Node*> nodes(0, tree->nnodes);
    getTreePostOrder(tree, &nodes);
    
    for (int l=0; l<nodes.size(); l++) {
        Node *node = nodes[l];
        int i = node->name;
        
        if (node->isLeaf()) {
            // initialize leaves from sequence
        
            // iterate over sites
            for (int j=0; j<seqlen; j++) {
                int base = dna2int[int(seqs[i][j])];

                if (base == -1) {
                    // handle gaps
                    lktable[i][matind(4, j, 0)] = 0.0;
                    lktable[i][matind(4, j, 1)] = 0.0;
                    lktable[i][matind(4, j, 2)] = 0.0;
                    lktable[i][matind(4, j, 3)] = 0.0;
                } else {
                    // initialize base
                    lktable[i][matind(4, j, 0)] = -INFINITY;
                    lktable[i][matind(4, j, 1)] = -INFINITY;
                    lktable[i][matind(4, j, 2)] = -INFINITY;
                    lktable[i][matind(4, j, 3)] = -INFINITY;

                    lktable[i][matind(4, j, base)] = 0.0;
                }
            }
        } else {
            // compute internal nodes from children
            Node *node1 = node->children[0];
            Node *node2 = node->children[1];
            
            calcLogLkTableRow(seqlen, model, 
                              lktable[node1->name], 
                              lktable[node2->name], 
                              lktable[node->name],
                              node1->dist, node2->dist);
        }
    }
}



// calculate log(P(D | T, B))
template <class Model>
floatlk getTotalLikelihood(floatlk** lktable, Tree *tree, 
                           int seqlen, Model &model, const float *bgfreq)
{
    // integrate over the background base frequency
    const floatlk *rootseq = lktable[tree->root->name];
    floatlk lk = 0.0;
    for (int k=0; k<seqlen; k++) {
        floatlk prob = 0.0;
        for (int x=0; x<4; x++)
            prob += bgfreq[x] * rootseq[matind(4, k, x)];
        lk += log(prob);
    }

    // return log likelihood
    return lk;
}


// calculate log(P(D | T, B))
template <class Model>
floatlk getTotalLogLikelihood(floatlk** lktable, Tree *tree, 
                              int seqlen, Model &model, const float *bgfreq)
{
    // integrate over the background base frequency
    const floatlk *rootseq = lktable[tree->root->name];
    floatlk lk = 0.0;
    for (int k=0; k<seqlen; k++) {
        floatlk prob = -INFINITY;
        for (int x=0; x<4; x++)
            prob = logadd(prob, log(bgfreq[x]) + rootseq[matind(4, k, x)]);
        lk += prob;
    }

    // return log likelihood
    return lk;
}


template <class Model>
floatlk calcSeqProb(Tree *tree, int nseqs, char **seqs, 
                    const float *bgfreq, Model &model)
{
    int seqlen = strlen(seqs[0]);
    
    LikelihoodTable table(tree->nnodes, seqlen);
    calcLogLkTable(table.lktable, tree, nseqs, seqlen, seqs, model);
    floatlk logl = getTotalLogLikelihood(table.lktable, tree, seqlen, 
                                         model, bgfreq);
    
    return logl;
}

extern "C" {

floatlk calcSeqProbHky(Tree *tree, int nseqs, char **seqs, 
                       const float *bgfreq, float ratio)
{
    HkyModel hky(bgfreq, ratio);
    return calcSeqProb(tree, nseqs, seqs, bgfreq, hky);
}

} // extern "C"


//=============================================================================
// find MLE branch lengths


// rootings are stored as edges (node1, node2).  A root order is 
// a concatenation of rooting pairs.
void getRootOrder(Tree *tree, ExtendArray<Node*> *nodes, Node *node=NULL)
{
    if (!node) {
        // start at tree root (but don't include root)
        getRootOrder(tree, nodes, tree->root->children[0]);
        getRootOrder(tree, nodes, tree->root->children[1]);
    } else {
        // record pre-process
        if (node->parent == tree->root) {
            nodes->append(tree->root->children[0]);
            nodes->append(tree->root->children[1]);
        } else {
            nodes->append(node);
            nodes->append(node->parent);
        }

        // recurse
        for (int i=0; i<node->nchildren; i++)
            getRootOrder(tree, nodes, node->children[i]);
    }
}

template <class Model>
class MLBranchAlgorithm
{
public:

    MLBranchAlgorithm(Tree *tree, int seqlen, Model *model) :
	table(tree->nnodes, seqlen),
        model(model),
        dmodel(model->deriv()),
        d2model(dmodel->deriv()),
        lk_deriv(seqlen, model, dmodel),
	lk_deriv2(seqlen, model, dmodel, d2model)
    {

        minx = 0.00001;
        maxx = 10;

        // setup branch length optimizer
        opt = gsl_root_fdfsolver_alloc(gsl_root_fdfsolver_newton);

        opt_func.f = &branch_f;
        opt_func.df = &branch_df;
        opt_func.fdf = &branch_fdf;
        opt_func.params = this;
    }

    ~MLBranchAlgorithm()
    {
        gsl_root_fdfsolver_free(opt);
        delete dmodel;
        delete d2model;
    }


    static double branch_f(double x, void *params)
    {       
        double y;
        MLBranchAlgorithm *p = (MLBranchAlgorithm*) params;
        if (x < p->minx || x > p->maxx)
            y = p->lk_deriv(p->minx);
        else
            y = p->lk_deriv(x);
        //printf("x %f y %f\n", x, y);
        return y;
    }

    static double branch_df(double x, void *params)
    {
        double dy;
        MLBranchAlgorithm *p = (MLBranchAlgorithm*) params;
        if (x < p->minx || x > p->maxx)
            dy = p->lk_deriv(p->minx) / (x - p->minx);
        else
            dy = p->lk_deriv2(x);
        
        //printf("x %f dy %f\n", x, dy);
        return dy;
    }

    static void branch_fdf(double x, void *params, 
                           double *f, double *df)
    {
        *f = branch_f(x, params);
        *df = branch_df(x, params);
    }


    // slower more stable branch fitting
    float fitBranch2(Tree *tree, const float *bgfreq, float initdist)
    {
        Node *node1 = tree->root->children[0];
        Node *node2 = tree->root->children[1];

        lk_deriv.set_params(table.lktable[node1->name], 
                            table.lktable[node2->name], 
                            bgfreq);

        return bisectRoot(lk_deriv, 0.0,
                          max(initdist*10.0, 0.01), .0001);

    }

    float fitBranch(Tree *tree, const float *bgfreq, float initdist)
    {
        double r = initdist, r0 = initdist;
        const double esp = 1e-3;

        if (initdist < minx)
            initdist = minx;

        Node *node1 = tree->root->children[0];
        Node *node2 = tree->root->children[1];
        
        lk_deriv.set_params(table.lktable[node1->name], 
                            table.lktable[node2->name], 
                            bgfreq);
        
        lk_deriv2.set_params(table.lktable[node1->name], 
                             table.lktable[node2->name], 
                             bgfreq);

        gsl_root_fdfsolver_set(opt, &opt_func, initdist);

        int status = GSL_CONTINUE;
        const int maxiter = 30;
        int iter;
        for (iter=0; iter<maxiter && status==GSL_CONTINUE; iter++) {
            //printf("root %f\n", r);

            // do one iteration
            status = gsl_root_fdfsolver_iterate(opt);
            r0 = r;
            r = gsl_root_fdfsolver_root(opt);
            if (status)
                break;
            status = gsl_root_test_delta(r, r0, 0, esp);
        }
        //printf("root done iters=%d r=%f\n", iter, r);

        if (iter == maxiter || status != GSL_SUCCESS) {
            r = fitBranch2(tree, bgfreq, initdist);
        }

        //printf("bisect: %f\n\n", r);

        // return final branch length
        return r;
    }
        



    floatlk fitBranches(Tree *tree, int nseqs, int seqlen, char **seqs, 
		      const float *bgfreq, 
		      ExtendArray<Node*> &rootingOrder)
    {
	float logl = -INFINITY;
	floatlk **lktable = table.lktable;

	for (int i=0; i<rootingOrder.size(); i+=2) {
	    // remembering old children of root
	    Node *oldnode1 = tree->root->children[0];
	    Node *oldnode2 = tree->root->children[1];

	    // choose new root
	    tree->reroot(rootingOrder[i], rootingOrder[i+1]);
	    if (tree->root->children[0] == oldnode1 &&
		tree->root->children[1] == oldnode2)
		continue;	    

	    // rebuild invaild likelihood values
	    // determine starting node to rebuild
	    Node *ptr;
	    if (oldnode1->parent == oldnode2)
		ptr = oldnode1;
	    else if (oldnode2->parent == oldnode1)
		ptr = oldnode2;
	    else
		assert(0);

	    // walk up to root of tree, rebuilding conditional likelihoods
	    for (; ptr; ptr = ptr->parent) {
		if (!ptr->isLeaf())
		    calcLkTableRow(seqlen, *model, 
				   lktable[ptr->children[0]->name], 
				   lktable[ptr->children[1]->name], 
				   lktable[ptr->name],
				   ptr->children[0]->dist, 
				   ptr->children[1]->dist);
	    }

	    // get total probability before branch length change
	    floatlk loglBefore = getTotalLikelihood(lktable, tree, 
                                                    seqlen, *model, bgfreq);

            Node *node1 = tree->root->children[0];
            Node *node2 = tree->root->children[1];
            float initdist = node1->dist + node2->dist;

            // find new MLE branch length for root branch
            float mle = fitBranch(tree, bgfreq, initdist);
	    node1->dist = mle / 2.0;
	    node2->dist = mle / 2.0;
	

	    // recompute the root node row in lktable
	    calcLkTableRow(seqlen, *model, 
			   lktable[node1->name], 
			   lktable[node2->name], 
			   lktable[tree->root->name],
			   node1->dist, 
			   node2->dist);
	
	    // get total probability after branch change    
	    logl = getTotalLikelihood(lktable, tree, 
				      seqlen, *model, bgfreq);
	
	    // don't accept a new branch length if it lowers total likelihood
	    if (logl < loglBefore) {
	         // revert
	    	 node1->dist = initdist / 2.0;
	    	 node2->dist = initdist / 2.0;
	    	 logl = loglBefore;
	    	 //assert(0);
	    }
        
	    printLog(LOG_HIGH, "hky: lk=%f\n", logl);        
	}

	return logl;
    }


    floatlk fitBranchesConverge(Tree *tree, int nseqs, char **seqs, 
                                const float *bgfreq,
                                ExtendArray<Node*> &rootingOrder,
                                int maxiter=10)
    {
        floatlk lastLogl = -INFINITY, logl = -INFINITY;
        const floatlk converge = logf(1.002);
    
        // initialize the condition likelihood table
        int seqlen = strlen(seqs[0]);
        calcLkTable(table.lktable, tree, nseqs, seqlen, seqs, *model);


        // remember original rooting for restoring later
        Node *origroot1 = tree->root->children[0];
        Node *origroot2 = tree->root->children[1];
    
        // iterate over branches improving each likelihood
        for (int j=0; j<maxiter; j++) {
            printLog(LOG_HIGH, "hky: iter %d\n", j);  

            logl = fitBranches(tree, nseqs, seqlen, seqs, bgfreq, 
                               rootingOrder);
        
            // determine whether logl has converged
            floatlk diff = fabs(logl - lastLogl);
            if (diff < converge) {
                //printf("conv %d %f %f\n", j, logl, lastLogl);
                printLog(LOG_HIGH, "hky: diff = %f < %f\n", diff, converge);
                break;
            } else {
                printLog(LOG_HIGH, "hky: diff = %f > %f\n", diff, converge);
            }
            lastLogl = logl;
        }
    
        // restore original rooting
        tree->reroot(origroot1, origroot2);

        return logl;
    }


    void setBranchRange(double _minx, double _maxx)
    {
        minx = _minx;
        maxx = _maxx;
    }


    double minx;
    double maxx;

    gsl_root_fdfsolver *opt;
    gsl_function_fdf opt_func;

    LikelihoodTable table;
    Model *model;
    typename Model::Deriv *dmodel;
    typename Model::Deriv::Deriv *d2model;
    
    DistLikelihoodDeriv<Model, typename Model::Deriv> lk_deriv;
    DistLikelihoodDeriv2<Model, typename Model::Deriv, 
                         typename Model::Deriv::Deriv> lk_deriv2;
};



// NOTE: assumes binary Tree
template <class Model>
floatlk findMLBranchLengths(Tree *tree, int nseqs, char **seqs, 
                            const float *bgfreq, Model &model,
                            int maxiter=10, 
                            double minlen=0.0, double maxlen=10.0)
{
    // timing
    Timer timer;
    

    int seqlen = strlen(seqs[0]);
    Timer timer2;
    MLBranchAlgorithm<Model> mlalg(tree, seqlen, &model);
    mlalg.setBranchRange(minlen, maxlen);
    printLog(LOG_MEDIUM, "mlalloc time: %f\n", timer2.time());
    
    
    // determine rooting order
    ExtendArray<Node*> rootingOrder(0, 2*tree->nnodes);
    getRootOrder(tree, &rootingOrder);
    
    // perform fitting
    floatlk logl = mlalg.fitBranchesConverge(tree, nseqs, seqs, bgfreq,
                                             rootingOrder, maxiter);
    
    
    printLog(LOG_MEDIUM, "mldist time: %f\n",  timer.time());
    
    return logl;
}


double findMLBranchLengthsHky(Tree *tree, int nseqs, char **seqs, 
                              const float *bgfreq, float kappa, int maxiter,
                              double minlen, double maxlen)
{
    HkyModel hky(bgfreq, kappa);
    return findMLBranchLengths(tree, nseqs, seqs, bgfreq, hky, maxiter,
                               minlen, maxlen);
}



double findMLKappaHky(Tree *tree, int nseqs, char **seqs, 
                      const float *bgfreq, float minkappa, float maxkappa,
                      float kappastep)
{
    const int maxiter = 1;
    floatlk maxlk = -INFINITY;
    float maxk = minkappa;

    // special case
    if (nseqs < 2)
        return 0.0;

    for (float k=minkappa; k<=maxkappa; k+=kappastep) {
        float l = findMLBranchLengthsHky(tree, nseqs, seqs, bgfreq, k, 
                                         maxiter);
        if (l > maxlk) {
            maxlk = l;
            maxk = k;
        }
    }

    return maxk;
}



extern "C" {

floatlk findMLBranchLengthsHky(int nnodes, int *ptree, int nseqs, char **seqs, 
                               float *dists, const float *bgfreq, float ratio, 
                               int maxiter, bool parsinit)
{
    //int seqlen = strlen(seqs[0]);
        
    gsl_set_error_handler_off();

    // create tree objects
    Tree tree(nnodes);
    ptree2tree(nnodes, ptree, &tree);
    tree.setDists(dists);

    if (parsinit)
        parsimony(&tree, nseqs, seqs);
    
    floatlk logl = findMLBranchLengthsHky(&tree, nseqs, seqs, bgfreq, 
                                          ratio, maxiter);
    tree.getDists(dists);
    
    return logl;
}


} // extern C


} // namespace spidir
