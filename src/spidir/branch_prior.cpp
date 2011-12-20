/*=============================================================================

  Matt Rasmussen
  Copyright 2010-2011

  Branch length prior of the SPIMAP model

=============================================================================*/


//#include <gsl/gsl_integration.h>

// c++ headers
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <time.h>

// spidir headers
#include "common.h"
#include "ExtendArray.h"
#include "birthdeath.h"
#include "branch_prior.h"
#include "gamma.h"
#include "logging.h"
#include "phylogeny.h"
#include "seq_likelihood.h"


namespace spidir {


BranchParams NULL_PARAM;




// TODO: remove
float estimateGeneRate(Tree *tree, SpeciesTree *stree, 
                       int *recon, int *events, SpidirParams *params)
{    
    return 0.0;
}




//=============================================================================
// branch prior functions


double approxGammaSum(int nparams, double x, float *gs_alpha, float *gs_beta,
                      bool approx)
{
    const double tol = .001;
    const double minfrac = .01;

    // filter for extreme parameters
    double mean = 0.0;
    for (int i=0; i<nparams; i++) {
	if (isinf(gs_beta[i]) || isnan(gs_beta[i])) {
	    gs_alpha[i] = gs_alpha[--nparams];
	    gs_beta[i] = gs_beta[nparams];
	    i--;
	} else {
	    mean += gs_alpha[i] / gs_beta[i];
	}
    }

    // remove params with effectively zero mean
    double mu = 0.0;
    double var = 0.0;
    for (int i=0; i<nparams; i++) {	
	if (gs_alpha[i] / gs_beta[i] < minfrac * mean) {
	    gs_alpha[i] = gs_alpha[--nparams];
	    gs_beta[i] = gs_beta[nparams];
	    i--;
	} else {
	    mu += gs_alpha[i] / gs_beta[i];
	    var += gs_alpha[i] / gs_beta[i] / gs_beta[i];
	}
    }    
    

    // there is nothing to do
    if (nparams == 0)
	return -INFINITY;
    
    //float dist = node->dist / generate;
    
    double logp;
    if (approx) {
	// approximation
	double a2 = mean*mean/var;
	double b2 = mean/var;
	logp = gammalog(x, a2, b2);
    } else {
	logp = log(gammaSumPdf(x, nparams, gs_alpha, gs_beta, tol));
    }
    
    if (isnan(logp))
	logp = -INFINITY;
    return logp;
}



// Generate a random sample of duplication points
void setRandomMidpoints(int root, Tree *tree, SpeciesTree *stree,
                        Node **subnodes, int nsubnodes, 
                        int *recon, int *events, 
                        ReconParams *reconparams,
			float birth, float death)
{
    const float esp = .001;

    // deal with pre-duplications
    if (root == tree->root->name) {
	do {
	    reconparams->pretime = 
		expovariate(reconparams->params->pretime_lambda);
	} while (reconparams->pretime <= esp);
    }

    for (int i=0; i<nsubnodes; i++) {
	Node *node = subnodes[i];
        
        if (events[node->name] == EVENT_DUP) {
            float lastpoint;
            
            if (node->parent != NULL && 
		recon[node->name] == recon[node->parent->name])
                // if I'm the same species branch as my parent 
                // then he is my last midpoint
                lastpoint = reconparams->midpoints[node->parent->name];
            else
                // i'm the first on this branch so the last midpoint is zero
                lastpoint = 0.0;
            
            // pick a midpoint based on a DB process
            float remain = 1.0 - lastpoint;
	    int snode = recon[node->name];
	    float time;

	    if (snode == stree->root->name) {
		// use preduplication time for pre-duplicates
		time = reconparams->pretime;
	    } else {
		time = stree->nodes[snode]->dist;
	    }

            const float k = sampleBirthWaitTime1(remain * time * 
                                                 (1.0 - esp), 
                                                 birth, death);
            reconparams->midpoints[node->name] = lastpoint +
                esp * remain + k / time;

            if (reconparams->midpoints[node->name] == 1.0) {
                reconparams->midpoints[node->name] = 
                    1.0 - ((1.0 - lastpoint) / 2.0);
            }

            // try again if samples lead to zero length branches
            if (reconparams->midpoints[node->name] == 1.0 ||
                reconparams->midpoints[node->name] == lastpoint)
            {
                return setRandomMidpoints(root, tree, stree,
                                          subnodes, nsubnodes, 
                                          recon, events, 
                                          reconparams,
                                          birth, death);
            }

        } else {
            // genes or speciations reconcile exactly to the end of the branch
            reconparams->midpoints[node->name] = 1.0;
        }
    }
}



// Reconcile a branch to the species tree
void reconBranch(int node, Tree *tree, SpeciesTree *stree,
		 int *recon, int *events, 
		 SpidirParams *params,
		 ReconParams *reconparams)
{
    Node** nodes = tree->nodes.get();
    Node** snodes = stree->nodes.get();
    
    ExtendArray<BranchPart> &parts = reconparams->parts[node];
    parts.clear();
    
    // set fractional branches
    if (recon[node] == recon[nodes[node]->parent->name]) {
        // start reconciles to a subportion of species branch
        if (events[node] == EVENT_DUP) {
            // only case k's are dependent
	    parts.append(BranchPart(recon[node], FRAC_DIFF));
        } else {
	    parts.append(BranchPart(recon[node], FRAC_PARENT));
	}

	// end reconcilies to nothing
    } else {
        if (events[nodes[node]->parent->name] == EVENT_DUP) {
            // start reconciles to last part of species branch
            int snode = recon[nodes[node]->parent->name];
	    parts.append(BranchPart(snode, FRAC_PARENT));
        }
	// else: start reconcilies to nothing

        if (events[node] == EVENT_DUP) {
            // end reconciles to first part of species branch
	    parts.append(BranchPart(recon[node], FRAC_NODE));
        }
	// else: end reconciles to nothing
    }
    
    // set midparams
    if (recon[node] != recon[nodes[node]->parent->name]) {
        // we begin and end on different branches
        int snode;

        // determine most recent species branch which we fully recon to
        if (events[node] == EVENT_DUP)
            snode = snodes[recon[node]]->parent->name;
        else
            snode = recon[node];

        // walk up species spath until starting species branch
        // starting species branch is either fractional or NULL
        int parent_snode;
        if (nodes[node]->parent != NULL)
            parent_snode = recon[nodes[node]->parent->name];
        else
            parent_snode = -1;
	
        while (snode != parent_snode && snodes[snode]->parent != NULL) {
	    parts.append(BranchPart(snode, FRAC_ONE));
            snode = snodes[snode]->parent->name;
        }
    } 

}



// get times vector from midpoints
void getReconTimes(Tree *tree, SpeciesTree *stree, 
		   Node *node, 
		   ReconParams *reconparams,
		   ExtendArray<float> &times)
{
    const float *k = reconparams->midpoints;
    const int name = node->name;
    const ExtendArray<BranchPart> &parts = reconparams->parts[name];


    times.clear();   

    for (int i=0; i<parts.size(); i++) {
	float time = stree->nodes[parts[i].species]->dist;
	
	if (parts[i].species == stree->root->name)
	    time = reconparams->pretime;

	switch (parts[i].frac) {
	case FRAC_DIFF:
	    times.append((k[name] - k[node->parent->name]) *
			 time);
	    break;

	case FRAC_PARENT: {
	    float kp = k[node->parent->name];
	    times.append((1.0 - kp) * time);
            } break;

	case FRAC_NODE:
	    times.append(k[name] * time);
	    break;

	case FRAC_ONE:
	    times.append(time);
	    break;
	}
    }
}


// get gamma sum parameters
void getReconParams(Tree *tree, Node *node, 
		    ReconParams *reconparams,
		    float generate,
		    float *times,
		    float *gs_alpha,
		    float *gs_beta, int nparams)
{
    SpidirParams *params = reconparams->params;
    ExtendArray<BranchPart> &parts = reconparams->parts[node->name];

    for (int j=0; j<parts.size(); j++) {
	int snode = parts[j].species;
	gs_alpha[j] = params->sp_alpha[snode];
	gs_beta[j] = (params->sp_beta[snode] / (generate * times[j]));

        // DEBUG
        if (isinf(gs_beta[j])) {
            printf("> %f %f %f\n", params->sp_beta[snode], 
                   generate, times[j]);
            assert(0);
        }
    }
}



// Calculate branch probability
double branchprob(Tree *tree, SpeciesTree *stree, Node *node,
                  float generate, ReconParams *reconparams,
                  bool approx=true)
{
    // get times
    ExtendArray<float> times(0, tree->nnodes);
    getReconTimes(tree, stree, node, reconparams, times);
    int nparams = times.size();

    // get gammaSum terms 
    float gs_alpha[nparams];
    float gs_beta[nparams];
    getReconParams(tree, node, reconparams,
		   generate, times, gs_alpha, gs_beta, nparams);
   
    // compute gamma sum
    return approxGammaSum(nparams, node->dist, gs_alpha, gs_beta, approx);
}


// Calculate branch probability
double branchprobUnfold(Tree *tree, SpeciesTree *stree, 
                        float generate, ReconParams *reconparams,
                        bool approx=true)
{
    Node *node0 = tree->root->children[0];
    Node *node1 = tree->root->children[1];

    // get times
    ExtendArray<float> times0(0, tree->nnodes);
    ExtendArray<float> times1(0, tree->nnodes);
    getReconTimes(tree, stree, node0, reconparams, times0);
    getReconTimes(tree, stree, node1, reconparams, times1);
    int nparams = times0.size() + times1.size();

    // get gamma sum terms 
    float gs_alpha[nparams];
    float gs_beta[nparams];

    getReconParams(tree, node0, reconparams,
		   generate, times0, gs_alpha, gs_beta, times0.size());
    getReconParams(tree, node1, reconparams,
		   generate, times1, 
		   gs_alpha + times0.size(), 
		   gs_beta + times0.size(), times1.size());

    // compute gamma sum
    return approxGammaSum(nparams, node0->dist + node1->dist, 
			  gs_alpha, gs_beta, approx);
}


// get roots of speciation subtrees
void getSpecSubtrees(Tree *tree, int *events, ExtendArray<Node*> *rootnodes)
{
    for (int i=0; i<tree->nnodes; i++) {
	if (i == tree->root->name) {
	    rootnodes->append(tree->nodes[i]);
	} else if (events[i] == EVENT_SPEC) {
	    for (int j=0; j<2; j++) {
		rootnodes->append(tree->nodes[i]->children[j]);
	    }
	}
    }
}


// Returns True if duplications encounters
bool getSubtree(Node *node, int *events, ExtendArray<Node*> *subnodes)
{
    subnodes->append(node);
    bool dupsPresent = false;

    // recurse
    if (events[node->name] == EVENT_DUP || node->parent == NULL) {
	if (events[node->name] == EVENT_DUP)
	    dupsPresent = true;
	
        for (int i=0; i<node->nchildren; i++) 
            dupsPresent |= getSubtree(node->children[i], events, subnodes);
    }

    return dupsPresent;
}


class BranchPriorCalculator
{
public:
    BranchPriorCalculator(Tree *_tree,
			  SpeciesTree *_stree,
			  int *_recon, int *_events, 
			  SpidirParams *_params,
			  float _birth, float _death,
			  int _nsamples=1000, bool _approx=true) :
        tree(_tree),
        stree(_stree),
        recon(_recon),
        events(_events),
        params(_params),
	birth(_birth),
	death(_death),
	reconparams(_tree->nnodes, _params),
	rootnodes(0, _tree->nnodes),
	nsamples(_nsamples),
	approx(_approx),
        subnodes(0, _tree->nnodes),
        times(0, tree->nnodes)
    {
	// determine speciation subtrees	
	getSpecSubtrees(tree, events, &rootnodes);
    }
    
    
    ~BranchPriorCalculator()
    {
    }
    


    // Calculate branch probability
    double branchprob(Tree *tree, SpeciesTree *stree, Node *node,
                      float generate, ReconParams *reconparams,
                      bool approx=true)
    {
        // get times
        times.clear();
        getReconTimes(tree, stree, node, reconparams, times);
        int nparams = times.size();

        // get gammaSum terms 
        float gs_alpha[nparams];
        float gs_beta[nparams];
        getReconParams(tree, node, reconparams,
                       generate, times, gs_alpha, gs_beta, nparams);
   
        // compute gamma sum
        double p = approxGammaSum(nparams, node->dist, gs_alpha, gs_beta, approx);

        return p;
    }



    // subtree prior conditioned on divergence times
    double subtreeprior_cond(Tree *tree, SpeciesTree *stree,
                             int *recon, float generate,
                             ReconParams *reconparams,
                             ExtendArray<Node*> &subnodes,
                             bool unfold)
    {
	double logp = 0.0;
        
	// loop through all branches in subtree
	for (int j=0; j<subnodes.size(); j++) {
	    Node *node = subnodes[j];
	
	    //if (recon[node->name] == sroot)
	    //	continue;

	    if (node == tree->root)
		continue;

	    if (unfold) {
		if (node == tree->root->children[1]) {
		    // skip right branch
		    continue;
		} else if (node == tree->root->children[0]) {
		    // unfold left branch
		    logp += branchprobUnfold(tree, stree,
                                             generate, reconparams,
                                             approx);
		    continue;
		}
	    }
	    logp += branchprob(tree, stree, node, generate, reconparams,
                               approx);
	}
        
	return logp;
    }


    // Calculate the likelihood of a subtree
    double subtreeprior(int root, float generate, ReconParams *reconparams)
    {
   
	// set reconparams by traversing subtree
	subnodes.clear();
	bool dupsPresent = getSubtree(tree->nodes[root], events, &subnodes);
    
	// reconcile each branch
	for (int i=0; i<subnodes.size(); i++)
	    if (subnodes[i] != tree->root)
		reconBranch(subnodes[i]->name, tree, stree, 
			    recon, events, params, reconparams);
	
	// root branch must be unfolded
	bool unfold = (root == tree->root->name &&
		       tree->root->nchildren == 2);
                
        //printf("unfold %d\n", unfold);

	if (!dupsPresent) {

            setRandomMidpoints(root, tree, stree,
                               subnodes, subnodes.size(),
                               recon, events, reconparams,
                               birth, death);

	    return subtreeprior_cond(tree, stree, recon, 
				     generate, reconparams, subnodes,
				     unfold);
	} else {
	    // perform integration by sampling
	    double prob = 1.0;
            //RunningStat stat;

	    for (int i=0; i<nsamples; i++) {
		// propose a setting of midpoints
		setRandomMidpoints(root, tree, stree,
				   subnodes, subnodes.size(),
				   recon, events, reconparams,
				   birth, death);
                
		double sampleLogl = subtreeprior_cond(tree, stree, recon, 
					       generate, reconparams, 
					       subnodes,
					       unfold);
                
                prob = logadd(prob, sampleLogl);
		//prob += exp(sampleLogl) / nsamples;
                //stat.push(exp(sampleLogl));
                //printf("%d %f %f %f\n", i, sampleLogl, prob -log(i+1), 
                //       log(stat.sdev()));
	    }

	    return prob - log(nsamples);
	}
    }



    
    double calc_cond(float generate)
    {
        double logp = 0.0;
	
	// multiple the probability of each subtree
	for (int i=0; i<rootnodes.size(); i++) {
	    logp += subtreeprior(rootnodes[i]->name, generate, &reconparams);
	}
        
        // gene rate probability
        if (params->gene_alpha > 0 && params->gene_beta > 0)
            logp += loginvgammaPdf(generate, params->gene_alpha, params->gene_beta);
        
        printLog(LOG_HIGH, "generate: %f %f\n", generate, exp(logp));
        return logp;
    }
    
    inline double calcprob_cond(float generate)
    {
        return exp(calc_cond(generate));
    }


    double calc()
    {
	// TODO: make this equal portions of gamma
        double logp = -INFINITY;
        float mid = params->gene_beta / (params->gene_alpha - 1.0);
        float gstart = mid * 0.01;
        float gend = mid * 3.0;
        float step = (gend - gstart) / 20.0;
        
	// integrate over gene rates
        for (float g=gstart; g<gend; g+=step) {
            float gi = g + step / 2.0;
            double p = calc_cond(gi);

            logp = logadd(logp, p);
            printLog(LOG_HIGH, "generate_int: %f %f\n", gi, p);
        }
        
        // multiply probabilty by integration step
        logp += log(step);
	return logp;
    }
    
    
protected:
    Tree *tree;    
    SpeciesTree *stree;
    int *recon;
    int *events;
    SpidirParams *params;  
    float birth;
    float death;    
    ReconParams reconparams;
    ExtendArray<Node*> rootnodes;
    int nsamples;
    bool approx;

    ExtendArray<Node*> subnodes;
    ExtendArray<float> times;
};

  

// TODO: move birth, death into param
double branchPrior(Tree *tree,
                   SpeciesTree *stree,
                   int *recon, int *events, SpidirParams *params,
                   float generate,
                   float pretime_lambda, float birth, float death,
                   int nsamples, bool approx)
{
    double logl = 0.0; // log likelihood
    Timer timer;

    
    BranchPriorCalculator priorcalc(tree, stree, recon, events, params,
				    birth, death, nsamples, approx);
    
    if (generate > 0) {
        // estimate with given generate
        logl = priorcalc.calc_cond(generate);
    } else {
        logl = priorcalc.calc();
    }
    
    //setLogLevel(LOG_MEDIUM);
    printLog(LOG_MEDIUM, "branch prior time: %f\n", timer.time());
    
    return logl;
}


extern "C" {

// Calculate the likelihood of a tree
double branchPrior(int nnodes, int *ptree, float *dists,
                   int nsnodes, int *pstree, float *sdists,
                   int *recon, int *events,
                   float *sp_alpha, float *sp_beta, float generate,
                   float pretime_lambda, float birth, float death,
                   float gene_alpha, float gene_beta,
                   int nsamples, bool approx)
{
    // create tree objects
    Tree tree(nnodes);
    ptree2tree(nnodes, ptree, &tree);
    tree.setDists(dists);
    
    SpeciesTree stree(nsnodes);
    ptree2tree(nsnodes, pstree, &stree);
    stree.setDepths();
    stree.setDists(sdists);
    
    SpidirParams params(nsnodes, NULL, sp_alpha, sp_beta, 
			gene_alpha, gene_beta, pretime_lambda);
    

    return branchPrior(&tree, &stree,
		       recon, events, &params, 
		       generate, 
		       pretime_lambda, birth, death,
		       nsamples, approx);
}

} // extern "C"









//=============================================================================
// Gene rate estimation


// Ignores sequence data, assumes given branch lengths are perfectly known
float maxPosteriorGeneRate(Tree *tree, SpeciesTree *stree,
                           int *recon, int *events, SpidirParams *params)
{
    return 0.0;
}


float maxPosteriorGeneRate(int nnodes, int *ptree, float *dists,
                           int nsnodes, int *pstree, 
                           int *recon, int *events,
                           float *mu, float *sigma,
                           float alpha, float beta)
{
    return 0.0;
}


void samplePosteriorGeneRate(Tree *tree, 
                             int nseqs, char **seqs, 
                             const float *bgfreq, float ratio,
                             SpeciesTree *stree, 
                             int *gene2species,
                             SpidirParams *params,
                             int nsamples,
                             geneRateCallback callback,
                             void *userdata)
{
}


// TODO: repare
// Uses MCMC to sample from P(B,G|T,D)
void samplePosteriorGeneRate(Tree *tree,
                             int nseqs, char **seqs, 
                             const float *bgfreq, float ratio,
                             SpeciesTree *stree,
                             int *recon, int *events, SpidirParams *params,
                             int nsamples,
                             geneRateCallback callback,
                             void *userdata)
{
    // state variables B, G
    ExtendArray<float> dists(tree->nnodes);
    float next_generate = 0, generate = 0;
    float logl = -INFINITY;
    float logl2 = 999;
    float next_logl, next_logl2;
    float alpha;

    // TODO: revive
    assert(0);

    const float generate_step = .2;
    const float min_generate = .0001;

    // initialize first state
    generate = gammavariate(params->gene_alpha, params->gene_beta);
    // TODO: replace
    //generateBranchLengths(tree, stree,
    //                      recon, events,
    //                      params, generate);

    
    // perform MCMC
    for (int i=0; i<nsamples; i++) {
        // generate a new state 
	//printf("sample %d\n", i);
     
	for (int j=0; j<tree->nnodes; j++) {
	    assert(events[j] <= 2 && events[j] >= 0);

	    // NAN distances
	    if (isnan(tree->nodes[j]->dist)) {
		tree->nodes[j]->dist = min_generate;
	    }	    
	}
 
        if (frand() < .5) {
	    //printf("sample G\n");
            printLog(LOG_HIGH, "sample G: ");

            // sample G_2
            next_generate = frand(max(generate - generate_step, min_generate),
                                  generate + generate_step);
            //next_generate = gammavariate(params->alpha, params->beta);

            // if P(B_1|T,G_1) not exist, make it
            if (logl2 > 0) {                
                logl2 = branchPrior(tree, stree, recon, events, params,
				    generate,
				    1, 1, 1);
            }
	    
            // set new branch lengths B_2
            for (int j=0; j<tree->nnodes; j++) {
                tree->nodes[j]->dist *= next_generate / generate;

		// sometimes when generate change is drastic we need to catch
		// NAN distances
		if (isnan(tree->nodes[j]->dist)) {
		    tree->nodes[j]->dist = min_generate;
		}
	    }

            
            // calculate P(D|B_2,T)
            next_logl = calcSeqProbHky(tree, nseqs, seqs, bgfreq, ratio);
	    
            // calculate P(B_2|T,G_2)
            next_logl2 = branchPrior(tree, stree, recon, events, params,
				     next_generate,
				     1, 1, 1);

            // q(G_1|G_2) = 1 /(next_generate + generate_step - 
            //                  max(next_generate - generate_step, 
            //                      min_generate))
            // q(G_2|G_1) = 1 /(generate + generate_step - 
            //                  max(generate - generate_step, 
            //                      min_generate))
            // qratio = (generate + generate_step - 
            //                  max(generate - generate_step, 
            //                      min_generate)) /
            //          (next_generate + generate_step - 
            //                  max(next_generate - generate_step, 
            //                      min_generate))


            float logl3 = gammalog(next_generate, 
                                   params->gene_alpha, params->gene_beta) - 
                          gammalog(generate, 
                                   params->gene_alpha, params->gene_beta) +
                          log((generate + generate_step - 
                               max(generate - generate_step, 
                                   min_generate)) /
                              (next_generate + generate_step - 
                               max(next_generate - generate_step, 
                                   min_generate)));

            alpha = exp(next_logl + next_logl2 - logl - logl2 + logl3);

        } else {
	    //printf("sample B\n");
            printLog(LOG_HIGH, "sample B: ");

            // keep same gene rate G_2 = G_1
            next_generate = generate;
            
            // sample B_2
	    // TODO: replace
            //generateBranchLengths(tree, stree,
            //                      recon, events,
            //                      params, next_generate, -2, -2);
	    
	    for (int j=0; j<tree->nnodes; j++) {
		// NAN distances
		if (isnan(tree->nodes[j]->dist)) {
		    tree->nodes[j]->dist = min_generate;
		}	    
	    }

            // calculate P(D|B_2,T)
            next_logl = calcSeqProbHky(tree, nseqs, seqs, bgfreq, ratio);
            next_logl2 = 999;
            
            alpha = exp(next_logl - logl);
        }
        

        if (frand() <= alpha) {
            // accept: logl, G, B
            printLog(LOG_HIGH, "accept\n");
            logl = next_logl;
            logl2 = next_logl2;
            generate = next_generate;
            tree->getDists(dists);
        } else {
            // reject
            printLog(LOG_HIGH, "reject\n");
            
            // restore previous B
            tree->setDists(dists);
        }
        
        // return sample
        callback(generate, tree, userdata);
    }
}




// Uses MCMC to sample from P(B,G|T,D)
void samplePosteriorGeneRate_old(Tree *tree,
                             int nseqs, char **seqs, 
                             const float *bgfreq, float ratio,
                             SpeciesTree *stree,
                             int *recon, int *events, SpidirParams *params,
                             int nsamples,
                             geneRateCallback callback,
                             void *userdata)
{
    // state variables B, G
    ExtendArray<float> dists(tree->nnodes);
    float next_generate = 0, generate = 0;
    float logl = -INFINITY;
    float next_logl;
    
    // TODO: revive
    assert(0);

    for (int i=0; i<nsamples; i++) {
        // generate a new state 

        // sample G
        next_generate = gammavariate(params->gene_alpha, params->gene_beta);

        // sample B
	// TODO: replace
        //generateBranchLengths(tree, stree,
        //                      recon, events,
        //                      params, next_generate);
        
        // calculate P(D|B,T)
        next_logl = calcSeqProbHky(tree, nseqs, seqs, bgfreq, ratio);

        //printf(">> %f %f\n", next_logl, logl);

        float alpha = exp(next_logl - logl);
        if (frand() <= alpha) {
            // accept: logl, G, B
            printLog(LOG_HIGH, "accept: %f, %f\n", next_logl, logl);
            logl = next_logl;
            generate = next_generate;
            tree->getDists(dists);
        } else {
            // reject
            printLog(LOG_HIGH, "reject: %f, %f\n", next_logl, logl);
            
            // restore previous B
            tree->setDists(dists);
        }
        
        // return sample
        callback(generate, tree, userdata);
    }
}





} // namespace spidir
