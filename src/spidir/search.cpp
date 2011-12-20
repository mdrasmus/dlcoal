/*=============================================================================

  Matt Rasmussen
  Copyright 2007-2011

  Gene tree search functions

=============================================================================*/


#include "common.h"
#include "distmatrix.h"
#include "logging.h"
#include "Matrix.h"
#include "model.h"
#include "nj.h"
#include "parsimony.h"
#include "phylogeny.h"
#include "search.h"
#include "top_change.h"
#include "top_prior.h"
#include "treevis.h"


namespace spidir {

//=============================================================================


TreeSet::~TreeSet()
{
    clear();
}

void TreeSet::clear()
{
    for (Set::iterator it=trees.begin();
         it != trees.end(); it++)
    {
        delete [] (int*) *it;
    }

    trees.clear();
}

bool TreeSet::insert(Tree *tree)
{
    int *key2 = new int [tree->nnodes+1];
    tree->hashkey(key2);
    key2[tree->nnodes] = -2; // cap key

    Set::iterator it = trees.find(key2);
    if (it == trees.end()) {
        trees.insert(key2);
        return true;
    } else {
        delete [] key2;
        return false;
    }
}

bool TreeSet::has(Tree *tree)
{
    key.ensureSize(tree->nnodes+1);
    tree->hashkey(key);
    key[tree->nnodes] = -2; // cap key
    return trees.find(key) != trees.end();
}



//=============================================================================
// NNI Proposer

NniProposer::NniProposer(int niter) :
    niter(niter),
    iter(0),
    nodea(NULL),
    nodeb(NULL)
{}

void NniProposer::propose(Tree *tree)
{
    // increase iteration
    iter++;
    
    // propose new tree
    proposeRandomNni(tree, &nodea, &nodeb);
    performNni(tree, nodea, nodeb);
}

void NniProposer::revert(Tree *tree)
{
    // undo topology change
    performNni(tree, nodea, nodeb);
}


//=============================================================================
// SPR Proposer

SprProposer::SprProposer(int niter) :
    NniProposer(niter)
{
}

void SprProposer::propose(Tree *tree)
{   
    // increase iteration
    iter++;
    
    // choose a SPR move
    proposeRandomSpr(tree, &nodea, &nodeb);
    
    // remember sibling of nodea
    const Node *p = nodea->parent;
    nodec = (p->children[0] == nodea) ? p->children[1] : p->children[0];
    
    // perform SPR move
    performSpr(tree, nodea, nodeb);
}

void SprProposer::revert(Tree *tree)
{
    performSpr(tree, nodea, nodec);
}


//=============================================================================
// Mixuture of Proposers

void MixProposer::addProposer(TopologyProposer *proposer, float weight)
{
    totalWeight += weight;
    methods.push_back(Method(proposer, weight));
}

void MixProposer::propose(Tree *tree)
{
    // increase iteration
    iter++;

    // randomly choose method
    float choice = frand() * totalWeight;
    float sum = methods[0].second;
    unsigned int i = 0;
    while (i < methods.size()-1 && sum < choice) {
        i++;
        sum += methods[i].second;
    }

    // make proposal
    lastPropose = i;
    methods[i].first->propose(tree);
}

void MixProposer::revert(Tree *tree)
{
    methods[lastPropose].first->revert(tree);
}



//=============================================================================
// SPR Neighborhood Proposer

SprNbrProposer::SprNbrProposer(int niter, int radius) :
    NniProposer(niter),
    radius(radius),
    basetree(NULL),
    reverted(false)
{
}

void SprNbrProposer::propose(Tree *tree)
{
    // ensure the same tree is used for each proposal
    if (!basetree)
        basetree = tree;
    else
        assert(basetree == tree);


    // start a new subtree
    if (iter == 0 || queue.size() == 0 || !reverted) {
        iter = 0;
        pickNewSubtree();
    }

    iter++; // increase iteration
    
    // go through skip
    while (queue.size() > 0) {
        // get new branch point
        nodea = queue.front();
        queue.pop_front();
    
        // remember sibling of subtree (nodeb)
        const Node *p = subtree->parent;
        nodeb = (p->children[0] == subtree) ? p->children[1] : p->children[0];
    
        // perform only valid SPR moves
        // NOTE: the tree may have changed, thus we need to double check
        // whether the Spr is valid.
        if (validSpr(tree, subtree, nodea)) {
            performSpr(tree, subtree, nodea);
            break;
        }
    }

    assert(tree->assertTree());
}

void SprNbrProposer::revert(Tree *tree)
{
    performSpr(tree, subtree, nodeb);
    assert(tree->assertTree());
}

void SprNbrProposer::pickNewSubtree()
{
    const Tree *tree = basetree;

    assert(basetree->nnodes >= 5);

    // find subtree (a) to cut off (any node that is not root or child of root)
    int choice;
    do {
        choice = irand(tree->nnodes);
    } while (tree->nodes[choice]->parent == NULL ||
             tree->nodes[choice]->parent->parent == NULL);
    Node *a = tree->nodes[choice];
    subtree = a;
    
    // find sibling (b) of a
    Node *c = a->parent;
    const int bi = (c->children[0] == a) ? 1 : 0;
    Node *b = c->children[bi];
    
    // uninitialize path distances
    pathdists.clear();
    for (int i=0; i<tree->nnodes; i++)
        pathdists.push_back(-1);
    
    // setup path distances and queue
    pathdists[a->name] = 0;
    pathdists[c->name] = 0;
    pathdists[b->name] = 0;
    queue.clear();
    list<Node*> tmpqueue;
    tmpqueue.push_back(c);
    tmpqueue.push_back(b);

    // traverse tree via depth first traversal
    while (tmpqueue.size() > 0) {
        Node *n = tmpqueue.front();
        tmpqueue.pop_front();
        
        // do not traverse beyond radius
        if (pathdists[n->name] >= radius)
            continue;

        // queue only valid new branch points:
        // n must not be root, a, descendant of a, c (parent of a), or  
        // b (sibling of a)
        if (n->parent && n != subtree && n != b && n != c) {
            queue.push_back(n);
        }

        // queue up unvisited neighboring edges
        Node *w = n->parent;

        if (w && pathdists[w->name] == -1) {
            pathdists[w->name] = pathdists[n->name] + 1;
            tmpqueue.push_back(w);
        }

        if (n->nchildren == 2) {
            Node *u = n->children[0];
            Node *v = n->children[1];

            if (pathdists[u->name] == -1) {
                pathdists[u->name] = pathdists[n->name] + 1;
                tmpqueue.push_back(u);
            }

            if (pathdists[v->name] == -1) {
                pathdists[v->name] = pathdists[n->name] + 1;
                tmpqueue.push_back(v);
            }
        }
    }
}


//=============================================================================
// Recon root proposer

void ReconRootProposer::propose(Tree *tree)
{
    const float rerootProb = 1.0;
    
    // propose new tree
    proposer->propose(tree);
    
    // reroot tree if stree is given
    if (frand() < rerootProb) {
        oldroot1 = tree->root->children[0];
        oldroot2 = tree->root->children[1];
        
        if (stree != NULL) {
            reconRoot(tree, stree, gene2species);
        }
    } else {
        oldroot1 = NULL;
        oldroot2 = NULL;
    }
}

void ReconRootProposer::revert(Tree *tree)
{
    // undo topology change
    if (oldroot1)
        tree->reroot(oldroot1, oldroot2);
    
    proposer->revert(tree);
}


//=============================================================================
// Dup/Loss proposer

DupLossProposer::DupLossProposer(TopologyProposer *proposer, 
                                 SpeciesTree *stree, int *gene2species,
                                 float dupprob, float lossprob,
                                 int quickiter, int niter) :
    proposer(proposer),
    quickiter(quickiter),
    niter(niter),
    iter(0),
    correctTree(NULL),
    correctSeen(false),
    stree(stree),
    gene2species(gene2species),
    dupprob(dupprob),
    lossprob(lossprob),
    recon(0),
    events(0),
    oldtop(NULL)
{
    doomtable = new double [stree->nnodes];
    calcDoomTable(stree, dupprob, lossprob, doomtable);

}


DupLossProposer::~DupLossProposer()
{
    delete [] doomtable;
}


void DupLossProposer::propose(Tree *tree)
{
    iter++;
    
    // do simple proposal if dup/loss probs are disabled
    if (dupprob < 0.0 || lossprob < 0.0 || quickiter <= 1) {
        proposer->propose(tree);
        return;
    }
    
    // save old topology
    oldtop = tree->copy();
    
    
    // recon tree to species tree
    recon.ensureSize(tree->nnodes);
    events.ensureSize(tree->nnodes);
    recon.setSize(tree->nnodes);
    events.setSize(tree->nnodes);
    
    ExtendArray<Tree*> trees(0, quickiter);
    ExtendArray<float> logls(0, quickiter);
    
    reconcile(tree, stree, gene2species, recon);
    labelEvents(tree, recon, events);
    double bestlogp = birthDeathTreePriorFull(tree, stree, recon, events, 
                                              dupprob, lossprob,
                                              doomtable);


    double sum = -INFINITY;

    // make many subproposals
    proposer->reset();
    for (int i=0; i<quickiter; i++) {
        proposer->propose(tree);
        
        // only allow unique proposals
        // but if I have been rejecting too much allow some non-uniques through
        if (uniques.has(tree) && trees.size() >= .1 * i) {
            proposer->revert(tree);
            continue;
        }

        reconcile(tree, stree, gene2species, recon);
        labelEvents(tree, recon, events);
        double logp = birthDeathTreePriorFull(tree, stree, recon, events, 
                                              dupprob, lossprob,
                                              doomtable);
        printLog(LOG_HIGH, "search: qiter %d %f %f\n", i, logp, bestlogp);
        
        Tree *tree2 = tree->copy();
        
        // save tree and logl
        trees.append(tree2);
        logls.append(logp);
        sum = logadd(sum, logp);
        
        if (logp > bestlogp)
            // make more proposals off this one
            bestlogp = logp;
        else
            proposer->revert(tree);
    }    
    
    // propose one of the subproposals 
    double choice = frand();
    double partsum = -INFINITY;
    
    for (int i=0; i<trees.size(); i++) {
        partsum = logadd(partsum, logls[i]);
    
        //printf("part %d %f (%f)\n", i, expf(partsum - sum), choice);
        
        if (choice < exp(partsum - sum)) {
            // propose tree i
            printLog(LOG_MEDIUM, "search: choose %d %f %f\n", i, 
                     logls[i], exp(logls[i] - sum));
            tree->setTopology(trees[i]);            
            break;
        }
        
    }

    // add tree to unqiues
    uniques.insert(tree);

    // clean up subproposals
    for (int i=0; i<trees.size(); i++)
        delete trees[i];

}

void DupLossProposer::revert(Tree *tree)
{
    // do simple proposal if dup/loss probs are disabled
    if (dupprob < 0.0 || lossprob < 0.0 || quickiter <= 1) {
        proposer->revert(tree);
        return;
    }
    
    //printf("set oldtop\n");
    tree->setTopology(oldtop);
    delete oldtop;
}



//=============================================================================

UniqueProposer::~UniqueProposer()
{
    seenTrees.clear();
}


void UniqueProposer::propose(Tree *tree)
{
    iter++;
    
    for (int i=0;; i++) {
        printLog(LOG_HIGH, "search: unique trees seen %d (tries %d)\n", 
                 seenTrees.size(), i+1);

        proposer->propose(tree);
        
        if (seenTrees.insert(tree)) {
            // return new tree
            break;
        } else {
            if (i < ntries) {
                // revert and loop again
                proposer->revert(tree);
            } else {
                // give up and return tree
                break;
            }
        }
    }
}


//=============================================================================
// get initial tree

// propose initial tree by Neighbor Joining
Tree *getInitialTree(string *genes, int nseqs, int seqlen, char **seqs)
{
    int nnodes = nseqs * 2 - 1;

    ExtendArray<int> ptree(nnodes);
    ExtendArray<float> dists(nnodes);
    Matrix<float> distmat(nseqs, nseqs);

    calcDistMatrix(nseqs, seqlen, seqs, distmat.getMatrix());
    neighborjoin(nseqs, distmat.getMatrix(), ptree, dists);

    Tree *tree = new Tree(nnodes);
    ptree2tree(nnodes, ptree, tree);
    tree->setLeafNames(genes);

    // special case
    if (nseqs == 2) {
        tree->nodes[0]->dist = dists[0];
        tree->nodes[1]->dist = dists[1];
    }

    return tree;
}


// propose initial tree and root by species
Tree *getInitialTree(string *genes, int nseqs, int seqlen, char **seqs,
                     SpeciesTree *stree, int *gene2species)
{
    Tree *tree = getInitialTree(genes, nseqs, seqlen, seqs);
    reconRoot(tree, stree, gene2species);

    return tree;
}


//=============================================================================
// search logging

void printSearchStatus(Tree *tree, SpeciesTree *stree, int *gene2species,
                       int *recon=NULL, int *events=NULL)
{
    if (stree) {
        bool cleanupRecon = false;
        if (!recon) {
            recon = new int [tree->nnodes];
            cleanupRecon = true;
        }

        bool cleanupEvents = false;    
        if (!events) {
            events = new int [tree->nnodes];
            cleanupEvents = true;
        }    

        //assert(tree->assertTree());
    
        reconcile(tree, stree, gene2species, recon);
        labelEvents(tree, recon, events);
        int losses = countLoss(tree, stree, recon);        
        
        // count dups
        int dups = 0;
        for (int i=0; i<tree->nnodes; i++)
            if (events[i] == EVENT_DUP)
                dups++;
        
        printLog(LOG_LOW, "search: dups = %d\n", dups);
        printLog(LOG_LOW, "search: loss = %d\n", losses);
        
        if (cleanupRecon)
            delete [] recon;
        if (cleanupEvents)
            delete [] events;        
    }                
    
    if (isLogLevel(LOG_LOW))
        displayTree(tree, getLogFile());
}


void printLogTree(int loglevel, Tree *tree)
{
    if (isLogLevel(loglevel)) {
        printLog(loglevel, "tree: ");
        writeNewickTree(getLogFile(), tree, 0, true);
        printLog(loglevel, "\n");
    }
}


class Prob
{
public:
    double seqlk;
    double branchp;
    double topp;
    double logp;

    double calcJoint(Model *model, Tree *tree)
    {
        model->setTree(tree);
        seqlk = model->likelihood();
        branchp = model->branchPrior();
        topp = model->topologyPrior();
        logp = seqlk + branchp + topp;
        return logp;
    }
};


void printLogProb(int loglevel, Prob *prob)
{
    printLog(loglevel, "search: lnl    = %f\n", prob->logp);
    printLog(loglevel, "search: seqlk  = %f\n", prob->seqlk);
    printLog(loglevel, "search: branch = %f\n", prob->branchp);          
    printLog(loglevel, "search: top    = %f\n", prob->topp);
}


//=============================================================================
// Search Climb

TreeSearchClimb::TreeSearchClimb(Model *model, TopologyProposer *proposer) :
    model(model),
    proposer(proposer)
{
}


TreeSearchClimb::~TreeSearchClimb()
{}


Tree *TreeSearchClimb::search(Tree *initTree, string *genes, 
			      int nseqs, int seqlen, char **seqs)
{
    Tree *toptree = NULL;
    double toplogp = -INFINITY, nextlogp;
    Tree *tree = initTree;
    Timer correctTimer;    
    Timer proposalTimer;

    Prob prob;
    
    // setup search debug: testing against know correct tree
    Tree *correct = proposer->getCorrect();
    double correctLogp = -INFINITY;
    if (correct) {
        // determine probability of correct tree
        parsimony(correct, nseqs, seqs); // get initial branch lengths
        correctLogp = prob.calcJoint(model, correct);
        printLog(LOG_LOW, "search: correct tree lnl = %f\n", correctLogp);
    }

    
    // determine initial tree topology
    if (initTree == NULL)
        tree = getInitialTree(genes, nseqs, seqlen, seqs,
                              model->getSpeciesTree(), 
                              model->getGene2species());


    // special cases (1 and 2 leaves)
    if (nseqs < 3) {
        return tree->copy();
    }
    
    
    ExtendArray<int> recon(tree->nnodes);
    ExtendArray<int> events(tree->nnodes);
    
    // calc probability of initial tree
    parsimony(tree, nseqs, seqs); // get initial branch lengths
    toplogp = prob.calcJoint(model, tree);
    toptree = tree->copy();
    
    
    // log initial tree
    printLog(LOG_LOW, "search: initial\n");
    printLogProb(LOG_LOW, &prob);
    printLogTree(LOG_LOW, tree);
    printSearchStatus(tree, model->getSpeciesTree(), 
                      model->getGene2species(), &recon[0], &events[0]);
    
    
    int naccept = 0;
    int nreject = 0;
    
    // search loop
    proposer->reset();
    for (int i=0; proposer->more(); i++) {
        printLog(LOG_LOW, "search: iter %d\n", i);
    
        // propose new tree 
        proposalTimer.start();
        proposer->propose(tree);
        proposer->testCorrect(tree);
        proposal_runtime += proposalTimer.time();
        
        // calculate probability of proposal
        nextlogp = prob.calcJoint(model, tree);
        bool accept = (nextlogp > toplogp);

        // log proposal
        if (accept)
            printLog(LOG_LOW, "search: accept\n");
        else
            printLog(LOG_LOW, "search: reject\n");
        printLogProb(LOG_LOW, &prob);
        printLogTree(LOG_LOW, tree);

        // log for search debug
        if (correct && (nextlogp >= correctLogp || 
                        tree->sameTopology(correct)))
        {
            printLog(LOG_LOW, "search: correct tree time = %f\n", 
                     correctTimer.time());
            printLog(LOG_LOW, "search: correct tree logl = %f\n", 
                     correctLogp);
            printLog(LOG_LOW, "search: correct tree is best = %d\n",
                     int(tree->sameTopology(correct)));
            printLog(LOG_LOW, "search: correct tree better logl = %f\n", 
                     nextlogp);

            correct = NULL;
        }


        // act on acceptance
        if (accept) {
            naccept++;
            proposer->accept(true);
            toplogp = nextlogp;
            delete toptree;
            toptree = tree->copy();

            printSearchStatus(tree, 
                              model->getSpeciesTree(), 
                              model->getGene2species(), &recon[0], &events[0]);

        } else {           
            // display rejected tree
            if (isLogLevel(LOG_MEDIUM))
                printSearchStatus(tree, model->getSpeciesTree(), 
				  model->getGene2species(), 
                                  &recon[0], &events[0]);
            
            // reject, undo topology change
            nreject++;
            proposer->accept(false);             
            proposer->revert(tree);
        }

        printLog(LOG_LOW, "\n");
    }
    
    // print final log messages
    printLog(LOG_LOW, "accept rate: %f\n", naccept / double(naccept+nreject));

    // call probability break down again of top tree
    prob.calcJoint(model, toptree);
    
    // log final tree
    printLog(LOG_LOW, "search: final\n");
    printLogProb(LOG_LOW, &prob);
    
    // clean up
    if (initTree == NULL)
        delete tree;

    return toptree;
}




extern "C" {


// Calculate the likelihood of a tree
Tree *searchClimb(int niter, int quickiter,
		  int nseqs, char **gene_names, char **seqs,
		  int nsnodes, int *pstree, float *sdists,
		  int *gene2species,
		  float *sp_alpha, float *sp_beta, float generate,
		  float pretime_lambda, float birth, float death,
		  float gene_alpha, float gene_beta,
		  float *bgfreq, float kappa,
		  int nsamples, bool approx)
{
    int nnodes = 2*nseqs - 1;

    setLogLevel(LOG_MEDIUM);

    // create stree
    SpeciesTree stree(nsnodes);
    ptree2tree(nsnodes, pstree, &stree);
    stree.setDepths();
    stree.setDists(sdists);
    
    // params
    SpidirParams params(nsnodes, NULL, sp_alpha, sp_beta, 
			gene_alpha, gene_beta, pretime_lambda);
    
    // make gene names
    string genes[nseqs];
    //char s[101];
    for (int i=0; i<nseqs; i++) {
	//snprintf(s, 100, "%d", i);
	genes[i] = gene_names[i];
    }

    // model
    const int maxiter = 2;
    int seqlen = strlen(seqs[0]);
    SpimapModel model(nnodes, &stree, &params, gene2species,
		      pretime_lambda, birth, death, nsamples, approx, false);
    model.setLikelihoodFunc(new HkySeqLikelihood(nseqs, seqlen, seqs, 
                                                 bgfreq, kappa, maxiter));
    // proposers
    NniProposer nni(niter);
    SprProposer spr(niter);
    MixProposer mix(niter);
    mix.addProposer(&nni, .5);
    mix.addProposer(&spr, .5);
    ReconRootProposer rooted(&mix, &stree, gene2species);
    UniqueProposer unique(&rooted, niter);
    DupLossProposer proposer(&mix, &stree, gene2species, 
			     birth, death,
                             quickiter, niter);
    
    TreeSearchClimb search(&model, &proposer);

    Tree *tree = search.search(NULL, genes, nseqs, seqlen, seqs);
    
    return tree;
}

} // extern "C"


} // namespace spidir
