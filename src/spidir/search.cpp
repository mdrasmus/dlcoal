/*=============================================================================

  Matt Rasmussen
  Copyright 2007-2011

  Gene tree search functions

=============================================================================*/


#include "common.h"
#include "branch_prior.h"
#include "distmatrix.h"
#include "logging.h"
#include "Matrix.h"
#include "nj.h"
#include "newick.h"
#include "parsimony.h"
#include "phylogeny.h"
#include "search.h"
#include "seq_likelihood.h"
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
// Nearest Neighbor Interchange Topology Proposal

/*

    Proposes a new tree using Nearest Neighbor Interchange
       
       Branch for NNI is specified by giving its two incident nodes (node1 and 
       node2).  Change specifies which  subtree of node1 will be swapped with
       the uncle.  See figure below.

         node2
        /     \
      nodeb    node1
               /  \
         nodea     * 

*/
void performNni(Tree *tree, Node *nodea, Node *nodeb)
{
    Node *node1 = nodea->parent;
    Node *node2 = nodeb->parent;
    
    // assert that node1 and node2 are incident to the same branch
    assert(node1->parent == node2 ||
           node2->parent == node1);
    
    // find child indexes
    int a = (node1->children[0] == nodea) ? 0 : 1;
    assert(node1->children[a] == nodea);

    int b = (node2->children[0] == nodeb) ? 0 : 1;
    assert(node2->children[b] == nodeb);
    
    // swap parent pointers
    nodea->parent = node2;
    nodeb->parent = node1;
    
    // swap child pointers
    node2->children[b] = nodea;
    node1->children[a] = nodeb;
}


void proposeRandomNni(Tree *tree, Node **a, Node **b)
{
    // find edges for NNI
    int choice;
    do {
        choice = irand(tree->nnodes);
    } while (tree->nodes[choice]->isLeaf() || 
             tree->nodes[choice]->parent == NULL);
    
    Node *node1 = tree->nodes[choice];
    Node *node2 = tree->nodes[choice]->parent;
    *a = node1->children[irand(2)];
    *b = (node2->children[0] == node1) ? node2->children[1] :
                                         node2->children[0];
    assert((*a)->parent->parent == (*b)->parent);
}


//=============================================================================
// Subtree pruning and regrafting (SPR)


/*
    a = subtree
    e = newpos
    
    BEFORE
            ....
        f         d
       /           \
      c             e
     / \           ...
    a   b
   ... ...

    AFTER

        f         d
       /           \
      b             c
     ...           / \
                  a   e
                 ... ...

    Requirements:
    1. a (subtree) is not root or children of root
    2. e (newpos) is not root, a, descendant of a, c (parent of a), or 
       b (sibling of a)
    3. tree is binary

*/
void performSpr(Tree *tree, Node *subtree, Node *newpos)
{
    Node *a = subtree;
    Node *e = newpos;

    Node *c = a->parent;
    Node *f = c->parent;
    const int bi = (c->children[0] == a) ? 1 : 0;
    Node *b = c->children[bi];
    const int ci = (f->children[0] == c) ? 0 : 1;
    Node *d = e->parent;
    const int ei = (d->children[0] == e) ? 0 : 1;

    d->children[ei] = c;
    c->children[bi] = e;
    f->children[ci] = b;
    b->parent = f;
    //    b->dist += c->dist;
    //e->dist /= 2.0;
    //    c->dist = e->dist;
    c->parent = d;
    e->parent = c;
}

/*
    What if e == f  (also equivalent to NNI) this is OK

    BEFORE
    
          d
         / \
        e  ...
       / \
      c  ...         
     / \           
    a   b
   ... ...

    AFTER
          d
         / \
        c
       / \
      a   e
     ... / \
        b  ...
       ...
       
  What if d == f  (also equivalent to NNI) this is OK
  
    BEFORE
          
        f
       / \
      c   e
     / \  ...
    a   b
   ... ...

    AFTER
          
        f
       / \
      b   c  
     ... / \ 
        a   e
       ... ...  
*/



/*
    Requirements:
    1. a (subtree) is not root or children of root
    2. e (newpos) is not root, a, descendant of a, c (parent of a), or 
       b (sibling of a)
    3. tree is binary
*/
void proposeRandomSpr(Tree *tree, Node **subtree, Node **newpos)
{
    assert(tree->nnodes >= 5);

    // find subtree (a) to cut off (any node that is not root or child of root)
    int choice;
        do {
        choice = irand(tree->nnodes);
    } while (tree->nodes[choice]->parent == NULL ||
             tree->nodes[choice]->parent->parent == NULL);
    Node *a = tree->nodes[choice];
    *subtree = a;
    
    // find sibling (b) of a
    const Node *c = a->parent;
    const int bi = (c->children[0] == a) ? 1 : 0;
    const Node *b = c->children[bi];
    
    // choose newpos (e)
    Node *e = NULL;
    do {
        choice = irand(tree->nnodes);
        e = tree->nodes[choice];
        
        // test if e is a valid choice
        if (e->parent == NULL || e == a || e == c || e == b)
            continue;
        
        // also test if e is a descendent of a
        bool under_a = false;
        for (Node *ptr = e->parent; ptr != NULL; ptr = ptr->parent) {
            if (ptr == a) {
                under_a = true;
                break;
            }
        }            
        
        if (under_a)
            continue;
        
        break;
    } while (true);
    *newpos = e;
}


bool validSpr(Tree *tree, const Node *subtree, const Node *newpos)
{
    const Node *a = subtree;
    const Node *e = newpos;
    
    // find sibling (b) of a
    const Node *c = a->parent;
    const int bi = (c->children[0] == a) ? 1 : 0;
    const Node *b = c->children[bi];

    // test if a is a valid choice
    if (a->parent == NULL || a->parent->parent == NULL)
        return false;

    // test if e is a valid choice
    if (e->parent == NULL || e == a || e == c || e == b) {
        //printf("ERROR: a=%d, e=%d, b=%d, c=%d\n", a->name, e->name, 
        //       b->name, c->name);
        return false;
    }
        
    // also test if e is a descendent of a
    for (Node *ptr = e->parent; ptr != NULL; ptr = ptr->parent)
        if (ptr == a)
            return false;
        
    return true;
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
    
    nodea = nodeb = NULL;

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
// Fitting branch lengths


HkyFitter::HkyFitter(int nseqs, int seqlen, char **seqs, 
                     float *bgfreq, float tsvratio, int maxiter,
                     double minlen, double maxlen) :
    nseqs(nseqs),
    seqlen(seqlen),
    seqs(seqs),
    bgfreq(bgfreq),
    tsvratio(tsvratio),
    maxiter(maxiter),
    minlen(minlen),
    maxlen(maxlen)
{}


double HkyFitter::findLengths(Tree *tree)
{ 
    Timer timer;
    double logl = findMLBranchLengthsHky(tree, nseqs, seqs, bgfreq, 
                                         tsvratio, maxiter, minlen, maxlen);
    runtime += timer.time();

    return logl;
}


//=============================================================================
// Prior function function

SpidirPrior::SpidirPrior(
    int nnodes, SpeciesTree *stree, 
    SpidirParams *params, 
    int *gene2species,
    float predupprob, float dupprob, float lossprob, 
    int nsamples, bool approx, bool useBranchPrior) :
    
    nnodes(nnodes),
    stree(stree),
    params(params),
    gene2species(gene2species),
    recon(nnodes),
    events(nnodes),
    predupprob(predupprob),
    dupprob(dupprob),
    lossprob(lossprob),
    nsamples(nsamples),
    approx(approx),
    useBranchPrior(useBranchPrior)
{
    doomtable = new double [stree->nnodes];
    calcDoomTable(stree, dupprob, lossprob, doomtable);
}


SpidirPrior::~SpidirPrior()
{
    delete [] doomtable;
}

double SpidirPrior::branchPrior(Tree *tree)
{
    Timer timer;

    if (!useBranchPrior)
        return 0.0;

    // reconcile tree to species tree
    reconcile(tree, stree, gene2species, recon);
    labelEvents(tree, recon, events);
    float generate = -99;
        
    double logp = spidir::branchPrior(tree, stree,
                                     recon, events, params,
                                     generate, predupprob, dupprob, lossprob,
                                     nsamples, approx);
    branch_runtime += timer.time();

    return logp;
}


double SpidirPrior::topologyPrior(Tree *tree)
{
    Timer timer;

    // DEBUG
    //assert(tree->assertTree());

    // reconcile tree to species tree
    reconcile(tree, stree, gene2species, recon);
    labelEvents(tree, recon, events);
    

    //reconAssert(tree, stree, recon);
    
    double logp = birthDeathTreePriorFull(tree, stree, recon, events, 
                                         dupprob, lossprob,
                                         doomtable);

    top_runtime += timer.time();

    return logp;
}



//=============================================================================


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


//=============================================================================
// Search Climb

TreeSearchClimb::TreeSearchClimb(Prior *prior,
				 TopologyProposer *proposer,
				 BranchLengthFitter *fitter) :
    prior(prior),
    proposer(proposer),
    fitter(fitter)
{
}


TreeSearchClimb::~TreeSearchClimb()
{}


void printLogTree(int loglevel, Tree *tree)
{
    if (isLogLevel(loglevel)) {
        printLog(loglevel, "tree: ");
        writeNewickTree(getLogFile(), tree, 0, true);
        printLog(loglevel, "\n\n");
    }
}


Tree *TreeSearchClimb::search(Tree *initTree, 
			      string *genes, 
			      int nseqs, int seqlen, char **seqs)
{
    Tree *toptree = NULL;
    double toplogp = -INFINITY, nextlogp;
    Tree *tree = initTree;
    Timer correctTimer;    
    Timer proposalTimer;
    
    // search testing
    Tree *correct = proposer->getCorrect();
    double correctLogp = -INFINITY;
    if (correct) {
        // determine probability of correct tree
        parsimony(correct, nseqs, seqs); // get initial branch lengths
        double seqlk = fitter->findLengths(correct);
        double branchp = prior->branchPrior(correct);
        double topp = prior->topologyPrior(correct);
        correctLogp = seqlk + branchp + topp;

        printLog(LOG_LOW, "search: correct tree lnl = %f\n", correctLogp);
    }

    
    // determine initial tree
    if (initTree == NULL)
        tree = getInitialTree(genes, nseqs, seqlen, seqs,
                              prior->getSpeciesTree(), 
                              prior->getGene2species());
    
    ExtendArray<int> recon(tree->nnodes);
    ExtendArray<int> events(tree->nnodes);
    
    // calc probability of initial tree
    parsimony(tree, nseqs, seqs); // get initial branch lengths
    double seqlk = fitter->findLengths(tree);
    double branchp = prior->branchPrior(tree);
    double topp = prior->topologyPrior(tree);
    toplogp = seqlk + branchp + topp;
    toptree = tree->copy();
    
    
    // log initial tree
    printLog(LOG_LOW, "search: initial\n");
    printLog(LOG_LOW, "search: lnl    = %f\n", toplogp);
    printLog(LOG_LOW, "search: seqlk  = %f\n", seqlk);
    printLog(LOG_LOW, "search: branch = %f\n", branchp);          
    printLog(LOG_LOW, "search: top    = %f\n", topp);
    printLogTree(LOG_LOW, tree);
    printSearchStatus(tree, prior->getSpeciesTree(), 
                      prior->getGene2species(), &recon[0], &events[0]);
    
    
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
        
	if (isLogLevel(LOG_MEDIUM)) {
	    displayTree(tree, getLogFile());
	}

        // calculate likelihood
        seqlk = fitter->findLengths(tree);
        branchp = prior->branchPrior(tree);
        topp = prior->topologyPrior(tree);
        nextlogp = seqlk + branchp + topp;
        

        // search test
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

        // acceptance rule
        bool accept = (nextlogp > toplogp);

        // print accept
        if (accept)
            printLog(LOG_LOW, "search: accept\n");
        else
            printLog(LOG_LOW, "search: reject\n");


        // print info
        printLog(LOG_LOW, "search: lnl    = %f\n", nextlogp);
        printLog(LOG_LOW, "search: seqlk  = %f\n", seqlk);
        printLog(LOG_LOW, "search: branch = %f\n", branchp);          
        printLog(LOG_LOW, "search: top   = %f\n", topp);
        printLogTree(LOG_LOW, tree);

        // act on acceptance
        if (accept) {
            // accept
            naccept++;
            proposer->accept(true);
            toplogp = nextlogp;
            delete toptree;
            toptree = tree->copy();

            printSearchStatus(tree, 
                              prior->getSpeciesTree(), 
                              prior->getGene2species(), &recon[0], &events[0]);

        } else {           
            // reject, undo topology change
            
            // display rejected tree
            if (isLogLevel(LOG_MEDIUM))
                printSearchStatus(tree, prior->getSpeciesTree(), 
				  prior->getGene2species(), 
                                  &recon[0], &events[0]);
            
            nreject++;
            proposer->accept(false);             
           proposer->revert(tree);
        }
    }
    
    // print final log messages
    printLog(LOG_LOW, "accept rate: %f\n", naccept / double(naccept+nreject));

    // call probability break down again of top tree
    seqlk = fitter->findLengths(toptree);
    branchp = prior->branchPrior(toptree);
    topp = prior->topologyPrior(toptree);
    double logp = seqlk + branchp + topp;
    
    // probability stats on final tree
    printLog(LOG_LOW, "search: final\n");
    printLog(LOG_LOW, "search: lnl    = %f\n", logp);
    printLog(LOG_LOW, "search: seqlk  = %f\n", seqlk);
    printLog(LOG_LOW, "search: branch = %f\n", branchp);          
    printLog(LOG_LOW, "search: top    = %f\n", topp);
    
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

    // prior
    SpidirPrior prior(nnodes, &stree, &params, gene2species,
		      pretime_lambda, birth, death, nsamples, approx, false);

    // seq likelihood
    const int maxiter = 2;
    int seqlen = strlen(seqs[0]);
    HkyFitter fitter(nseqs, seqlen, seqs, 
		     bgfreq, kappa, maxiter);

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
    
    TreeSearchClimb search(&prior, &proposer, &fitter);

    Tree *tree = search.search(NULL, genes, nseqs, seqlen, seqs);
    
    return tree;
}

} // extern "C"






//=============================================================================
// MCMC search


Tree *searchMCMC(Tree *initTree, 
                 string *genes, int nseqs, int seqlen, char **seqs,
                 SampleFunc *samples,
                 Prior *prior,
                 TopologyProposer *proposer,
                 BranchLengthFitter *fitter)
{
    Tree *toptree = NULL;
    double toplogl = -INFINITY, logl=-INFINITY, nextlogl;
    Tree *tree = initTree;
    
    
    // determine initial tree
    if (initTree == NULL)
        tree = getInitialTree(genes, nseqs, seqlen, seqs);
    
    // used for printSearchStatus
    ExtendArray<int> recon(tree->nnodes);
    ExtendArray<int> events(tree->nnodes);
    
    // init likelihood score
    parsimony(tree, nseqs, seqs); // get initial branch lengths
    logl = fitter->findLengths(tree);
    logl += prior->branchPrior(tree);
    toplogl = logl;
    toptree = tree->copy();
    
    int accept = 0;
    int reject = 0;
    
        
    // MCMC loop
    for (int i=0; proposer->more(); i++) {
        printLog(LOG_LOW, "search: iter %d\n", i);
    
        // propose new tree 
        proposer->propose(tree);
        //assert(tree->assertTree());
        
        double seqlk = fitter->findLengths(tree);
        double branchlk = prior->branchPrior(tree);
        nextlogl = seqlk + branchlk;       
        
        
        // acceptance rule
        if (nextlogl > logl ||
            nextlogl - logl > log(frand()))
        {
            // accept
            printLog(LOG_MEDIUM, "search: accept %f (%f)\n", nextlogl, logl);
            if (isLogLevel(LOG_MEDIUM))
                printSearchStatus(tree, prior->getSpeciesTree(), 
                                  prior->getGene2species(), recon, events);            
            
            accept++;
            logl = nextlogl;

            // keep track of toptree            
            if (logl > toplogl) {
                delete toptree;
                toptree = tree->copy();
                toplogl = logl;
            }
            
        } else {                
            // reject, undo topology change
            printLog(LOG_MEDIUM, "search: reject %f (%f)\n", nextlogl, logl);
            if (isLogLevel(LOG_MEDIUM))
                printSearchStatus(tree, prior->getSpeciesTree(), 
                                  prior->getGene2species(), recon, events);
            
            reject++;
            proposer->revert(tree);
        }
        
        // return a tree sample
        (*samples)(tree);
    }
    
    printLog(LOG_LOW, "accept rate: %f\n", accept / double(accept+reject));
    
    return toptree;
}



Tree *searchClimb(Tree *initTree, 
                  string *genes, int nseqs, int seqlen, char **seqs,
                  Prior *prior,
                  TopologyProposer *proposer,
                  BranchLengthFitter *fitter)
{
    return NULL;
}

} // namespace spidir





/*=============================================================================
OLD DUP LOSS PROPOSER

DupLossProposer::DupLossProposer(TopologyProposer *proposer, 
                                 SpeciesTree *stree, int *gene2species,
                                 float dupprob, float lossprob,
                                 int quickiter, int niter,
                                 int nsamples, bool extend) :
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
    oldtop(NULL),
    nsamples(nsamples),
    samplei(0),
    treesize(0)
{
    doomtable = new float [stree->nnodes];
    calcDoomTable(stree, dupprob, lossprob, maxdoom, doomtable);

}


DupLossProposer::~DupLossProposer()
{
    if (oldtop)
        delete oldtop;
    delete [] doomtable;

    clearQueue();
}


bool treePropCmp(const DupLossProposer::TreeProp &a, 
                 const DupLossProposer::TreeProp &b)
{
    return a.second > b.second;
}


void DupLossProposer::queueTrees(Tree *tree)
{
    
    // recon tree to species tree
    recon.ensureSize(tree->nnodes);
    events.ensureSize(tree->nnodes);
    recon.setSize(tree->nnodes);
    events.setSize(tree->nnodes);
    
    reconcile(tree, stree, gene2species, recon);
    labelEvents(tree, recon, events);

    // init queue for subproposals
    queue.clear();
    
    sum = -INFINITY;
    
    // make many subproposals
    proposer->reset();
    for (int i=0; i<quickiter; i++) {
        proposer->propose(tree);
        
        Node *oldroot1 = tree->root->children[0];
        Node *oldroot2 = tree->root->children[1];

        reconcile(tree, stree, gene2species, recon);
        labelEvents(tree, recon, events);
        float logl = birthDeathTreePriorFull(tree, stree, recon, events, 
                                             dupprob, lossprob,
                                             doomtable);
        sum = logadd(sum, logl);
        
        // save tree and logl
        Tree *tree2 = tree->copy();
        queue.push_back(TreeProp(tree2, logl));

        // restore tree
        tree->reroot(oldroot1, oldroot2);
        proposer->revert(tree);
    }
    
    // makes the random sample slightly faster
    std::sort(queue.begin(), queue.end(), treePropCmp);
    
    // start a new sample count
    samplei = 0;
}


void DupLossProposer::propose(Tree *tree)
{
    iter++;
    treesize = tree->nnodes;

    // do simple proposal if dup/loss probs are disabled
    if (dupprob < 0.0 || lossprob < 0.0 || quickiter <= 1) {
        proposer->propose(tree);
        return;
    }
    
    if (queue.size() == 0 || samplei >= nsamples)
        queueTrees(tree);
    samplei++;
    
    // save old topology
    if (oldtop)
        delete oldtop;
    oldtop = tree->copy();
    
    // randomly sample a tree from the queue
    float choice = log(frand()) + sum;
    float partsum = -INFINITY;
    
    for (unsigned int i=0; i<queue.size(); i++) {
        partsum = logadd(partsum, queue[i].second);
        
        if (choice < partsum) {
            // propose tree i
            tree->setTopology(queue[i].first);

            // remove tree i from queue
            sum -= queue[i].second;
            queue[i] = queue.back();
            queue.pop_back();            
            break;
        }
    }
}


void DupLossProposer::accept(bool accepted)
{
    if (accepted) {
        clearQueue();
        
        // extend more iterations
        if (extend)
            niter = max(niter, iter + treesize * 2);
    }
}


void DupLossProposer::revert(Tree *tree)
{
    // do simple proposal if dup/loss probs are disabled
    if (dupprob < 0.0 || lossprob < 0.0 || quickiter <= 1) {
        proposer->revert(tree);
        return;
    }
    
    tree->setTopology(oldtop);
    delete oldtop;
    oldtop = NULL;
}

void DupLossProposer::reset()
{
    iter = 0;
    clearQueue();
}


void DupLossProposer::clearQueue()
{
    // clear queue
    while (queue.size() > 0) {
        Tree *tree = queue.back().first;
        queue.pop_back();
        delete tree;
    }
}




=============================================================================
*/
