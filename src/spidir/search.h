/*=============================================================================

  Matt Rasmussen
  Copyright 2007-2011

  Gene tree search functions

=============================================================================*/


#ifndef SPIDIR_SEARCH_H
#define SPIDIR_SEARCH_H

#include <set>

#include "model_params.h"


namespace spidir {

using namespace std;


struct lttree
{
    bool operator()(const int* k1, const int* k2) const
    {
        for (int i=0; k1[i] != -2; i++) {
            if (k1[i] < k2[i])
                return true;
            if (k1[i] > k2[i])
                return false;
        }
        return false;
    }
};

class TreeSet
{
public:
    TreeSet() : key(100) {}
    ~TreeSet();

    void clear();
    bool insert(Tree *tree);
    bool has(Tree *tree);
    int size() { return trees.size(); }

    ExtendArray<int> key;
    typedef set<int*, lttree> Set;
    Set trees;
};


class TopologyProposer
{
public:
    TopologyProposer() :
        correctTree(NULL),
        correctSeen(false)
    {}

    virtual ~TopologyProposer() {}
    virtual void propose(Tree *tree) {}
    virtual void revert(Tree *tree) {}
    virtual bool more() { return false; }
    virtual void reset() {}
    virtual void accept(bool accepted) {}

    virtual void setCorrect(Tree *tree) { correctTree = tree; }
    virtual Tree *getCorrect() { return correctTree; }
    virtual bool seenCorrect() { return correctSeen; }
    virtual void testCorrect(Tree *tree)
    {
        // debug: keep track of correct tree in search
        if (correctTree) {
            if (tree->sameTopology(correctTree))
                correctSeen = true;
        }
    }

protected:
    Tree *correctTree;
    bool correctSeen;
};



class NniProposer: public TopologyProposer
{
public:
    NniProposer(int niter=500);

    virtual void propose(Tree *tree);
    virtual void revert(Tree *tree);
    virtual bool more() { return iter < niter; }
    virtual void reset() { iter = 0; }


protected:    
    int niter;
    int iter;

    Node *nodea;
    Node *nodeb;
};


class SprProposer: public NniProposer
{
public:
    SprProposer(int niter=500);

    virtual void propose(Tree *tree);
    virtual void revert(Tree *tree);

protected:    
    int niter;
    int iter;

    Node *nodea;
    Node *nodeb;
    Node *nodec;
};


class SprNbrProposer: public NniProposer
{
public:
    SprNbrProposer(int niter=500, int radius=4);

    virtual void propose(Tree *tree);
    virtual void revert(Tree *tree);
    void reset() { 
        iter = 0; 
        reverted = false;
        basetree = NULL;
    }

    void pickNewSubtree();

protected:
    int radius;
    Tree *basetree;
    Node *subtree;
    list<Node*> queue;
    vector<int> pathdists;
    bool reverted;
};


class MixProposer: public TopologyProposer
{
public:
    MixProposer(int niter=500) : totalWeight(0), niter(niter), iter(0) {}

    virtual void propose(Tree *tree);
    virtual void revert(Tree *tree);    
    virtual bool more() { return iter < niter; }
    virtual void reset() { 
        iter = 0; 

        // propagate reset
        for (unsigned int i=0; i<methods.size(); i++)
            methods[i].first->reset();
    }

    void addProposer(TopologyProposer *proposer, float weight);

protected:

    float totalWeight;
    typedef pair<TopologyProposer*,float> Method;
    vector<Method> methods;
    int lastPropose;
    int niter;
    int iter;
};


class ReconRootProposer: public TopologyProposer
{
public:
    ReconRootProposer(TopologyProposer *proposer,
                      SpeciesTree *stree=NULL, int *gene2species=NULL) :
        proposer(proposer),
        stree(stree),
        gene2species(gene2species),
        oldroot1(NULL),
        oldroot2(NULL)
    {}

    virtual void propose(Tree *tree);
    virtual void revert(Tree *tree);
    virtual bool more() { return proposer->more(); }
    virtual void reset() { return proposer->reset(); }
    virtual void accept(bool accepted) { proposer->accept(accepted); }

protected:
    TopologyProposer *proposer;
    SpeciesTree *stree;
    int *gene2species;
    Node *oldroot1;
    Node *oldroot2;
};


class DupLossProposer: public TopologyProposer
{
public:
    DupLossProposer(TopologyProposer *proposer, 
                    SpeciesTree *stree,
                    int *gene2species,
                    float dupprob,
                    float lossprob,
                    int quickiter=100, int niter=500);

    virtual ~DupLossProposer();

    virtual void propose(Tree *tree);
    virtual void revert(Tree *tree);
    virtual bool more() { return iter < niter; }
    virtual void reset() {
        iter = 0;
        uniques.clear();
        proposer->reset();
    }
    
protected:
    TopologyProposer *proposer;
    int quickiter;
    int niter;
    int iter;
    Tree *correctTree;
    bool correctSeen;
    SpeciesTree *stree;
    int *gene2species;
    float dupprob;
    float lossprob;
    double *doomtable;
    TreeSet uniques;

    ExtendArray<int> recon;
    ExtendArray<int> events;
    Tree *oldtop;
};


class UniqueProposer: public TopologyProposer
{
public:
    UniqueProposer(TopologyProposer *proposer, int niter=-1, int ntries=10) : 
        proposer(proposer), niter(niter), iter(0), ntries(ntries) {}
    virtual ~UniqueProposer();

    virtual void propose(Tree *tree);
    virtual void revert(Tree *tree) { return proposer->revert(tree); }
    virtual bool more() { 
        if (niter == -1)
            return proposer->more(); 
        else 
            return iter < niter;
    }
    virtual void reset() { 
        seenTrees.clear();
        proposer->reset(); 
    }
    
    TopologyProposer *proposer;
    TreeSet seenTrees;
    int niter;
    int iter;
    int ntries;
};


class DefaultSearch
{
public:
    DefaultSearch(int niter, int quickiter,
                  SpeciesTree *stree, int *gene2species,
                  float duprate, float lossrate, float sprrate=.5,
                  int radius=3) :
        stree(stree),
        gene2species(gene2species),
        radius(radius),

        nni(niter),
        spr(niter),
        sprnbr(niter, radius),

        mix(niter),

        rooted(&mix, stree, gene2species),
        unique(&rooted, niter),
        dl(&unique, stree, gene2species, 
           duprate, lossrate, quickiter, niter),

        mix2(niter)
        
    {
        mix.addProposer(&nni, (1 - sprrate) / 2.0);
        mix.addProposer(&sprnbr, sprrate);
        mix.addProposer(&spr, (1 - sprrate) / 2.0);

        mix2.addProposer(&unique, .2);
        mix2.addProposer(&dl, .8);
    }
        
    int niter;
    int quickiter;
    SpeciesTree *stree;
    int *gene2species;
    float duprate;
    float lossrate;
    int radius;
    
    NniProposer nni;
    SprProposer spr;
    SprNbrProposer sprnbr;
    
    MixProposer mix;

    ReconRootProposer rooted;
    UniqueProposer unique;
    DupLossProposer dl;

    MixProposer mix2;
};



//=============================================================================



class SampleFunc
{
public:
    SampleFunc(FILE *output) :
        output(output)
    {
    }
    
    virtual ~SampleFunc()
    {
        fclose(output);
    }

    void operator()(Tree *tree)
    {
        writeNewickTree(output, tree, 0, true);
        fprintf(output, "\n");
    }
    
protected:
    FILE *output;
};


class TreeSearch
{
public:

    TreeSearch() :
        proposal_runtime(0.0)
    {}

    virtual ~TreeSearch()
    {}

    virtual Tree *search(Tree *initTree, 
			 string *genes, 
			 int nseqs, int seqlen, char **seqs)
    { return NULL; }

    double proposal_runtime;
};


class TreeSearchClimb : public TreeSearch
{
public:

    TreeSearchClimb(Model *model, TopologyProposer *proposer);
    virtual ~TreeSearchClimb();

    virtual Tree *search(Tree *initTree, 
			 string *genes, 
			 int nseqs, int seqlen, char **seqs);

protected:
    Model *model;
    TopologyProposer *proposer;

};




Tree *getInitialTree(string *genes, int nseqs, int seqlen, char **seqs,
                     SpeciesTree *stree, int *gene2species);
Tree *getInitialTree(string *genes, int nseqs, int seqlen, char **seqs);



} // namespace spidir


#endif // SPIDIR_SEARCH_H




