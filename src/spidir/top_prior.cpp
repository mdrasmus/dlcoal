/*=============================================================================

  Matt Rasmussen
  Copyright 2007-2011

  Gene tree topology prior of the SPIMAP model

=============================================================================*/


// c/c++ includes
#include <math.h>


// spidir includes
#include "birthdeath.h"
#include "common.h"
#include "phylogeny.h"
#include "Tree.h"



namespace spidir
{

extern "C" {


// TODO: make this into a class, so that allocations can happen once

// returns the number of labeled histories exist for the given tree topology
// uses the subtree starting at root and going until leaves.
// NOTE: assumes binary tree
double subtreeCorrection(Tree *tree, Node *root, ExtendArray<Node*> &leaves)
{
    double n = 1.0;
    
    // get nodes in post order
    ExtendArray<Node*> queue(0, tree->nnodes);

    // helper array for ensuring postorder traversal
    ExtendArray<int> visited(tree->nnodes);
    for (int i=0; i<visited.size(); i++)
        visited[i] = 0;

    // count number of descendant internal nodes
    ExtendArray<int> ninternals(tree->nnodes);

    // process leaves
    for (int i=0; i<leaves.size(); i++) {
        Node *node = leaves[i];
        ninternals[node->name] = 0;
        
        // queue parent
        visited[node->parent->name]++;
        queue.append(node->parent);
    }

    
    // go up tree until root
    for (int i=0; i<queue.size(); i++) {
        Node *node = queue[i];
        
        // do not process a node until both children are processed
        if (visited[node->name] != 2)
            continue;
        
        // count internal children
        const int right = ninternals[node->children[1]->name];
        const int left = ninternals[node->children[0]->name];
        const int x = ninternals[node->name] = 1 + right + left;
        n /= x;

        if (node == root)
            return n;

        visited[node->name]++;
        visited[node->parent->name]++;
        queue.append(node->parent);
    }

    // Note: this will occur for subtree like
    //    root
    //     |
    //     +
    //    / \  
    //   /   \  
    //  A     B
    // 
    return n;
}


// N_n(T) = prod_u (c_u!)^{-1}
double fullTreeCorrection(Tree *tree, Tree *stree, int *recon)
{
    double val = 0.0;
    ExtendArray<int> spcount(stree->nnodes);

    // count genes per species
    for (int i=0; i<stree->nnodes; i++)
        spcount[i] = 0;
    for (int i=0; i<tree->nnodes; i++) {
        if (tree->nodes[i]->isLeaf()) {
            double x = ++spcount[recon[tree->nodes[i]->name]];
            val -= logf(x);
        }
    }

    return val;
}


// computes the entries of the doom probabilty table
// uses analytic solution to infinite summation
void calcDoomTable(Tree *tree, float birthRate, float deathRate, 
                   double *doomtable)
{
    const double l = birthRate;
    const double u = deathRate;
    const double r = l - u;
    const double lu = l / u;
    double p0, p1;

    // get nodes in post order
    ExtendArray<Node*> nodes(0, tree->nnodes);
    getTreePostOrder(tree, &nodes);    

    
    for (int i=0; i<tree->nnodes; i++) {
        Node *node = nodes[i];
        
        if (node->isLeaf()) {
            doomtable[node->name] = -INFINITY;
        } else {          
            double prod = 0.0;            
            for (int j=0; j<node->nchildren; j++) {
                Node *child = node->children[j];

                // compute u_t and P(t(c))
                const double t = child->dist;
                const double dc = exp(doomtable[child->name]);

                if (birthRate == deathRate) {
                    const double lt = l * t;
                    const double lt1 = 1.0 + lt;
                    p0 = lt / lt1;
                    p1 = 1.0 / lt1 / lt1;
                } else {
                    const double ert = exp(-r * t);
                    const double luert = l - u*ert;
                    p0 = (u - u * ert) / luert;
                    p1 = r*r * ert / luert / luert;
                }

                prod += log(p0 + dc * p1 / (1.0 - lu * p0 * dc));
            }
	    
            doomtable[node->name] = prod;
        }
    }
}


// stores the leaves of the subtree below node that reconciles to snode
void getSpecSubtree(Node *node, Node *snode, int *recon, int *events,
                    ExtendArray<Node*> &nodes)
{
    for (int i=0; i<node->nchildren; i++) {
        Node *child = node->children[i];

        // only consider nodes that reconcile to snode
        if (recon[child->name] == snode->name) {
            if (events[child->name] == EVENT_SPEC ||
                events[child->name] == EVENT_GENE)
            {
                nodes.append(child);
            } else { // EVENT_DUP
                getSpecSubtree(child, snode, recon, events, nodes);
            }
        }
    }
}


// TODO: does not handle branches above the species tree root yet
// NOTE: assumes binary species tree
double birthDeathTreePrior(Tree *tree, Tree *stree, int *recon, 
                           int *events, float birthRate, float deathRate,
                           double *doomtable)
{
    const double l = birthRate;
    const double u = deathRate;
    const double r = l - u;
    const double lu = l / u;
    double p0, p1;
    double prob = 0.0;
    ExtendArray<Node*> subleaves(0, tree->nnodes);    

    // preroot duplications
    //if (events[tree->root->name] == EVENT_DUP)

    // loop through speciation nodes in tree
    for (int i=0; i<tree->nnodes; i++) {
        Node *node = tree->nodes[i];        
        if (events[node->name] == EVENT_SPEC) {

            // loop through nodes u \in child(R(v))
            Node *snode = stree->nodes[recon[node->name]];
            for (int j=0; j<snode->nchildren; j++) {
                Node *schild = snode->children[j];

                // get subtree that reconciles to snode
                subleaves.clear();
                getSpecSubtree(node, schild, recon, events, subleaves);

                // compute u_t and P(t(c))
                const double t = schild->dist;
                const double dc = exp(doomtable[schild->name]);
                const int s = subleaves.size();
                
                if (birthRate == deathRate) {
                    const double lt = l * t;
                    const double lt1 = 1.0 + lt;
                    p0 = lt / lt1;
                    p1 = 1.0 / lt1 / lt1;
                } else {
                    const double ert = exp(-r * t);
                    const double luert = l - u*ert;
                    p0 = (u - u * ert) / luert;
                    p1 = r*r * ert / luert / luert;
                }

                if (s > 1) {
                    prob += log(subtreeCorrection(tree, node, subleaves));
                    prob += log(pow(lu * p0, s-1) * p1 / 
                            (pow(1.0 - lu * p0 * dc, s+1)));
		} else if (s == 1) {
                    const double a = 1.0 - lu * p0 * dc;
                    prob += log(p1 / a / a);
                } else { // s == 0
                    prob += log(p0 + p1 * dc / (1.0 - lu * p0 *dc));
                }
                
                // sum over ndoom; non-analytic solution
                /*
                double sum = 0.0;
                for (int ndoom=0;  ndoom<=maxdoom; ndoom++) {
                    int totleaves = subleaves.size() + ndoom;
                    sum += fchoose(totleaves, ndoom) *
                        birthDeathCount(totleaves, 
                                        schild->dist,
                                        birthRate, deathRate) *
                        ipow(exp(doomtable[schild->name]), ndoom);
                }
                prob += log(sum);
                */
            }
        } else if (events[node->name] == EVENT_DUP) {
            // 2^{|dup(T, R, S)|}
            prob += log(2.0);
        }
    }

    
    prob += fullTreeCorrection(tree, stree, recon);
    return prob;
}


// Convenience function
// Adds and removes implied species nodes to the gene tree
// NOTE: assumes binary species tree
double birthDeathTreePriorFull(Tree *tree, Tree *stree, int *recon, 
                               int *events, float birthRate, float deathRate,
                               double *doomtable)
{
    ExtendArray<int> recon2(0, 2 * tree->nnodes);
    recon2.extend(recon, tree->nnodes);

    ExtendArray<int> events2(0, 2 * tree->nnodes);
    events2.extend(events, tree->nnodes);


    int addedNodes = addImpliedSpecNodes(tree, stree, recon2, events2);
    double p = birthDeathTreePrior(tree, stree, recon2, events2, 
                                   birthRate, deathRate,
                                   doomtable);
    removeImpliedSpecNodes(tree, addedNodes);

    return p;
}


}

} // namespace spidir

