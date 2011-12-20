/*=============================================================================

  Matt Rasmussen
  Copyright 2007-2011

  SPIMAP model parameters

=============================================================================*/



#include "common.h"
#include "model_params.h"
#include "logging.h"
#include "phylogeny.h"

namespace spidir {


//=============================================================================
// Spidir Parameters

SpidirParams *readSpidirParams(const char* filename)
{
    FILE *infile = NULL;
    
    if ((infile = fopen(filename, "r")) == NULL) {
        fprintf(stderr, "cannot read file '%s'\n", filename);
        return NULL;
    }
    
    const int MAX_NAME = 51;
    float param1, param2;
    float alpha = -1, beta = -1;
    ExtendArray<float> mu(0, 40);
    ExtendArray<float> sigma(0, 40);
    ExtendArray<string> names(0, 40);
    
    char name[MAX_NAME];
    
    while (!feof(infile)) {
        int ntokens = fscanf(infile, "%50s\t%f\t%f", name, &param1, &param2);
        if (ntokens <= 0)
            break;
        if (ntokens != 3) {
            return NULL;
        }
        
        if (!strcmp(name, "baserate")) {
            alpha = param1;
            beta = param2;
        } else {
            names.append(name);
            mu.append(param1);
            sigma.append(param2);
        }
    }
    fclose(infile);    
    
    return new SpidirParams(names.size(), names, mu, sigma, alpha, beta);
}



bool SpidirParams::order(SpeciesTree *stree)
{
    if (stree->nnodes != nsnodes) {
        printError("wrong number of parameters: %d %d\n", stree->nnodes, nsnodes);
        return false;
    }
    
    ExtendArray<Node*> nodeorder(0, stree->nnodes);
    getTreePreOrder(stree, &nodeorder);
        

    // make interior node names
    ExtendArray<int> inodes(stree->nnodes);
    
    int inodename = 1;
    for (int i=0; i<stree->nnodes; i++) {
        Node *node = nodeorder[i];
        if (node->isLeaf()) {
            inodes[node->name] = -1;
        } else {
            inodes[node->name] = inodename++;
        }
    }
    
    
    // loop through preordered nodes to construct permutation
    ExtendArray<int> invperm(0, stree->nnodes);
        
    for (int j=0; j<nsnodes; j++) {
        if (invperm.size() != j) {
            printError("unable to match '%s' to the species tree", 
                       names[j].c_str());
            return false;
        }

        // try to parse node id as an int
        int id;
        bool isint = (sscanf(names[j].c_str(), "%d", &id) == 1);
        
        for (int i=0; i<stree->nnodes; i++) {
            if (stree->nodes[i]->isLeaf()) {
                // if leaf, check if names match
                if (names[j] == stree->nodes[i]->longname) {
                    invperm.append(i);
                    break;
                }
            } else {
                if (isint && id == inodes[i]) {
                    invperm.append(i);
                    break;
                }
            }
        }
    }

    
    // apply permutation
    ExtendArray<int> perm(stree->nnodes);
    invertPerm(invperm, perm, nsnodes);    
      
    permute(names, perm, nsnodes);
    permute(sp_alpha, perm, nsnodes);
    permute(sp_beta, perm, nsnodes);
    
    return true;
}


} // namespace spidir
