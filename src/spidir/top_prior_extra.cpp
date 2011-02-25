
// c/c++ includes
#include <math.h>


// spidir includes
#include "birthdeath.h"
#include "common.h"
#include "phylogeny.h"
#include "top_prior.h"
#include "Tree.h"



namespace spidir
{

extern "C" {


// returns the number of labeled histories with 'ngenes' surviving lineages
int inumHistories(int ngenes)
{
    // gaurd against overflow
    assert(ngenes <= 9);

    int n = 1;
    for (int i=2; i<=ngenes; i++) {
        n *= i*(i-1) / 2;
    }
    return n;
}


// returns the number of labeled histories with 'ngenes' surviving lineages
double numHistories(int ngenes)
{
    double n = 1;
    for (int i=2; i<=ngenes; i++) {
        n *= i*(i-1) / 2;
    }
    return n;
}

// returns the number of labeled histories exist for the given tree topology
// NOTE: assumes binary tree
int inumTopologyHistories(Tree *tree)
{
    int n = 1;

    // get nodes in post order
    ExtendArray<Node*> nodes(0, tree->nnodes);
    getTreePostOrder(tree, &nodes);

    // count number of descendant internal nodes
    ExtendArray<int> ninternals(tree->nnodes);

    for (int i=0; i<tree->nnodes; i++) {
        Node *node = nodes[i];

        if (node->isLeaf()) {
            ninternals[node->name] = 0;
        } else {
            // count internal children
            const int right = ninternals[node->children[0]->name];
            const int left = ninternals[node->children[1]->name];
            ninternals[node->name] = 1 + right + left;
            n *= choose(right + left, right);
        }
    }

    return n;
}




// returns the number of labeled histories exist for the given tree topology
// NOTE: assumes binary tree
double numTopologyHistories(Tree *tree)
{
    double n = 1;

    // get nodes in post order
    ExtendArray<Node*> nodes(0, tree->nnodes);
    getTreePostOrder(tree, &nodes);

    // count number of descendant internal nodes
    ExtendArray<int> ninternals(tree->nnodes);

    for (int i=0; i<tree->nnodes; i++) {
        Node *node = nodes[i];

        if (node->isLeaf()) {
            ninternals[node->name] = 0;
        } else {
            // count internal children
            const int right = ninternals[node->children[0]->name];
            const int left = ninternals[node->children[1]->name];
            ninternals[node->name] = 1 + right + left;
            n *= fchoose(right + left, right);
        }
    }

    return n;
}


// computes the entries of the doom probabilty table more slowly
void calcDoomTableSlow(Tree *tree, float birthRate, float deathRate, 
                       int maxdoom,
                       double *doomtable)
{
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
                double sum = 0.0;

                for (int ndoom=0; ndoom<=maxdoom; ndoom++) {
                    sum += (birthDeathCount(ndoom, child->dist, 
                                            birthRate, deathRate) *
                            ipow(exp(doomtable[child->name]), ndoom));
                }

                prod += log(sum);
            }
	    
            doomtable[node->name] = prod;
        }
    }
}



// returns the number of labeled histories exist for the given tree topology
// uses the subtree starting at root and going until leaves.
// NOTE: assumes binary tree
double numSubtopologyHistories(Tree *tree, Node *root, ExtendArray<Node*> &leaves)
{
    double n = 1;
    
    // get nodes in post order
    ExtendArray<Node*> queue(0, tree->nnodes);

    // count number of descendant internal nodes
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
        ninternals[node->name] = 1 + right + left;
        n *= fchoose(right + left, right);

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


int intCmp(const void *_a, const void *_b)
{
    return *((int*) _a) - *((int*) _b);
}


// N_n(T, allLeaves=false) = 2^{-M(T)} prod_u c_u!
// N_n(T, allLeaves=true) = 2^{-M(T)} L(T)!
double numRedundantTopologies(Tree *tree, Node *root, 
                              ExtendArray<Node*> &leaves, 
                              int *hashids, bool allLeaves)
{

    // get nodes in post order
    ExtendArray<Node*> queue(0, 2 * tree->nnodes);

    // initialize visited array to zeros
    ExtendArray<int> visited(tree->nnodes);
    for (int i=0; i<visited.size(); i++)
        visited[i] = 0;
    
    // process leaves
    for (int i=0; i<leaves.size(); i++) {
        Node *node = leaves[i];        
        // queue parent
        visited[node->parent->name]++;
        queue.append(node->parent);
    }

    int nmirrors = 0;
    
    // go up tree until root
    for (int i=0; i<queue.size(); i++) {
        Node *node = queue[i];
        
        // do not process a node until all children are processed
        if (visited[node->name] != node->nchildren)
            continue;
        
        // count internal children
        if (node->nchildren == 2) {
            nmirrors += int(hashids[node->children[0]->name] == 
                            hashids[node->children[1]->name]);
        } else {
            // we do not handle multifurcating nodes
            assert(node->nchildren < 2);
        }

        if (node == root)
            break;

        visited[node->name]++;
        visited[node->parent->name]++;
        queue.append(node->parent);
    }

    //printf("queue.size = %d\n", queue.size());

    double val = 0.0;

    if (allLeaves) {
        // val = log(factorial(leaves.size()))
        for (double i=2.0; i<=leaves.size(); i+=1.0)
            val += logf(i);
    } else {
        // get hashes
        ExtendArray<int> leafhashes(0, leaves.size());
        for (int i=0; i<leaves.size(); i++)
            leafhashes.append(hashids[leaves[i]->name]);
    
        qsort((void*) leafhashes.get(), leafhashes.size(), sizeof(int), intCmp);

        double colorsize = 1;
        for (int i=1; i<leaves.size(); i++) {
            if (leafhashes[i] != leafhashes[i-1]) {
                // val *= factorial(colorsize)
                for (double j=2; j<=colorsize; j+=1.0)
                    val += logf(j);
                colorsize = 1.0;
            } else {
                colorsize += 1.0;
            }
        }
        for (double j=2; j<=colorsize; j+=1.0)
            val += logf(j);
    }

    // divide by 2^M
    val -= nmirrors * logf(2.0);

    return val;
}



class KeyList {
public:
    KeyList(int key) : 
        key(key),
        next(NULL)
    {}

    int key;
    KeyList *next;
};

class HashNode {
public:
    HashNode(int nodeid, KeyList *start, KeyList *end, int len) :
        nodeid(nodeid),
	start(start),
	end(end),
	len(len)
    {}

    int nodeid;
    KeyList *start;
    KeyList *end;
    int len;
};


int hashNodeCmp(const void *_a, const void *_b)
{
    HashNode *a = *((HashNode**) _a);
    HashNode *b = *((HashNode**) _b);
    
    // first compare diffs
    int diff = a->len - b->len;
    if (diff) {
    	return diff;
    } else {
        // compare keys
	KeyList *keya = a->start;
	KeyList *keyb = b->start;

	while (true) {
	    diff = keya->key - keyb->key;
	    if (diff)
		return diff;
            if (keya == a->end || keyb == b->end)
                break;
            keya = keya->next;
            keyb = keyb->next;
	}

	// both hashes are the same
	// 1. same length
	// 2. same key subsequence
	return 0;
    }
}



void getHashIds(Tree *tree, int *recon, int *hashids)
{

    ExtendArray<HashNode*> hashnodes(tree->nnodes);
    ExtendArray<KeyList*> keylist(tree->nnodes);

    // get post order of nodes
    ExtendArray<Node*> postnodes(0, tree->nnodes);
    getTreePostOrder(tree, &postnodes);    

    // build hash nodes
    for (int i=0; i<postnodes.size(); i++)
    {
        Node *node=postnodes[i];
        
        if (node->isLeaf()) {
            KeyList *key = new KeyList(recon[node->name]);
            keylist[node->name] = key;
            hashnodes[node->name] = new HashNode(node->name, key, key, 1);
        } else {
            if (node->nchildren == 1) {                
                KeyList *key = new KeyList(-1);
                keylist[node->name] = key;
                HashNode *hnode1 = hashnodes[node->children[0]->name];

                // join lists: list1 -> key
                hashnodes[node->name] = 
                    new HashNode(node->name, hnode1->start, key, 
                                 hnode1->len + 1);
                hnode1->end->next = key;
                
            } else if (node->nchildren == 2) {
                KeyList *key = new KeyList(-2);
                keylist[node->name] = key;
                HashNode *hnode1 = hashnodes[node->children[0]->name];
                HashNode *hnode2 = hashnodes[node->children[1]->name];
                int len = hnode1->len + hnode2->len + 1;
                int cmp = hashNodeCmp(&hnode1, &hnode2);

                if (cmp <= 0) {
                    // join lists: list1 -> list2 -> key
                    hashnodes[node->name] = new HashNode(node->name,
                                                         hnode1->start, 
                                                         key,
                                                         len);
                    hnode1->end->next = hnode2->start;
                    hnode2->end->next = key;
                } else {
                    // join lists: list2 -> list1 -> key
                    hashnodes[node->name] = new HashNode(node->name,
                                                         hnode2->start, 
                                                         key,
                                                         len);
                    hnode2->end->next = hnode1->start;
                    hnode1->end->next = key;
                }
            } else {
                // cannot handle multifurcating nodes
                assert(0);
            }
        }
    }

    // sort hashnodes
    qsort((void*) hashnodes.get(), hashnodes.size(), sizeof(HashNode*),
	  hashNodeCmp);
    
    int hashid = 0;
    hashids[hashnodes[0]->nodeid] = hashid;
    for (int i=1; i<hashnodes.size(); i++) {
	// use new hashid if nodes differ
	if (hashNodeCmp(&hashnodes[i], &hashnodes[i-1]))
	    hashid++;
	hashids[hashnodes[i]->nodeid] = hashid;
    }
    

    // clean up
    for (int i=0; i<tree->nnodes; i++) {
	delete hashnodes[i];
        delete keylist[i];
    }

    //printf("hashids = ");
    //printIntArray(hashids, tree->nnodes);
    //printf("\n");
}



// TODO: does not handle branches above the species tree root yet
// NOTE: assumes binary species tree
double birthDeathTreePrior2(Tree *tree, Tree *stree, int *recon, 
                            int *events, float birthRate, float deathRate,
                            double *doomtable, int maxdoom)
{

    double prob = 0.0;
    ExtendArray<Node*> subleaves(0, tree->nnodes);
    
    // catch undefined params
    if (birthRate == deathRate)
        deathRate = .99 * birthRate;
    
    ExtendArray<int> hashids(tree->nnodes);
    getHashIds(tree, recon, hashids);

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

		double nhist, thist;

		if (subleaves.size() == 0) {
		    nhist = 1.0;
		    thist = 1.0;
		} else {
		    nhist = numSubtopologyHistories(tree, node, subleaves);
		    thist = numHistories(subleaves.size());		    

		    if (subleaves[0]->isLeaf()) {
                        // correct subtrees that have leaves
			nhist *= exp(numRedundantTopologies(tree, node, 
                                                            subleaves, 
                                                            hashids,
                                                            false));
		    } else {
                        // correct subtrees that have leaves
			double a = exp(numRedundantTopologies(tree, node, 
                                                              subleaves, 
                                                              hashids,
                                                              true));
                        nhist *= a;
                    }
		}
		
                // sum over ndoom
                double sum = 0.0;
                for (int ndoom=0;  ndoom<=maxdoom; ndoom++) {
                    int totleaves = subleaves.size() + ndoom;
                    sum += fchoose(totleaves, ndoom) *
                        birthDeathCount(totleaves, 
                                        schild->dist,
                                        birthRate, deathRate) *
                        ipow(exp(doomtable[schild->name]), ndoom);
                }

                prob += log(nhist) - log(thist) + log(sum);
            }
        }
    }

    ExtendArray<Node*> leaves(0, tree->nnodes);
    for (int i=0; i<tree->nnodes; i++) {
	if (tree->nodes[i]->isLeaf())
	    leaves.append(tree->nodes[i]);
    }
    double x = numRedundantTopologies(tree, tree->root, leaves, 
                                      hashids, false);
    prob -= x;
    
    return prob;
}



} // extern C

} // namespace spidir
