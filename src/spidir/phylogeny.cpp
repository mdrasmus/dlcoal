/*=============================================================================

  Matt Rasmussen
  Copyright 2007-2011

  Gene trees, species trees, reconciliations

=============================================================================*/


#include "logging.h"
#include "phylogeny.h"
#include "parsing.h"

namespace spidir {


//=============================================================================
// reconciliation functions


typedef pair<Node*, Node*> Edge;
void getReconRootOrder(Node *node, ExtendArray<Edge> *edges)
{
    edges->append(Edge(node, node->parent));
    
    if (!node->isLeaf()) {
        for (int i=0; i<node->nchildren; i++)
            getReconRootOrder(node->children[i], edges);
        edges->append(Edge(node, node->parent));
    }
}



// NOTE: assumes binary tree
void reconRoot(Tree *tree, SpeciesTree *stree, int *gene2species)
{

    // special cases
    if (tree->nnodes < 3)
        return;
    
    // determine rooting order
    ExtendArray<Edge> edges(0, tree->nnodes);
    edges.append(Edge(tree->root->children[0],
                      tree->root->children[1]));
        
    for (int i=0; i<tree->root->nchildren; i++) {
        Node *node = tree->root->children[i];
        for (int j=0; j<node->nchildren; j++) {
            getReconRootOrder(node->children[j], &edges);
        }
        edges.append(Edge(tree->root->children[0],
                          tree->root->children[1]));
    }
    
    // try initial root and recon
    tree->reroot(edges[0].first, edges[0].second);
    ExtendArray<int> recon(tree->nnodes);
    ExtendArray<int> events(tree->nnodes);
    reconcile(tree, stree, gene2species, recon);
    labelEvents(tree, recon, events);

    int minroot = 0;
    int mincost = countDuplications(events.size(), events);
    int cost = mincost;
    
    // try other roots
    for (int i=0; i<edges.size(); i++) {
        // get new edge
        Edge edge = edges[i];
        if (edge.first->parent != edge.second)
            swap(edge.first, edge.second);
    
        // uncount cost
        if (events[tree->root->name] == EVENT_DUP)
            cost--;
        if (events[edge.second->name] == EVENT_DUP)
            cost--;
        
        // reroot
        tree->reroot(edge.first, edge.second);
        
        // Recompute recon and events
        recon[edge.second->name] = reconcileNode(edge.second, stree, recon);
        recon[tree->root->name] = reconcileNode(tree->root, stree, recon);
        events[edge.second->name] = labelEventsNode(edge.second, recon);
        events[tree->root->name] = labelEventsNode(tree->root, recon);

        // count any new duplications
        if (events[tree->root->name] ==  EVENT_DUP)
            cost ++;        
        if (events[edge.second->name] ==  EVENT_DUP)
            cost++;
        
        // record mincost root
        if (cost < mincost) {
            mincost = cost;
            minroot = i;
        }
    }
    
    // root tree by minroot
    tree->reroot(edges[minroot].first, edges[minroot].second);
}



// Find Last Common Ancestor
Node *treeLca(SpeciesTree *stree, Node *node1, Node *node2)
{
    int index1 = stree->preorder[node1->name];
    int index2 = stree->preorder[node2->name];
    
    while (index1 != index2) {
        if (index1 > index2) {
            node1 = node1->parent;
            index1 = stree->preorder[node1->name];
        } else {
            node2 = node2->parent;
            index2 = stree->preorder[node2->name];
        }
    }
    
    return node1;
}


// NOTE: assumes binary species tree
void reconcile_recurse(Tree *tree, Node *node, SpeciesTree *stree, int *recon)
{
    // recurse
    for (int i=0; i<node->nchildren; i++)
        reconcile_recurse(tree, node->children[i], stree, recon);
    
    // post process
    if (node->nchildren > 0) {
        int sname1 = recon[node->children[0]->name];
        int sname2 = recon[node->children[1]->name];
    
        // this node's species is lca of children species
        recon[node->name] = treeLca(stree, 
                                    stree->nodes[sname1], 
                                    stree->nodes[sname2])->name;
    }
}


// TODO: implement more efficiently with post order traversal
// reconcile a gene tree with a species tree
void reconcile(Tree *tree, SpeciesTree *stree,
               int *gene2species, int *recon)
{  
    // label gene leaves with their species
    for (int i=0; i<tree->nnodes; i++)
        if (tree->nodes[i]->isLeaf())
            recon[i] = gene2species[i];
    
    reconcile_recurse(tree, tree->root, stree, recon);    
}


// test whether a reconciliation is valid
void reconAssert(Tree *tree, SpeciesTree *stree, int *recon)
{
    for (int i=0; i<tree->nnodes; i++) {
        Node *node = tree->nodes[i];
        Node *snode = stree->nodes[recon[node->name]];

        if (node->parent == NULL)
            continue;
        
        Node *snode2 = stree->nodes[recon[node->parent->name]];
        
        // ensure every parent's species is above its child's species
        for(Node *ptr = snode; ptr != snode2; ptr = ptr->parent) {
            assert(ptr != NULL);
        }
    }
}



// label events for each node in tree
// NOTE: assumes binary gene tree
void labelEvents(Tree *tree, int *recon, int *events)
{
    Node **nodes = tree->nodes;

    for (int i=0; i<tree->nnodes; i++) {
        if (nodes[i]->nchildren == 0)
            events[i] = EVENT_GENE;
        else 
        if (recon[i] == recon[nodes[i]->children[0]->name] ||
            recon[i] == recon[nodes[i]->children[1]->name])
            events[i] = EVENT_DUP;
        else
            events[i] = EVENT_SPEC;
    }
}


int countLoss_recurse(Node *node, SpeciesTree *stree, int *recon)
{
    int loss = countLossNode(node, stree, recon);

    // recurse
    for (int i=0; i<node->nchildren; i++)
        loss += countLoss_recurse(node->children[i], stree, recon);

    return loss;
}


int countLoss(Tree *tree, SpeciesTree *stree, int *recon)
{
    return countLoss_recurse(tree->root, stree, recon);
}


// assumes binary tree
int countLossNode(Node *node, SpeciesTree *stree, int *recon)
{
    int loss = 0;
    
    // if not parent, then no losses
    if (node->parent == NULL)
        return 0;
    
    // determine starting and ending species
    Node *sstart = stree->nodes[recon[node->name]];
    Node *send = stree->nodes[recon[node->parent->name]];
    
    // the species path is too short to have losses
    if (sstart == send)
        return 0;
    
    // determine species path of this gene branch (node, node->parent)
    Node *ptr = sstart->parent;
    while (ptr != send) {
        // process ptr
        loss += ptr->nchildren - 1;
        
        // go up species tree
        ptr = ptr->parent;
    }
    
    // determine whether node->parent is a dup
    // if so, send (a.k.a. species end) is part of species path
    if (send->name == recon[node->parent->children[0]->name] ||
        send->name == recon[node->parent->children[1]->name])
         loss += send->nchildren - 1;
        
    return loss;
}


// insert new speciation node above node
void addSpecNode(Node *node, Node *snode, Tree *tree, 
                 ExtendArray<int> &recon, ExtendArray<int> &events)
{
    //printf("node: %d %d\n", node->name, snode->name);
    
    Node *newnode = tree->addNode(new Node(1));
    Node *parent = node->parent;
    
    // find index of node in parent's children
    int nodei=0;
    for (; nodei<parent->nchildren; nodei++)
        if (parent->children[nodei] == node)
            break;
    assert(nodei != parent->nchildren);
    
    // insert new node into tree
    parent->children[nodei] = newnode;
    newnode->parent = parent;
    newnode->children[0] = node;
    node->parent = newnode;
    
    // add recon and events info
    recon.append(snode->name);
    events.append(EVENT_SPEC);
}


int addImpliedSpecNodes(Tree *tree, Tree *stree, 
    ExtendArray<int> &recon, ExtendArray<int> &events)
{
    int addedNodes = 0;
    
    // recurse
    int nnodes = tree->nnodes;
    for (int i=0; i<nnodes; i++) {
        Node *node = tree->nodes[i];
        // process this node and the branch above it
        
        // if no parent, then no implied speciation nodes above us
        //if (node->parent == NULL) 
        //   continue;

        // handle root node specially
        if (node->parent == NULL) {
            // ensure root of gene tree properly reconciles to
            // root of species tree
            if (recon[node->name] == stree->root->name)
                continue;
            assert(tree->root == node);
            // NOTE: root is not last node (may need to relax this condition)
            tree->root = tree->addNode(new Node(1));
            tree->root->children[0] = node;
            node->parent = tree->root;
            recon.append(stree->root->name);
            events.append(EVENT_SPEC);
            addedNodes++;
        }
        
        // determine starting and ending species
        Node *sstart = stree->nodes[recon[node->name]];
        Node *send = stree->nodes[recon[node->parent->name]];
    
        // the species path is too short to have implied speciations
        if (sstart == send)
            continue;
    
        Node *parent = node->parent;
    
        // determine species path of this gene branch (node, node->parent)    
        for (Node *ptr = sstart->parent; ptr != send; ptr = ptr->parent) {
            // process ptr
            addSpecNode(node, ptr, tree, recon, events);
            addedNodes++;
            node = node->parent;
        }
    
        // determine whether node->parent is a dup
        // if so, send (a.k.a. species end) is part of species path
        if (events[parent->name] == EVENT_DUP) {
            addSpecNode(node, send, tree, recon, events);
            addedNodes++;
        }
    }

    return addedNodes;
}


// NOTE: this function adds back distance
void removeImpliedSpecNodes(Tree *tree, int addedNodes)
{
    int nnodes = tree->nnodes;

    // remove nodes from end of node list
    for (int i=nnodes-1; i>=nnodes-addedNodes; i--)
    {
        Node *oldnode = tree->nodes[i];
        Node *parent = oldnode->parent;
        Node *child = oldnode->children[0];

        assert(oldnode->nchildren == 1);

        // was node added as new root?
        if (parent == NULL) {
            tree->root = oldnode->children[0];
            tree->root->parent = NULL;
            tree->root->dist += oldnode->dist;
            tree->nodes.pop();
            tree->nnodes--;
            delete oldnode;
            continue;
        }
        
        // find index of oldnode in parent's children
        int nodei = 0;
        for (; nodei<parent->nchildren; nodei++)
            if (parent->children[nodei] == oldnode)
                break;
        assert(nodei != parent->nchildren);
    
        // remove old node from tree
        parent->children[nodei] = child;
        child->parent = parent;
        child->dist += oldnode->dist;
        tree->nodes.pop();
        tree->nnodes--;
        delete oldnode;
    }
}


void writeRecon(FILE *out, Tree *tree, SpeciesTree *stree,
                int *recon, int *events)
{
    const char* eventstr[] = { "gene", "spec", "dup" };

    for (int i=0; i<tree->nnodes; i++) {
        Node *node = tree->nodes[i];
        fprintf(out, "%s\t%s\t%s\n", node->longname.c_str(), 
                stree->nodes[recon[i]]->longname.c_str(), 
                eventstr[events[i]]);
    }
}


bool writeRecon(const char *filename, Tree *tree, SpeciesTree *stree,
                int *recon, int *events)
{
    FILE *out = NULL;
    
    if ((out = fopen(filename, "w")) == NULL) {
        printError("cannot write file '%s'\n", filename);
        return false;
    }

    writeRecon(out, tree, stree, recon, events);
    fclose(out);
    return true;
}



// sets the longnames of the internal nodes of a tree that do not already 
// have longnames.
// uses preorder traversal, by default root is "n1"
// returns last name used in subtree
int setInternalNames(Tree *tree, Node *node, int name)
{
    const int maxsize = 20;
    char numstr[maxsize + 1];

    if (node == NULL) {
        node = tree->root;
    }


    if (node->longname == "") {
        snprintf(numstr, maxsize, "n%d", name);
        node->longname = string(numstr);
    }

    for (int i=0; i<node->nchildren; i++) {
        if (!node->children[i]->isLeaf()) {
            name = setInternalNames(tree, node->children[i], name+1);
        }
    }

    return name;
}


//=============================================================================
// Gene2species

const string Gene2species::NULL_SPECIES;

bool Gene2species::read(const char *filename)
{
    BufferedReader reader;
    if (!reader.open(filename, "r"))
        return false;
    
    char *line;
    string expr, species;
    char *ptr;
    
    // process each line of the file    
    while ((line = reader.readLine())) {
        // skip blank lines
        if (strlen(line) < 2)
            continue;
    
        char *e = strtok_r(line, "\t", &ptr);
        char *s = strtok_r(NULL, "\n", &ptr);
        
        // if bad format, quit
        if (e == NULL || s == NULL)
            return false;
        
        expr = e;
        species = s;
        
        if (expr.size() == 0) {
            // bad gene name expression
            return false;
        } else if (expr[0] == '*') {
            // suffix
            m_rules.append(Gene2speciesRule(Gene2speciesRule::SUFFIX,
                                            expr.substr(1, expr.size()-1), 
                                            species));
        } else if (expr[expr.size() - 1] == '*') {
            // prefix
            m_rules.append(Gene2speciesRule(Gene2speciesRule::PREFIX,
                                            expr.substr(0, expr.size()-1), 
                                            species));
        } else {
            // exact match
            m_exactLookup[expr] = species;
        }
    }

    return true;
}

string Gene2species::getSpecies(string gene)
{
    for (int i=0; i<m_rules.size(); i++) {
        switch (m_rules[i].rule) {
            case Gene2speciesRule::PREFIX:
                if (gene.find(m_rules[i].expr, 0) == 0)
                    return m_rules[i].species;
                break;

            case Gene2speciesRule::SUFFIX:
                if (gene.rfind(m_rules[i].expr, gene.size()-1) == 
                    gene.size() - m_rules[i].expr.size())
                    return m_rules[i].species;
                break;
        }
    }

    // try to find gene in exact match hashtable
    return m_exactLookup[gene];
}

// Returns the gene to species mapping in the 'map' output array
// Only leaf nodes will have defined values for gene2species
bool Gene2species::getMap(string *genes, int ngenes, 
                          string *species, int nspecies, int *map)
{
    // process each gene
    for (int i=0; i<ngenes; i++) {
        // get the species name for the gene
        string sp = getSpecies(genes[i]);

        if (sp.size() == 0) {
            // no species mapping
            // map[i] = -1;
            return false;
        } else {
            // find the index of the species name
            map[i] = -1;
            for (int j=0; j<nspecies; j++) {
                if (sp == species[j]) {
                    map[i] = j;
                    break;
                }
            }
            
            // no species found
            if (map[i] == -1)
                return false;
        }
    }

    return true;
}



} // namespace spidir
