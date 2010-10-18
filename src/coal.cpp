
// c/c++ includes
#include <math.h>
#include <stdio.h>

#include "common.h"
//#include "Tree.h"

using namespace spidir;

namespace dlcoal
{

extern "C" {

//=============================================================================
// primitive tree format conversion functions

struct intnode
{
    int parent;
    int child[2];
};

// creates a int tree from a parent tree
// Note: assumes binary tree
intnode *make_itree(int nnodes, int *ptree)
{
    intnode *itree = new intnode [nnodes];
    
    // initialize
    for (int i=0; i<nnodes; i++) {
        itree[i].parent = ptree[i];
        itree[i].child[0] = -1;
        itree[i].child[1] = -1;
    }
    
    // populate
    for (int i=0; i<nnodes; i++) {
        int parent = ptree[i];
        
        if (parent != -1) {
            if (itree[parent].child[0] == -1)
                itree[parent].child[0] = i;
            else
                itree[parent].child[1] = i;
        }
    }

    return itree;
}


void free_itree(intnode *itree)
{
    delete [] itree;
}



class LineageCounts
{
public:
    LineageCounts(int nsnodes) {
        starts = new int [nsnodes];
        ends = new int [nsnodes];
    }

    ~LineageCounts() {
        delete [] starts;
        delete [] ends;
    }
    
    int *starts;
    int *ends;
};


class TopStats
{
public:
    TopStats(int nnodes, int nsnodes) {    
        nodes_per_species = new int [nsnodes];
        descend_nodes = new int [nnodes];
    }

    ~TopStats() {
        delete [] nodes_per_species;
        delete [] descend_nodes;
    }

    // How many gene nodes per species
    int *nodes_per_species;
     // How many descendent nodes recon to the same species
    int *descend_nodes;
};


class ProbCounts
{
public:
    ProbCounts(int nsnodes) :
        nsnodes(nsnodes)
    {    
        starts = new double* [nsnodes];
        for (int i=0; i<nsnodes; i++) starts[i] = NULL;
        ends = new double* [nsnodes];
        for (int i=0; i<nsnodes; i++) ends[i] = NULL;
    }

    ~ProbCounts() {
        for (int i=0; i<nsnodes; i++) {
            if (starts[i])
                delete [] starts[i];
            if (ends[i])
                delete [] ends[i];
        }
    }

    int nsnodes;
    double **starts;
    double **ends;
};



//=============================================================================

// The probabiluty of going from 'a' lineages to 'b' lineages in time 't'
// with population size 'n'
double prob_coal_counts(int a, int b, double t, double n)
{
    double C = 1.0;
    
    for (int y=0; y<b; y++)
        C *= (b+y)*(a-y)/double(a+y);

    double s = exp(-b*(b-1)*t/2.0/n) * C;

    for (int k=b+1; k<a+1; k++) {
        const double k1 = double(k - 1);
        C *= double(b+k1)*(a-k1)/(a+k1)/(b-k);
        s += exp(-k*k1*t/2.0/n) * (2*k-1) / double(k1+b) * C;
    }
    
    for (int i=1; i<=b; i++)
        s /= i;

    return s;
}


double num_labeled_histories(int nleaves, int nroots)
{
    double n = 1.0;
    for (int i=nroots + 1; i<=nleaves; i++)
        n *= i * (i - 1) / 2.0;
    return n;
}


// Returns the count of gene lineages present at each node in the species
// tree 'tree' given a gene tree 'tree' and reconciliation 'recon'
void count_lineages_per_branch(LineageCounts *counts,
    int nnodes, int *recon, int *pstree, int nsnodes)
{
    const int nleaves = (nnodes + 1) / 2;
    //const int nsleaves = (nsnodes + 1) / 2;
    int *lineages_a = counts->starts;
    int *lineages_b = counts->ends;

    // init lineage counts
    for (int i=0; i<nsnodes; i++) {
        lineages_a[i] = 0;
        lineages_b[i] = 0;
    }

    for (int node=0; node<nleaves; node++)
        lineages_a[recon[node]]++;; // leaf lineage

    for (int node=nleaves; node<nnodes; node++)
        lineages_b[recon[node]]--; // coal

    for (int snode=0; snode<nsnodes-1; snode++) {
        lineages_b[snode] += lineages_a[snode]; // propogate up branch
        lineages_a[pstree[snode]] += lineages_b[snode]; // collect children
    }
    lineages_b[nsnodes-1] += lineages_a[nsnodes-1]; // set sroot
}


// The function computes terms necessary for many topology calculations
void get_topology_stats(TopStats *top_stats, 
   int *ptree, int nnodes, int *recon, int *pstree, int nsnodes)
{
    const int nleaves = (nnodes + 1) / 2;
    int *nodes_per_species = top_stats->nodes_per_species;
    int *descend_nodes = top_stats->descend_nodes;

    // clear stats data structure
    for (int snode=0; snode<nsnodes; snode++)
        nodes_per_species[snode] = 0;

    for (int node=0; node<nnodes; node++)
        descend_nodes[node] = 0;

    // iterate through tree
    for (int node=nleaves; node<nnodes; node++) {
        nodes_per_species[recon[node]]++;
        descend_nodes[node]++;

        const int parent = ptree[node];
        if (recon[node] == recon[parent])
            descend_nodes[parent] += descend_nodes[node];
    }
}


// The function computes terms necessary for many topology calculations
void get_topology_stats2(TopStats *top_stats, 
   intnode *itree, int nnodes, int *recon, int *pstree, int nsnodes)
{
    const int nleaves = (nnodes + 1) / 2;
    int *nodes_per_species = top_stats->nodes_per_species;
    int *descend_nodes = top_stats->descend_nodes;

    // clear stats data structure
    for (int snode=0; snode<nsnodes; snode++)
        nodes_per_species[snode] = 0;

    for (int node=0; node<nnodes; node++)
        descend_nodes[node] = 0;

    // iterate through tree
    for (int node=nleaves; node<nnodes; node++) {
        nodes_per_species[recon[node]]++;
        descend_nodes[node]++;

        const int parent = itree[node].parent;
        if (recon[node] == recon[parent])
            descend_nodes[parent] += descend_nodes[node];
    }
}



// Returns the log probability of a reconciled gene tree ('tree', 'recon')
// from the coalescent model given a species tree 'stree' and
// population sizes 'n'
double prob_multicoal_recon_topology(int *ptree, int nnodes, int *recon, 
                                     int *pstree, int nsnodes, 
                                     double *sdists, double *popsizes)
{
    LineageCounts counts(nsnodes);
    count_lineages_per_branch(&counts, nnodes, recon, pstree, nsnodes);

    TopStats top_stats(nnodes, nsnodes);
    get_topology_stats(&top_stats, ptree, nnodes, recon, pstree, nsnodes);
    
    // iterate through species tree branches
    double lnp = 0.0; // log probability
    for (int snode=0; snode<nsnodes-1; snode++) {
        const int a = counts.starts[snode];
        const int b = counts.ends[snode];
        //if (a == 1)
        //    continue;

        const int n = top_stats.nodes_per_species[snode];
        double fact=1.0; for (int i=2; i<=n; i++) fact *= i; // fact(n)
            
        double p = log(prob_coal_counts(a, b, sdists[snode], popsizes[snode])
                   * fact / num_labeled_histories(a, b));
        lnp += p;
    }

    // root branch
    const int snode = nsnodes - 1;
    const int a = counts.starts[snode];
    const int n = top_stats.nodes_per_species[snode];
    double fact=1.0; for (int i=2; i<=n; i++) fact *= i; // fact(n)
    lnp += log(fact / num_labeled_histories(a, 1));

    const int nleaves = (nnodes + 1) / 2;
    for (int node=nleaves; node<nnodes; node++)
        lnp -= log(top_stats.descend_nodes[node]);

    return lnp;
}


// Returns the log probability of a reconciled gene tree ('tree', 'recon')
// from the coalescent model given a species tree 'stree' and
// population sizes 'n'
double prob_multicoal_recon_topology2(intnode *itree, int nnodes, int *recon, 
                                      int *pstree, int nsnodes, 
                                      double *stimes, double *popsizes)
{
    LineageCounts counts(nsnodes);
    count_lineages_per_branch(&counts, nnodes, recon, pstree, nsnodes);

    TopStats top_stats(nnodes, nsnodes);
    get_topology_stats2(&top_stats, itree, nnodes, recon, pstree, nsnodes);
    
    // iterate through species tree branches
    double lnp = 0.0; // log probability
    for (int snode=0; snode<nsnodes-1; snode++) {
        const int a = counts.starts[snode];
        const int b = counts.ends[snode];
        if (a == 1)
            continue;

        const int n = top_stats.nodes_per_species[snode];
        double fact=1.0; for (int i=2; i<=n; i++) fact *= i; // fact(n)
            
        const double t = stimes[pstree[snode]] - stimes[snode];
        lnp += log(prob_coal_counts(a, b, t, popsizes[snode])
                   * fact / num_labeled_histories(a, b));
    }

    // root branch
    const int snode = nsnodes - 1;
    const int a = counts.starts[snode];
    const int n = top_stats.nodes_per_species[snode];
    double fact=1.0; for (int i=2; i<=n; i++) fact *= i; // fact(n)
    lnp += log(fact / num_labeled_histories(a, 1));

    const int nleaves = (nnodes + 1) / 2;
    for (int node=nleaves; node<nnodes; node++)
        lnp -= log(top_stats.descend_nodes[node]);

    return lnp;
}



// use dynamic programming to calc prob of lineage counts
void calc_prob_counts_table(ProbCounts *prob_counts,
                            int *gene_counts, double T, 
                            intnode *istree, int nsnodes, 
                            double *popsizes,
                            int sroot, int *sleaves, int nsleaves,
                            double *stimes)
{
    // array of max number of lineages per snode
    int* sizes = new int [nsnodes];

    // stack of snodes to visit
    int *stack = new int [nsnodes];
    int *stack_pushes = new int [nsnodes];
    for (int i=0; i<nsnodes; i++) stack_pushes[i] = 0;

    // push sleaves onto stack
    for (int i=0; i<nsleaves; i++)
        stack[i] = sleaves[i];
    int stack_len = nsleaves;

    // post-order traversal of species tree
    for (int stack_i=0; stack_i<stack_len; stack_i++) {
        const int snode = stack[stack_i];

        // push parent
        // reached root
        const int sparent = istree[snode].parent;
        if (snode != sroot && sparent >= 0)
            if (++stack_pushes[sparent] == 2)
                stack[stack_len++] = sparent;

        double *start;
        int M;

        if (stack_i < nsleaves) {
            // leaf case
            M = gene_counts[snode];
            
            // populate starting lineage counts
            start = new double [M+1];
            for (int i=0; i<M; i++) start[i] = 0.0;
            start[M] = 1.0;
        } else {
            // internal node case
            const int c1 = istree[snode].child[0];
            const int c2 = istree[snode].child[1];
            const int M1 = sizes[c1];
            const int M2 = sizes[c2];
            M = M1 + M2; // max lineage counts in this snode
            double *end1 = prob_counts->ends[c1];
            double *end2 = prob_counts->ends[c2];

            // populate starting lineage counts
            start = new double [M+1];
            start[0] = 0.0; start[1] = 0.0;
            for (int k=2; k<=M; k++) {
                start[k] = 0.0;
                for (int i=0; i<k; i++)
                    if (i <= M1 && k-i <= M2)
                        start[k] += end1[i] * end2[k-i];
            }
        }
        sizes[snode] = M;

        // populate ending lineage counts
        const double n = popsizes[snode];
        double ptime = (sparent >= 0) ? stimes[istree[snode].parent] : T;
        double *end = new double [M+1];
        if (ptime < 0) {
            // unbounded end time, i.e. complete coalescence
            for (int i=0; i<=M; i++) end[i] = 0.0;
            end[1] = 1.0;
        } else {
            // fixed end time
            const double t = ptime - stimes[snode];

            end[0] = 0.0;
            for (int k=1; k<=M; k++) {
                end[k] = 0.0;
                for (int i=k; i<=M; i++) 
                    end[k] += prob_coal_counts(i, k, t, n) * start[i];
            }
        }

        prob_counts->starts[snode] = start;
        prob_counts->ends[snode] = end;
    }

    delete [] sizes;
}



double prob_locus_coal_recon_topology(int *ptree, int nnodes, int *recon, 
                                      int *plocus_tree, intnode *iltree,
                                      int nlocus_nodes, 
                                      double *popsizes, double *ltimes,
                                      int *daughters, int ndaughters)
{
    bool own_iltree = false;
    if (!iltree) {
        own_iltree = true;
        iltree = make_itree(nlocus_nodes, plocus_tree);
    }
    intnode *itree = make_itree(nnodes, ptree);

    double lnp = prob_multicoal_recon_topology2(itree, nnodes, recon, 
                                                plocus_tree, nlocus_nodes, 
                                                ltimes, popsizes);

    const int nleaves = (nnodes + 1) / 2;

    ProbCounts prob_counts(nlocus_nodes);
    int *stack = new int [nnodes];
    int stack_len = 0;

    int *gene_counts = new int [nlocus_nodes];
    int *subleaves = new int [nlocus_nodes];
    int nsubleaves;

    LineageCounts counts(nlocus_nodes);
    count_lineages_per_branch(&counts, nnodes, recon, 
                              plocus_tree, nlocus_nodes);

    // make daughter set
    bool *daughters_set = new bool [nlocus_nodes];
    for (int i=0; i<nlocus_nodes; i++) daughters_set[i] = false;
    for (int i=0; i<ndaughters; i++) daughters_set[daughters[i]] = true;

    // find relevant subtree
    for (int i=0; i<ndaughters; i++) {
        const int daughter = daughters[i];

        // determine leaves of the coal subtree
        stack_len = 1;
        stack[0] = daughter;
        nsubleaves = 0;
        for (int stack_i=0; stack_i<stack_len; stack_i++) {
            const int lnode = stack[stack_i];
            if (lnode < nleaves) {
                // leaf
                gene_counts[lnode] = counts.starts[lnode];
                subleaves[nsubleaves++] = lnode;
            } else {
                for (int i=0; i<2; i++) {
                    const int child = iltree[lnode].child[i];
                    if (daughters_set[child]) {
                        gene_counts[child] = 1;
                        subleaves[nsubleaves++] = child;
                    } else {
                        // push child
                        stack[stack_len++] = child;
                    }
                }
            }
        }

        const double T = ltimes[plocus_tree[daughter]];
        calc_prob_counts_table(&prob_counts,
                               gene_counts, T, 
                               iltree, nlocus_nodes, 
                               popsizes,
                               daughter, subleaves, nsubleaves,
                               ltimes);
        lnp -= log(prob_counts.ends[daughter][1]);

        if (lnp == -INFINITY)
            break;
    }


     // deallocates 'forward tree'
    if (own_iltree)
        free_itree(iltree);
    free_itree(itree);

    delete [] stack;
    delete [] gene_counts;
    delete [] subleaves;


    return lnp;
}


    /*
def prob_locus_coal_recon_topology(tree, recon, locus_tree, n, daughters):
    """
    Returns the log probability of a reconciled gene tree ('tree', 'recon')
    from the coalescent model given a locus tree 'locus_tree',
    population sizes 'n', and daughters set 'daughters'
    """

    # initialize popsizes, lineage counts, and divergence times
    popsizes = coal.init_popsizes(locus_tree, n)
    lineages = coal.count_lineages_per_branch(tree, recon, locus_tree)
    locus_times = treelib.get_tree_timestamps(locus_tree)


    # calc log probability
    lnp = coal.prob_multicoal_recon_topology(
        tree, recon, locus_tree, popsizes, lineages=lineages)
    

    def walk(node, gene_counts, leaves):
        if node.is_leaf():
            gene_counts[node.name] = lineages[node][0]
            leaves.add(node)
        else:
            for child in node.children:
                if child in daughters:
                    gene_counts[child.name] = 1
                    leaves.add(child)
                else:
                    walk(child, gene_counts, leaves)

    for daughter in daughters:
        # determine leaves of the coal subtree
        gene_counts = {}
        leaves = set()
        walk(daughter, gene_counts, leaves)

        p = coal.cdf_mrca_bounded_multicoal(
            gene_counts, locus_times[daughter.parent], locus_tree, popsizes,
            sroot=daughter, sleaves=leaves, stimes=locus_times)

        if p == -util.INF:
            return -util.INF

        lnp -= p
    
    return lnp



    */

//=============================================================================
// events
enum {
    EVENT_GENE = 0,
    EVENT_SPEC = 1,
    EVENT_DUP = 2
};


//  Probability density for for next birth at time 't' given
//  'n'=1 lineages starting at time 0, evolving until time 'T' with a
//  'birth' and 'death' rates for a reconstructed process.
double birthWaitTime1(float t, float T, float birth, float death,
		    float denom)
{
    const double r = birth - death;

    return r * exp(-r*t) / denom;
}

//  Probability density for for next birth at time 't' given
//  'n'=1 lineages starting at time 0, evolving until time 'T' with a
//  'birth' and 'death' rates for a reconstructed process.
double birthWaitTimeDenom1(float T, float birth, float death)
{
    const double r = birth - death;
    const double a = death / birth;

    return 1.0 - a * exp(-r * T);
}


// Sample the next birth event from a reconstructed birthdeath process.
// Let there be 'n'=1 lineages at time 0 that evolve until time 'T' with
// 'birth' and 'death' rates.
// Conditioned that a birth will occur
double sampleBirthWaitTime1(float T, float birth, float death)
{
    
    // TODO: could make this much more efficient (use straight line instead of
    // flat line).
    
    // uses rejection sampling
    double denom = birthWaitTimeDenom1(T, birth, death);
    double start_y = birthWaitTime1(0, T, birth, death, denom);
    double end_y = birthWaitTime1(T, T, birth, death, denom);
    double M = max(start_y, end_y);
    
    while (true) {
        double t = frand(T);
        double f = birthWaitTime1(t, T, birth, death, denom);

        if (frand() <= f / M)
            return t;
    }
}



// Sample duplication times for only a subtree
void sample_dup_times_subtree(double *times, intnode *itree, 
                              double start_time, double time_span, 
                              int duproot, 
                              int *recon, int *events,
                              double birth, double death,
                              int *stack)
{
    // init stack
    int stack_len = 1;
    stack[0] = duproot;

    // walk duplication subtree
    for (int stack_i=0; stack_i<stack_len; stack_i++) {
        const int node = stack[stack_i];
        const double parent_time = (stack_i == 0) ? 
            start_time : times[itree[node].parent];
        const double remain = time_span - (start_time - parent_time);

        if (remain <= 0.0) {
            printf("node %d %e %e %e\n", node, 
                   time_span, start_time, parent_time);
            assert(false);
        }
        double t;
        do {
            t = sampleBirthWaitTime1(remain, birth, death);
            times[node] = parent_time - t;
        } while (times[node] == parent_time);

        const int snode = recon[node];
        for (int i=0; i<2; i++) {
            const int child = itree[node].child[i];
            if (events[child] == EVENT_DUP && recon[child] == snode) {
                // push child on to stack
                stack[stack_len++] = child;
            }
        }
    }
}



void sample_dup_times(double *times,
                      intnode *itree, int nnodes, int *pstree, int nsnodes,
                      double *stimes,
                      int *recon, int *events, double birth, double death,
                      double pretime, double premean, int *stack)
{
    

    const int root = nnodes - 1;
    const int sroot = nsnodes - 1;


    // set pretimes
    if (events[root] != EVENT_SPEC) {
        int snode;
        double start_time, time_span;
        
        if (recon[root] != sroot) {
            // tree root is a dup within species tree
            snode = recon[root];
            start_time = stimes[pstree[snode]];
            time_span = stimes[pstree[snode]] - stimes[snode];
        } else {
            // tree root is a pre-spec dup
            if (pretime < 0) {
                assert (premean >= 0);

                pretime = 0.0;
                do {
                    pretime = expovariate(1.0/premean);
                } while (pretime == 0.0);
            }
            start_time = stimes[sroot] + pretime;
            time_span = pretime;
        }
        
        sample_dup_times_subtree(times, itree, 
                                 start_time, time_span, 
                                 root, 
                                 recon, events,
                                 birth, death, stack);
    }

    // set times
    for (int node=nnodes-1; node>=0; node--) {
        const int parent = itree[node].parent;
        const int snode = recon[node];

        if (events[node] == EVENT_SPEC) {
            // set speciation time
            times[node] = stimes[recon[node]];

        } else if (events[node] == EVENT_DUP && parent >= 0 &&
                   snode != recon[parent])
        {
            // set duplication times within duplication subtree
            // node is duproot
            const double start_time = stimes[pstree[snode]];
            const double time_span = start_time - stimes[snode];
            
            sample_dup_times_subtree(times, itree, 
                                     start_time, time_span, 
                                     node, 
                                     recon, events,
                                     birth, death, stack);

        } else if (events[node] == EVENT_GENE) {
            times[node] = 0.0;
        }
    }
}


double prob_locus_coal_recon_topology_samples(
    int *ptree, int nnodes, int *recon, 
    int *plocus_tree, int nlocus_nodes, 
    int *locus_recon, int *locus_events,
    double *popsizes, 
    int *pstree, int nsnodes, double *stimes,
    int *daughters, int ndaughters, 
    double birth, double death,
    int nsamples, double pretime, double premean)
{
    // alloc datastructures
    double *ltimes = new double [nlocus_nodes];
    intnode *iltree = make_itree(nlocus_nodes, plocus_tree);
    int *stack = new int [nnodes];
    
    // integrate over duplication times using sampling
    double prob = -INFINITY;
    for (int i=0; i<nsamples; i++) {
        // sample duplication times
        sample_dup_times(ltimes,
                         iltree, nlocus_nodes, pstree, nsnodes,
                         stimes,
                         locus_recon, locus_events, birth, death,
                         pretime, premean, stack);
        
        // coal topology probability
        double const coal_prob = prob_locus_coal_recon_topology(
            ptree, nnodes, recon, 
            plocus_tree, iltree, nlocus_nodes, 
            popsizes, ltimes,
            daughters, ndaughters);

        prob = logadd(prob, coal_prob);
    }

    // clean up
    delete [] ltimes;
    free_itree(iltree);
    delete [] stack;

    return prob - log(nsamples);
}




}

} // namespace dlcoal
