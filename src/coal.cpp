
// c/c++ includes
#include <math.h>
#include <stdio.h>

namespace dlcoal
{

extern "C" {



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


// Returns the count of gene lineages present at each node in the species
// tree 'tree' given a gene tree 'tree' and reconciliation 'recon'
void count_lineages_per_branch(LineageCounts *counts,
    int *ptree, int nnodes, int *recon, int *pstree, int nsnodes)
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




// Returns the log probability of a reconciled gene tree ('tree', 'recon')
// from the coalescent model given a species tree 'stree' and
// population sizes 'n'
double prob_multicoal_recon_topology(int *ptree, int nnodes, int *recon, 
                                     int *pstree, int nsnodes, 
                                     double *sdists, double *popsizes)
{
    LineageCounts counts(nsnodes);
    count_lineages_per_branch(&counts, ptree, nnodes, recon, pstree, nsnodes);

    TopStats top_stats(nnodes, nsnodes);
    get_topology_stats(&top_stats, ptree, nnodes, recon, pstree, nsnodes);
    
    // iterate through species tree branches
    double lnp = 0.0; // log probability
    for (int snode=0; snode<nsnodes-1; snode++) {
        const int a = counts.starts[snode];
        const int b = counts.ends[snode];
        if (a == b)
            continue;

        const int n = top_stats.nodes_per_species[snode];
        double fact=1.0; for (int i=2; i<=n; i++) fact *= i; // fact(n)
            
        lnp += log(prob_coal_counts(a, b, sdists[snode], popsizes[snode])
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


    /*

def calc_prob_counts_table(gene_counts, T, stree, popsizes,
                            sroot, sleaves, stimes):

    # use dynamic programming to calc prob of lineage counts
    # format: prob_counts[node] = [a, b]
    prob_counts = {}
    def walk(node):
        if node in sleaves:
            # leaf case
            M = gene_counts[node.name]

            # populate starting lineage counts
            start = [0.0] * (M+1)
            start[M] = 1.0
        
        else:
            # internal node case
            assert len(node.children) == 2
                        
            c1 = node.children[0]
            c2 = node.children[1]
            M1 = walk(c1)
            M2 = walk(c2)
            M = M1 + M2 # max lineage counts in this snode
            end1 = prob_counts[c1][1]
            end2 = prob_counts[c2][1]

            # populate starting lineage counts
            start = [0.0, 0.0]
            for k in xrange(2, M+1):
                start.append(sum(end1[i] * end2[k-i]
                                 for i in xrange(1, k)
                                 if i <= M1 and k-i <= M2))

        # populate ending lineage counts
        n = popsizes[node.name]
        ptime = stimes[node.parent] if node.parent else T
        if ptime is None:
            # unbounded end time, i.e. complete coalescence
            end = [0.0, 1.0] + [0.0] * (M-1)
        else:
            # fixed end time
            t = ptime - stimes[node]

            end = [0.0]
            for k in xrange(1, M+1):
                end.append(
                    sum(prob_coal_counts(i, k, t, n) * start[i]
                        for i in xrange(k, M+1)))

        prob_counts[node] = [start, end]

        assert abs(sum(start) - 1.0) < .001
            
        return M
    M = walk(sroot)

    return prob_counts



    */


}

} // namespace dlcoal
