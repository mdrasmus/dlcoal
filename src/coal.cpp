
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






/*

def count_lineages_per_branch(tree, recon, stree):
    """
    Returns the count of gene lineages present at each node in the species
    tree 'tree' given a gene tree 'tree' and reconciliation 'recon'
    """

    # init lineage counts
    lineages = {}
    for snode in stree:
        lineages[snode] = [0, 0]

    for node in tree.postorder():
        snode = recon[node]
        if node.is_leaf():
            lineages[snode][0] += 1 # leaf lineage
        else:
            lineages[snode][1] -= 1 # coal

    for snode in stree.postorder():
        if not snode.is_leaf():
            lineages[snode][0] = sum(lineages[x][1] for x in snode.children)
        lineages[snode][1] += lineages[snode][0]
        
    return lineages


// The function computes terms necessary for many topology calculations
def get_topology_stats(tree, recon, stree):
    """
    The function computes terms necessary for many topology calculations
    """
    
    # How many gene nodes per species
    nodes_per_species = dict.fromkeys(stree, 0)

    # How many descendent nodes recon to the same species
    descend_nodes = {} 

    # iterate through tree
    for node in tree.postorder():
        if len(node.children) > 1:
            nodes_per_species[recon[node]] += 1
            if not node.is_leaf():
                descend_nodes[node] = 1 + sum(descend_nodes.get(child, 0)
                                              for child in node.children
                                              if recon[child] == recon[node])

    return nodes_per_species, descend_nodes



// Returns the log probability of a reconciled gene tree ('tree', 'recon')
// from the coalescent model given a species tree 'stree' and
// population sizes 'n'
double prob_multicoal_recon_topology(
     Tree *tree, int *recon, Tree *stree, float n,
     root=None, leaves=None,
     lineages=None, top_stats=None):
 
    
    popsizes = init_popsizes(stree, n)
    rev_recon = None
    if lineages is None:
        rev_recon = get_rev_recon(tree, recon, stree)
        lineages = count_lineages_per_branch(tree, recon, stree,
                                             rev_recon=rev_recon)
    if top_stats is None:
        top_stats = get_topology_stats(tree, recon, stree,
                                       rev_recon=rev_recon)

    # iterate through species tree branches
    lnp = 0.0 # log probability
    for snode in stree.postorder():
        if snode.parent:
            # non root branch
            a, b = lineages[snode]
            
            lnp += (log(prob_coal_counts(a, b, snode.dist,
                                         popsizes[snode.name]))
                    + stats.logfactorial(top_stats[0].get(snode, 0))
                    - log(num_labeled_histories(a, b)))
        else:
            a = lineages[snode][0]
            lnp += (stats.logfactorial(top_stats[0].get(snode, 0)) -
                    log(num_labeled_histories(a, 1)))

    for node, cnt in top_stats[1].iteritems():
        lnp -= log(cnt)
    
    return lnp;
}
    */


}

} // namespace dlcoal
