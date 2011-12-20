
import sys, copy

import dlcoal
from dlcoal import duploss

from rasmus import util, stats, treelib

from compbio import phylo, coal


#=============================================================================
# reconciliation

def dlcoal_recon(tree, stree, gene2species,
                 n, duprate, lossrate,
                 pretime=None, premean=None,
                 nsearch=1000,
                 nsamples=100, nprescreen=20,
                 search=None,
                 init_locus_tree=None,
                 log=sys.stdout):
    """
    Perform reconciliation using the DLCoal model

    Returns maxrecon defined as

    maxrecon = {'coal_recon': coal_recon,
                'locus_tree': locus_tree,
                'locus_recon': locus_recon,
                'locus_events': locus_events,
                'daughters': daughters,
                'data': extra_information }
    
    """

    if search is None:
        search = lambda tree: DLCoalTreeSearch(tree, stree, gene2species,
                                               duprate, lossrate,
                                               nprescreen=nprescreen)

    reconer = DLCoalRecon(tree, stree, gene2species,
                          n, duprate, lossrate,
                          pretime=pretime, premean=premean,
                          nsamples=nsamples, log=log,
                          init_locus_tree=init_locus_tree)
    reconer.set_proposer(DLCoalReconProposer(
        tree, stree, gene2species, search=search))
    return reconer.recon(nsearch).get_dict()



class DLCoalRecon (object):

    def __init__(self, tree, stree, gene2species,
                 n, duprate, lossrate,
                 pretime=None, premean=None,
                 nsamples=100,
                 init_locus_tree=None,
                 name_internal="n", log=sys.stdout):

        # init coal tree
        self.coal_tree = tree
        self.stree = stree
        self.gene2species = gene2species
        self.n = n
        self.duprate = duprate
        self.lossrate = lossrate
        self.pretime = pretime
        self.premean = premean
        self.nsamples = nsamples
        self.name_internal = name_internal
        self.log_stream = log
        self.init_locus_tree = init_locus_tree \
                               if init_locus_tree else tree.copy()

        self.proposer = DLCoalReconProposer(tree, stree, gene2species)


    def set_proposer(self, proposer):
        """Set the proposal algorithm"""
        self.proposer = proposer


    def set_log(self, log):
        self.log_stream = log
        

    def recon(self, nsearch=1000):
        """Perform reconciliation"""
        
        self.init_search()
        proposal = self.proposer.init_proposal()
        self.maxrecon = proposal.copy()
        for i in xrange(nsearch):
            if i % 10 == 0:
                print "search", i

            util.tic("eval")
            p = self.eval_proposal(proposal)
            util.toc()

            util.tic("prop")
            self.eval_search(p, proposal)
            proposal = self.proposer.next_proposal()
            util.toc()
        
        # rename locus tree nodes
        dlcoal.rename_nodes(self.maxrecon.locus_tree, self.name_internal)
        
        return self.maxrecon


    def init_search(self):
        """Initialize new search"""

        # init locus tree as congruent to coal tree
        # equivalent to assuming no ILS
        #self.proposer.set_locus_tree(self.coal_tree.copy())
        self.proposer.set_locus_tree(self.init_locus_tree.copy())
        
        self.maxp = - util.INF
        self.maxrecon = None


    def next_proposal(self):
        """Returns next proposal"""
        self.proposal.next_proposal()


    def eval_proposal(self, proposal):
        """Compute probability of proposal"""

        #util.tic("eval")
        # compute recon probability
        info = {}

        # DEBUG
        counts = coal.count_lineages_per_branch(self.coal_tree,
                                                proposal.coal_recon,
                                                proposal.locus_tree)
        maxcount = max(x[0] for x in counts.values())
        #util.logger("max lineage count %d" % maxcount)
        if maxcount > 10:
            return -util.INF
        
        p = dlcoal.prob_dlcoal_recon_topology(self.coal_tree,
                                              proposal.coal_recon,
                                              proposal.locus_tree,
                                              proposal.locus_recon,
                                              proposal.locus_events,
                                              proposal.daughters,
                                              self.stree, self.n,
                                              self.duprate, self.lossrate,
                                              self.pretime, self.premean,
                                              nsamples=self.nsamples,
                                              add_spec=False,
                                              info=info)
        proposal.data = info
        #util.toc()
        
        return p


    def eval_search(self, p, proposal):
        """Evaluate a proposal for search"""
        
        self.log_proposal(proposal)

        if p > self.maxp:
            self.maxp = p
            self.maxrecon = proposal.copy()

            # search with a new copy
            self.proposer.accept()
        else:
            self.proposer.reject()


    def log_proposal(self, proposal):
        self.log_stream.write(repr(proposal) + "\n")
        self.log_stream.flush()



class DLCoalReconProposer (object):

    def __init__(self, coal_tree, stree, gene2species,
                 search=phylo.TreeSearchNni,
                 num_coal_recons=1): # DEBUG
        self._coal_tree = coal_tree
        self._stree = stree
        self._gene2species = gene2species
        self._locus_search = search(None)

        # coal recon search
        self._num_coal_recons = num_coal_recons
        self._i_coal_recons = 0
        self._coal_recon_enum = None
        self._coal_recon_depth = 2
        self._accept_locus = False

        self._recon = None
        

    def set_locus_tree(self, locus_tree):
        self._locus_search.set_tree(locus_tree)

    def init_proposal(self):
        """Get first proposal"""

        if self._locus_search.get_tree() is None:
            self._locus_search.set_tree(self._coal_tree.copy())
        self._i_coal_recons = 0
        self._recon = self._recon_lca(self._locus_search.get_tree().copy())
        
        return self._recon


    def next_proposal(self):

        if len(self._locus_search.get_tree().leaves()) <= 2:
            return self._recon
        
        if self._i_coal_recons >= self._num_coal_recons:
            # propose new locus_tree
            
            # if locus_tree has not yet been accepted, then revert it
            if not self._accept_locus:
                self._locus_search.revert()
                
            self._locus_search.propose()
            self._accept_locus = False
            self._i_coal_recons = 0
            locus_tree = self._locus_search.get_tree().copy()
            
            # TODO: make recon root optional
            phylo.recon_root(locus_tree, self._stree,
                             self._gene2species,
                             newCopy=False)
            dlcoal.rename_nodes(locus_tree)

            # propose remaining parts of dlcoal recon
            self._recon = self._recon_lca(locus_tree)
        else:
            # modify coal_recon
            
            try:
                self._i_coal_recons += 1
                self._coal_recon_enum.next()
            except StopIteration:
                self._i_coal_recon = self._num_coal_recons
                return self.next_proposal()

        return self._recon


    def _recon_lca(self, locus_tree):
        # get locus tree, and LCA locus_recon
        locus_recon = phylo.reconcile(locus_tree, self._stree,
                                      self._gene2species)
        locus_events = phylo.label_events(locus_tree, locus_recon)

        # propose LCA coal_recon
        coal_recon = phylo.reconcile(self._coal_tree,
                                     locus_tree, lambda x: x)

        # propose daughters (TODO)
        daughters = self._propose_daughters(
            self._coal_tree, coal_recon,
            locus_tree, locus_recon, locus_events)


        self._coal_recon_enum = phylo.enum_recon(
            self._coal_tree, locus_tree,
            recon=coal_recon,
            depth=self._coal_recon_depth)


        return Recon(coal_recon, locus_tree, locus_recon, locus_events,
                     daughters)


    def _propose_daughters(self, coal_tree, coal_recon,
                           locus_tree, locus_recon, locus_events):
        return propose_daughters(coal_tree, coal_recon,
                                 locus_tree, locus_events)

    def accept(self):
        self._accept_locus = True
    

    def reject(self):
        pass


def propose_daughters(coal_tree, coal_recon, locus_tree, locus_events):

    lineages = coal.count_lineages_per_branch(
        coal_tree, coal_recon, locus_tree)
    daughters = set()

    for node, event in locus_events.iteritems():
        if event == "dup":
            # choose one of the children of node to be a daughter
            children = [child for child in node.children
                        if lineages[child][1] == 1]
            if len(children) > 0:
                daughters.add(children[stats.sample([1] * len(children))])

    return daughters



class Recon (object):
    """
    The reconciliation datastructure for the DLCoal model
    """
    
    def __init__(self, coal_recon, locus_tree, locus_recon, locus_events,
                 daughters, data=None):
        self.coal_recon = coal_recon
        self.locus_tree = locus_tree
        self.locus_recon = locus_recon
        self.locus_events = locus_events
        self.daughters = daughters

        if data is None:
            self.data = {}
        else:
            self.data = data

    def copy(self):
        return Recon(self.coal_recon,
                     self.locus_tree, self.locus_recon, self.locus_events,
                     self.daughters, data=copy.deepcopy(self.data))


    def get_dict(self):
        return {"coal_recon": self.coal_recon,
                "locus_tree": self.locus_tree,
                "locus_recon": self.locus_recon,
                "locus_events": self.locus_events,
                "daughters": self.daughters}
    
    
    def __repr__(self):
        return repr({"coal_recon": [(x.name, y.name) for x,y in
                                    self.coal_recon.iteritems()],
                     "locus_tree": self.locus_tree.get_one_line_newick(
                         root_data=True),
                     "locus_top": phylo.hash_tree(self.locus_tree),
                     "locus_recon": [(x.name, y.name) for x,y in
                                     self.locus_recon.iteritems()],
                     "locus_events": [(x.name, y) for x,y in
                                      self.locus_events.iteritems()],
                     "daughters": [x.name for x in self.daughters],
                     "data": self.data})


#=============================================================================
# tree search

class DLCoalTreeSearch (phylo.TreeSearch):

    def __init__(self, tree, stree, gene2species, duprate, lossrate,
                 tree_hash=None, nprescreen=20, weight=.2):
        phylo.TreeSearch.__init__(self, tree)

        self.stree = stree
        self.gene2species = gene2species
        self.duprate = duprate
        self.lossrate = lossrate
        
        #self.search = UniqueTreeSearch(tree, phylo.TreeSearchNni(tree),
        #                               tree_hash)
        #self.search = UniqueTreeSearch(tree, phylo.TreeSearchSpr(tree),
        #                               tree_hash)

        mix = phylo.TreeSearchMix(tree)
        mix.add_proposer(phylo.TreeSearchNni(tree), .4)
        mix.add_proposer(phylo.TreeSearchSpr(tree), .6)        
        #self.search = phylo.TreeSearchUnique(tree, mix, tree_hash)

        prescreen = phylo.TreeSearchPrescreen(tree, mix,
                                              self.prescreen,
                                              poolsize=nprescreen)

        mix2 = phylo.TreeSearchMix(tree)
        mix2.add_proposer(prescreen, 1.0-weight)
        mix2.add_proposer(mix, weight)

        self.search = mix2
        



    def set_tree(self, tree):
        self.tree = tree
        self.search.set_tree(tree)

    def reset(self):
        self.search.reset()

    def propose(self):
        self.search.propose()
        return self.tree
        
    def revert(self):
        self.search.revert()
        return self.tree


    def prescreen(self, tree):

        recon = phylo.reconcile(tree, self.stree, self.gene2species)
        events = phylo.label_events(tree, recon)

        #print tree.root.name
        #treelib.draw_tree_names(tree, maxlen=8)
        
        return duploss.prob_dup_loss(
            tree, self.stree, recon, events,
            self.duprate, self.lossrate)

    
