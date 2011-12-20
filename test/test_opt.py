# test dlcoal
# test optimizations

import dlcoal

#from rasmus.common import *
from rasmus.testing import *

from compbio import coal
import dlcoal

try:
    import scipy.optimize
except ImportError:
    pass

#=============================================================================

def sample_coal_cond_counts_cdf(a, b, t, n, p=None):
    """
    Samples the next coalescent between 'a' lineages in a population size of
    'n', conditioned on there being 'b' lineages at time 't'.
    """
    
    # this code solves this equation for t
    #   cdf(t) - p = 0
    # where p ~ U(0, 1)

    import scipy.optimize

    if p is None:
        p = random.random()

    # compute constants
    lama = -a*(a-1)/2.0/n
    C0 = stats.prod((b+y)*(a-1-y)/(a-1+y) for y in xrange(b))
    c = -b*(b-1)/2.0/n
    d = 1.0/stats.factorial(b) * (-lama) / coal.prob_coal_counts(a, b, t, n)

    # CDF(t) - p
    def f(x):
        if x <= 0:
            return x - p
        if x >= t:
            return 1.0 - p + (x - t)

        C = C0
        s = exp(c*t) * (exp((lama-c)*x)-1.0) / (lama-c) * C
        for k in xrange(b+1, a):
            k1 = k - 1
            lam = -k*k1/2.0/n
            C = (b+k1)*(a-1-k1)/(a-1+k1)/(b-k) * C
            s += exp(lam*t) * (exp((lama-lam)*x) - 1.0) / (lama - lam) \
                 * (2*k-1) / (k1+b) * C
        print "f", p, x, s * d - p

        return s * d - p
    
    return f


class Opt (unittest.TestCase):

    def test_brent(self):
        def f(x):
            return (x - 2.0)
        
        print scipy.optimize.brentq(f, 0.0, 5.0, disp=False)
        print stats.bisect_root(f, 0.0, 5.0)
        
        fequal(scipy.optimize.brentq(f, 0.0, 5.0, disp=False),
               stats.bisect_root(f, 0.0, 5.0))


    def test_brent2(self):

        a = 5
        b = 2
        t = 2.0
        n = 1.0

        for p in frange(0, 1, .02):
            f = sample_coal_cond_counts_cdf(a, b, t, n, p)
            print "v1"
            v1 = scipy.optimize.brentq(f, 0.0, t, disp=False)
            print "v2"
            v2 = stats.bisect_root(f, 0.0, t, 1e-4)
            print p, v1, v2
            fequal(v1, v2, 1e-2)


    def test_coal(self):

        a = 5
        b = 2
        t = 2.0
        n = 1.0

        for i in range(10):
            x = coal.sample_coal_cond_counts(a, b, t, n)
            


#=============================================================================

if __name__ == "__main__":
    test_main()



