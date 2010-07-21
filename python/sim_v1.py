#!/usr/bin/env python

## first attempt at removing the immediate fixation assumption from the D/L sim
##  should work for a single branch of the 'species tree'


from math import *
import random
from rasmus import stats
from compbio import coal

DEBUGPRINT = True

def debugprint(s):
    if DEBUGPRINT:
        print s

def sim_branch(N=1e6, T=1e4, S=1e1, dr=.0012, lr=.0011, Fd=None, Fl=None, Fi=[1e0]):
    """ 
    Run a Poisson-based D/L simulation on a single branch.
    N is population size, 
    T is total time to run the sim, 
    S is the number of steps in the time period T,
    dr is the duplication rate (in events/individual/myr), 
    lr is the loss rate (in events/individual/myr),
    Fd is the frequency for dups,
    Fl is the frequency for losses, 
    Fi is the initial frequency. 
    """
    if Fd == None:
        Fd = 1.0/N
    if Fl == None:
        Fl = 1.0/N
    
    if type(Fi) != list: # ensure Fi is a list or an int
        Fi = [Fi] if type(Fi) == int or type(Fi) == float else []
    
    if len(Fi) == 0:  # should occur only if the branch goes extinct (and all sub-branches)
        return None # TODO: may want a different return value
    
    if T <= 0.0 or S <= 0.0: # simulation has finished its time
        return flatten(Fi) # TODO: ensure this is the proper return value
    
    for j in xrange(len(Fi)):
        Fi[j] = min(1.0,max(0.0, Fi[j])) # ensure each frequency falls in [0.0,1.0]
        Fi[j] = round(N * Fi[j]) / N # kill off unrealistic frequencies (esp. near 0.0)
    
    ## kill dead branches
    j = 0
    while j < len(Fi):
        if Fi[j] <= 0.0:
            Fi = Fi[:j] + Fi[j+1:]
        else:
            j += 1
    if len(Fi) == 0:
        return None # TODO: may want a different return value
    
    debugprint('Time left: ' + str(T))
    debugprint(' Before step: ' + str(Fi))
    
    newFi = sim_step(N,T/S,dr,lr,Fd,Fl,Fi) # run simulation for a time step
    
    return sim_branch(N,T-T/S,S-1,dr,lr,Fd,Fl,newFi) # iterate to next step, decreasing total time left


def sim_step(N, time, dr=.0012, lr=.0011, Fd=1e-5, Fl=1e-5, Fi=[1e0]):
    """
    Run a Poisson-based D/L simulation for a single time step.
    Note that all inputs are assumed to be sanitized (no checks, unlike sim_branch).
    """
    for j in xrange(len(Fi)):
        Fi[j] = sim_step_indiv(N, time, dr, lr, Fd, Fl, Fi[j]) # get frequency after D/L events and WF transforms
    
    # TODO: determine if the following steps are needed here...
    Fi = flatten(Fi)
    
    debugprint(" After step: " + str(Fi))
    
    j = 0
    while j < len(Fi):
        if Fi[j] <= 0.0:
            Fi = Fi[:j] + Fi[j+1:]
            debugprint("  Removed a dead branch")
        else:
            j += 1
    
    return flatten(Fi) # flatten to avoid [foo, [bar, baz]] nesting


def sim_step_indiv(N, time, dr=.0012, lr=.0011, Fd=1e-5, Fl=1e-5, pi=1e0):
    """
    Run a Poisson-based D/L simulation on an individual frequency.
    All inputs are assumed sanitized
    """
    def chooseevent(subrate, fullrate):
        return random.random() <= subrate / fullrate

    drpyr = dr / 1e6 # dr in events/myr/indiv; drpyr in events/yr/indiv
    lrpyr = lr / 1e6 # lr in events/myr/indiv; lrpyr in events/yr/indiv
    edrpyr = drpyr * pi * N # get population dup rate per year   TODO: may need *N?
    elrpyr = lrpyr * pi * N # get population loss rate per year  TODO: see note above
    eventrate = edrpyr + elrpyr # get total event rate per year

    numlosses = 0

    clock = time
    newp = [pi]

    # simulate Poisson process for D/L events
    while clock > 0.0:
        eventtime = stats.exponentialvariate(eventrate) # time to next event
        if eventtime > clock:
            clock = 0.0
            break
        clock -= eventtime
        if chooseevent(edrpyr, eventrate): # determine whether event is D or L
            newp.append(sim_step_indiv(N, clock, dr, lr, Fd, Fl, Fd))
            debugprint("  Duplication")
        else:
            numlosses += 1 # all losses in period grouped together (approx)

    # poor approx?
    if numlosses > 0:
        newp[0] = max(newp[0] - numlosses * Fl, 0.0)
        debugprint("  Losses: " + str(numlosses))
    newp[0] = (coal.sample_freq_CDF(newp[0], N, time))

    return flatten(newp) # needs flattening?





## taken from http://kogs-www.informatik.uni-hamburg.de/~meine/python_tricks
## TODO: rewrite without borrowing code
def flatten(x):
    """flatten(sequence) -> list

    Returns a single, flat list which contains all elements retrieved
    from the sequence and all recursively contained sub-sequences
    (iterables).

    Examples:
    >>> [1, 2, [3,4], (5,6)]
    [1, 2, [3, 4], (5, 6)]
    >>> flatten([[[1,2,3], (42,None)], [4,5], [6], 7, MyVector(8,9,10)])
    [1, 2, 3, 42, None, 4, 5, 6, 7, 8, 9, 10]
    """

    result = []
    for el in x:
        #if isinstance(el, (list, tuple)):
        if hasattr(el, "__iter__") and not isinstance(el, basestring):
            result.extend(flatten(el))
        else:
            result.append(el)
    return result



if __name__ == "__main__":
    N = 1e4
    T = 1e5
    S = 1e1
    dr = .0012
    lr = .0011
    Fd = .05 #1.0/N
    Fl = .05 #1.0/N
    Fi = .5
    ans = sim_branch(N=N, T=T, S=S, dr=dr, lr=lr, Fd=Fd, Fl=Fl, Fi=Fi)
    print ans
    
