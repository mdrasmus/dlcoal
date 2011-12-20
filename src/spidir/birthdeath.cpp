/*=============================================================================

  Matt Rasmussen
  Copyright 2010-2011

  Birth-death models

=============================================================================*/

// c/c++ includes
#include <math.h>

// spidir includes
#include "birthdeath.h"
#include "common.h"


namespace spidir
{

extern "C" {


// returns the probability of 1 gene giving rise to ngenes after time 'time'
double birthDeathCount(int ngenes, float time, 
                       float birthRate, float deathRate)
{
    const double l = birthRate;
    const double u = deathRate;

    if (birthRate == deathRate) {
        const double ut = time / (1.0 / birthRate + time);
        if (ngenes == 0)
            return ut;
        return (1.0 - ut)*(1.0 - ut) * ipow(ut, ngenes-1);
    }

    const double r = l - u;
    const double a = u / l;

    const double ut = (1.0 - exp(-r*time)) / (1.0 - a * exp(-r*time));
    const double p0 = a*ut;

    if (ngenes == 0)
        return p0;
    
    // p(0, t) = ut
    // p(1, t) = ...
    // (1.0 - p0)*(1.0 - ut) * ut^{ngenes-1}
    return (1.0 - p0)*(1.0 - ut) * ipow(ut, ngenes-1);
}


// returns the probability of 'start' genes giving rise to 'end' genes after 
// time 'time'
// slower more stable computation
double birthDeathCountsSlow(int start, int end, float time, 
                            float birth, float death)
{
    const double l = birth;
    const double u = death;
    const double r = l - u;
    const double a = u / l;

    const double ertime = exp(-r*time);
    const double ut = (1.0 - ertime) / (1.0 - a * ertime);
    const double p0 = a*ut;
    
    // all 'start' genes die out
    if (end == 0) {
        return ipow(p0, start);
    }
    
    const int iter = (start < end) ? start : end;
    double p = 0.0;

    for (int j=0; j<=iter; j++) {
        p += fchoose(start, j) *
             fchoose(start + end - j - 1, start - 1) *
             ipow(p0, start-j) *
             ipow(ut, end-j) *
             ipow(1 - p0 - ut, j);
    }
    
    // do not allow invalid values to propogate
    if (isnan(p) || isinf(p) || p > 1.0) {
        printf("p=%e genes=(%d, %d) b=%f d=%f t=%f\n", p, start, end,
               birth, death, time);
        fflush(stdout);
        assert(0);
    }
    return p;
}


// returns the probability of 'start' genes giving rise to 'end' genes after 
// time 'time'
// much faster computation than birthDeathCountsSlow
double birthDeathCountsLog(int start, int end, float time, 
                           float birth, float death)
{
    if (start == 0) {
        if (end == 0)
            return 0.0;
        else
            return -INFINITY;
    }

    const double ertime = exp((birth-death)*time);
    const double tmp = (ertime-1.0) / (birth*ertime - death);
    const double a = death * tmp;
    const double b = birth * tmp;
    
    // all 'start' genes die out
    if (end == 0) {
        //return ipow(a, start);
        return log(a) * start;
    }
    
    // log scale
    // compute base case
    double f = log(a) * start + log(b) * end;
    int sf = 1;
    if (start > 1)
        f += log(start + end - 1);
    for (int k=2; k<start; k++)
        f += log((start + end - k) / double(k));
    
    double p = f;
    int sp = 1;
    double x = start;
    double y = end;
    double z = start + end - 1;
    const double oneab = (1.0 - a - b) / (a * b);
    const int iter = (start < end) ? start : end;
    for (int j=1; j<=iter; j++) {
        double c = oneab * x * y / (j * z);        
        sf *= (2 * int(c > 0.0) - 1);
        f += log(fabs(c));
        logadd_sign(sp, p, sf, f, &sp, &p);
        x--;
        y--;
        z--;
        //printf("p=%f, s=%f, f=%f, c=%f\n", p, s, f, c);
    }

    // round to zero prob
    if (sp < 0)
        return -INFINITY;
    
    return p;
}


// returns the probability of 'start' genes giving rise to 'end' genes after 
// time 'time'
// much faster computation than birthDeathCountsSlow
double birthDeathCounts(int start, int end, float time, 
                           float birth, float death)
{
    if (start == 0) {
        if (end == 0)
            return 1.0;
        else
            return 0.0;
    }

    const double ertime = exp((birth-death)*time);
    const double tmp = (ertime-1.0) / (birth*ertime - death);
    const double a = death * tmp;
    const double b = birth * tmp;
    
    // all 'start' genes die out
    if (end == 0) {
        return ipow(a, start);
    }
    
    // compute base case
    double f = ipow(a, start) * ipow(b, end);
    if (start > 1)
        f *= (start + end - 1);
    for (int k=2; k<start; k++)
        f *= (start + end - k) / double(k);


    double p = f;
    double x = start;
    double y = end;
    double z = start + end - 1;
    const double oneab = 1.0 - a - b;
    const int iter = (start < end) ? start : end;
    for (int j=1; j<=iter; j++) {
        f *= (oneab * x * y / (j * a * b * z));
        p += f;
        x--;
        y--;
        z--;
    }

    if (p < 0.0)
        p = 0.0;
    

    if (p > 1.0) {
        // resort to a slower more stable function
        return birthDeathCountsSlow(start, end, time, birth, death);
    }

    return p;
}


// Probability of no birth from 'n' lineages starting at time 0, 
// evolving until time 'T' with 'birth' and 'death' rates
// for a reconstructed process.
double probNoBirth(int n, float T, float birth, float death) 
{
    if (birth == 0.0)
        return 1.0;
    else if (birth == death) {
        return 1.0 / ipow(1.0 + birth * T, n);
    }

    const double r = birth - death;
    const double exprt = exp(-r * T);

    return pow(1.0 - (birth*(1.0 - exprt)) / (birth - death * exprt), n);
}





//=============================================================================
// sampling


//  Probability density for for next birth at time 't' given
//  'n' lineages starting at time 0, evolving until time 'T' with a
//  'birth' and 'death' rates for a reconstructed process.
double birthWaitTime(float t, int n, float T, float birth, float death)
{    
    if (birth == death) {
        const double t2 = t - T;
        const double nl = 1.0 / birth;
        return birth * n * ipow(-nl+t2, n) / ipow(-nl-T, n) / (1.0-birth*t2);
    }

    const double r = birth - death;
    const double a = death / birth;

    return n * r * exp(-n*r*t) * \
           pow(1.0 - a * exp(-r * (T - t)), n-1) / \
	   pow(1.0 - a * exp(-r * T), n);
}

//  Probability density for for next birth at time 't' given
//  'n'=1 lineages starting at time 0, evolving until time 'T' with a
//  'birth' and 'death' rates for a reconstructed process.
double birthWaitTime1(float t, float T, float birth, float death)
{
    // special case
    if (birth == death) {
        const double t2 = t - T;
        const double nl = 1.0 / birth;
        return birth * (-nl+t2) / (-nl-T) / (1.0-birth*t2);
    }

    const double r = birth - death;
    const double a = death / birth;

    return r * exp(-r*t) / (1.0 - a * exp(-r * T));
}


// numerator for birthWaitTime
double birthWaitTimeNumer(float t, int n, float T, 
                              float birth, float death,
                              double denom)
{    
    const double r = birth - death;
    const double a = death / birth;

    return n * r * exp(-n*r*t) * 
           pow(1.0 - a * exp(-r * (T - t)), n-1) / denom;
}

// denominator for birthWaitTime
double birthWaitTimeDenom(int n, float T, float birth, float death)
{    
    const double r = birth - death;
    const double a = death / birth;

    return pow(1.0 - a * exp(-r * T), n);
}



//  Probability density for for next birth at time 't' given
//  'n'=1 lineages starting at time 0, evolving until time 'T' with a
//  'birth' and 'death' rates for a reconstructed process.
//   The denominator 'denom' must be precomputed
double birthWaitTimeNumer1(float t, float T, float birth, float death, 
                         float denom)
{
    const double r = birth - death;
    return r * exp(-r*t) / denom;
}

//  Computes the denominator for birthWaitTime1
double birthWaitTimeDenom1(float T, float birth, float death)
{
    const double r = birth - death;
    const double a = death / birth;
    return 1.0 - a * exp(-r * T);
}




// Sample the next birth event from a reconstructed birthdeath process.
// Let there be 'n' lineages at time 0 that evolve until time 'T' with
// 'birth' and 'death' rates.
// Conditioned that a birth will occur
double sampleBirthWaitTime(int n, float T, float birth, float death)
{
    
    // TODO: could make this more efficient

    if (birth == death) {
        double start_y = birthWaitTime(0, n, T, birth, death);
        double end_y = birthWaitTime(T, n, T, birth, death);
        double M = max(start_y, end_y);
    
        while (true) {
            double t = frand(T);
            double f = birthWaitTime(t, n, T, birth, death);
            
            if (frand() <= f / M)
                return t;
        }
    } else {
        // uses rejection sampling
        double denom = birthWaitTimeDenom(n, T, birth, death);
        double start_y = birthWaitTimeNumer(0, n, T, birth, death, denom);
        double end_y = birthWaitTimeNumer(T, n, T, birth, death, denom);
        double M = max(start_y, end_y);
    
        while (true) {
            double t = frand(T);
            double f = birthWaitTimeNumer(t, n, T, birth, death, denom);

            if (frand() <= f / M)
                return t;
        }
    }
}



// Sample the next birth event from a reconstructed birthdeath process.
// Let there be 'n'=1 lineages at time 0 that evolve until time 'T' with
// 'birth' and 'death' rates.
// Conditioned that a birth will occur
double sampleBirthWaitTime1(float T, float birth, float death)
{    
    // TODO: could make this much more efficient

    if (birth == death) {
        // uses rejection sampling
        double start_y = birthWaitTime1(0, T, birth, death);
        double end_y = birthWaitTime1(T, T, birth, death);
        double M = max(start_y, end_y);
    
        while (true) {
            double t = frand(T);
            double f = birthWaitTime1(t, T, birth, death);
            
            if (frand() <= f / M)
                return t;
        }

    } else {
        // uses rejection sampling
        double denom = birthWaitTimeDenom1(T, birth, death);
        double start_y = birthWaitTimeNumer1(0, T, birth, death, denom);
        double end_y = birthWaitTimeNumer1(T, T, birth, death, denom);
        double M = max(start_y, end_y);
    
        while (true) {
            double t = frand(T);
            double f = birthWaitTimeNumer1(t, T, birth, death, denom);

            if (frand() <= f / M)
                return t;
        }
    }
}





}

} // namespace spidir
