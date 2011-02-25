/*=============================================================================

  Matt Rasmussen
  Copyright 2007-2011

  Common miscellaneous functions

=============================================================================*/

// c++ headers
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>


// spidir headers
#include "common.h"

namespace spidir {


//=============================================================================
// Math

int choose(int n, int k)
{       
    return int(fchoose(n, k) + .5);
}


double fchoose(int n, int k)
{   
    if (n < 0 || k < 0 || k > n)
        return 0;
    
    // optimization for speed
    if (k > n/2)
        k = n - k;
    
    double t = 1.0;
    double m = n;
    for (double i=1; i<=k; i++)
        t *= (m - i + 1) / i;
    return t;
}


//=============================================================================
// distributions

extern "C" {

// probability density distribution of the Poisson
float poisson(int x, float lambda)
{
    if (x < 0 || lambda <= 0)
        return 0.0;
    
    float a = 0.0;
    for (float i=1.0; i<x+1; i+=1.0)
        a += log(lambda / i);
    return exp(-lambda + a);
}


// Normal distribution.
//
// mu is the mean, and sigma is the standard deviation.
//
float normalvariate(float mu, float sigma)
{
    // Uses Kinderman and Monahan method. Reference: Kinderman,
    // A.J. and Monahan, J.F., "Computer generation of random
    // variables using the ratio of uniform deviates", ACM Trans
    // Math Software, 3, (1977), pp257-260.

    const static float NV_MAGICCONST = 4 * exp(-0.5)/sqrt(2.0);
    float u1, u2, z, zz;

    do {
        u1 = frand();
        u2 = 1.0 - frand();
        z = NV_MAGICCONST*(u1-0.5)/u2;
        zz = z*z/4.0;
    } while (zz > -log(u2));
    
    return mu + z*sigma;
}


} // extern "C"


//=============================================================================
// sorting

// Invert a permutation
void invertPerm(int *perm, int *inv, int size)
{
    for (int i=0; i<size; i++)
        inv[perm[i]] = i;
}



//=============================================================================
// input/output

void printIntArray(int *array, int size)
{
    for (int i=0; i<size; i++)
        printf("%d ", array[i]);
    printf("\n");
}

void printFloatArray(float *array, int size)
{
    for (int i=0; i<size; i++)
        printf("%f ", array[i]);
    printf("\n");
}



} // namespace spidir
