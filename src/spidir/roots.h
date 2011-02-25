/*=============================================================================

  Matt Rasmussen
  Copyright 2007-2011

  Simple root finding methods

=============================================================================*/


#ifndef SPIDIR_ROOTS_H
#define SPIDIR_ROOTS_H

// headers c++ 
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string>
#include <vector>
#include <string.h>


using namespace std;

namespace spidir {

// Find a root of a function func(x) using the secant method
// x0 and x1 are initial estimates of the root
template <class Func>
float secantRoot(Func &f, float x0, float x1, int maxiter, 
                 float minx=.000001, float esp=.002)
{
    float f0 = f(x0);
    for (int i=0; i<maxiter; i++) {
        if (fabs((x1 - x0)*2.0 / (x0+x1)) < esp)
            return x0;
        float f1 = f(x1);
        float x2 = x1 - (x1 - x0) * f1 / (f1 - f0);
        
        x0 = x1;
        x1 = (x2 > minx) ? x2 : minx;
        f0 = f1;
    }

    return x1;
}


// Find a root of a function func(x) using the bisection method
// This is less efficient but is more robust than Newton's or Secant
// x0 and x1 are initial estimates of the root
template <class Func>
float bisectRoot(Func &f, float x0, float x1, const float err=.001)
{
    float f0 = f(x0);
    float f1 = f(x1);
    
    while (fabs(x1 - x0) > 2 * err) {
        //printf("in:  %f %f; %f %f\n", x0, x1, f0, f1);

        float x2 = (x0 + x1) / 2.0;
        float f2 = f(x2);
        
        if (f0 * f2 > 0) {
            x0 = x2;
            f0 = f2;
        } else {
            x1 = x2;
            f1 = f2;
        }
    }

    return (x0 + x1) / 2.0;
}


} // namespace spidir

#endif // SPIDIR_ROOTS_H
