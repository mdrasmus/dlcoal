/*=============================================================================

  Matt Rasmussen
  Copyright 2007-2011

  Gamma distribution

=============================================================================*/


#ifndef SPIDIR_GAMMA_H
#define SPIDIR_GAMMA_H


// headers c++ 
#include <math.h>
#include <assert.h>

namespace spidir {

extern "C" {

double gammln(double xx);
double gamm(double x);
double gammalog(double x, double a, double b);
double gammaPdf(double x, double a, double b);
double invgammaPdf(double x, double a, double b);
double loginvgammaPdf(double x, double a, double b);
double invgammaDerivA(double x, double a, double b);
double invgammaDerivB(double x, double a, double b);
double invgammaDerivG(double x, double g);
double invgammaDerivG2(double x, double g);
//float invgammaCdf(float x, float a, float b);
//double quantInvgamma(double p, double a, double b);

// Derivative of Gamma distribution with respect to x
double gammaDerivX(double x, double a, double b);
    
// Derivative of Gamma distribution with respect to a
double gammaDerivA(double x, double a, double b);

// Derivative of Gamma distribution with respect to b
double gammaDerivB(double x, double a, double b);

// Derivative of Gamma distribution with respect to nu (its variance)
double gammaDerivV(double x, double v);

// Second Derivative of Gamma distribution with respect to nu (its variance)
double gammaDerivV2(double x, double v);

float gammavariate(float alpha, float beta);

inline float invgammavariate(float alpha, float beta)
{ return 1.0 / gammavariate(alpha, beta); }

double gammaSumPdf(double y, int n, float *alpha, float *beta, 
		   float tol);


} // extern "C"

} // namespace spidir

#endif // SPIDIR_GAMMA_H
