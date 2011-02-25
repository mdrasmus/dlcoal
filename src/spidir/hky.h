/*=============================================================================

  Matt Rasmussen
  Copyright 2007-2011

  HKY sequence model

=============================================================================*/


#ifndef SPIDIR_HKY_H
#define SPIDIR_HKY_H

#include "assert.h"
#include "math.h"

namespace spidir {


class HkyModelDeriv2
{
public:
    HkyModelDeriv2(const float *bgfreq, float ratio_kappa);
    void getMatrix(float t, float *matrix);
    
    // parameters
    float ratio;
    float pi[4];
    float pi_r;
    float pi_y;
    float rho;
    float b;
    float a_y;
    float a_r;
};


class HkyModelDeriv
{
public:
    HkyModelDeriv(const float *bgfreq, float kappa);
    void getMatrix(float t, float *matrix);
    typedef HkyModelDeriv2 Deriv;
    Deriv *deriv();

    // parameters
    float kappa;
    float ratio;
    float pi[4];
    float pi_r;
    float pi_y;
    float rho;
    float b;
    float a_y;
    float a_r;
};


class HkyModel
{
public:
    HkyModel(const float *bgfreq, float kappa);
    void getMatrix(float t, float *matrix);
    typedef HkyModelDeriv Deriv;
    Deriv *deriv();

    // parameters
    float kappa;
    float ratio;
    float pi[4];
    float pi_r;
    float pi_y;
    float rho;
    float b;
    float a_y;
    float a_r;
};

extern "C" {

void makeHkyMatrix(const float *bgfreq, float ratio, float t, float *matrix);
void makeHkyDerivMatrix(const float *bgfreq, float ratio, float t, float *matrix);
void makeHkyDeriv2Matrix(const float *bgfreq, float ratio, float t, float *matrix);

} // extern "C"

} // namespace spidir

#endif // SPIDIR_HKY_H

