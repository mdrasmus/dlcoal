/*=============================================================================

  Matt Rasmussen
  Copyright 2007-2011

  HKY sequence model

=============================================================================*/


#include "hky.h"
#include "common.h"
#include "seq.h"

namespace spidir {


/*=============================================================================
    From: Felsenstein. Inferring Phylogenies. p 202.

    NOTE: HKY is a special case of Tamura-Nei where 
        alpha_r / alpha_y = rho = pi_r / pi_y 

    NOTE: parameters are chosen such that 
        P(j != i | i, t=1, pi, ratio) = 1
        thus, time also corresponds to sub/site

    Definitions:
        i     = destination base
        j     = source base
        t     = time
        pi_i  = background/prior distribution of base i
        R     = Transition/Transversion ratio
        kappa = Transition/Transversion ratio (easier to specify)

    Parameterization:
        R    = (pi_t*pi_c + pi_a*pi_g) * kappa / (pi_y pi_r)
        beta = 1 / (2.0 * pi_r * pi_y * (1+R))
        rho  = pi_r / pi_y

        alpha_y = (pi_r * pi_y * R - pi_a*pi_g - pi_c*pi_t) / 
                  (2.0*(1+R)*(pi_y*pi_a*pi_g*rho + pi_r*pi_c*pi_t))
        alpha_r = rho * alpha_y    

    Convenience variables:
        if (dnatype[i] == DNA_PURINE) {
            alpha_i = alpha_r;
            pi_ry = pi_r;
        } else {
            alpha_i = alpha_y;
            pi_ry = pi_y;
        }
        int delta_ij =  int(i == j);
        int e_ij = int(dnatype[i] == dnatype[j]);
        
    Formula:
        prob(j | i, t, pi, R) = 
            exp(-(alpha_i + beta)*t) * delta_ij + 
            exp(-beta*t) * (1 - exp(-alpha_i*t)) * (pi_j*e_ij/pi_ry) + 
            (1 - exp(-beta*t)) * pi_j


  d/dt P(j|i,t,pi,R) = delta_ij(alpha_i + beta) exp(-(alpha_i+beta)t) +
                       (pi_j e_ij /pi_ry) 
		       [-beta exp(-beta t) + 
		        (alpha_i+beta)exp(-(alpha_i + beta)t)] +
		       pi_j beta exp(-beta t)


  d^2/dt^2 P(j|i,t,pi,R) = 
  \delta_{ij} (\alpha_i + \beta)^2 \exp(-(\alpha_i + \beta) t)  + 
  (\pi_j e_{ij}/\pi_{ry}) [\beta^2 \exp(-\beta t) -
          (\alpha_i + \beta)^2 \exp(-(\alpha_i + \beta) t) ] - 
   \pi_j \beta^2 \exp(-\beta t))

*/


HkyModelDeriv2::HkyModelDeriv2(const float *bgfreq, float ratio_kappa)
{
    // set background base frequencies
    for (int i=0; i<4; i++)
        pi[i] = bgfreq[i];

    pi_r = pi[DNA_A] + pi[DNA_G];
    pi_y = pi[DNA_C] + pi[DNA_T];
    rho = pi_r / pi_y;

    // convert the usual ratio definition (kappa) to Felsenstein's 
    // definition (R)
    ratio = (pi[DNA_T]*pi[DNA_C] + pi[DNA_A]*pi[DNA_G]) * ratio_kappa / \
        (pi_y * pi_r);
        
    // determine HKY parameters alpha_r, alpha_y, and beta
    b = 1.0 / (2.0 * pi_r * pi_y * (1.0+ratio));
    a_y = (pi_r*pi_y*ratio - 
           pi[DNA_A]*pi[DNA_G] - 
           pi[DNA_C]*pi[DNA_T]) / 
        (2.0*(1+ratio)*(pi_y*pi[DNA_A]*pi[DNA_G]*rho + 
                        pi_r*pi[DNA_C]*pi[DNA_T]));
    a_r = rho * a_y;
}

// transition probability P(j | i, t)
void HkyModelDeriv2::getMatrix(float t, float *matrix)
{
    for (int i=0; i<4; i++) {
        for (int j=0; j<4; j++) {
            // convenience variables
            // NOTE: it is ok to assign pi_ry, because it is only used when
            // dnatype[i] == dnatype[j]
            float a_i, pi_ry;
            switch (dnatype[i]) {
            case DNA_PURINE:
                a_i = a_r;
                pi_ry = pi_r;  
                break;
            case DNA_PRYMIDINE:
                a_i = a_y;
                pi_ry = pi_y;
                break;
            default:
                assert(0);
            }
            int delta_ij = int(i == j);
            int e_ij = int(dnatype[i] == dnatype[j]);
        
            // return transition probability
            float ab = a_i + b;
            float eabt = expf(-ab*t);
            float ebt = expf(-b*t);
        
            matrix[matind(4, i, j)] = delta_ij * ab*ab * eabt +
                (pi[j] * e_ij / pi_ry) * (b*b*ebt - ab*ab*eabt) -
                pi[j]*b*b*ebt;
        }
    }
}

//=============================================================================

HkyModelDeriv::HkyModelDeriv(const float *bgfreq, float kappa) :
    kappa(kappa)
{
    // set background base frequencies
    for (int i=0; i<4; i++)
        pi[i] = bgfreq[i];

    pi_r = pi[DNA_A] + pi[DNA_G];
    pi_y = pi[DNA_C] + pi[DNA_T];
    rho = pi_r / pi_y;

    // convert the usual ratio definition (kappa) to Felsenstein's 
    // definition (R)
    ratio = (pi[DNA_T]*pi[DNA_C] + pi[DNA_A]*pi[DNA_G]) * kappa / \
        (pi_y * pi_r);
        
    // determine HKY parameters alpha_r, alpha_y, and beta
    b = 1.0 / (2.0 * pi_r * pi_y * (1.0+ratio));
    a_y = (pi_r*pi_y*ratio - 
           pi[DNA_A]*pi[DNA_G] - 
           pi[DNA_C]*pi[DNA_T]) / 
        (2.0*(1+ratio)*(pi_y*pi[DNA_A]*pi[DNA_G]*rho + 
                        pi_r*pi[DNA_C]*pi[DNA_T]));
    a_r = rho * a_y;
}

// transition probability P(j | i, t)
void HkyModelDeriv::getMatrix(float t, float *matrix)
{
    for (int i=0; i<4; i++) {
        for (int j=0; j<4; j++) {
            // convenience variables
            // NOTE: it is ok to assign pi_ry, because it is only used when
            // dnatype[i] == dnatype[j]
            float a_i, pi_ry;
            switch (dnatype[i]) {
            case DNA_PURINE:
                a_i = a_r;
                pi_ry = pi_r;  
                break;
            case DNA_PRYMIDINE:
                a_i = a_y;
                pi_ry = pi_y;
                break;
            default:
                assert(0);
            }
            int delta_ij = int(i == j);
            int e_ij = int(dnatype[i] == dnatype[j]);
        
            // return transition probability
            float ab = a_i + b;
            float eabt = expf(-ab*t);
            float ebt = expf(-b*t);
        
            matrix[matind(4, i, j)] = - delta_ij * ab * eabt +
                (pi[j] * e_ij / pi_ry) * (- b * ebt + ab * eabt) +
                pi[j] * b * ebt;
        }
    }
}

HkyModelDeriv::Deriv *HkyModelDeriv::deriv()
{
    return new Deriv(pi, kappa);
}


//=============================================================================

HkyModel::HkyModel(const float *bgfreq, float kappa) :
    kappa(kappa)
{
    // set background base frequencies
    for (int i=0; i<4; i++)
        pi[i] = bgfreq[i];        

    pi_r = pi[DNA_A] + pi[DNA_G];
    pi_y = pi[DNA_C] + pi[DNA_T];
    rho = pi_r / pi_y;

    // convert the usual ratio definition (kappa) to Felsenstein's 
    // definition (R)
    ratio = (pi[DNA_T]*pi[DNA_C] + pi[DNA_A]*pi[DNA_G]) * kappa / \
        (pi_y * pi_r);
        
    // determine HKY parameters alpha_r, alpha_y, and beta
    b = 1.0 / (2.0 * pi_r * pi_y * (1.0+ratio));
    a_y = (pi_r*pi_y*ratio - 
           pi[DNA_A]*pi[DNA_G] - 
           pi[DNA_C]*pi[DNA_T]) / 
        (2.0*(1+ratio)*(pi_y*pi[DNA_A]*pi[DNA_G]*rho + 
                        pi_r*pi[DNA_C]*pi[DNA_T]));
    a_r = rho * a_y;
}


/* For Jukes Cantor we have: 

    pi_r = pi[DNA_A] + pi[DNA_G] = .5;
    pi_y = pi[DNA_C] + pi[DNA_T] = .5;
    rho = pi_r / pi_y = 1.0;

    // convert the usual ratio definition (kappa) to Felsenstein's 
    // definition (R)
    ratio = (pi[DNA_T]*pi[DNA_C] + pi[DNA_A]*pi[DNA_G]) * kappa / \
        (pi_y * pi_r) =
        (.25*.25 + .25*.25) * 1.0 / (.5 * .5) =
        (1/16 + 1/16) / .25 = 1/8 * 4 = .5;
        
    // determine HKY parameters alpha_r, alpha_y, and beta
    b = 1.0 / (2.0 * pi_r * pi_y * (1.0+ratio))
      = 1.0 / (2.0 * .5 * .5 * (1.0+.5))
      = 1.0 / (3/4) = 4/3

    a_y = (pi_r*pi_y*ratio - 
           pi[DNA_A]*pi[DNA_G] - 
           pi[DNA_C]*pi[DNA_T]) / 
        (2.0*(1+ratio)*(pi_y*pi[DNA_A]*pi[DNA_G]*rho + 
                        pi_r*pi[DNA_C]*pi[DNA_T]))
        = (.5*.5*.5 - 
           .25*.25 - 
           .25*.25) / 
           (2.0*(1+.5)*(.5*.25*.25*1 + .5*.25*.25))
        = (1/8 - 1/16 - 1/16) / (2.0*(1+.5)*(.5*.25*.25*1 + .5*.25*.25))
        = 0;

    a_r = rho * a_y = 0;

*/

// transition probability P(j | i, t)
void HkyModel::getMatrix(float t, float *matrix)
{
    for (int i=0; i<4; i++) {
        for (int j=0; j<4; j++) {
            // convenience variables
            // NOTE: it is ok to assign pi_ry, because it is only used when
            // dnatype[i] == dnatype[j]
            double a_i, pi_ry;
            switch (dnatype[i]) {
            case DNA_PURINE:
                a_i = a_r;
                pi_ry = pi_r;  
                break;
            case DNA_PRYMIDINE:
                a_i = a_y;
                pi_ry = pi_y;
                break;
            default:
                assert(0);
            }
            int delta_ij = int(i == j);
            int e_ij = int(dnatype[i] == dnatype[j]);
        
            // return transition probability
            double ait = exp(-a_i*t);
            double ebt = exp(-b*t);

            matrix[matind(4,i,j)] = ait*ebt * delta_ij + 
                ebt * (1.0 - ait) * (pi[j]*e_ij/pi_ry) + 
                (1.0 - ebt) * pi[j];
        }
    }
}


HkyModel::Deriv *HkyModel::deriv()
{
    return new Deriv(pi, kappa);
}



extern "C" {

void makeHkyMatrix(const float *bgfreq, float ratio, float t, float *matrix)
{
    HkyModel model(bgfreq, ratio);
    model.getMatrix(t, matrix);
}

void makeHkyDerivMatrix(const float *bgfreq, float ratio, float t, float *matrix)
{
    HkyModelDeriv model(bgfreq, ratio);
    model.getMatrix(t, matrix);
}

void makeHkyDeriv2Matrix(const float *bgfreq, float ratio, float t, float *matrix)
{
    HkyModelDeriv2 model(bgfreq, ratio);
    model.getMatrix(t, matrix);
}

} // namespace spidir

} // extern "C"
