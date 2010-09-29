/*=============================================================================

    SPIDIR - SPecies Informed DIstanced-based Reconstruction
    
    Matt Rasmussen
    Wed Jun 13 22:09:24 EDT 2007


=============================================================================*/

// c++ headers
#include <assert.h>
#include <math.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

// gsl
#include <gsl/gsl_sf.h>


// spidir headers
#include "common.h"

namespace spidir {


// stream for logging
static FILE *g_logstream = stderr;
static int g_loglevel = LOG_QUIET;


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


/* natural log of the gamma function
   gammln as implemented in the
 * first edition of Numerical Recipes in C */
double gammln(double xx)
{
    double x, tmp, ser;
    const static double cof[6]={76.18009172947146,    -86.50532032941677,
                                24.01409824083091,    -1.231739572450155,
                                0.1208650973866179e-2,-0.5395239384953e-5};
    int j;

    x=xx-1.0;
    tmp=x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser=1.000000000190015;
    for (j=0;j<=5;j++) {
        x += 1.0;
        ser += cof[j]/x;
    }
    return -tmp+log(2.5066282746310005*ser);
}

//
// Lanczos approximation to the gamma function. 
// 
// found on http://www.rskey.org/gamma.htm   
//
double gamm(double x) 
{
    double ret = (1.000000000190015 + 
                 76.18009172947146 / (x + 1) +  
                 -86.50532032941677 / (x + 2) + 
                 24.01409824083091 / (x + 3) +  
                 -1.231739572450155 / (x + 4) + 
                 1.208650973866179e-3 / (x + 5) + 
                 -5.395239384953e-6 / (x + 6));
    
    return ret * sqrt(2*M_PI)/x * pow(x + 5.5, x+.5) * exp(-x-5.5);
}




// natural log of the gamma distribution PDF
double gammalog(double x, double a, double b)
{
    if (x <= 0 || a <= 0 || b <= 0)
        return 0.0;
    else
        return -x * b + (a - 1.0) * log(x) + a * log(b) - gammln(a);
}

// the gamma distribution PDF
double gammaPdf(double x, double a, double b)
{
    if (x <= 0 || a <= 0 || b <= 0)
        return 0.0;
    else
        return exp(-x * b) * pow(x, a - 1.0) * pow(b, a) / gamm(a);
}


// Inverse gamma distribution PDF
double invgammaPdf(double x, double a, double b)
{
    return exp(loginvgammaPdf(x, a, b));
    //return pow(b, a) / gamm(a) * pow(1/x, a+1) * exp(-b/x);
}


// Log Inverse gamma distribution PDF
double loginvgammaPdf(double x, double a, double b)
{
    return a*log(b) - gammln(a) + (a+1)*log(1/x) - b/x;
}


double invgammaDerivA(double x, double a, double b)
{
    return gammaDerivA(1/x, a, b) / x / x;
}


double invgammaDerivB(double x, double a, double b)
{
    return gammaDerivB(1/x, a, b) / x / x;
}


double invgammaDerivG(double x, double d)
{
    double logd = log(d);
    double logix = log(1/x);
    double p = - gammln(1+d) + d*logd - d/x + (3+d)*logix;
    return exp(p) * (d * (-1 + x) + x + d*x*(logd + logix) - 
                     d*x*gsl_sf_psi_n(0, 1+d));
    /*
    return (1/gamm(1+d) * pow(d,d) * exp(-d/x) * pow(1/x, 3+d)) * 
            (d * (-1 + x) + x + d*x*(log(d) + log(1/x)) - 
             d*x*gsl_sf_psi_n(0, 1+d));
    */
}

double invgammaDerivG2(double x, double d)
{

    double logd = log(d);
    double logix = log(1/x);
    double dxlogdix = d*x*(logd + logix);
    double psi1 = gsl_sf_psi_n(0, 1 + d);
    double psi2 = gsl_sf_psi_n(1, 1 + d);

    double p = - gammln(1 + d) + d*log(d) -d/x + (4+d)*log(1/x);

    return exp(p) * (d - 2*(1 + d)*x + (3 + d)*x*x + 
         x*(logd + logix)*(2*(d*(-1 + x) + x) + dxlogdix) + 
         x*(-2*(d*(-1 + x) + x + dxlogdix)*
            psi1 + d*x*psi1*psi1 - d*x*psi2));
    return exp(p);
    
    /*
    return (1/gamm(1 + d))*pow(d,d)*exp(-d/x)*pow(1/x, 4 + d) *
        (d - 2*(1 + d)*x + (3 + d)*x*x + 
         x*(log(d) + log(1/x))*(2*(d*(-1 + x) + x) + 
                                d*x*(log(d) + log(1/x))) + 
         x*(-2*(d*(-1 + x) + x + d*x*(log(d) + log(1/x)))*
            gsl_sf_psi_n(0, 1 + d) + 
            d*x*pow(gsl_sf_psi_n(0, 1 + d), 2) - 
            d*x*gsl_sf_psi_n(1, 1 + d)));
    */
}


/*
float invgammaCdf(float x, float a, float b)
{
    return incompletegammac(a, b / x) / gamm(a);
}

double quantInvgamma(double p, double a, double b)
{
    return b / invincompletegammac(a, gamm(a) * p);
}
*/

// Derivative of Gamma distribution with respect to x
double gammaDerivX(double x, double a, double b)
{ 
    return pow(b, a) / gamm(a) * exp(-b*x) * pow(x, a-2) * (a - b*x - 1);
}
    
// Derivative of Gamma distribution with respect to a
double gammaDerivA(double x, double a, double b)
{
    return exp(-b*x) * pow(x,a-1) * pow(b, a) / gamm(a) * 
        (log(b) + log(x) - gsl_sf_psi_n(0, a));
}

// Derivative of Gamma distribution with respect to b
double gammaDerivB(double x, double a, double b)
{
    return pow(x, a-1) / gamm(a) * exp(-b*x) * pow(b, a-1) * (a - x*b);
}


// Derivative of Gamma distribution with respect to nu (its variance)
double gammaDerivV(double x, double v)
{
    float q = 1.0 / v;
    
    double lnq = log(q);
    double lnx = log(x);
    double w = lnq*(q+2) + lnx*(q-1) - x * q - gammln(q);
    double z = 1 - x + lnq + lnx - gsl_sf_psi_n(0, q);

    if (z > 0)
	return -exp(w + log(z));
    else
	return exp(w + log(-z));
    

    //return - pow(q, q+2) * pow(x, q-1) * exp(-x * q) / gamm(q) * 
    //	(1 - x + log(q) + log(x) - gsl_sf_psi_n(0, q));
}

// Second Derivative of Gamma distribution with respect to nu (its variance)
double gammaDerivV2(double x, double v)
{
    double q = 1.0 / v;
    double lnx = log(x);
    double lnq = log(q);
    double A = 1.0 + v - x + lnx;
    double psi0 = gsl_sf_psi_n(0, q);
    double psi1 = gsl_sf_psi_n(1, q);

    double w = -x*q + lnq*(q+4) + lnx*(q-1) - gammln(q);
    double z = 1 + 3*v - 2*x - 2*v*x + x*x + lnq*lnq - lnx*lnx +
	2*A*lnx + 2*A*lnq -
	2*(A+lnq)*psi0 + psi0*psi0 - psi1;
    double y;

    if (z > 0.0)
	y = exp(w + log(z));
    else
	y = -exp(w + log(-z));
    
    //printf("x=%f, v=%f, y=%f %f; %f %f\n", x, v, y, gammln(q), w, z);

    //assert(y < INFINITY);
    //assert(y > -INFINITY);
    //assert(!isnan(y));

    return y;

    /*
    return exp(-x*q) * pow(q,q+4)*pow(x,q-1)/gamm(q) *
	(1 + 3*v - 2*x - 2*v*x + x*x + lnq*lnq - lnx*lnx +
	 2*A*lnx + 2*A*lnq -
	 2*(A+lnq)*psi0 + psi0*psi0 - psi1);
    */
}



// PDF of a sum of n gamma variables
double gammaSumPdf(double y, int n, float *alpha, float *beta, 
		   float tol)
{
    const int nterms = 100; // maximum number of terms to compute

    // y must be positive
    if (y <= 0.0)
        return 0.0;

    // handle case where a single gamma is present
    if (n == 1)
	return gammaPdf(y, alpha[0], beta[0]);

    double b1 = 1.0 / beta[0];
    double beta2[n];
    for (int i=0; i<n; i++) {

	// convert betas into their reciprocals
	beta2[i] = 1.0 / beta[i];

	// find min beta (b1)
	if (beta2[i] < b1)
	    b1 = beta2[i];
    }

    // compute terms C and p
    double C = 1.0;
    double p = 0.0;
    for (int i=0; i<n; i++) {
	C *= pow(b1 / beta2[i], alpha[i]);
	p += alpha[i];	
    }


    double prob = 0.0, prob2 = 0.0;
    double gammas[nterms]; // originally 1-indexed
    double deltas[nterms]; // zero-indexed
    deltas[0] = 1.0;
    for (int i=0; i<nterms-1; i++) {
        int k = i + 1;
        gammas[i] = 0.0;
	for (int j=0; j<n; j++)
            gammas[i] += alpha[j] * pow(1.0 - b1 / beta2[j], k) / k;
	
        k = i;
        deltas[k+1] = 0.0;
	for (int j=1; j<k+2; j++)
	    deltas[k+1] += j * gammas[j-1] * deltas[k+1-j];
	deltas[k+1] /= (k+1);

        prob2 = deltas[k] * pow(y, p+k-1) * exp(-y/b1) /
	    (gamm(p+k) * pow(b1, p+k));

        prob += prob2;

	// test for convergence
	if (prob2 / prob < tol)
	    break;
    }

    return C * prob;
}




double negbinomPdf(int k, double r, double p)
{
    return gsl_ran_negative_binomial_pdf(k, p, r);
}



// Log Derivative of Negative Binomial distribution with respect to r
double negbinomDerivR(int k, double r, double p)
{
    // (1-p)^k p^r G(k+r) / (k! G(r)) ( psi^(0)(k+r) - psi^(0)(r) + log(p) )
    double lnp = log(p);
    double A = k * log(1-p) + r*lnp + gammln(k+r) - gammln(k+1) - gammln(r);
    double B = gsl_sf_psi_n(0, k+r) - gsl_sf_psi_n(0, r) + lnp;
    return exp(A) * B;
}


// Derivative of Negative Binomial distribution with respect to p
double negbinomDerivP(int k, double r, double p)
{
    // (1-p)^k p^r G(k+r) / (k! G(r)) (k/(p-1) + r/p)
    double lnp = log(p);
    double A = k * log(1-p) + r*lnp + gammln(k+r) - gammln(k+1) - gammln(r);
    double B = k/(p-1) + r/p;
    return exp(A) * B;
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


float gammavariate(float alpha, float beta)
{
    const static float LOG4 = 1.3862943611198906;
    const static float SG_MAGICCONST = 1.0 + log(4.5);
    
    assert(alpha > 0.0 && beta > 0.0);
    
    // convert beta
    beta = 1.0 / beta;
    
    if (alpha > 1.0) {
        // Uses R.C.H. Cheng, "The generation of Gamma
        // variables with non-integral shape parameters",
        // Applied Statistics, (1977), 26, No. 1, p71-74

        float ainv = sqrt(2.0 * alpha - 1.0);
        float bbb = alpha - LOG4;
        float ccc = alpha + ainv;

        while (1) {
            float u1 = frand();
            if (u1 < 1e-7 || u1 > .9999999)
                continue;
            float u2 = 1.0 - frand();
            float v = log(u1 / (1.0-u1)) / ainv;
            float x = alpha * exp(v);
            float z = u1*u1*u2;
            float r = bbb+ccc*v-x;
            if (r + SG_MAGICCONST - 4.5*z >= 0.0 || r >= log(z))
                return x * beta;
        }
        
    } else if (alpha == 1.0) {
        // expovariate(1)
        float u = 0;
        while (u <= 1e-7)
            u = frand();
        return -log(u) * beta;
        
    } else { 
        // alpha in (0, 1)
        // Uses ALGORITHM GS of Statistical Computing - Kennedy & Gentle
        float x;

        while (1) {
            float u = frand();
            float b = (M_E + alpha)/M_E;
            float p = b*u;
            if (p <= 1.0)
                x = pow(p, (1.0/alpha));
            else
                x = -log((b-p)/alpha);
            float u1 = frand();
            if (p > 1.0)
                if (u1 <= pow(x, (alpha - 1.0)))
                    break;
                else if (u1 <= exp(-x))
                    break;
        }
        return x * beta;
    }
}




} // extern "C"

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


bool inChars(char c, const char *chars)
{
   if (!chars)
      return false;
   for (;*chars; chars++)
      if (c == *chars) return true;
   return false;
}


bool chomp(char *str)
{
   int len = strlen(str);
   if (str[len-1] == '\n') {
      str[len-1] = '\0';
      return true;
   } else
      return false;
}


vector<string> split(const char *str, const char *delim, bool multiDelim)
{
    vector<string> tokens;   
    int i=0, j=0;
   
    while (str[i]) {
        // walk to end of next token
        for (; str[j] && !inChars(str[j], delim); j++);
        
        if (i == j)
            break;

        // save token
        tokens.push_back(string(&str[i], j-i));
        
        if (!str[j])
            break;
        j++;
        i = j;
    }
    
    return tokens;
}


string trim(const char *word)
{
    char buf[101];
    if (sscanf(word, "%100s", buf) == 1)
        return string(buf);
    else
        return string("");
}


//=============================================================================
// Errors and Logging

void printError(const char *fmt, ...)
{
    va_list ap;   
    va_start(ap, fmt);
   
    fprintf(stderr, "error: ");
    vfprintf(stderr, fmt, ap);
    fprintf(stderr, "\n");
   
    va_end(ap);
}


void printLog(int level, const char *fmt, ...)
{
    if (level <= g_loglevel) {
        va_list ap;   
        va_start(ap, fmt);
        vfprintf(g_logstream, fmt, ap);
        va_end(ap);
	fflush(g_logstream);
    }
}


bool openLogFile(const char *filename)
{
    FILE *stream = fopen(filename, "w");
    
    if (stream != NULL) {
        openLogFile(stream);
        return true;
    } else {
        return false;
    }
}

void openLogFile(FILE *stream)
{
    g_logstream = stream;
}


void setLogLevel(int level)
{
    g_loglevel = level;
}

bool isLogLevel(int level)
{
    return level <= g_loglevel;
}

void closeLogFile()
{
    fclose(g_logstream);
}

FILE *getLogFile()
{
    return g_logstream;
}



} // namespace spidir
