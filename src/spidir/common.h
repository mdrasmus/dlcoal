#ifndef SPIDIR_COMMON_H
#define SPIDIR_COMMON_H

/*=============================================================================

    SPIDIR - SPecies Informed DIstanced-based Reconstruction
    
    Matt Rasmussen
    Wed Jun 13 22:09:24 EDT 2007


=============================================================================*/

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

// constants
#ifndef INFINITY
#   define INFINITY 1e1000
#endif



//=============================================================================
// Math

// indexing a matrix stored as a single array (row-major)
// m: number of columns
// i: row index
// j: column index
#define matind(m, i, j) ((m)*(i) + (j))


inline float frand(float max=1.0)
{ return rand() / float(RAND_MAX) * max; }

inline float frand(float min, float max)
{ return min + (rand() / float(RAND_MAX) * (max-min)); }

inline int irand(int max)
{
    const int i = int(rand() / float(RAND_MAX) * max); 
    return (i == max) ? max - 1 : i;
}

inline int irand(int min, int max)
{
    const int i = min + int(rand() / float(RAND_MAX) * (max - min)); 
    return (i == max) ? max - 1 : i;
}


// computes the log(normalPdf(x | u, s^2))
inline float normallog(float x, float u, float s)
{
    const float log_sqrt_2_pi = 0.91893853320467267;
    if (s == 0.0)
        return -INFINITY;
    else
        return - logf(s) - log_sqrt_2_pi - (x-u)*(x-u) / (2.0*s*s);
    
    //return - log(s) - log(sqrt(2.0*M_PI)) - (x-u)*(x-u) / (2.0*s*s);
    //return log(1.0/(s * sqrt(2.0*M_PI)) * exp(- (x-u)*(x-u) / (2.0 * s*s)));
}

extern "C" {

float poisson(int x, float lambda);

float normalvariate(float mu, float sigma);

inline float expovariate(float lambda)
{ return -log(frand()) / lambda; }



} // extern "C"

template <class T>
double variance(T *vals, int size)
{
    double mean = 0.0;
    for (int i=0; i<size; i++)
        mean += vals[i];
    mean /= size;
    
    double tot = 0.0;
    for (int i=0; i<size; i++)
        tot += (vals[i] - mean) * (vals[i] - mean);
    
    return tot / (size - 1);
}

template <class T>
double variance(T *vals, int size, T mean)
{
    double tot = 0.0;
    for (int i=0; i<size; i++) {
        T d = vals[i] - mean;
        tot += d * d;
    }
    
    return tot / (size - 1);
}


template <class T>
double stdev(T* vals, int size)
{
    return sqrt(variance(vals, size));
}


template <class T>
double stdev(T* vals, int size, T mean)
{
    return sqrt(variance(vals, size, mean));
}


class RunningStat
{
public:
    RunningStat() : n(0) {}

    void clear()
    {
        n = 0;
    }

    void push(double x)
    {
        n++;
        
        // See Knuth TAOCP vol 2, 3rd edition, page 232
        if (n == 1) {
            oldM = newM = x;
            oldS = 0.0;
        } else {
            newM = oldM + (x - oldM)/n;
            newS = oldS + (x - oldM)*(x - newM);
    
            // set up for next iteration
            oldM = newM; 
            oldS = newS;
        }
    }

    int numValues() const
    {
        return n;
    }

    double mean() const
    {
        return (n > 0) ? newM : 0.0;
    }

    double variance() const
    {
        return ((n > 1) ? newS/(n - 1) : 0.0 );
    }

    double sdev() const
    {
        return sqrt(variance());
    }

private:
    int n;
    double oldM, newM, oldS, newS;
};
 



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


// computes log(a + b) given log(a) and log(b)
inline double logadd(double lna, double lnb)
{
    double diff = lna - lnb;
    if (lna == 1.0)
        return lnb;
    if (lnb == 1.0)
        return lna;
    if (diff < 500.0)
        return log(exp(diff) + 1.0) + lnb;
    else
        return lna;
}

void invertPerm(int *perm, int *inv, int size);

template <class T>
void permute(T* array, int *perm, int size)
{
    T *tmp = new T [size];
    
    // transfer permutation to temp array
    for (int i=0; i<size; i++)
        tmp[i] = array[perm[i]];
    
    // copy permutation back to original array
    for (int i=0; i<size; i++)
        array[i] = tmp[i];

    delete [] tmp;
}

template <class T>
T ipow(T val, int expo)
{
    T result = 1.0;
    unsigned int e = expo;

    if ((int)e < 0) {
	e = -e;
	val = 1.0 / val;
    }

    while (true) {
	if (e & 1)
	    result *= val;
	if ((e >>= 1) == 0)
	    break;
	val *= val;
    }

    return result;
}


int choose(int n, int k);

double fchoose(int n, int k);


template <class T>
int findval(T *array, int size, const T &val)
{
    for (int i=0; i<size; i++)
        if (array[i] == val)
            return i;
    return -1;
}


//=============================================================================
// sorting

template <class KeyType, class ValueType>
struct RankSortCmp
{
    RankSortCmp(ValueType *values): values(values) {}
    
    bool operator()(KeyType i, KeyType j)
    { return values[i] < values[j]; }
    
    ValueType *values;
};

template <class KeyType, class ValueType>
void ranksort(KeyType *keys, ValueType *values, int size)
{
    RankSortCmp<KeyType, ValueType> cmp(values);
    sort(keys, keys + size, cmp);
}


//=============================================================================
// input/output

void printIntArray(int *array, int size);
void printFloatArray(float *array, int size);



} // namespace spidir

#endif // SPIDIR_COMMON_H
