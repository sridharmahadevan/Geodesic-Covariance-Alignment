#ifndef NNETWORK_H
#define NNETWORK_H

#include <cmath>
#include "mextools/sigmoids.h"
#include "mextools/random.h"

static Random<MT_RNG> rng;

// Structure state stores the current weight of one entry of a network, together with gradient information
inline double sign(const double x) {
    if (x > 0) 
        return 1.0;
    else if (x < 0)
        return -1.0;
    else
        return 0.0;
}

inline double squash(const double x) 
{
    if (x > 0) 
        return (x < 1.0 ? x : log(x)+1.0);
    else
        return (x > -1.0 ? x : -(log(-x)+1.0));       
}


// for this weight, stepsize etc.
class state {
    public:
    double weight;
    
#ifndef ONLY_CALC    
    double grad;
    double lastgrad;
    double stepsize;
    double laststep;
    
    void initialize(const double initial_weightspan, const double initial_stepsize) {
        if (initial_weightspan >= 0) {
            weight = initial_weightspan * (rng()+rng()-1.0);
            grad = 0.0;
            lastgrad = 0.0;
            stepsize = initial_stepsize;
            laststep = 0.0;
        } else {
            grad = 0.0;
            lastgrad = 0.0;
            stepsize *= initial_stepsize;
            laststep = 0.0;
        }
    }

    double rprop_step(const double err, const double lasterr, const double weightdecay = 0.0, const double weightlimit = 4.0) {
            grad -= weightdecay * weight;    
        
            const double signchange = grad * lastgrad;
            
            if (signchange > 0) {
                const double newstepsize = min(stepsize * 1.3, 0.4);
                lastgrad = grad;
                laststep = sign(grad) * newstepsize;
                weight += laststep;
                
                if (weight > weightlimit) {
                    const double d = weight - weightlimit;
                    weight -= d;
                    laststep -= d;
                } else if (weight < -weightlimit) {
                    const double d = weight+weightlimit;
                    weight -= d;
                    laststep -= d;
                } else {
                    stepsize = newstepsize; 
                }
            } else if (signchange < 0) {
                stepsize = max(stepsize * 0.25, 0.0000004);
                lastgrad = 0;
                laststep = 0.0;
                //if (err > lasterr) {
                    weight -= 0.61803 * laststep;
                //}    
            } else {
                lastgrad = grad;
                laststep = sign(grad) * stepsize;
                weight += laststep;
                
                if (weight > weightlimit) {
                    const double d = weight - weightlimit;
                    weight -= d;
                    laststep -= d;
                } else if (weight < -weightlimit) {
                    const double d = weight+weightlimit;
                    weight -= d;
                    laststep -= d;
                }
            }
            
            grad = 0.0;     // clear gradient information for next round
            return weight * weight;
    }
    
    // Idea : Each sample will be treated twice, we compare error and sign of gradients between two
    // trials to see if stepwidth is right, otherwise stepwidth is decreased
    // Still a mu is given, stepsize is supposed to change around 1.0
    void stochgrad_step(const double mu, const int trial, const double e, const double olde) {
           const double signchange = grad * lastgrad;
            if (trial) {
                if (signchange > 0) {
                    stepsize = min(stepsize * 1.1, 10.0);
                    weight += mu * sign(grad) * stepsize;         
                } else if (signchange < 0) {
                    stepsize = max(stepsize * 0.8, 0.0001);
                    if (e > olde) {
                        weight -= 0.61803 * laststep;
                    }    
                } else {
                    weight += mu * sign(grad) * stepsize;
                }
            } else {
                lastgrad = grad;
                laststep = mu * sign(grad) * stepsize;
                weight += laststep;
            }
           
            grad = 0.0;     // clear gradient information for next round
    }
#endif    
    
    state() : weight(0.0) {};
};

        

// Square and absolute epsilon-insensitive loss functions


inline pair<double, double> epsquadloss(const double x, const double y, const double eps)
{
    double d = x - y;

    if (d > eps) {
        d -= eps;
    } else if (d < -eps)
    {
        d += eps;
    } else {
        d = 0.0;
    }
    
    return make_pair(d*d,  2.0*d);
}


inline pair<double, double> epsabsloss(const double x, const double y, const double eps)
{
    const double d = x - y;

    if (d > eps) {
        return make_pair(d-eps,  1.0);
    } else if (d < -eps)
    {
        return make_pair(-d-eps,  -1.0);
    } else {
        return make_pair(0.0,  0.0);
    }
}


// loss for binary classification (0/1 coding), errors only linear with eps margin
inline pair<double, double> epsclassloss(const double x, const double y, const double eps)
{
    if (x) {
        const double d = x - y;  // when x (desired) is greater than y, we have a nonzero loss
        
        if (d > eps) 
            return make_pair(d-eps,  1.0);
        else    
            return make_pair(0.0,  0.0);
    } else {
          const double d = y - x;  // when x (desired) is lower than y, we have a nonzero loss
        
            if (d > eps) 
                return make_pair(d-eps,  -1.0);
            else    
                return make_pair(0.0,  0.0);
    }
}


// take negative log-likelihood function for binary classification, x must be 0 for class 0 and 1 for class 1
inline pair<double, double> logitloss(const double x, const double y, const double eps)
{
    const double q = exp(y);
    const double p = q / (1.0 + q);
    
    
    if (x) {
        if ((p+eps) >= 1.0) {
            return make_pair(0.0, 0.0);         
        } else {
            return make_pair(-log(p+eps), p/((1+q)*(p+eps)));         
        }
    } else {
        if ((1.0-p+eps) >= 1.0) {
            return make_pair(0.0, 0.0);         
        } else {
            return make_pair(-log(1.0-p+eps), -p/((1+q)*(1.0-p+eps)));         
        }
    }
}

// Fast activation function (piecewise linear activation function)

inline double pl(const double x)
{
    if (x < 0.5)
        return x;
    else if (x < 1.5)
        return 0.5 + 0.5*(x-0.5);
    else if (x < 4.5)
        return 1.0 + 0.05 * (x-1.5);
    else 
        return 1.15 + 0.01 * (x-4.5);    
}

inline double pld(const double x)
{
    if (x < 0.5)
        return 1.0;
    else if (x < 1.5)
        return 0.5;
    else if (x < 4.5)
        return 0.05;
    else
        return 0.01;    
}

inline double piecewice_linear(const double x) 
{

        if (x >= 0)
            return pl(x);
        else
            return -pl(-x);
}

inline double piecewice_linear_derivative(const double x) 
{
    if (x >= 0)
        return pld(x);
    else
        return pld(-x);
}

// tanh like basis function activation function, return function value and derivative
inline pair<double, double> sigmoid(const double x) 
{
#ifdef LINEARIZE	
	return make_pair(x, 1.0);
#else
	const double q = tanh(x);
	return make_pair(q + 0.15 * x, 1.15 - q * q);
#endif
}

// tanh like basis function activation function, return function value and derivative
inline pair<double, double> linear(const double x) 
{	
	return make_pair(x, 1.0);
}



// Radial basis function activation function, return function value and derivative
inline pair<double, double> rbf(const double x) 
{
#ifdef LINEARIZE	
	return make_pair(x, 1.0);
#else
	const double p = exp(-x*x);
    return make_pair(p + 0.02 * x, -2*x*p + 0.02);
#endif   
}

inline pair<double, double> myloss(const double x, const double y, const double eps, const int mode)
{
    if (mode == 1) {
        // absolute loss
        return(epsabsloss(x, y, eps));
    } else if (mode == 2) {
        // logit loss, but we need y to be 0 for class 0 and 1 for class 1
        return(logitloss(x, y, eps));
    } else if (mode == 3) {
        // binary classification loss (eps, absolute)
        return(epsclassloss(x, y, eps));
    }
    else {
        // quadratic loss
        return(epsquadloss(x, y, eps));
    }    
}

#endif
