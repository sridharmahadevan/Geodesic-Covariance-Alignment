#ifndef SIGMOIDS_H
#define SIGMOIDS_H

// Header file containing some useful activation functions for neural networks applications

#include <cmath>
#include <utility>

/* fast, but crude  replacement for libm's exp */
/* proposed by Nicol Schraudolph */


#ifndef M_LN2
#define M_LN2 0.69314718055994530942 
#endif

using namespace std;

static union {
	double d;
	struct {
		int j,i;
	} n;
} _ns_eco;

#define EXP_A (1048576/M_LN2)
#define EXP_C 60801

#define EXP( y ) (_ns_eco.n.i = (int) ( EXP_A * ( y ) + (1072693248 - EXP_C)), _ns_eco.d)

/* fast replacement for libm's tanh */

inline double sigmoid_tanh(const double x) {
	const pair<double, double> y = (x >= 0) ? make_pair(x, 1.0) : make_pair(-x, -1.0);
	
	if (y.first < 22.0) {
		const double ex = EXP(-(x+x)) - 1.0;
		return -ex/(ex + 2.0);	
	} else 
		return y.second;
}

inline double sigmoid_exp(const double x) {
	return EXP( x );
}


#undef EXP_A
#undef EXP_C
#undef EXP 

#endif
