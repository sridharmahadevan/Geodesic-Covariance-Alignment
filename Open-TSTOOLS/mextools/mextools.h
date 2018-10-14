#ifndef MEXTOOLS_H
#define MEXTOOLS_H

#include <cmath>
#include <cstdio>
#include <iostream>
#include <vector>
#include <algorithm>
#include <iterator>
#include <functional>
#include <numeric>

// Determine if this C++ code is compiled to a mex-file
// Otherwise, attach code that will produce a stand-alone executable (using mex2main.*)
#ifdef MATLAB_MEX_FILE		// compile to standard matlab mex-dll
	#include "mex.h"
	#undef malloc
	#undef realloc
	#undef free
	#define malloc mxMalloc
	#define realloc mxRealloc
	#define free mxFree
	#define printf mexPrintf
#elif defined OCTAVE_OCT_FILE	// compile to octave oct-dll
	#include "mexwrap.h"
	#include "mex2oct.h"
	#include "mex2oct.cpp"
	#include "mexwrap.cpp"
#else
	#include "mexwrap.h"		// compile to a.out binary
	#include "mex2main.h"
	#ifndef DONT_INCLUDE_MEX_2_MAIN_CPP 
		#include "mex2main.cpp"
		#include "mexwrap.cpp"
	#endif
#endif

// declare main class names 
class mvector;
class mmatrix;
// This file gives an example class for the implementation of a point_set which can
// be used by the nearest neighbor algorithm
// This particular implementation can be used for use with Matlab mex-files, where
// a point set is given as a Fortran style matrix (column major). The coordinates of
// one point is given by one row of this matrix

// interleaved_pointer	: a smart pointer that iterates with an interleave of #inc elements over an array
// May be used to iterate over rows or columns of a dense matrix of type T
template<class T>
class interleaved_pointer // : public forward_iterator<T, int>
{
	private:
		T* ptr;
		const long increment;	
	public:
		typedef interleaved_pointer self;
		 
		inline interleaved_pointer(T* const p, const long inc) : ptr(p), increment(inc) {};
		inline T operator*() const { return *ptr; }
		inline T& operator*()  { return *ptr; }

		inline T operator[](const long i) const {
			return ptr[i*increment];
		}		
		inline self& operator++() {
    		ptr+=increment;
    		return *this;
  		}
  		inline self& operator--() {
    		ptr-=increment;
    		return *this;
  		}		
		bool operator==(const interleaved_pointer& x) {
			return (ptr == x.ptr); 
		}
		bool operator!=(const interleaved_pointer& x) {
			return (ptr != x.ptr); 
		}		
};	


// Simple double vector class
// The base used for indexing (zero in C, one in Fortran) is adjustable
// No automatic resizing, no bound checking, just a lightweight vector class
template<class T>
class nbased_vector {
	protected:
		const long L;				// length
		T* v;
		const bool owner;			// flag to indicate ownership of allocated memory
		const long base_index;	// index at which the first element of this array can be addressed (e.g 1 for Matlab)
			
	public:
		nbased_vector(T* const V, const long Base_index = 1, const long l = 1) : v(V-Base_index), owner(false), base_index(Base_index), L(l) {};
		nbased_vector(const long l, const long Base_index = 1) : v(0), owner(true), base_index(Base_index), L(l)
		{
			if (L > 0) {
				v = new T[L];
				v -= base_index;
			}
		}	
		~nbased_vector() { 
			if (owner) {
				v += base_index; 
				delete[] v; 
			}
		}
		
		unsigned long length() const { return L; }		

		T operator()(const long i) const { return v[i]; }
		T& operator()(const long i) { return v[i]; }
		
		typedef T* iterator;
		iterator begin() { return &(v[base_index]); }
		iterator end() { return &(v[L+base_index]); }	// past the end (STL style)
};

// Simple vector class for use in Matlab mex-files, one based indexing, so the C++ code will look very
// similar to the m-file code
//
// Can be used to access input arguments : mvector vec(prhs[0]);
// or to create output arguments : mvector out(plhs[0] = mxCreateDoubleMatrix(24, 1, mxREAL));
//
// for (long j=1; j <= 24; j++) out(j) = 8.05 + 0.23 * j;

class mvector : public nbased_vector<double> {
	protected:	
		mvector(const mvector& a); 	 // give dummy copy constructor to prevent copying 
		mvector(double* V, const long length) :
			nbased_vector<double>((double*) V, length) {};
				
	public:
		mvector(const mxArray* const a) : 
			nbased_vector<double>((double*) mxGetPr(a), 1, mxGetM(a) * mxGetN(a)) {};
			
		// create a vector from a column (1..N) of a matlab M by N matrix
		mvector(const mxArray* a, const long col) : 
			nbased_vector<double>(((double*) mxGetPr(a))+ (col-1)*mxGetM(a), 1, mxGetM(a)) {};

		mvector(const long l) : nbased_vector<double>(l) {};
				
		friend class mmatrix;
};

// Simple N by M double matrix with one based indexing for rows, columns
class mmatrix : public nbased_vector<double> {
	protected:
		const long M;
		const long N;
		
		mmatrix(const mmatrix& a);		// give dummy copy constructor to prevent copying 
		
	public:
		mmatrix(const mxArray* const a) :
			nbased_vector<double>((double*) mxGetPr(a), mxGetM(a)+1, mxGetM(a)*mxGetN(a)), M(mxGetM(a)), N(mxGetN(a)) {
				if (mxGetNumberOfDimensions(a) > 2) {
					mexErrMsgTxt("Class mmatrix cannot handle multidimensional Matlab arrays");
					return;
				}
			};
		mmatrix(const long m, const long n) : 
			nbased_vector<double>(m*n, m+1), M(m), N(n) {};	
			
		long getM() const { return M; }
		long getN() const { return N; }
		
		double operator()(const long i, const long j) const { return v[i+j*M]; }
		double& operator()(const long i, const long j) { return v[i+j*M]; }
		
		typedef double* column_iterator;
		typedef interleaved_pointer<double> row_iterator;
		
		column_iterator column_begin(int col) { 
			return (v + col * M + 1);
		}
		column_iterator column_end(int col) { 
			return (v + (col+1) * M + 1);
		}	
		row_iterator row_begin(int row) { 
			return interleaved_pointer<double>(v + M + row, M);
		}
		row_iterator row_end(int row) { 
			return interleaved_pointer<double>(v + M*(N+1) + row, M);
		}			
};

#endif
