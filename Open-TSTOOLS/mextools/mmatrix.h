#ifndef MMATRIX_H
#define MMATRIX_H


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
class interleaved_pointer	
{
	private:
		const T* ptr;
		const long increment;	
	public:
		typedef interleaved_pointer self;
		inline interleaved_pointer(const T* const p, const long inc) : ptr(p), increment(inc) {};
		inline T operator*() const { return *ptr; }
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
			return ptr == x.ptr; 
		}
		bool operator!=(const interleaved_pointer& x) {
			return ptr != x.ptr; 
		}		
};	


// Simple double vector class
// The base used for indexing (zero in C, one in Fortran) is adjustable
// No automatic resizing, no bound checking, just a lightweight vector class
template<class T>
class onebased_vector {
	protected:
		T* v;
		const int owner;			// flag to indicate ownership of allocated memory
		const long base_index;	// index at which the first element of this array can be addressed (e.g 1 for Matlab)
			
	public:
		onebased_vector(T* const V, const long Base_index = 1) : v(V-Base_index), owner(0), base_index(Base_index) {};
		onebased_vector(const long L, const long Base_index = 1) : v(0), owner(1), base_index(Base_index)
		{
			if (L > 0) {
				v = new T[L];
				v -= base_index;
			}
		}	
		~onebased_vector() { 
			if (owner) {
				v += base_index; 
				delete[] v; 
			}
		}
		
		T operator()(const long i) const { return v[i]; }
		T& operator()(const long i) { return v[i]; }
};

// Simple vector class for use in Matlab mex-files, one based indexing, so the C++ code will look very
// similar to the m-file code
//
// Can be used to access input arguments : mvector vec(prhs[0]);
// or to create output arguments : mvector out(plhs[0] = mxCreateDoubleMatrix(24, 1, mxREAL));
//
// for (long i=1; i <= 24; i++) out(i) = 8.05 + 0.23 * i;

class mvector : public onebased_vector<double> {
	protected:
		const unsigned long Length;
		
		mvector(const mvector& a); 	 // give dummy copy constructor to prevent copying 
		mvector(double* V, const long length) :
			onebased_vector<double>((double*) V), Length(length) {};
				
	public:
		mvector(const mxArray* const a) : 
			onebased_vector<double>((double*) mxGetPr(a)), Length(mxGetM(a) * mxGetN(a)) {};
			
		// create a vector from a column (1..N) of a matlab M by N matrix
		mvector(const mxArray* a, const long col) : 
			onebased_vector<double>(((double*) mxGetPr(a))+ (col-1)*mxGetM(a)), Length(mxGetM(a)) {};

		mvector(const long L) : onebased_vector<double>(L), Length(L) {};
		unsigned long length() const { return Length; }
		
		typedef double* iterator;
		iterator begin() { return &(v[1]); }
		iterator end() { return &(v[Length+1]); }	// past the end (STL style)
		
		friend class mmatrix;
};

// Simple N by M double matrix with one based indexing for rows, columns
class mmatrix : protected onebased_vector<double> {
	protected:
		const long M;
		const long N;
		
		mmatrix(const mmatrix& a);		// give dummy copy constructor to prevent copying 
		
	public:
		mmatrix(const mxArray* const a) :
			onebased_vector<double>((double*) mxGetPr(a), mxGetM(a)+1), M(mxGetM(a)), N(mxGetN(a)) {
				if (mxGetNumberOfDimensions(a) > 2) {
					mexErrMsgTxt("Class mmatrix cannot handle multidimensional Matlab arrays");
					return;
				}
			};
		mmatrix(const long m, const long n) : 
			onebased_vector<double>(m*n, m+1), M(m), N(n) {};	
			
		long getM() const { return M; }
		long getN() const { return N; }
		
		double operator()(const long i, const long j) const { return v[i+j*M]; }
		double& operator()(const long i, const long j) { return v[i+j*M]; }
		
		typedef interleaved_pointer<double> column_iterator;
		typedef interleaved_pointer<double> row_iterator;
		
		column_iterator column_begin(int col) { 
			return interleaved_pointer<double>(v + col * M + 1, 1);
		}
		column_iterator column_end(int col) { 
			return interleaved_pointer<double>(v + (col+1) * M + 1, 1);
		}	
		row_iterator row_begin(int row) { 
			return interleaved_pointer<double>(v + M + row, M);
		}
		row_iterator row_end(int row) { 
			return interleaved_pointer<double>(v + M*(N+1) + row, M);
		}			
};

#endif
