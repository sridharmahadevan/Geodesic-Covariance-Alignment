#include "mexwrap.h"
#include "mex2oct.h"

#define mexPrintf printf
#define NAME_TO_STRING( s ) #s

void mexErrMsgTxt(const char* text)
{
	error(text);
}

mxArray* Matrix2mxArray(Matrix& m) {
	const long M = m.rows();
	const long N = m.columns(); 
	
	mxArray* x = mxCreateDoubleMatrix(M, N, mxREAL);
	
	const double* const source = m.fortran_vec();
	double* const drain = mxGetPr(x);
	
	memcpy(drain, source, M*N*sizeof(double));
		
// 	mxSetPr(x, source);	
// 		
// 	if (mxGetPr(x) != source) {
// 		mexErrMsgTxt("cannot set fortran_vec pointer in mxArray");
// 	}	
	
	return x;
}

DEFUN_DLD ( MEXNAME , args, nargout, "")
{
	int nlhs;

  	octave_value_list retval;	// The list of values to return.  See the declaration in oct-obj.h
   
  	int nrhs = args.length();
	
	if (nargout < 1)
		nlhs = 1;	// require at least one output argument
	else
		nlhs = nargout;
		
	mxArray** plhs = new mxArray*[nlhs+16];
	mxArray** prhs = new mxArray*[nrhs+16];
	
	for (int i=0; i < nlhs+16; i++) {
		plhs[i] = 0;
	}
	
 	// Eingabeargumente in mxArrays verwandeln
	
	for (int i = 0; i < nrhs; i++)  
    {
		octave_value arg = args(i); 	// args wird von 0 bis nrhs-1 adressiert
		if (arg.is_real_type ()) {	
			Matrix m = arg.matrix_value();

			if (error_state)
				return retval;
				
			prhs[i] = Matrix2mxArray(m);	
		}
//         else if (arg.) {
//         
//         }
		else {
			gripe_wrong_type_arg( NAME_TO_STRING( MEXNAME ), arg);
			return retval;
		}
    }
	
	mexFunction(nlhs, plhs, nrhs, (const mxArray**) prhs); 

	// Ausgabeargumente in Octave Matrizen verwandeln
	
	long counter = 0;
 	for (int i = 0; i < nlhs; i++)
    {	
		mxArray* arr = plhs[i];
		
		if (arr) {
			const long M = mxGetM(arr);
			const long N = mxGetN(arr); 
			
			Matrix m(M, N);

			const double* const source = mxGetPr(arr);
			double* const drain = m.fortran_vec();

			memcpy(drain, source, M*N*sizeof(double));
						 
			retval(counter) = m;
			counter++;
		}
    }	
	
	delete[] plhs;
	delete[] prhs;
	
	// free all mxArrays that have not been freed previously
	mxArray::cleanup();

	return retval;
}



	
