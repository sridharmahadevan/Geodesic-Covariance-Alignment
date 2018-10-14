// mex2main.h
// Create a standalone application from an mex file by linking the mex-code
// with this wrapper (and mex2main.cpp) which gives the basic matlab functionality that 
// the user's mex-code expects to have

#ifndef MEX2MAIN_H
#define MEX2MAIN_H

#define mexPrintf printf 

namespace mex2main {
	// prototypes for argument passing fucntions
	long parse_input_list(const char* filename, mxArray** prhs);	// substitute a matlab command call by an file based arguement list 
	mxArray* readmatrix(FILE* fp);	// arguments that are given from command line are now read in from ASCII file
	void save_array(const mxArray* x, const char* filename); 	// output arguments are written as ASCII file
}
#endif




