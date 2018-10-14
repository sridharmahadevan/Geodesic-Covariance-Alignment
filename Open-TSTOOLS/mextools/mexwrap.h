// mex2wrap.h
//
// Emulate mex-environment so that a mex-file thinks it is linked against
// Matlab

#ifndef MEXWRAP_H
#define MEXWRAP_H

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <set>

using namespace std;

#define mxMAXNAM  32    /* maximum name length */

typedef unsigned int mxChar;
typedef char mxName[mxMAXNAM];

typedef enum {
	mxCELL_CLASS = 1,
	mxSTRUCT_CLASS,
	mxOBJECT_CLASS,
	mxCHAR_CLASS,
	mxSPARSE_CLASS,
	mxDOUBLE_CLASS,
	mxSINGLE_CLASS,
	mxINT8_CLASS,
	mxUINT8_CLASS,
	mxINT16_CLASS,
	mxUINT16_CLASS,
	mxINT32_CLASS,
	mxUINT32_CLASS,
	mxINT64_CLASS,		/* place holder - future enhancements */
	mxUINT64_CLASS,		/* place holder - future enhancements */
    mxOPAQUE_CLASS,
	mxUNKNOWN_CLASS = -1
} mxClassID;

typedef enum {
    mxREAL,
    mxCOMPLEX
} mxComplexity;

typedef enum {
    mxDOUBLE_ARRAY      = 2,    /* start here to align with dispatch */
    mxSPARSE_ARRAY,
    mxCHARACTER_ARRAY,
    mxCELL_ARRAY,
    mxSTRUCTURE_ARRAY,
    mxFLOAT_ARRAY,
    mxINT8_ARRAY,
    mxUINT8_ARRAY,
    mxINT16_ARRAY,
    mxUINT16_ARRAY,
    mxINT32_ARRAY,
    mxUINT32_ARRAY,
    mxINT64_ARRAY,      /* place holder - future enhancements */
    mxUINT64_ARRAY,     /* place holder - future enhancements */
    mxOBJECT_ARRAY,
    mxUNKNOWN_ARRAY,
} mxArrayType;

typedef enum {
    mxLOCAL_SCOPE,
    mxLOCAL_STATIC_SCOPE,
    mxGLOBAL_SCOPE,
    mxMEMBER_SCOPE,
    mxTEMPORARY_SCOPE,
    mxPERSISTENT_SCOPE,
    mxCONSTANT_SCOPE,
    mxAUTO_SCOPE
} mxVariableScope;


class mxArray;
typedef mxArray* mxArray_ptr;

// Give prototypes for functions in mexwrap.cpp

/*
 * allocate memory, notifying registered listener
 */
extern void *mxMalloc(size_t n);


/*
 * allocate cleared memory, notifying registered listener.
 */
extern void *mxCalloc(size_t n, size_t	size);


/*
 * free memory, notifying registered listener.
 */
extern void mxFree(void *ptr);	/* pointer to memory to be freed */


/*
 * reallocate memory, notifying registered listener.
 *
 * Note: If realloc() fails to allocate a block, it is the end-user's
 * responsibility to free the block because the ANSI definition of realloc
 * states that the block will remain allocated.  realloc() returns a NULL in
 * this case.  This means that calls to realloc of the form:
 *
 * x = mxRealloc(x, size)
 *
 * will cause memory leaks if realloc fails (and returns a NULL).
 *
 */
extern void *mxRealloc(void *ptr, size_t size);

/*
 * Return the class (catergory) of data that the array holds.
 */
extern mxClassID mxGetClassID(const mxArray* pa);


/* 
 * Get pointer to array name.
 */
extern const char *mxGetName(const mxArray* pa);

/* 
 * Set array name.  This routine copies the string pointed to by s
 * into the mxMAXNAM length character name field.
 */
extern void mxSetName(mxArray    *pa, const char *s);


/*
 * Get pointer to data
 */
extern void *mxGetData(const mxArray* pa);

/*
 * Set pointer to data
 */
extern void mxSetData(
    mxArray* pa,
    void  *newdata		/* pointer to data */
    );

/*
 * Get real data pointer for numeric array
 */
extern double *mxGetPr(const mxArray* pa);


/*
 * Set real data pointer for numeric array
 */
extern void mxSetPr(mxArray* pa, double  *pr);


/*
 * Get imaginary data pointer for numeric array
 */
extern void *mxGetImagData(const mxArray* pa);

/*
 * Set imaginary data pointer for numeric array
 */
extern void mxSetImagData(
    mxArray* pa,
    void    *newdata		/* imaginary data array pointer */
    );


/*
 * Get imaginary data pointer for numeric array
 */
extern double *mxGetPi(const mxArray* pa);


/*
 * Set imaginary data pointer for numeric array
 */
extern void mxSetPi(mxArray* pa, double  *pi);


/* 
 * Determine whether the specified array contains numeric (as opposed 
 * to cell or struct) data.
 */
extern bool mxIsNumeric(const mxArray* pa);


/* 
 * Determine whether the given array is a cell array.
 */
extern bool mxIsCell(const mxArray* pa);


/*  
 * Determine whether the given array contains character data. 
 */
extern bool mxIsChar(const mxArray* pa);


/*
 * Determine whether the given array is a sparse (as opposed to full). 
 */
extern bool mxIsSparse(const mxArray* pa);


/*
 * Determine whether the given array is a structure array.
 */
extern bool mxIsStruct(const mxArray* pa);


/*
 * Determine whether the given array contains complex data.
 */
extern bool mxIsComplex(const mxArray* pa);


/*
 * Determine whether the specified array represents its data as 
 * double-precision floating-point numbers.
 */
extern bool mxIsDouble(const mxArray* pa);


/*
 * Determine whether the specified array represents its data as 
 * single-precision floating-point numbers.
 */
extern bool mxIsSingle(const mxArray* pa);


/*
 * Determine whether the given array's logical flag is on.
 */ 
extern bool mxIsLogical(const mxArray* pa);


/*
 * Determine whether the specified array represents its data as 
 * signed 8-bit integers.
 */
extern bool mxIsInt8(const mxArray* pa);


/*
 * Determine whether the specified array represents its data as 
 * unsigned 8-bit integers.
 */
extern bool mxIsUint8(const mxArray* pa);


/*
 * Determine whether the specified array represents its data as 
 * signed 16-bit integers.
 */
extern bool mxIsInt16(const mxArray* pa);


/*
 * Determine whether the specified array represents its data as 
 * unsigned 16-bit integers.
 */
extern bool mxIsUint16(const mxArray* pa);


/*
 * Determine whether the specified array represents its data as 
 * signed 32-bit integers.
 */
extern bool mxIsInt32(const mxArray* pa);


/*
 * Determine whether the specified array represents its data as 
 * unsigned 32-bit integers.
 */
extern bool mxIsUint32(const mxArray* pa);


/*
 * Get 8 bits of user data stored in the mxArray header.  NOTE: This state
 * of these bits are not guaranteed to be preserved after API function
 * calls.
 */
extern int mxGetUserBits(
    const mxArray	*pa
    );


/*
 * Set 8 bits of user data stored in the mxArray header. NOTE: This state
 * of these bits are not guaranteed to be preserved after API function
 * calls.
 */ 
extern void mxSetUserBits(mxArray	*pa, int	value);

/*
 * Get the real component of the specified array's first data element.
 */
extern double mxGetScalar(const mxArray* pa);


/*
 * Specify that the data in an array is to be treated as Boolean data.
 */
extern void mxSetLogical(mxArray* pa);


/*
 * Specify that the data in an array is to be treated as numerical
 * (as opposed to Boolean) data. 
 */
extern void mxClearLogical(mxArray* pa);


/*
 * Is the isFromGlobalWorkspace bit set?
 */
extern bool mxIsFromGlobalWS(const mxArray* pa);


/*
 * Set the isFromGlobalWorkspace bit.
 */
extern void mxSetFromGlobalWS(mxArray* pa, bool global);


/*
 * Get number of dimensions in array
 */
extern int mxGetNumberOfDimensions(const mxArray* pa);


/* 
 * Get row dimension
 */
extern int mxGetM(const mxArray* pa);	


/* 
 * Set row dimension
 */
extern void mxSetM(mxArray* pa, int m /* row dimension */);


/* 
 * Get column dimension
 */
extern int mxGetN(const mxArray* pa);

/*
 * Get pointer to dimension array
 */
extern const int *mxGetDimensions(const mxArray* pa);


/*
 * Is array empty
 */
extern bool mxIsEmpty(const mxArray* pa);


/*
 * Get row data pointer for sparse numeric array
 */
extern int *mxGetIr(const mxArray* pa);

/*
 * Set row data pointer for numeric array
 */
extern void mxSetIr(
    mxArray* pa,
    int     *newir		/* row data array pointer */
    );


/*
 * Get column data pointer for sparse numeric array
 */
extern int *mxGetJc(
    const mxArray* pa
    );


/*
 * Set column data pointer for numeric array
 */
extern void mxSetJc(mxArray* pa, int *newjc		/* column data array pointer */);

/*
 * Get maximum nonzero elements for sparse numeric array
 */
extern int mxGetNzmax(const mxArray* pa);

/*
 * Set maximum nonzero elements for numeric array
 */
extern void mxSetNzmax(mxArray* pa,	int nzmax);

/* 
 * Get number of elements in array
 */
extern int mxGetNumberOfElements(const mxArray* pa);

/*
 * Get array data element size
 */
extern int mxGetElementSize(const mxArray* pa);

/* 
 * Return the offset (in number of elements) from the beginning of 
 * the array to a given subscript.
 */
extern int mxCalcSingleSubscript(const mxArray* pa, int nsubs, const int *subs);

/*
 * Get number of structure fields in array
 */
extern int mxGetNumberOfFields(const mxArray* pa);

/*
 * Get a pointer to the specified cell element. 
 */ 
extern mxArray *mxGetCell(const mxArray* pa, int i);

/*
 * Set an element in a cell array to the specified value.
 */
extern void mxSetCell(mxArray* pa, int i, mxArray *value);

/*
 * Get the index to the named field.
 */ 
extern int mxGetFieldNumber(const mxArray* pa, const char *name);

/*
 * Return a pointer to the contents of the named field for 
 * the ith element (zero based).
 */ 
extern mxArray *mxGetFieldByNumber(const mxArray* pa, int i, int fieldnum);

/*
 * Set pa[i][fieldnum] = value 
 */
extern void mxSetFieldByNumber(mxArray* pa, int i, int fieldnum, mxArray *value);

/*
 * Return a pointer to the contents of the named field for the ith 
 * element (zero based).  Returns NULL on no such field or if the
 * field itself is NULL
 */
extern mxArray *mxGetField(const mxArray* pa, int i, const char *fieldname);

/*
 * Set pa[i]->fieldname = value  
 */
extern void mxSetField(mxArray* pa, int i, const char *fieldname, mxArray *value);

/*
 * Return pointer to the nth field name
 */
extern const char *mxGetFieldNameByNumber(const mxArray* pa, int n);

 /* 
 * Return the name of an array's class.  
 */
extern const char *mxGetClassName(const mxArray* pa);

/*
 * Determine whether an array is a member of the specified class. 
 */
extern bool mxIsClass(const mxArray* pa, const char *name);

 /* 
 * Set column dimension
 */
extern void mxSetN(mxArray* pa, int n);

/*
 * Set dimension array and number of dimensions.  Returns 0 on success and 1
 * if there was not enough memory available to reallocate the dimensions array.
 */
extern int mxSetDimensions(mxArray* pa, const int *size, int ndims);

/*
 * Deallocate (free) the heap memory held by the specified array.
 */
extern void mxDestroyArray(mxArray* pa);

/*
 * Create a numeric array and initialize all its data elements to 0.
 */ 
extern mxArray *mxCreateNumericArray(int ndim, const int *dims, mxClassID classid, mxComplexity flag);

/*
 * Create a two-dimensional array to hold double-precision 
 * floating-point data; initialize each data element to 0.
 */
extern mxArray *mxCreateDoubleMatrix(int m, int n, mxComplexity flag);

/*
 * Create a 2-Dimensional sparse array.
 */
extern mxArray *mxCreateSparse(int m, int n, int nzmax, mxComplexity flag);

/*
 * Copies characters from a MATLAB array to a char array
 * nChars is the number of bytes in the output buffer
 */
extern void mxGetNChars(const mxArray* pa, char *buf, int nChars);

/* 
 * Converts a string array to a C-style string.
 */
extern int mxGetString(const mxArray* pa, char *buf, int buflen);

/*
 * Create a NULL terminated C string from an mxArray of type mxCHAR_CLASS
 * Supports multibyte character sets.  The resulting string must be freed
 * with mxFree.  Returns NULL on out of memory.
 */
extern char *mxArrayToString(const mxArray* pa);

/*
 * Create a 1-by-n string array initialized to str.
 */
extern mxArray *mxCreateStringFromNChars(const char *str, int n);

/*
 * Create a 1-by-n string array initialized to null terminated string
 * where n is the length of the string.
 */
extern mxArray *mxCreateString(const char *str);

/*
 * Create an N-Dimensional array to hold string data;
 * initialize all elements to 0.
 */
extern mxArray *mxCreateCharArray(int ndim, const int *dims);

/*
 * Create a string array initialized to the strings in str. 
 */
extern mxArray *mxCreateCharMatrixFromStrings(int m, const char **str);

/*
 * Create a 2-Dimensional cell array, with each cell initialized
 * to NULL.
 */
extern mxArray *mxCreateCellMatrix(int m, int n);

/*
 * Create an N-Dimensional cell array, with each cell initialized
 * to NULL. 
 */
extern mxArray *mxCreateCellArray(int ndim, const int *dims);

/*
 * Create a 2-Dimensional structure array having the specified fields;
 * initialize all values to NULL.
 */ 
extern mxArray *mxCreateStructMatrix(int m, int n, int nfields, const char **fieldnames);

/*
 * Create an N-Dimensional structure array having the specified fields;
 * initialize all values to NULL.
 */ 
extern mxArray *mxCreateStructArray(int ndim, const int *dims, int nfields, const char **fieldnames);

/*
 * Make a deep copy of an array, return a pointer to the copy. 
 */
extern mxArray *mxDuplicateArray(const mxArray *in);


/*
 * Set classname of an unvalidated object array.  It is illegal to 
 * call this function on a previously validated object array. 
 * Return 0 for success, 1 for failure.
 */
extern int mxSetClassName(mxArray* pa, const char *classname);


/* 
 * Add a field to a structure array. Returns field number on success or -1 if inputs
 * are invalid or an out of memory condition occurs.
 */
extern int mxAddField(mxArray* pa, const char *fieldname);

/*
 * Remove a field from a structure array.  Does nothing if no such field exists.
 * Does not destroy the field itself.
*/
extern void mxRemoveField(mxArray* pa, int field);


/*
 * Function for obtaining MATLAB's concept of EPS
 */
extern double mxGetEps(void);


/*
 * Function for obtaining MATLAB's concept of INF (Used in MEX-File callback).
 */
extern double mxGetInf(void);


/*
 * Function for obtaining MATLAB's concept of NaN (Used in MEX-File callback).
 */
extern double mxGetNaN(void);

/*
 * test for finiteness in a machine-independent manner
 */
extern bool mxIsFinite(double x);


/*
 * test for infinity in a machine-independent manner
 */
extern bool mxIsInf(double x);


/*
 * test for NaN in a machine-independent manner
 */
extern bool mxIsNaN(double x);

// create a minimal matlab compatible mxArray class (up to now, only real (== double)
// matrix are supported)

class mxArray {
	private:
		static set<mxArray_ptr, less<mxArray_ptr> > array_table;	// used to keep track of all mxArrays
																	// that are created	
		
		static set<void*, less<void*> > memory_table;															
																	// keep track of all memory that is allocated
                                                                    // mx*** calls
	public: 
    	mxName            name;
    	mxArrayType       type;
    	int               scope;
    	mxArray*          link;
    	int               number_of_dims;
    	int               nelements_allocated;
    	struct {
        	unsigned int    scalar_flag : 1;
        	unsigned int    logical_flag : 1;
        	unsigned int    empty_flag : 1;
        	unsigned int    global_flag : 1;
        	unsigned int    on_arraylist_flag : 1;
        	unsigned int    zero_imag_flag : 1;
        	unsigned int    static_flag : 1;
        	unsigned int    colon_flag : 1;
        	unsigned int    private_data_flag : 1;
        	unsigned int    reference_count : 7;
        	unsigned int    kernel_bits : 8;
        	unsigned int    user_bits : 7;
        	unsigned int    string_flag : 1;
    	}   flags;

    	union {
        	struct {
            	int  m;
            	int  n;
        	}   matrix_dims;
        	int  *dim_array;
    	}   size;

    	union {
        	struct {
            	void  *pdata;
            	void  *pimag_data;
            	int   *ir;
            	int   *jc;
        	}   number_array;

        	struct {
            	mxArray  **cells;
        	}   cell_array;

        	struct {
            	mxArray      **fields;
            	mxName        *field_names;
            	char          *object_classname;
            	int            object_tag;  /* if > 0, structure is object */
            	unsigned int   object_chksum;
            	int            number_of_fields;
        	}   structure_array;

        	struct {
            	mxArray  *pglobal;
        	}   global_array;
    	} data;	
																												
		mxArray();	
		~mxArray();

		friend mxArray *mxCreateDoubleMatrix(int m, int n, mxComplexity flag);
		friend void mxDestroyArray(mxArray *x);
		friend void* mxMalloc(size_t n);
		friend void* mxCalloc(size_t n, size_t	size);
		friend void mxFree(void *ptr);
		friend void* mxRealloc(void *ptr, size_t size);
		
		static void cleanup();	// clear all mxArrays in array_table and all memory
};


// These functions to this prototypes need to be implemented by the user 

void mexErrMsgTxt(const char* text);

// Prototype for THE "mexFunction" call
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]);    
#endif


// struct mxArray_tag {
//     mxName            name;
//     mxArrayType       type;
//     int               scope;
//     mxArray          *link;
//     int               number_of_dims;
//     int               nelements_allocated;
//     struct {
//         unsigned int    scalar_flag : 1;
//         unsigned int    logical_flag : 1;
//         unsigned int    empty_flag : 1;
//         unsigned int    global_flag : 1;
//         unsigned int    on_arraylist_flag : 1;
//         unsigned int    zero_imag_flag : 1;
//         unsigned int    static_flag : 1;
//         unsigned int    colon_flag : 1;
//         unsigned int    private_data_flag : 1;
//         unsigned int    reference_count : 7;
//         unsigned int    kernel_bits : 8;
//         unsigned int    user_bits : 7;
//         unsigned int    string_flag : 1;
//     }   flags;
//   
//     union {
//         struct {
//             int  m;
//             int  n;
//         }   matrix_dims;
//         int  *dim_array;
//     }   size;
//   
//     union {
//         struct {
//             void  *pdata;
//             void  *pimag_data;
//             int   *ir;
//             int   *jc;
//         }   number_array;
//     
//         struct {
//             mxArray  **cells;
//         }   cell_array;
//     
//         struct {
//             mxArray      **fields;
//             mxName        *field_names;
//             char          *object_classname;
//             int            object_tag;  /* if > 0, structure is object */
//             unsigned int   object_chksum;
//             int            number_of_fields;
//         }   structure_array;
// 
//         struct {
//             mxArray  *pglobal;
//         }   global_array;
//     }   data;
// };

