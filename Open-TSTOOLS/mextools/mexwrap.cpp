#include "mexwrap.h"

using namespace std;

mxArray::mxArray() 
	: type(mxUNKNOWN_ARRAY), number_of_dims(2), scope(mxLOCAL_SCOPE), link(0),
	nelements_allocated(0)	
{
	name[0] = 0;

	pair<set<mxArray_ptr, less<mxArray_ptr> >::const_iterator, bool > p;
	
	p = array_table.insert(this);
	
	if (!p.second) {
		mexErrMsgTxt("Error creating new array : pointer is already in use");
	}

#ifdef MEXWRAP_DEBUG
	fprintf(stderr, "Created mxArray : %p\n", this);
	fprintf(stderr, "Table size : %d\n", array_table.size());
#endif
}

mxArray::~mxArray() 
{
	const long count = array_table.erase(this);
	
	if (count != 1) {
		mexErrMsgTxt("Error deleting array : pointer is already cleared");
	}
	
#ifdef MEXWRAP_DEBUG	
	fprintf(stderr, "%d arrays erased\n", count);
	fprintf(stderr, "Table size : %d\n", array_table.size());
#endif
}

// Clear all mxArrays that are still in array_table.
// All mxArrays that are freed by the user by calling mxDestroyArray
// are already erased from the array_table, so we only free left over
// mxArrays
void mxArray::cleanup() {
	set<mxArray_ptr, less<mxArray_ptr> >::iterator i = array_table.begin();
	
	while (i != array_table.end()) {	
		#ifdef MEXWRAP_DEBUG
			fprintf(stderr, "Deleting mxArray : %p\n", *i);
		#endif		
		
		mxDestroyArray(*(i++));
	}
	
	array_table.erase(array_table.begin(), array_table.end());
	
	set<void*, less<void*> >::iterator j = memory_table.begin();
	
	while (j != memory_table.end()) {	
		#ifdef MEXWRAP_DEBUG
			fprintf(stderr, "Deleting memory : %p\n", *j);
		#endif		
		
		mxFree(*(j++));
	}
		
	#ifdef MEXWRAP_DEBUG
		fprintf(stderr, "Array table size : %d\n", array_table.size());
		fprintf(stderr, "Memory table size : %d\n", memory_table.size());
	#endif		
}

mxArray* mxCreateDoubleMatrix(int m, int n, mxComplexity flag)
{	
	mxArray* x = new mxArray;
	
	x->type = mxDOUBLE_ARRAY;
	
	if (n*m > 0) {
		x->data.number_array.pdata = mxMalloc(m * n * sizeof(double));	
	}
	else {
		x->data.number_array.pdata = 0;
	}
	
	if ((flag == mxCOMPLEX) && (n*m > 0)) 
		x->data.number_array.pimag_data = mxMalloc(m * n * sizeof(double));
	else
		x->data.number_array.pimag_data = 0;
	
	x->number_of_dims = 2;

	x->size.matrix_dims.m = m;
	x->size.matrix_dims.n = n;
	
	return x;
}


mxArray *mxCreateString(const char *str) {
	const long length = strlen(str);

	mxArray* x = new mxArray;;
	x->type = mxCHARACTER_ARRAY;
	x->number_of_dims = 2;
	
	x->size.matrix_dims.m = 0;
	x->size.matrix_dims.n = 0;
	
	if (length > 0) {
		x->data.number_array.pdata = (void*) mxMalloc(length * sizeof(mxChar));	
	
		x->size.matrix_dims.m = 1;
		x->size.matrix_dims.n = length;
	
		mxChar* y = (mxChar*) x->data.number_array.pdata; 
		
		for (long i=0; i<length;i++) {
			y[i] = (mxChar) str[i];
		}
	}

	return x;
}


// destroy mxArray structure
void mxDestroyArray(mxArray *x)
{
	if (x) {
		if (x->type == mxDOUBLE_ARRAY) {
			mxFree(x->data.number_array.pdata);
			mxFree(x->data.number_array.pimag_data);	
		}

		delete x;	
	}
}

void *mxMalloc(size_t n) {
	void* const mem = malloc(n);
	
	if (mem) {	
		pair<set<void*, less<void*> >::const_iterator, bool > p;

		p = mxArray::memory_table.insert(mem);

		#ifdef MEXWRAP_DEBUG
			fprintf(stderr, "Allocating %d bytes at %p\n", n, mem);
			fprintf(stderr, "Memory table size : %d\n", mxArray::memory_table.size());
		#endif			

		if (!p.second) {
			mexErrMsgTxt("Error allocating new memory : pointer is already in use");
		}
	}
	return mem; 
}

void *mxCalloc(size_t n, size_t size) {
	void* const mem = calloc(n, size);
	
	if (mem) {	
		pair<set<void*, less<void*> >::const_iterator, bool > p;

		p = mxArray::memory_table.insert(mem);

		if (!p.second) {
			mexErrMsgTxt("Error allocating new memory : pointer is already in use");
		}	
	}
	
	return mem;
}

void mxFree(void *ptr) {
	if (ptr) {
		const long count = mxArray::memory_table.erase(ptr);

		if (count == 1) {	// only free pointer when it was allocated by one of the mx*** functions
			free(ptr);
		} 
	}
}

void *mxRealloc(void *ptr, size_t size) {
	void* const mem = realloc(ptr, size);
	
	if (mem) {	
		const long count = mxArray::memory_table.erase(ptr);

		if (count != 1) {
			mexErrMsgTxt("Error deleting memory : pointer is already freed");
		}

		pair<set<void*, less<void*> >::const_iterator, bool > p;

		p = mxArray::memory_table.insert(mem);
		if (!p.second) {
			mexErrMsgTxt("Error allocating new memory : pointer is already in use");
		}	
	}	
	
	return mem;
}


/*
 * Return the class (catergory) of data that the array holds.
 */
mxClassID mxGetClassID(const mxArray *x) {
	switch(x->type) {
		case mxDOUBLE_ARRAY:   
			return mxDOUBLE_CLASS;
		break;  
    	case mxSPARSE_ARRAY:
    		return mxSPARSE_CLASS;
		break;
		case mxCHARACTER_ARRAY:
			return mxCHAR_CLASS;
		break; 
    	case mxCELL_ARRAY:
			return mxCELL_CLASS;
		break;
    	case mxSTRUCTURE_ARRAY:
			return mxSTRUCT_CLASS;
		break;
    	case mxFLOAT_ARRAY:
			return mxSINGLE_CLASS;
		break;
    	case mxINT8_ARRAY:
			return mxINT8_CLASS;
		break;
    	case mxUINT8_ARRAY:
			return mxUINT8_CLASS;
		break;
    	case mxINT16_ARRAY:
			return mxINT16_CLASS;
		break;
    	case mxUINT16_ARRAY:
			return mxUINT16_CLASS;
		break;
    	case mxINT32_ARRAY:
			return mxINT32_CLASS;
		break;
    	case mxUINT32_ARRAY:
			return mxUINT32_CLASS;
		break;
    	case mxINT64_ARRAY:
			return mxINT64_CLASS;
		break;    
    	case mxUINT64_ARRAY:
			return mxUINT64_CLASS;
		break;   
    	case mxOBJECT_ARRAY:
			return mxOBJECT_CLASS;
		break;
    	case mxUNKNOWN_ARRAY:
			return mxUNKNOWN_CLASS;
		break;
		default:
			return mxUNKNOWN_CLASS;
		break;
	}
}

int mxGetString(const mxArray *x, char *buf, int buflen) {
	if (mxIsChar(x)) {
		long maxlen = (buflen-1) > mxGetN(x) ? mxGetN(x) : (buflen-1);
		
		mxChar* y = (mxChar*) x->data.number_array.pdata; 
		
		for (long i=0; i < maxlen; i++) {
			buf[i] = (char) y[mxGetM(x)*i];
		}
		
		buf[maxlen] = '\0';		

		return 0;	// SUCCESS
	}
	
	return 1;	// FAILURE
}

/* 
 * Get pointer to array name.
 */
const char *mxGetName(const mxArray *x) {
	return x->name;
}

/* 
 * Set array name.  This routine copies the string pointed to by s
 * into the mxMAXNAM length character name field.
 */
void mxSetName(mxArray *x, const char *s) {
	strncpy(x->name, s, mxMAXNAM);
}

int mxGetM(const mxArray* x)
{
	return x->size.matrix_dims.m;
}

int mxGetN(const mxArray* x)
{
	return x->size.matrix_dims.n;
}

double* mxGetPr(const mxArray* x) 
{
	return (double*) x->data.number_array.pdata;
}

// Get imaginary data pointer for numeric array
double *mxGetPi(const mxArray *x) {
	return (double*) x->data.number_array.pimag_data;
}

/*
 * Get pointer to data
 */
void* mxGetData(const mxArray *x) {
	return (void*) &(x->data.number_array.pdata);
}

/*
 * Set pointer to data
 */
void mxSetData(mxArray *x, void  *newdata) {
	x->data.number_array.pdata = newdata;
}

/*
 * Set imaginary data pointer for numeric array
 */
void mxSetImagData(mxArray *x, void *newdata) {
	x->data.number_array.pimag_data = newdata;
}


void mxSetPr(mxArray *x, double* pr) {
	x->data.number_array.pdata = (void*) pr;
}

/*
 * Set imaginary data pointer for numeric array
 */
void mxSetPi(mxArray *x, double  *pi) {
	x->data.number_array.pimag_data = (void*) pi;
}

int mxGetNumberOfDimensions(const mxArray *x) {
	return x->number_of_dims;
}

const int* mxGetDimensions(const mxArray *x) {
	int* dims = (int*) mxMalloc(x->number_of_dims * sizeof(int));
	
	if (dims) {
		if (x->number_of_dims == 2) {
			dims[0] = x->size.matrix_dims.m;
			dims[1] = x->size.matrix_dims.n;
		} else {
			for (int i=0; i < x->number_of_dims; i++)
				dims[i] = x->size.dim_array[i];
		}
	}
	
	return dims;
}

int mxSetDimensions(mxArray *x, const int *dims, int ndims) {
	x->number_of_dims = ndims;
	
	if (x->number_of_dims == 2) {
		x->size.matrix_dims.m = dims[0]; 
		x->size.matrix_dims.n = dims[1];
	} else {
		x->size.dim_array = (int*) mxMalloc(x->number_of_dims * sizeof(int));
		
		if (x->size.dim_array == 0)
			return 1;	// failure
		
		for (int i=0; i < x->number_of_dims; i++)
			x->size.dim_array[i] = dims[i];
	}
	
	return 0;	// success
}

// Determine whether the specified array contains numeric (as opposed 
// to cell or struct) data.
bool mxIsNumeric(const mxArray *x) {
	switch(mxGetClassID(x)) {
		case mxDOUBLE_CLASS:   
			return true;
		break;  
    	case mxSPARSE_CLASS:
    		return true;
		break;
		case mxCHAR_CLASS:
			return false;
		break; 
    	case mxCELL_CLASS:
			return false;
		break;
    	case mxSTRUCT_CLASS:
			return false;
		break;
    	case mxSINGLE_CLASS:
			return true;
		break;
    	case mxINT8_CLASS:
			return true;
		break;
    	case mxUINT8_CLASS:
			return true;
		break;
    	case mxINT16_CLASS:
			return true;
		break;
    	case mxUINT16_CLASS:
			return true;
		break;
    	case mxINT32_CLASS:
			return true;
		break;
    	case mxUINT32_CLASS:
			return true;
		break;
    	case mxINT64_CLASS:
			return true;
		break;    
    	case mxUINT64_CLASS:
			return true;
		break;   
    	case mxOBJECT_CLASS:
			return false;
		break;
    	case mxUNKNOWN_CLASS:
			return false;
		break;
		default:
			return false;
		break;
	}
}

// Determine whether the given array is a cell array.
bool mxIsCell(const mxArray *x) {
	return mxGetClassID(x) == mxCELL_CLASS;
}

// Determine whether the given array contains character data. 
bool mxIsChar(const mxArray *x) {
	return mxGetClassID(x) == mxCHAR_CLASS;
}

// Determine whether the given array is a sparse (as opposed to full). 
bool mxIsSparse(const mxArray *x) {
	return mxGetClassID(x) == mxSPARSE_CLASS;
}

// Determine whether the given array is a structure array.
bool mxIsStruct(const mxArray *x) {
	return mxGetClassID(x) == mxSTRUCT_CLASS;
}

// Determine whether the given array contains complex data.
bool mxIsComplex(const mxArray *x) {
	if ((mxIsNumeric(x)) && (x->data.number_array.pimag_data != 0))
		return true;
	else
		return false;
}

// Determine whether the specified array represents its data as 
// double-precision floating-point numbers.
bool mxIsDouble(const mxArray *x) {
	return mxGetClassID(x) == mxDOUBLE_CLASS;
}

// Determine whether the specified array represents its data as 
// single-precision floating-point numbers.
bool mxIsSingle(const mxArray *x) {
	return mxGetClassID(x) == mxSINGLE_CLASS;
}

// Determine whether the given array's logical flag is on.
bool mxIsLogical(const mxArray *x) {
	return false;
}

// Determine whether the specified array represents its data as 
// signed 8-bit integers.
bool mxIsInt8(const mxArray *x) {
	return mxGetClassID(x) == mxINT8_CLASS;
}

// Determine whether the specified array represents its data as 
// unsigned 8-bit integers.
bool mxIsUint8(const mxArray *x) {
	return mxGetClassID(x) == mxUINT8_CLASS;
}

// Determine whether the specified array represents its data as 
// signed 16-bit integers.
bool mxIsInt16(const mxArray *x) {
	return mxGetClassID(x) == mxINT16_CLASS;
}

// Determine whether the specified array represents its data as 
// unsigned 16-bit integers.
bool mxIsUint16(const mxArray *x) {
	return mxGetClassID(x) == mxUINT16_CLASS;
}

// Determine whether the specified array represents its data as 
// signed 32-bit integers.
bool mxIsInt32(const mxArray *x) {
	return mxGetClassID(x) == mxINT32_CLASS;
}

// Determine whether the specified array represents its data as 
// unsigned 32-bit integers.
bool mxIsUint32(const mxArray *x) {
	return mxGetClassID(x) == mxUINT32_CLASS;
}

/*
 * Get 8 bits of user data stored in the mxArray header.  NOTE: This state
 * of these bits are not guaranteed to be preserved after API function
 * calls.
 */
int mxGetUserBits(const mxArray	*x) {
	return 0;
}

/*
 * Set 8 bits of user data stored in the mxArray header. NOTE: This state
 * of these bits are not guaranteed to be preserved after API function
 * calls.
 */ 
void mxSetUserBits(mxArray *x, int value) {
	
}
	
/*
 * Get the real component of the specified array's first data element.
 */
double mxGetScalar(const mxArray *x) {
	return *((double*) x->data.number_array.pdata);
}

/*
 * Specify that the data in an array is to be treated as Boolean data.
 */
void mxSetLogical(mxArray *x) {

}

/*
 * Specify that the data in an array is to be treated as numerical
 * (as opposed to Boolean) data. 
 */
void mxClearLogical(mxArray *x) {

}

// define static tables of mxArrays and memory
set<mxArray_ptr, less<mxArray_ptr> > mxArray::array_table;
set<void*, less<void*> > mxArray::memory_table;	

