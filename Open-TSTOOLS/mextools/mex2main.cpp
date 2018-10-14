#include "mexwrap.h"
#include "mex2main.h"

using namespace std;

// First define mexErrMsgTxt
// In this setup, just terminate the application in case mexErrMsgTxt is invoked
// (in Matlab, of course, only the current command is stopped)
void mexErrMsgTxt(const char* text)
{
	fprintf(stderr, "An error ocurred : ");
	fprintf(stderr, text);
	fprintf(stderr, "\n");
	exit(-1);
}

// maximal number of input and output arrays for the mex-function call
#define MAX_ARGUMENTS 256

int main(int argc, char** argv)
{	
	using namespace mex2main;
	long i;
	
	mxArray_ptr plhs[MAX_ARGUMENTS];
	mxArray_ptr prhs[MAX_ARGUMENTS];
	
	// clear arguments
	for (i=0 ; i < MAX_ARGUMENTS; i++) { plhs[i] = prhs[i] = NULL; }
	
	long nlhs, nrhs;
	
	if (argc < 2) {
		fprintf(stderr, "Mex2main-Driver : Two arguments are needed : filename of input_list and desired number of output arrays (1)\n");
		exit(-1);
	}
	if (argc > 2)
		nlhs = atol(argv[2]);
	else
		nlhs = 1;
	
	if ((nlhs < 0) || (nlhs > MAX_ARGUMENTS)) {
		fprintf(stderr, "Bad number of output arrays requested\n");
		exit(-1);
	}
	
	// parse input_list and load all input arrays

	nrhs = parse_input_list(argv[1], prhs);

	// call the mexFunction
	
	mexFunction(nlhs, plhs, nrhs, (const mxArray**) prhs);
	
	// save output arrays
	
	for (i = 0; i < nlhs; i++) {
		char number[64];
		char* outname = (char*) strdup(argv[0]);
		outname = (char*) realloc(outname, strlen(outname) + 256);
		
		strcat(outname, "_out");
		
		sprintf(number, "%ld", i);
		strcat(outname, number);
		
		strcat(outname, ".dat");
		
		// Now find the last slash or backslash to separate 
		
		int lastslash = -1;
		for (int n=0; n < strlen(outname)-1; n++) 
			if ((outname[n] == '\\') || (outname[n] == '/'))
				lastslash = n; 
		
		save_array(plhs[i], outname+lastslash+1);
		
		free(outname);
	}
	
	// free all mxArrays that have not been freed previously
	mxArray::cleanup();
	
	return 0;
}

namespace mex2main {
// Now simulate arguments passing : While Matlab just gives mxArray pointers to objects
// already stored somewhere in memory (the so called Matlab workspace), here we load
// arguments from ASCII files and again store output arguments as such. Therefore, we
// have to write some functions to read matrices etc.

#define MAX_LINE_SIZE 512000

mxArray* readmatrix(FILE* fp)
{
	long m, n;
	char buffer[MAX_LINE_SIZE];
	char* token;
	char minus[] = ", \t\n\r";		// tokens that divide numbers
	long rows = 0;					// counter for number of rows
	long columns = 0;				// counter for number of columns
	
	mxArray* x = NULL;
	
	// first make a dry run to see how many rows and columns we have
	while(fgets(buffer, MAX_LINE_SIZE, fp) != NULL) {
		if (buffer[0] != '#') { 		/* comments start with # */
			if (rows <= 100) {			// check first 100 lines of file to see if number of columns are constant	
				long newcolumns = 0;
				token = strtok(buffer, minus);
				while (token != NULL) {
					newcolumns++;
					token = strtok(NULL, minus);
				}
				if ((columns != 0) && (newcolumns != columns)) {
					fprintf(stderr, "Number of columns must be the same for all rows of input data file\n");
					return NULL;
				}
				columns = newcolumns;	
			}
			rows++;
		}
	} 
	
	if (rows * columns == 0) {
		fprintf(stderr, "There seem to be no data in input file\n");
		return NULL;
	}
	
	// now do real read
	
	rewind(fp);
	
	x = mxCreateDoubleMatrix(rows, columns, mxREAL);
	
	//printf("%d %d   %d %d \n", rows, columns, mxGetM(x), mxGetN(x));
	
	m = 0; 
	while(fgets(buffer, MAX_LINE_SIZE, fp) != NULL) {
		n = 0;
		if (buffer[0] != '#') { 		/* Kommentarzeile mit # einleiten */
			token = strtok(buffer, minus);
			if (m >= rows) {
				fprintf(stderr, "Unexpected change of number of rows while reading input data file\n");
				return NULL;
			}
			if (token != NULL) {
				do 			
				{
					if (n >= columns) {
						fprintf(stderr, "Number of columns must be the same for all rows of input data file (row %ld)\n", m+1);
						return NULL;
					}					
					(mxGetPr(x))[m + n * mxGetM(x)] = atof(token);
					n++;
				} while ((token = strtok(NULL, minus)) != NULL);
				m++;
			}
		}
	} 	
	
	return x;
}	

#define INPUT_LIST_MAX_LINESIZE 32768

long parse_input_list(const char* filename, mxArray** prhs)
{
	long arguments_found = 0;
	
	FILE* fid = fopen(filename, "r");

	if (fid != NULL) {
		char buffer[INPUT_LIST_MAX_LINESIZE];
		char buffer2[INPUT_LIST_MAX_LINESIZE];
		
		char* token;
		char minus[] = "[],; \n\r";
		char* endchar;
		
		while(fgets(buffer, INPUT_LIST_MAX_LINESIZE-2, fid) != NULL) {
			buffer[INPUT_LIST_MAX_LINESIZE-1] = '\0';
			
			if (buffer[0] != '#') { 		/* Kommentarzeile wird mit # einleiten und ignoriert */
				if (buffer[0] == '[') { 		/* Vektor darf in Form [3.234 -324.34 5] gegeben werden */
					long count = 0;
					
#ifdef MEXWRAP_DEBUG
					//printf("%d : %s\n", strlen(buffer), buffer);
#endif 
					strncpy (buffer2, buffer, INPUT_LIST_MAX_LINESIZE-1);

					token = strtok(buffer+1, minus);

					// first count number of elements		
					while (token != NULL) {
#ifdef MEXWRAP_DEBUG					
						printf("%s ", token);  
#endif
						if (!strcmp(token, "]")) {
							break;
						} else {
							count++; 
						}

						token = strtok(NULL, minus); 
					}

					if (count > 0) {
						if (buffer2[strlen(buffer2)-1] == '\'') {	// see if vector is transposed
							prhs[arguments_found] = mxCreateDoubleMatrix(count, 1, mxREAL);	
						} else {
							prhs[arguments_found] = mxCreateDoubleMatrix(1, count, mxREAL);	
						}
#ifdef MEXWRAP_DEBUG	
						printf("Found vector argument of length %ld \n", count);
#endif
						// parse this line again and fill values into vector
				
						count = 0;
						token = strtok(buffer2+1, minus);
						while (token != NULL) {
							//printf("%s ", token); 
							double value = strtod(token, &endchar);	// try to see if we got a number value
							if ((value == 0) && (endchar == token)) {
								continue; //break;
							} else {
								(mxGetPr(prhs[arguments_found]))[count++] = value;
#ifdef MEXWRAP_DEBUG									 
								printf("%lf ", value);
#endif								
							}

							token = strtok(NULL, minus); 	
						}
#ifdef MEXWRAP_DEBUG						
						printf("\n");
#endif
						arguments_found++;
					} else {
						fprintf(stderr, "Error parsing vector argument %s \n", buffer2);
						exit(-1);
					}
				} else {
					token = strtok(buffer, minus);

					while (token != NULL) {

						double value = strtod(token, &endchar);	// try to see if we got a number value

						if ((value == 0) && (endchar == token)) {	
							// we got a string value instead of a number

							const long token_length = strlen(token);

							// first see if this string is enclosed by ''
							if ((token_length > 2) && (token[0] == '\'') && (token[token_length-1] == '\'')) {
								char* tmpstring = (char*) malloc(token_length+4);

								for (long j=0; j < token_length-2; j++) {
									tmpstring[j] = token[j+1];
								}
								
								tmpstring[token_length-2] = '\0';

								prhs[arguments_found] = mxCreateString(tmpstring);

								free((void*) tmpstring);	
#ifdef MEXWRAP_DEBUG									
								printf("Found string argument %s\n", token);
#endif
							
							} 
							else {

								// this string is interpreted as a file name of a data file

								FILE* fp = fopen(token, "r");
								if (fp != NULL) {
									prhs[arguments_found] = readmatrix(fp);
									if (prhs[arguments_found] == NULL) {
										fprintf(stderr, "Error reading data file %s \n", token);
										exit(-1);
									}
#ifdef MEXWRAP_DEBUG									
									printf("Found matrix argument\n");
#endif
									fclose(fp);	
								} else {
									fprintf(stderr, "Error opening data file %s \n", token);
									exit(-1);
								}
							}
						} 
						else {	// we got a single number
							prhs[arguments_found] = mxCreateDoubleMatrix(1, 1, mxREAL);		
							(mxGetPr(prhs[arguments_found]))[0] = value;
#ifdef MEXWRAP_DEBUG								
							printf("Found scalar argument\n");
#endif							
						}
						
						arguments_found++;
						token = strtok(NULL, minus);
					}
				}
			}	
		}	  

		fclose(fid);
	} else {
		fprintf(stderr, "Error opening argument list file %s \n", filename);
		exit(-1);
	}
	
	return arguments_found;
}

void save_array(const mxArray* x, const char* filename) {
	if (x != NULL) {	
		FILE* fid = fopen(filename, "w");

		for (long m = 0; m < mxGetM(x); m++) {
			for (long n = 0; n < mxGetN(x); n++) {
				fprintf(fid, "%f ", (mxGetPr(x))[m + n * mxGetM(x)]);
			}
			fprintf(fid, "\n");
		}

		fclose(fid);
	} else {
		fprintf(stderr, "Warning: One or more output arguments not assigned\n");
	}
}

}
