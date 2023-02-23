/*mex -lmwlapack -lmwblas -largeArrayDims tangent_orth_omega.c 

TTeMPS Toolbox. 
Michael Steinlechner, 2013-2016
Questions and contact: michael.steinlechner@epfl.ch
BSD 2-clause license, see LICENSE.txt
*/
#define U_SLICE(i,j) &U[i][(ind[d*j+i]-1)*r[i]*r[i+1]]
#define V_SLICE(i,j) &V[i][(ind[d*j+i]-1)*r[i]*r[i+1]]
#define RES_SLICE(i,j) &result[i][(ind[d*j+i]-1)*r[i]*r[i+1]]

#include "mex.h"
#include "blas.h"

/* calling: 
	TTeMPS_tangent_omega_orth( n, r, CU, CV, ind, vals)
*/
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] ) {

	/* input variables */
	double* n_raw;
	double* r_raw;
	double** U;
	double** V;
	double* ind_raw;
	double* vals;
	
	/* output variables */
	double** result;
	mxArray** result_cells;
	
	/* internal variables */
	double* L;
	double* current;
	double* tmp;
	
	mwSignedIndex* n;
	mwSignedIndex* r;
	mwSignedIndex* ind;

	mwSignedIndex numSubsref;
	mwSignedIndex d;
	mwSignedIndex i;
	mwSignedIndex j;
	mwSignedIndex maxrank = 1;
					

	/* get sizes */
	n_raw = mxGetPr( prhs[0] );
	/* get ranks */
	r_raw = mxGetPr( prhs[1] );
	/* get indices */
	ind_raw = mxGetPr( prhs[4] );
	d = mxGetM( prhs[4] );
	numSubsref = mxGetN( prhs[4] );
	vals = mxGetPr( prhs[5] );
	
	n = mxMalloc( d*sizeof(mwSignedIndex) );
	r = mxMalloc( (d+1)*sizeof(mwSignedIndex) );
	ind = mxMalloc( d*numSubsref*sizeof(mwSignedIndex) );
	
	/* Convert index arrays to integer arrays as they get converted
	 * to double arrays when passing to mex.
	 * Converting beforehand allows to avoid multiple typecasts inside the inner loop */
	for( i = 0; i < d; ++i ) {
		n[i] = (mwSignedIndex) n_raw[i];
		r[i] = (mwSignedIndex) r_raw[i];
		if( r[i] > maxrank )
			maxrank = r[i];
	}
	r[d] = (mwSize) r_raw[d];
	
	for( i = 0; i < numSubsref*d; ++i ) {
		ind[i] = (mwSignedIndex) ind_raw[i];
	}
	

	/* Get pointers to the matrices within the cell array */
	U = mxMalloc( d*sizeof(double*) );
	V = mxMalloc( d*sizeof(double*) );
	
	for( i = 0; i < d; ++i ) {
		U[i] = mxGetPr( mxGetCell( prhs[2], i ) );
    	V[i] = mxGetPr( mxGetCell( prhs[3], i ) );
	}
	
	/* Allocate space for output */
	plhs[0] = mxCreateCellMatrix( 1, d );
	result_cells = mxMalloc( d*sizeof(mxArray*) );
	result = mxMalloc( d*sizeof(double*) );
	
	for( i=0; i < d; i++){
		result_cells[i] = mxCreateDoubleMatrix( r[i]*r[i+1]*n[i], 1, mxREAL);
		result[i] = mxGetPr( result_cells[i] );
		mxSetCell( plhs[0], i, result_cells[i] );
	}
	
	/* Allocate enough space for internal intermediate results */
	L = mxMalloc( maxrank*(d-1)*sizeof(double) );
	current = mxMalloc( maxrank*sizeof(double) );
	tmp = mxMalloc( maxrank*sizeof(double) );
	
	/* helper variables for dgemv call */
	char transa = 'T';
	char no_transa = 'N';
	mwSignedIndex ONE_i = 1;
	double ONE_d = 1.0;
	double ZERO_d = 0.0;

	for( j = 0; j < numSubsref; ++j ) {
		
		/* LEFT TO RIGHT FIRST (PRECOMPUTE)*/
		/* ... copy first core to L: */
		dcopy( &r[1], U_SLICE(0,j), &ONE_i, &L[0], &ONE_i );
		/* ... and then multiply with the other cores and store results in columns of L: */
		for( i = 1; i < d-1; ++i ) {
			dgemv( &transa, &r[i], &r[i+1], &ONE_d, 
					U_SLICE(i,j), 
					&r[i],   
				 	&L[maxrank*(i-1)], 
					&ONE_i, &ZERO_d, &L[maxrank*i], &ONE_i);
		}
		
		/* RIGHT TO LEFT PRODUCTS NOW -- USING PRECOMPUTED LEFT SIDES FROM ABOVE */
		/* last dU is without any contributions from the right */
		daxpy( &r[d-1], &vals[j], &L[maxrank*(d-2)], &ONE_i, RES_SLICE(d-1,j), &ONE_i );
	
		/* copy rightmost slice to current variable */
		dcopy( &r[d-1], V_SLICE(d-1,j), &ONE_i, current, &ONE_i );
		
		/* sweep right-left to form dU{i-1} to dU{1} */
		for( i = d-2; i > 0; --i ) {
			/* Outer product update: 
			 * result(:,:,idx) = result(:,:,idx) + L(1:r(i), i-1)*current' */
			dger( &r[i], &r[i+1], &vals[j], 
				  &L[maxrank*(i-1)], &ONE_i, 
				  current, &ONE_i, 
				  RES_SLICE(i,j), &r[i] );
			
			/* update current */	
			dgemv( &no_transa, &r[i], &r[i+1], &ONE_d, 
					V_SLICE(i,j), 
					&r[i],   
				 	current, 
					&ONE_i, &ZERO_d, tmp, &ONE_i);
			/* ... and copy result back to current */
			dcopy( &r[i], tmp, &ONE_i, current, &ONE_i );
		}
		
		/* last core */
		daxpy( &r[1], &vals[j], current, &ONE_i, RES_SLICE(0,j), &ONE_i );
		
	}

	mxFree( n );
	mxFree( r );
	mxFree( ind );
	mxFree( current );
	mxFree( tmp );
	mxFree( U );
	mxFree( V ); 
	mxFree( L );
}
