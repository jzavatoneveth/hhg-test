/*==========================================================
 * inversions.c
 *
 * Computes the inversion vector of a permutation
 *
 * The calling syntax is:
 *
 *		[ invs ] = inversions(p)
 *
 * This is a MEX-file for MATLAB.
 *
 *
 * Copyright (c) 2018 Jacob Zavatone-Veth, MIT License
 *
 *========================================================*/

#include "mex.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

/*
O(n*log(n)) inversion counting without merge sort, thanks to
https://stackoverflow.com/questions/48977484/inversion-count-algorithm-without-using-merge-sort-c
*/
void Inversions(uint64_t *a, uint64_t *b, uint64_t dim)
{
  /* Initialize a counter of the number of zero values in the array a */
  uint64_t z = 0;

  /* Allocate a temporary container */
  size_t n = (dim+2) * sizeof(uint64_t);
  uint64_t *c = (uint64_t*) malloc(n);

  /* Iterate until all values of a are zero */
  do {
    /* Set all elements of c to zero at each iteration */
    memset(c, 0, n);

    /* Iterate over elements */
    for (int i = 0; i < dim; i++) {
      if (a[i] % 2 == 0) {
        b[i] += c[a[i]+1];
      } else{
        c[a[i]] += 1;
      }
    }

    z = 0; /* Reset the zero counter */

    /* Iterate over elements */
    for (int k = 0; k < dim; k++){
      /* Divide each element of a by 2 (as integers) */
      a[k] /= 2;
      /* Accumulate the counter of zeros */
      if (a[k] == 0) { z++; }
    }
  } while (z != dim);

  free(c); /* Free the temporary container */
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  uint64_t *source;
  uint64_t *permutation;
  uint64_t dim;
  uint64_t *inversion_count;

  /* Validate the number of inputs and outputs */
  if (nrhs != 1) {mexErrMsgIdAndTxt("inversions:nrhs", "One input required.");}
  if (nlhs != 1) {mexErrMsgIdAndTxt("inversions:nlhs", "One output required.");}

  /* Validate the types of the inputs */
  if (!mxIsUint64(prhs[0])) {
    mexErrMsgIdAndTxt("nversions:intype", "Input must be of type uint64.");
  }

  /* Validate the dimensions of the inputs */
  if (mxGetN(prhs[0]) != 1) {
    mexErrMsgIdAndTxt("inversions:incolumn", "Input must be a column vector.");
  }
  if (mxGetM(prhs[0]) <= 1) {
    mexErrMsgIdAndTxt("inversions:inshape", "Input must have at least two elements.");
  }

  /* Get a pointer to a deep copy of the data in the input matrix */
  permutation = (uint64_t*) mxGetData(mxDuplicateArray(prhs[0]));

  /* Get the dimensions of the input vectors */
  dim = (uint64_t) mxGetM(prhs[0]);

  /* Create the output matrix */
  plhs[0] = mxCreateNumericMatrix((mwSize)dim, 1, mxUINT64_CLASS, mxREAL);

  /* Get a pointer to the output matrix */
  inversion_count = (uint64_t*) mxGetData(plhs[0]);

  /* Call the computational routine */
  Inversions(permutation, inversion_count, dim);
}
