/*
 *   This is an implementation of RGBM algorithm for Gene Regulatory Network
 *   inference from any type of expression data built on top of ennet package, in form of an R package.
 *   Copyright (C) 2016  Raghvendra Mall 
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program, see LICENSE.
 */
// The random number generation part of this code was handled by Khalid Kunji
#include <stdio.h>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <cfloat>
#include "regression_stump.h"
#include "pcg_basic.h"
#include <math.h>
#include <limits>

//#include <R.h>

inline double signnum_c(double x) {
  if (x > 0.0) return 1.0;
  if (x < 0.0) return -1.0;
  return x;
}

inline int compvar(const void *a, const void *b)
{
  double result = *(double *)a - *(double *)b;
  if (result > 0 )
    return 1;
  if (result < 0)
    return -1;
  return 0;
}

inline double calculate_median(double *y_temp, int start, int end)
{
  double med_y = 0.0;
  int mid,left,right;
  int index=0,N_t;
  
  // No of elements in y_temp
  //N = sizeof(y_temp)/sizeof(double);
  
  // No of elements to calculate median over
  N_t = end-start;
  double *y = new double[N_t];
  for (int i=start;i<end;i++)
  {
    y[index] = y_temp[i];
    index++;
  }
  
  //Sort and select median
  qsort(y,N_t,sizeof(double),compvar);
  if (N_t%2==1)
  {
    mid = N_t/2;
    med_y = y[mid];
  }
  else{
    right = N_t/2;
    left = right-1;
    med_y = 1.0*(y[left] + y[right])/2;
  }
  delete[] y;
  return(med_y);
}

inline int compare(const void *a, const void *b) {
	double result = **(double **) a - **(double**) b;
	if (result > 0)
		return 1;
	if (result < 0)
		return -1;
	return 0;
}

static pcg32_random_t pcg32_global = PCG32_INITIALIZER;

const Model train_regression_stump(const int N, const int P, const double *x,
		const double *y, const double col_sampling_rate,
		const double row_sampling_rate, const int lf, const int M, const double nu) 
{
	/*
	 * reset random seed
	 */
	//srand(time(NULL));
	// In this version of the code, we'll use the global rng, rather than a
	// local one.
	
	// You should *always* seed the RNG.  The usual time to do it is the
	// point in time when you create RNG (typically at the beginning of the
	// program).
	//
	// pcg32_srandom_r takes two 64-bit constants (the initial state, and the
	// rng sequence selector; rngs with different sequence selectors will
	// *never* have random sequences that coincide, at all) - the code below
	// shows three possible ways to do so.
	bool seed_given = false;
  uint64_t seed = 123;
  if (seed_given) {
    pcg32_srandom(42u, seed);
  } else {
    // Seed with a fixed constant
    pcg32_srandom(42u, 54u);
  }
	/*
	 * allocate memory for arrays ...
	 */
	const double **x_to_sort = (const double **) calloc(N,
			sizeof(const double*));
	int *x_sorted_index = (int *) calloc(N * P, sizeof(int));
	double *F = (double*) calloc(N, sizeof(double));
	double *h = (double *) calloc(N, sizeof(double));
	double *r = (double *) calloc(N, sizeof(double));
	bool *col_w = (bool *) calloc(P, sizeof(bool));
	double *row_w = (double *) calloc(N, sizeof(double));
	int *active_sorted_rows = (int *) calloc(N, sizeof(int));
	double *I = (double *) calloc(P, sizeof(double));
  double *y_temp = (double *)calloc(N,sizeof(double));
  double *r_temp = (double *)calloc(N,sizeof(double));
	double *left_temp = (double *)calloc(N,sizeof(double));
	double *right_temp = (double *)calloc(N,sizeof(double));
	int *index_map = (int*)calloc(N,sizeof(int));
	
	/*
	 * ... and scalars ...
	 */
	int P_unique_inbag;
	int P_already_bagged;
	int N_of_selections;
	int N_active_rows;
	double total_w;
	int row;
	int next_row;
	int col;
	int iteration;
	int s;
	int s_i;
	int s_left;
	int s_right;

	double local_importance;
	double local_threshold;
	double best_inbag_importance;
	double best_inbag_threshold = 0.0;
	int best_inbag_column = 0;
	double sum_r;
	double ib_w_l;
	double ib_w_r;
	double ib_sum_l;
	double ib_sum_r;
	double ib_gamma_l;
	double ib_gamma_r;
	double best_gamma_l;
	double best_gamma_r;
	double best_rho_l;
	double best_rho_r;
	const double * origin;
	double total_importance;
  int t,k_t;
  
	/*
	 * ... and vectors for returned object
	 */
	
	Model result = Model(M, P);
  for (int i=0; i<N; i++)
  {
    y_temp[i] = y[i];
  }
	/*
	 * Initialize importance I with zeros
	 */
	for (col = 0; col < P; col++) {
		I[col] = 0;
	}

	/*
	 * get indices of sorted columns of x
	 */
	for (col = 0; col < P; col++) {
		for (row = 0; row < N; row++) {
			x_to_sort[row] = x + col * N + row;
		}
		origin = x_to_sort[0];
		qsort(x_to_sort, N, sizeof(double*), compare);
		for (row = 0; row < N; row++) {
			x_sorted_index[col * N + row] = x_to_sort[row] - origin;
		}
	}

	double avg_y=0.0,med_y=0.0, overall_r=0.0, prev_r = 0.0, avg_r = 0.0;
  double precision = std::numeric_limits<float>::denorm_min();
	if (lf==1){
	  /*
	    * calculate initial prediction
	    * arithmetic mean of y for LS
	  */
	  for (row = 0; row < N; row++) {
		  avg_y += y[row];
	  }
	  avg_y /= N;
	}
	else if (lf==2){
	  /*
	   *  Calculate initial prediction
	   *  median of y value for LAD
	   */
	  med_y = calculate_median(y_temp,0,(int)N);
	}
	
	for (row = 0; row < N; row++) {
	  if (lf==1){
	    F[row] = avg_y;
	    }
	  else if (lf==2){
	    F[row] = med_y;
	  }
	}
	
	
	/*
	 * Start the main part of the algorithm
	 */
	for (iteration = 0; iteration < M; iteration++) 
	{
	  /*
	   * calculate working response for all observations
	   * according to squared-error loss or least absolute deviation loss function
	   */
	  overall_r = 0.0;
	  for (row = 0; row < N; row++) 
	  {
	    /*
	     * For Least-Squares
	     * L(y,f) = 0.5*(f-y)^2
	     * -dL/df = y-f
	     */
	    if (lf==1) { 
	      r[row] = y[row] - F[row]; 
	      overall_r = overall_r + abs(r[row]);
	    }
	    else if (lf==2) { 
	      r[row] = y[row]-F[row];
	      r_temp[row] = signnum_c(r[row]);  //Keep signed residuals in a temp variable for LAD
	      overall_r = overall_r + abs(r[row]);
	    }
	  }
	  
	  /*
	   * select sample of columns without replacement
	   * no replicates
	   */
	  if (col_sampling_rate > 0 && col_sampling_rate <= 1) {
	    P_unique_inbag = ceil(col_sampling_rate * P);
	    P_already_bagged = 0;
	    for (s = 0; s < P; s++) {
	      if (((float)ldexp((double)pcg32_random_r(&pcg32_global), -32)) * (P - s)
              < P_unique_inbag - P_already_bagged) {
        //   if ((1.0 * rand() / RAND_MAX) * (P - s)
        //     < P_unique_inbag - P_already_bagged) {
         
	        col_w[s] = true;
	        P_already_bagged++;
	      } else {
	        col_w[s] = false;
	      }
	      if (P_already_bagged >= P_unique_inbag) {
	        s++;
	        break;
	      }
	    }
	    // the remainder is not in the bag
	    for (; s < P; s++) {
	      col_w[s] = false;
	    }
	  }
	  
	  /*
	   * select sample of rows with replacement
	   * replicates are marked with row_w = 2, 3, etc
	   */
	  if (row_sampling_rate > 0) 
	  {
	    N_of_selections = ceil(row_sampling_rate * N);
	    for (s = 0; s < N; s++) {
	      row_w[s] = 0.0;
	    }
	    for (s = 0; s < N_of_selections; s++) {
	      row_w[(int)pcg32_boundedrand(N - .00000000001)] += 1.0;
	      //row_w[rand() % N] += 1.0;
	    }
	  }
	  
	  /*
	   * Calculate INBAG sum of all pseudo-residuals
	   */
	  sum_r = 0.0;
	  total_w = 0.0;
	  for (row = 0; row < N; row++) 
	  {
	    if (lf==1) 
	    {
	      // For LS calculate sum of contributions of residuals
	      sum_r += row_w[row] * r[row];
	    }
	    else if (lf==2) {
	      // For LAD calculate sum of contributions of signed residuals
	      sum_r += row_w[row] * r_temp[row];
	    }
	    total_w += row_w[row];
	  }
	  
	  best_inbag_importance = 0.0;
	  if (lf==1)
	  {
	    best_gamma_l = 0.0;
	    best_gamma_r = 0.0;
	  }
	  if (lf==2)
	  {
	    best_gamma_l = 0.0;
	    best_gamma_r = 0.0;
	    best_rho_l = 0.0;
	    best_rho_r = 0.0;
	  }
	  
	  /*
	   * go through all the columns in sample
	   * begin the search
	   */
	  
	  for (col = 0; col < P; col++) 
	  {
	    if (col_w[col])
	    {
	      /*
	       * Which rows are active
	       * with respect to the order
	       * of x
	       *
	       */
	      s_i = 0;
	      for (s = 0; s < N; s++)
	      {
	        if (row_w[x_sorted_index[col * N + s]] > 0.0) {
	          active_sorted_rows[s_i] = s;
	          s_i++;
	        }
	      }
	      N_active_rows = s_i;
	      
	      /*
	       * Initialize variables
	       */
	      ib_sum_l = 0.0;
	      ib_sum_r = sum_r;
	      ib_w_l = 0.0;
	      ib_w_r = total_w;
	      

	      /*
	       * go through all the thresholds
	       */
	      for (s = 0; s < N_active_rows - 1; s++) 
	      {
	        row = x_sorted_index[col * N + active_sorted_rows[s]];
	        next_row = x_sorted_index[col * N
	          + active_sorted_rows[s + 1]];
	        
	        if (lf==1){
	          //Residuals put in binary tree for LS
	          ib_sum_l += row_w[row] * r[row];
	          ib_sum_r -= row_w[row] * r[row];
	          ib_w_l += row_w[row];
	          ib_w_r -= row_w[row];
	        }
	        else if (lf==2)
	        {
	          ib_w_l += row_w[row];
	          ib_w_r -= row_w[row];
	          //Signed Residuals put in binary tree for LAD
  	        ib_sum_l += row_w[row] * r_temp[row];
	          ib_sum_r -= row_w[row] * r_temp[row];
	        }
	        
	        /*
	         * test if it's a valid split point
	         * OK: 1 2 | 3 4
	         * WRONG: 1 1 | 1 3
	         */
	        if (x[col * N + row] < x[col * N + next_row]) 
	        {
	          /*
	           * it's a valid split
	           */
	          local_threshold = (x[col * N + row]
                                + x[col * N + next_row]) / 2.0;
	          
	            ib_gamma_l = ib_sum_l / ib_w_l;
	            ib_gamma_r = ib_sum_r / ib_w_r;
	          
	          local_importance = ((ib_w_l * ib_w_r) / total_w )*(ib_gamma_l - ib_gamma_r)*(ib_gamma_l - ib_gamma_r);
	          
	          if (local_importance > best_inbag_importance) 
	          {
	            best_inbag_column = col;
	            best_inbag_importance = local_importance;
	            best_inbag_threshold = local_threshold;
	            // Got the best gamma_l and gamma_r for LS
	            if (lf==1){
	              best_gamma_l = ib_gamma_l;
	              best_gamma_r = ib_gamma_r;
	            }
	            else if (lf==2){
	              best_rho_l = ib_gamma_l;
	              best_rho_r = ib_gamma_r;
	            }
            }
	        }
	      }
	    }
	  }
	  if (best_inbag_importance > 0.0) {
	    /*
	     * Update relative influence
	     */
	    I[best_inbag_column] += best_inbag_importance;
	  }
	  
	  if (lf==2)
	  {
	    
	    // Find the list of active_sorted_rows w.r.t. sorted_index of x
	    s_i = 0;
	    for (s = 0; s < N; s++) {
	       if (row_w[x_sorted_index[best_inbag_column * N + s]] > 0.0) {
	        active_sorted_rows[s_i] = s;
	        s_i++;
	      }
	    }
	    N_active_rows = s_i;
	   
	    s_left = 0;
	    s_right = 0;
	    for (s = 0; s < N_active_rows - 1; s++)
	    {
	      row = x_sorted_index[best_inbag_column * N + active_sorted_rows[s]];
	      if (x[best_inbag_column * N + row] < best_inbag_threshold)
	      {
	        for (t=1; t <= (int)row_w[row]; t++)
	        {
	          left_temp[s_left] = r[row];
	          s_left++;
          }
	      }
	      else
	      {
	        for (t=1; t <= (int)row_w[row]; t++)
	        {
	          right_temp[s_right] = r[row];
	          s_right++;
	        }
	      }
	    }
	    best_gamma_l = calculate_median(left_temp,0,s_left);
	    best_gamma_r = calculate_median(right_temp,0,s_right);
	  }
	  
	  /*
	   * update f using all observations
	   * nu * h(x) for LS
	   */
	  for (row = 0; row < N; row++) {
	    if (x[best_inbag_column * N + row] < best_inbag_threshold) {
	      F[row] += nu * best_gamma_l;
	    } else {
	      F[row] += nu * best_gamma_r;
	    }
	  }
		  
	 /*
	  * Prepare result
	  * part 1
	  */
		result.setFeatSplitI(iteration, best_inbag_column);
		result.setFeatSplitT(iteration, best_inbag_threshold);
		result.setGammaL(iteration, best_gamma_l);
		result.setGammaR(iteration, best_gamma_r);
		
		avg_r = overall_r/N;
    if (avg_r <= precision && abs(avg_r-prev_r) <= precision)
		{
		  //printf("Print Number of iterations performed M = %d\n",iteration);
		  break;
		}
		prev_r = avg_r;
	}

	/*
	 * Scale importance
	 */
	total_importance = 0.0;
	for (col = 0; col < P; col++) {
		total_importance += I[col];
	}
	for (col = 0; col < P; col++) {
	  if(total_importance>0.0)
	  {
		  I[col] = I[col] / total_importance;
	  }
	}

	/*
	 * Preapre result
	 * part 2
	 */
	if (lf==1)
	{
	  result.setF0(avg_y);
	}
	else if (lf==2)
	{
	  result.setF0(med_y);
	}
	result.setNu(nu);
	for (col = 0; col < P; col++) {
		result.setImportance(col, I[col]);
	}

	/*
	 * Deallocate all variables allocated dynamically
	 */
	free(x_to_sort);
	free(x_sorted_index);
	free(F);
	free(h);
	free(r);
	free(col_w);
	free(row_w);
	free(active_sorted_rows);
	free(I);
	free(r_temp);
	free(y_temp);
	/*
	 * Return a copy of model
	 */
	return result;
}

