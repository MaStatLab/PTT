#include <R.h>
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <R_ext/Rdynload.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <R_ext/Linpack.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h> 
#include <string.h>
#include <float.h>
#include <list>
#include <map>
#define RLIB
#include "modeltree_node.cpp"
#include "gbt.cpp"

using namespace std;

extern "C" {

  SEXP markov_apt_kD_C(SEXP x_mat, SEXP n_pred, SEXP max_dim, SEXP Rrho0, SEXP Rrho0_mode, SEXP Rtran_mode, SEXP lognu_lowerbound, SEXP lognu_upperbound, SEXP n_grid, SEXP n_s, SEXP beta,SEXP nsamples, SEXP x_new_mat);

};

SEXP markov_apt_kD_C(SEXP x_mat, SEXP n_pred, SEXP max_dim, SEXP Rrho0, SEXP Rrho0_mode, SEXP Rtran_mode, SEXP lognu_lowerbound, SEXP lognu_upperbound, SEXP n_grid, SEXP n_s, SEXP beta, SEXP nsamples, SEXP x_new_mat) {
  
  SEXP RXwork = PROTECT(duplicate(x_mat));
  int nProtected = 1;
  SEXP RXnew = PROTECT(duplicate(x_new_mat)); ++nProtected;


  int *r_part_points;
  double *r_predictive_densities;

  double *Xwork; // observed data
  double *Xnew; // x values at which densities are to be estimated
  uint *Xwork_uint;
  uint *Xnew_uint;
  int nobs, nobs_new, p, k, n_post_sample;
  int n_part, n_part_sample;
  vector< vector<ushort> > part_points_vec;
  

  vector<double> predictive_densities_vec;

  SEXP logrho;
  SEXP logphi;
  SEXP part_points;    
  SEXP predictive_densities;
  SEXP ans;
  

  n_post_sample = INTEGER(nsamples)[0];

  vector< vector< vector<ushort> > > part_points_vec_post_samples(n_post_sample);
  vector< vector< vector<double> > > nu_and_prob_vec_post_samples(n_post_sample);
  

  SEXP part_points_post_samples;
  SEXP part_points_post_sample;  
  double *r_part_points_sample;

 
  if (n_post_sample > 0) {
    PROTECT(ans = allocVector(VECSXP,5));
  }

  else {
    PROTECT(ans = allocVector(VECSXP,4)); 
  }
  
  ++nProtected;

  p = INTEGER(n_pred)[0];
  k = INTEGER(max_dim)[0];
  nobs = LENGTH(x_mat)/p;
  nobs_new = LENGTH(x_new_mat)/p;
  
  Xwork = REAL(RXwork);
  Xwork_uint = new uint[LENGTH(x_mat)];

  Xnew = REAL(RXnew);
  Xnew_uint = new uint[LENGTH(x_new_mat)];

  for (int i = 0; i < LENGTH(x_mat); i++) {
    Xwork_uint[i] = convert_to_inverse_base_2(Xwork[i],k);
  }

  for (int i = 0; i < LENGTH(x_new_mat); i++) {
    Xnew_uint[i] = convert_to_inverse_base_2(Xnew[i],k);
  }
  

  GBT my_gbt(Xwork_uint,nobs,k,p,REAL(Rrho0)[0],INTEGER(Rrho0_mode)[0],INTEGER(Rtran_mode)[0],REAL(lognu_lowerbound)[0], REAL(lognu_upperbound)[0], INTEGER(n_grid)[0],INTEGER(n_s)[0],REAL(beta)[0]);

  my_gbt.update();


  part_points_vec = my_gbt.find_part();
  PROTECT(part_points = allocMatrix(INTSXP, part_points_vec.size(),2*p+2)); ++nProtected;
  r_part_points = INTEGER(part_points);

  n_part = part_points_vec.size();

  for(int i = 0; i < n_part; i++) {    
    for (int j = 0; j < 2*p+2; j++) {
	r_part_points[i+j*n_part] = part_points_vec[i][j];
    }  
  } 

  PROTECT(logphi = allocVector(REALSXP,1)); ++nProtected;
  PROTECT(logrho = allocVector(REALSXP,1)); ++nProtected;
 
  REAL(logphi)[0] = my_gbt.get_root_logphi();
  REAL(logrho)[0] = my_gbt.get_root_logrho();



  predictive_densities_vec = my_gbt.compute_predictive_density(Xnew_uint,nobs_new);

  PROTECT(predictive_densities = allocVector(REALSXP, nobs_new)); ++nProtected;
  r_predictive_densities = REAL(predictive_densities);

  for (int i = 0; i < nobs_new; i++) {
    r_predictive_densities[i] = predictive_densities_vec[i];
  }

  SET_VECTOR_ELT(ans, 0, logrho);
  SET_VECTOR_ELT(ans, 1, logphi);
  SET_VECTOR_ELT(ans, 2, part_points);
  SET_VECTOR_ELT(ans, 3, predictive_densities);
  
  if (n_post_sample > 0) {
  
 
    PROTECT(part_points_post_samples = allocVector(VECSXP,n_post_sample)); nProtected++; 
   
    for (int i = 0; i < n_post_sample; i++) {

      my_gbt.sample();
      n_part_sample = my_gbt.find_sample_part(part_points_vec_post_samples[i],nu_and_prob_vec_post_samples[i]);


      PROTECT(part_points_post_sample = allocMatrix(REALSXP, n_part_sample,2*p+1+2)); ++nProtected;
      r_part_points_sample = REAL(part_points_post_sample);

 
      for(int l = 0; l < n_part_sample; l++) {    

	for (int j = 0; j < 2*p+1; j++) {
	  r_part_points_sample[l+j*n_part_sample] = part_points_vec_post_samples[i][l][j];
	}  
	for (int j = 0; j < 2; j++) {
	  r_part_points_sample[l+(2*p+1+j)*n_part_sample] = nu_and_prob_vec_post_samples[i][l][j];

	}
      } 

      SET_VECTOR_ELT(part_points_post_samples,i,part_points_post_sample);

    }
    SET_VECTOR_ELT(ans,4,part_points_post_samples);

  }


  UNPROTECT(nProtected);
  
  delete [] Xwork_uint;
  delete [] Xnew_uint;

  return ans;

};


