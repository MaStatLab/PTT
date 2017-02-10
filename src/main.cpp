#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "gbt.h"
#include "helpers.h"


using namespace Rcpp;
using namespace arma;
using namespace std;


Rcpp::List fitPTTcpp(
    arma:mat X,
    int k = 10, // the maximum resolution level (or max_dim in R)
    double rho0 = 0.5, // baseline prior stopping probability
    int rho0_mode = 0, // how the prior stopping probabilities are specified
    int tran_mode,
    double lognu_lowerbound = -1, // lowerbound for log shrinkage parameter
    double lognu_upperbound = 4, // upperbound for log shrinkage parameter
    int n_grid = 5, // number of grid points for computing the marginal likelihood
    int n_s = 5, // number of non-stopping shrinkages states
    double beta, // hyparameter for the transition probabilities
    int n_post_samples, // number of posterior samples
    arma:mat Xnew, // points at which predictive densities are calculated
    arma:mat Omega // the joint sample space as a p-dimensional hypercube
    )
{


// number of dimensions
    int p = X.n_cols;
// number of observations
    int nobs = X.n_rows;
    int nobs_new = Xnew.n_rows;

    // constants for rescaling the observations according to the sample space
    vec a = 1.0 / (Omega.col(1) - Omega.col(0));
    vec b = - Omega.col(0) % a;

// Recode the the observations using 0/1 binary indicators
    vector<unsigned int> X_binary(nobs*p);
    for (int i = 0; i < nobs; i++) {
      for (int d = 0; d < p; d++) {
        X_binary(i,d) = convert_to_inverse_base_2(a(d)*X(i,d)+b(d), k+1);
        // X_binary[i*p+d] = convert_to_inverse_base_2(X(i,d), k);
      }
    }

    if (nobs_new > 0) {
      vector<unsigned int> Xnew_binary(nobs_new*p);
      for (int i = 0; i < nobs_new; i++) {
        for (int d = 0; d < p; d++) {
          Xnew_binary(i,d) = convert_to_inverse_base_2(a(d)*Xnew(i,d)+b(d), k+1);
          // Xnew_binary[i*p+d] = convert_to_inverse_base_2(X(i,d), k);
        }
      }
    }

    GBT my_gbt(X_binary,
               nobs,
               k,
               p,
               rho0,
               rho0_mode,
               tran_mode,
               lognu_lowerbound,
               lognu_upperbound,
               n_grid,
               n_s,
               beta);

    my_gbt.update();

// outputs to return
    vector< vector<unsigned integer> > part_points_vec = my_gbt.find_part(); // the HMAP partition
    double logphi = my_gbt.get_root_logphi(); // log marginal likelihood
    double logrho = my_gbt.get_root_logrho(); // log posterior stopping probability on the root

    if (nobs_new > 0) {
        vector<double> predictive_densities_vec = my_gbt.compute_predictive_density(Xnew_binary);
    }

    else {
        predictive_densities_vec = 0;
    }

    if (n_post_sample > 0) {
      vector< vector< vector<unsigned int> > > part_points_post_samples(n_post_sample);
      vector< vector< vector<double> > > nu_and_prob_post_samples(n_post_sample);


      for (int i = 0; i < n_post_sample; i++) {

        my_gbt.sample();
        n_part_sample = my_gbt.find_sample_part(part_points_post_samples[i],nu_and_prob_post_samples[i]);

      }

      Rcpp::List part_points_post_samples_R = wrap(part_points_post_samples);
      Rcpp::List nu_and_prob_post_samples_R = wrap(nu_and_prob_post_samples);
    }

    else {
      Rcpp:List part_points_post_samples_R = create();
      Rcpp::List nu_and_prob_post_samples_R = create();
    }

    return Rcpp::List::create(
      Rcpp::Named("logrho") = logrho,
      Rcpp::Named("logphi") = logphi,
      Rcpp::Named("predictive_densities") = predictive_densities_vec,
      Rcpp::Named("part_points_post_samples") = part_points_post_samples_R,
      Rcpp::Named("nu_and_prob_post_samples") = nu_and_prob_post_samples_R
    )
}
