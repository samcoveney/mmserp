# Stan Model for ERP
# -----------------------------------------------
# Author: Sam Coveney
# Date:   31-10-2021
#
# -----------------------------------------------

import pystan
import pickle

print('Using pystan {0}'.format(pystan.__version__))

stancode = \
"""
functions {
  vector spectral_density(vector Q, real rho, real alpha) {

    // NOTE: instead of supplying w = sqrt(Q) to SD function and using w^2 in SD expression, I just supply Q and use w^2 = Q
    // NOTE: including the GP amplitude in this function
    // NOTE: rbf spectral density included power Dimension / 2.0, but our dimension is 2, so ignoring this term
  
    int M = rows(Q);
    vector[M] SD = alpha^2 * (2.0 * pi() * rho^2) * exp( -2 * pi()^2 * rho^2 * Q );
    return SD;
  }
  real top_hat(vector y, vector psi, int N, real deltaS2, int num) {
    // use approximation of top-hat for log-likelihood function
    // y: lower bracket of the interval containing ERP
    // deltaS2: width of ERP bracket
    // num: number of sub-intervals for approximation of top-hat

    real llh;
    vector [2*N] llh_terms;
    real v;
    real normFac;

    vector [num] tmp;

    v = (deltaS2 / num)^2;
    normFac = 1.0 / (num*sqrt(2.0 * pi() * v));

    // loop over observations
    for (n in 1:(2*N)) {
      
      // centres of sub-intervals for the top-hat
      for (idx in 1:num) {
        tmp[idx] = -0.5 * ( psi[n] - (y[n] + (deltaS2/num)*(idx - 0.5) ) )^2 / v;
      }

      llh_terms[n] = log_sum_exp( tmp );

    }

    llh = 2*N*log(normFac) + sum(llh_terms);
    return llh;

  }
}
data {
    // measuremments
    // -------------
    int<lower=0> N;      // num. of observations locations
    vector[2*N] y;       // observations; two types of ERP per observation location

    real<lower=0> deltaS2; // S2 resolution for ERP)
    int<lower=0> num;      // number of subintervals for top-hat approx.

    // theta (parameter) eigen
    // -----------------------
    int<lower=0> M;    // number of eigen functions
    matrix[N,M] eigen; // Laplacian eigenvectors at measurement locations
    vector[M] Q;       // Laplacian eigenvectors
    
    // surrogate for ERP(theta1, theta2)
    // ------------------------------------
    int<lower=0> S;       // num. coefficients for surrogate models
    vector[S] surrogate1; // coefficients for ERPS1 surrogate model
    vector[S] surrogate2; // coefficients for ERPS2 surrogate model
}
parameters {
    vector[M] beta1; // for tau_out GP
    vector[M] beta2; // for APD_max GP

    // hyperparameters
    real mean1;
    real mean2;
    real<lower=0> rho1;
    real<lower=0> alpha1;
    real<lower=0> rho2;
    real<lower=0> alpha2;
}
model {
  vector[N] theta1;
  vector[N] theta2;
  vector[2*N] psi;
  matrix[N,S] sbasis;
  vector[N] tmp; // helper vector
  {

    theta1 = eigen * (beta1 .* sqrt(spectral_density(Q, rho1, alpha1))) + mean1;
    theta2 = eigen * (beta2 .* sqrt(spectral_density(Q, rho2, alpha2))) + mean2;

     // surrogate polynomial eigen i.e. (1, theta1, theta2, theta1*2, ... etc)
    for (n in 1:N) sbasis[n,1] = 1.0;
    tmp = theta1; for (n in 1:N) sbasis[n,2] = tmp[n];
    tmp = theta2; for (n in 1:N) sbasis[n,3] = tmp[n];
    tmp = theta1 .* theta1; for (n in 1:N) sbasis[n,4] = tmp[n];
    tmp = theta2 .* theta2; for (n in 1:N) sbasis[n,5] = tmp[n];
    tmp = theta1 .* theta2; for (n in 1:N) sbasis[n,6] = tmp[n];
    tmp = theta1 .* theta1 .* theta1; for (n in 1:N) sbasis[n,7] = tmp[n];
    tmp = theta2 .* theta2 .* theta2; for (n in 1:N) sbasis[n,8] = tmp[n];
    tmp = theta1 .* theta1 .* theta2; for (n in 1:N) sbasis[n,9] = tmp[n];
    tmp = theta1 .* theta2 .* theta2; for (n in 1:N) sbasis[n,10] = tmp[n];
    
    // create measurements
    psi[1:N]         = sbasis * surrogate1;
    psi[(1+N):(2*N)] = sbasis * surrogate2; 
  }

  // MUST use these priors for beta in order to specify GP model
  beta1 ~ std_normal();
  beta2 ~ std_normal();
  
  // hyperparameter priors
  rho1 ~ inv_gamma(1, 5);
  rho2 ~ inv_gamma(1, 5);

  //alpha1 ~ std_normal();
  alpha1 ~ inv_gamma(1, 5);
  //alpha2 ~ std_normal();
  alpha2 ~ inv_gamma(1, 5);

  // * implicitly use improper prior on mean1 and mean2 *

  // use custom likelihood function
  target += top_hat(y, psi, N, deltaS2, num);

}
"""

# compile stan model
sm = pystan.StanModel(model_code=stancode)

with open('stanmodel_tophat.pkl', 'wb') as f:
    pickle.dump(sm, f)


