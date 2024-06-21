/**-----------------**/
/**  BP_SCM  
/**-----------------**/
/** Created: 03 may  2024
/**-----------------**/
// Commentary:
//
// Stan code to estimate BP-SCM
//
/**-----------------**/

functions {
#include ./join_fun.stan
}
data{
  int prior_only;  // should the likelihood be ignored?
  int<lower=1> n_obs;  // total number of observations
  // Geographical data
  int<lower=0> n_geo;
  int<lower=1> id_geo[n_obs];  // grouping indicator per observation
  int<lower=1> n_comp; // number of connected subgraphs
  vector[n_comp] scales; // scale of subgraphs
  int comp_size[n_comp]; // observational units per group
  int comp_idx[n_geo]; // index of geo, ordered by group
  int<lower=1> n_edges; 
  int<lower=1, upper=n_obs> node1[n_edges];
  int<lower=1, upper=n_obs> node2[n_edges];
  // Proxy incidence data
  int inc[n_obs,2]; // response
  matrix[n_obs,2] off; // offsets
  // Covar : splines
  int n_age;
  int d_z;
  matrix[n_obs, d_z] Z;  
  // Type of random effect for b1, b2 and c (0 = none, 1 = iid, 2 = icar, 3 = bym2)
  int<lower=0,upper=3> re_typ[3];
}

transformed data {
  int has_re[3];
  for (j in 1:3)
  has_re[j] = (re_typ[j]>0);    
}

parameters{
  /* Int */
  real Int[3];
  /* Param RE */
  real<lower=0, upper = 2> s_t[3]; // total sd for the x_i
  real<lower=0, upper = 1> p_b[2]; // proportion due to the b_j 
  matrix[n_geo,3] re_un; /* Un-structured */
  matrix[n_geo,3] re_ic; /* Structured */
  real<lower=0,upper=1> phi[3]; /*propotion of structured variance*/
  /*Params splines*/
  matrix[d_z,3] b_s; // Coef spline std
  real<lower=0, upper = 2> sds_nl;
}

transformed parameters{
  real sd_b[2]; // total sd for the x_i
  real sd_c[3]; // total sd for the x_i
  matrix[n_obs,3] llbd ;
  matrix[n_geo,3] bc ;
  matrix[n_geo,3] lnu ;
  matrix[n_obs,3] s_a ; // s effects

  // sd for the b_j and the c_j
  for (j in 1:2) {
    sd_b[j] = p_b[j] * s_t[j];
    sd_c[j] = (1-p_b[j]) * s_t[j];
  }
  sd_c[3] = s_t[3];

  // Setup random intercepts
  for (j in 1:3) {
    bc[,j] = set_re(re_typ[j], re_un[,j], re_ic[,j], 1, phi[j],n_geo, n_comp, comp_size, comp_idx, scales);  
  }

  for (j in 1:2){
    lnu[,j]  = has_re[j] * sd_b[j] * bc[,j] + has_re[3] * sd_c[j]* bc[,3]  ;
  }
  lnu[,3]  = sd_c[3] * bc[,3] ;

  // Setup splines
  for (j in 1:3){
    s_a[,j] = pen_spl(b_s[,j],n_obs,d_z,Z,sds_nl);
  }
  
  /* log(lambda) */
  for (j in 1:3){
    llbd[,j]  = off[,1] + Int[j] + s_a[,j] +
       assign_geo(lnu[,j],n_obs,n_geo,id_geo)  ;
  }

}
model{
  // Priors
  // Constants
   target += normal_lpdf(Int | 0,5);
  // re
  for (j in 1:3){
    target += std_normal_lpdf(re_un[,j]);
    if (re_typ[j] < 2 ){
      target += std_normal_lpdf(re_ic[,j]);
    }
    if (re_typ[j] == 2){
      re_ic[,j] ~ icar_normal(1, node1, node2, n_comp, comp_size, comp_idx, 0);
    }
    if (re_typ[j] == 3){
      re_ic[,j] ~ icar_normal(1, node1, node2, n_comp, comp_size, comp_idx, 1);
    }
    target += beta_lpdf(phi[j]|1,1);
    target += student_t_lpdf(s_t[j] | 3, 0, 2.5) - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  }
  // Splines
  for (j in 1:3){
    target += std_normal_lpdf(b_s[,j]);
  }
  target += student_t_lpdf(sds_nl | 3, 0, 2.5) - 1 * student_t_lccdf(0 | 3, 0, 2.5);

  // likelihood
  if(!prior_only){
    for (n in 1:n_obs){
      target += bp_lpmf( {inc[n,1], inc[n,2]} |
			 llbd[n,1],
			 llbd[n,2],
			 llbd[n,3] );
    }
  }
}

generated quantities{
  matrix[n_geo,3] lbd;
  matrix[n_geo,2] mu;
  matrix[n_geo,3] nu;
  matrix[n_geo,2] theta;
  matrix[n_age,3] rho; 
  matrix[n_age,2] r; 
  vector[2] w;
  vector[5] tau;
  vector[2] tau_sh;
  vector[2] tau_sp;
  vector[4] p_sh;
  vector[n_obs] log_lik;

  //Log - lik
  for (n in 1:n_obs){
    log_lik[n]=bp_lpmf( {inc[n,1], inc[n,2]} |
			llbd[n,1],
			llbd[n,2],
			llbd[n,3] );
  }
  
  // Age rates
  for (i in 1:n_age){
    for (j in 1:3){
      rho[i,j] = exp( Int[j] + s_a[i,j] );
    }
    for (j in 1:2){
      r[i,j] = rho[i,j]+rho[i,3];
    }
  }

  // Expected relative risks (nu) and total number of cases (lbd) per geo for the X_i
  for (j in 1:3){
    nu[,j] = exp(lnu[,j]);
    for (i in 1:(n_geo)){
        lbd[i,j] = sum(exp( llbd[((i-1)*n_age + 1):(i*n_age),j] ));
    }
  }

  // Expected relative risks (theta) and total number of cases per geo (mu) for the Y_i
  for (j in 1:2){
    mu[,j] = lbd[,j]  + lbd[,3] ;
    for (i in 1:(n_geo)){
      theta[i,j]= sum(exp( off[((i-1)*n_age + 1):(i*n_age),1] + log( r[,j] ))) ; // "Expected" cases
    }
    theta[,j] = mu[,j]./theta[,j];
  }

  // Mean proportion of X3 in Y1 and Y2
  w[1] = sum(lbd[,3])/sum(mu[,1]);
  w[2] = sum(lbd[,3])/sum(mu[,2]);

  // sd due to specific component at the Y_i scale
  tau_sp[1] = has_re[1] * ( (1-w[1])*sd_b[1] );
  tau_sp[2] = has_re[2] * ( (1-w[2])*sd_b[2] );
  
  // sd due to shared components at the Y_i scale
  tau_sh[1] = has_re[3] * ( (1-w[1])*sd_c[1] + w[1]*sd_c[3] );
  tau_sh[2] = has_re[3] * ( (1-w[2])*sd_c[2] + w[2]*sd_c[3] );

  // total sd for X_i and Y_i
  tau[1] = sqrt(tau_sh[1]^2+tau_sp[1]^2);
  tau[2] = sqrt(tau_sh[2]^2+tau_sp[2]^2);
  tau[3] = sqrt(sd_c[1]^2+sd_b[1]^2);
  tau[4] = sqrt(sd_c[2]^2+sd_b[2]^2);
  tau[5] = sd_c[3]^2;  

  // Proportion of variance due to shared random effect
  p_sh[1] = tau_sh[1]^2/tau[1]^2;
  p_sh[2] = tau_sh[2]^2/tau[2]^2;
  p_sh[3] = sd_c[1]^2/tau[3]^2;
  p_sh[4] = sd_c[2]^2/tau[4]^2;

}


