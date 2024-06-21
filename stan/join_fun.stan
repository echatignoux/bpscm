/**-----------------**/
/**  JOIN_FUN
/**-----------------**/
/** Created: 03 may 2024
/**-----------------**/
// Commentary:
//
// Stan functions needed to estimate BP-SCM:
//  - spatial priors for latent variables
//  - set up of penalized spline
//  - BP likelihood
//
/**-----------------**/


/**
 * Log probability of the intrinsic conditional autoregressive (ICAR) prior,
 * excluding additive constants. 
 *
 * @param phi Vector of parameters for spatial smoothing (on unit scale)
 * @param spatial_scale Scale parameter for the ICAR model
 * @param node1 
 * @param node2
 * @param k number of groups
 * @param comp_size number of observational units in each group
 * @param comp_idx index of observations in order of their group membership
 * @param has_theta If the model contains an independent partial pooling term, phi for singletons can be zeroed out; otherwise, they require a standard normal prior. Both BYM and BYM2 have theta.
 *
 * @return Log probability density of ICAR prior up to additive constant
 **/
real icar_normal_lpdf(vector phi, real spatial_scale,
              int[] node1, int[] node2, 
              int n_comp, int[] comp_size, int[] comp_idx,
              int has_theta) {
  real lp;
  int pos=1;
  lp = -0.5 * dot_self(phi[node1] - phi[node2]);
  if (has_theta) {
    for (j in 1:n_comp) {
      /* sum to zero constraint for each connected group; singletons zero out */
      lp += normal_lpdf(sum(phi[segment(comp_idx, pos, comp_size[j])]) | 0, 0.001 * comp_size[j]);
      pos += comp_size[j];
    }
  } else {
    /* does not have theta */
    for (j in 1:n_comp) {
      if (comp_size[j] > 1) {
    /* same as above for non-singletons: sum to zero constraint */
    lp += normal_lpdf(sum(phi[segment(comp_idx, pos, comp_size[j])]) | 0, 0.001 * comp_size[j]);
      } else {
    /* its a singleton: independent Gaussian prior on phi */
    lp += normal_lpdf(phi[ segment(comp_idx, pos, comp_size[j]) ] | 0, spatial_scale);
      }      
      pos += comp_size[j];
    }
  }
  return lp;
}


/**
 * Combine local and global partial-pooling components into the convolved BYM2 term.
 *
 * @param phi_tilde local (spatially autocorrelated) component
 * @param theta_tilde global component
 * @param spatial_scale scale parameter for the convolution term
 * @param n number of spatial units
 * @param n_comp number of connected groups
 * @param comp_size number of observational units in each group
 * @param comp_idx index of observations in order of their group membership
 * @param inv_sqrt_scale_factor The scaling factor for the ICAR variance (see scale_c R function, using R-INLA); 
 *                              transformed from 1/scale^2 --> scale. Or, a vector of ones.
 * @param rho proportion of convolution that is spatially autocorrelated
 *
 * @return BYM2 convolution vector
 */
vector convolve_bym2(vector phi_tilde, vector theta_tilde,
		     real spatial_scale,
		     int n, int n_comp,
		     int[] comp_size, int[] comp_idx,
		     real rho, vector scale_factor
		     ) {
  vector[n] convolution;
  int pos=1;
  for (j in 1:n_comp) {
    if (comp_size[j] == 1) {
        convolution[ segment(comp_idx, pos, comp_size[j]) ] = spatial_scale * theta_tilde[ segment(comp_idx, pos, comp_size[j]) ];
    } else {
    convolution[ segment(comp_idx, pos, comp_size[j]) ] =
      spatial_scale * (sqrt(rho) * inv(sqrt(scale_factor[j])) * phi_tilde[ segment(comp_idx, pos, comp_size[j]) ] +
		       sqrt(1 - rho) * theta_tilde[ segment(comp_idx, pos, comp_size[j]) ]
		       );
  }
  pos += comp_size[j];
  }
  return convolution;
}


/**
 * Scale underliing icar graphe
 *
 * @param phi_tilde unscaled icar
 * @param n Number of factor levels
 * @param n_comp number of connected groups
 * @param comp_size number of observational units in each group
 * @param comp_idx index of observations in order of their group membership
 * @param scale_factor Scaling factor for the ICAR variance (see scale_c R function, using R-INLA); 
 * @return Vector of non-linear predictor
 **/
vector scale_icar(vector phi_tilde,
		  int n, int n_comp,
		  int[] comp_size, int[] comp_idx,
		  vector scale_factor
		  ) {
  vector[n] sc_icar;
  int pos=1;
  for (j in 1:n_comp) {
    if (comp_size[j] == 1) {
        sc_icar[ segment(comp_idx, pos, comp_size[j]) ] = phi_tilde[ segment(comp_idx, pos, comp_size[j]) ];
    } else {
    sc_icar[ segment(comp_idx, pos, comp_size[j]) ] =
      inv(sqrt(scale_factor[j])) * phi_tilde[ segment(comp_idx, pos, comp_size[j]) ] ;
  }
  pos += comp_size[j];
  }
  return sc_icar;
}


/**
 * Expand vector of values for each geo to vector of values for each obs
 *
 * @param x_geo Vector of value for each geo
 * @param n_obs number of obs
 * @param n_geo number of geo
 * @param id_geo vector of geo for each obs
 *
 * @return Vector of re assigned to each obs
 */
vector assign_geo(vector x_geo,
		  int n_obs,
		  int n_geo,
		  int[] id_geo) {

  vector[n_obs] x_obs;
  for (i in 1:n_obs) {
    x_obs[i] = x_geo[id_geo[i]];
  }
  return x_obs;

}

/**
 * Setting up random effects
 *
 * @param re_typ type of re (1 = iid, 2 = icar, 3 = bym)
 * @param theta Unstructures component
 * @param phi local (spatially autocorrelated) component
 * @param sds sd 
 * @param rho part of spatially autocorrelated coef in the BYM formulation
 * @param n_geo Number of factor levels
 * @param n_comp number of connected groups
 * @param comp_size number of observational units in each group
 * @param comp_idx index of observations in order of their group membership
 * @param scale_factor Scaling factor for the ICAR variance (see scale_c R function, using R-INLA); 
 * @return Vector of non-linear predictor
 **/
vector set_re(
	      int re_typ,
	      vector theta,
	      vector phi,
	      real sds,
	      real rho,
	      int n_geo,
	      int n_comp,
	      int[] comp_size,
	      int[] comp_idx,
	      vector scale_factor
	      ) {

  vector[n_geo] re ;
  /* re_k */
  if (re_typ == 0 ){
    re =  rep_vector(0,n_geo);
  }
  if (re_typ == 1 ){
    re =  theta*sds;
  }
  if (re_typ == 2 ){
    re = scale_icar( phi,  n_geo, n_comp, comp_size, comp_idx, scale_factor)  * sds;
  }
  if (re_typ == 3 ){
    re = convolve_bym2(phi, theta, sds,
		       n_geo, n_comp, comp_size, comp_idx,
		       rho, scale_factor);
  }
  return(re);
}


/**
 * Create spline function from spline basis and coefs
 *
 * @param bs standardized coef for nl part
 * @param n_obs number of obs
 * @param n_col spline basis dimension
 * @param Z spline basis
 * @param sig_nl spline sd (i.e. penalty) for nl part
 *
 * @return Spline basis * coef as a vector
 */
vector pen_spl(vector bs,
	       int n_obs , 
	       int n_col,  
	       matrix Z,
	       real sig_nl
	       ) {

  vector[n_obs] s ;
  vector[n_col] b ;
  b[1] = bs[1] ;
  b[2:n_col] = bs[2:n_col] * sig_nl;
  s = Z * b ;
  return s;
}

/**
* Log probability of bivariate poisson distribution
*
* @param r vector of observations (Y1,Y2)
* @param lmu_x1 log intensity for x1
* @param lmu_x2 log intensity for x2
* @param lmu_x3 log intensity for x3 (common cases)
*
* @return Log probability density of bivariate poisson distribution
**/
real bipois_lpmf(int[] r , real lmu_x1, real lmu_x2, real lmu_x3) {
  real ss;
  real log_s;
  real lmus;
  int  miny;
  
  miny = min(r[1], r[2]);
  
  ss = poisson_log_lpmf(r[1] | lmu_x1) + poisson_log_lpmf(r[2] | lmu_x2) - exp(lmu_x3);
  if(miny > 0) {
    lmus = -lmu_x1-lmu_x2+lmu_x3;
    log_s = ss; 
    /* log_s = 0; */
    
    for(k in 1:miny) {
      log_s +=
      lmus +
      log(r[1] - k + 1) +
      log(r[2] - k + 1) +
      - log(k);
      ss = log_sum_exp(ss, log_s);
    }
  }
  return(ss);
}

/**
* Laplace approximation to the log probability of bivariate poisson distribution
*
* @param r vector of observations (Y1,Y2)
* @param lmu_x1 log intensity for x1
* @param lmu_x2 log intensity for x2
* @param lmu_x3 log intensity for x3 (common cases)
*
* @return Log probability density of bivariate poisson distribution
**/
real bipoislap_lpmf(int[] r , real lmu_x1, real lmu_x2, real lmu_x3) {
  real ss;
  real z;
  real mc;
  real ms;
  real logbp;
  real b;
  real delta;
  real xx;
  real yy;
  int  miny;
  miny = min(r[1], r[2]);

  logbp=0;
  z = exp(-lmu_x1-lmu_x2+lmu_x3) ;
  ss = poisson_log_lpmf(r[1] | lmu_x1) + poisson_log_lpmf(r[2] | lmu_x2) - exp(lmu_x3);
  xx = r[1]+.5;
  yy = r[2]+.5;
  b = -(xx+yy)*z-1 ;
  delta = b^2-4*z*(xx*yy*z-0.5);
  mc = (-b-sqrt(delta))/(2*z);
  ms = abs(trigamma(r[1]-mc+1)+trigamma(r[2]-mc+1)+trigamma(mc+1));
  
  logbp = 0.5*log(2*pi()/ms) + 
  - lgamma(r[1] - mc + 1) +
  - lgamma(r[2] - mc + 1) - lgamma(mc + 1) +
  mc*log(z)  + 
  lgamma(r[1]+1) + lgamma(r[2]+1); 
  
  logbp +=  ss + log(normal_cdf( miny , mc, sqrt(ms)) - normal_cdf( 0 , mc, sqrt(ms))) ;
  
  return(logbp);
}


/**
* Log probability of bivariate poisson distribution
* combinaison of exact (min(Y1,Y2)<50) and laplace
* approximation (min(Y1,Y2)>=50)
*
* @param r vector of observations (Y1,Y2)
* @param lmu_x1 log intensity for x1
* @param lmu_x2 log intensity for x2
* @param lmu_x3 log intensity for x3 (common cases)
*
* @return Log probability density of bivariate poisson distribution
**/

real bp_lpmf(int[] r , real lmu_x1, real lmu_x2, real lmu_x3) {
  real ll;
  int  miny;
  miny = min(r[1], r[2]);
  if ((miny > 50)) {
    ll = bipoislap_lpmf(r|lmu_x1,lmu_x2,lmu_x3);
  } else {
    ll = bipois_lpmf(r|lmu_x1,lmu_x2,lmu_x3);
  }
  return(ll);
}

