//
// Author: Anita Nandi
// Date: 2019-02-14

// Data: Spatial field mesh and matrices, polygon data, covariate pixel data


#define TMB_LIB_INIT R_init_disaggregation
#include <TMB.hpp>

template <class Type>
Type objective_function<Type>::operator()()
{
  
  using namespace R_inla;
  using namespace density;
  using namespace Eigen;
  
  // ------------------------------------------------------------------------ //
  // Spatial field data
  // ------------------------------------------------------------------------ //
  
  // The A matrices are for projecting the mesh to a point for the pixel and point data respectively.
  DATA_SPARSE_MATRIX(Apixel);
  DATA_STRUCT(spde, spde_t);
  
  // ------------------------------------------------------------------------ //
  // Polygon level data
  // ------------------------------------------------------------------------ //
  
  // Covariate pixel data
  DATA_MATRIX(x);
  
  // two col matrix with start end indices for each shape case.
  DATA_IARRAY(startendindex);
  
  // Shape data. Cases and region id.
  DATA_VECTOR(polygon_response_data);
  DATA_VECTOR(response_sample_size);
  
  // Use to aggreagte pixel response values to polygon level
  DATA_VECTOR(aggregation_values);
  
  // ------------------------------------------------------------------------ //
  // Likelihood and link functions
  // ------------------------------------------------------------------------ //
  
  DATA_INTEGER(family);
  DATA_INTEGER(link);
  
  // ------------------------------------------------------------------------ //
  // Parameters
  // ------------------------------------------------------------------------ //
  
  PARAMETER(intercept);
  PARAMETER_VECTOR(slope);
  
  DATA_SCALAR(priormean_intercept);
  DATA_SCALAR(priorsd_intercept);
  DATA_SCALAR(priormean_slope);
  DATA_SCALAR(priorsd_slope);
  
  // Priors for likelihood
  PARAMETER(log_tau_gaussian);
  Type tau_gaussian = exp(log_tau_gaussian);
  Type gaussian_sd = 1 / sqrt(tau_gaussian);
  
  // INLA defines a loggamma prior on log tau. 
  // We evaluate a gamma prior on tau, but the parameters are 
  // therefore the same.
  Type prior_gamma_shape = 1;
  Type prior_gamma_rate = 5e-05;
  
  PARAMETER_VECTOR(iideffect);
  PARAMETER(iideffect_log_tau);
  Type iideffect_tau = exp(iideffect_log_tau);
  Type iideffect_sd = 1 / sqrt(iideffect_tau);
  
  Type iideffect_mean = 0.0;
  
  // Priors on iid random effect for polygons
  DATA_SCALAR(prior_iideffect_sd_max);
  DATA_SCALAR(prior_iideffect_sd_prob);
  
  // spde hyperparameters
  PARAMETER(log_sigma);
  PARAMETER(log_rho);
  Type sigma = exp(log_sigma);
  Type rho = exp(log_rho);
  
  // Priors on spde hyperparameters
  DATA_SCALAR(prior_rho_min);
  DATA_SCALAR(prior_rho_prob);
  DATA_SCALAR(prior_sigma_max);
  DATA_SCALAR(prior_sigma_prob);
  
  // Convert hyperparameters to natural scale
  // todo
  Type kappa = sqrt(8) / rho;
  Type nu = 1;
  Type tau = sigma * pow(kappa, nu) * sqrt(4 * M_PI);
  
  // Random effect parameters
  PARAMETER_VECTOR(nodemean);
  
  // Model component flags
  DATA_INTEGER(field);
  DATA_INTEGER(iid);
  
  // Number of polygons
  int n_polygons = polygon_response_data.size();
  // Number of pixels
  int n_pixels = x.rows();
  
  Type nll = 0.0;
  
  // ------------------------------------------------------------------------ //
  // Likelihood from priors
  // ------------------------------------------------------------------------ //
  
  nll -= dnorm(intercept, priormean_intercept, priorsd_intercept, true);
  for (int s = 0; s < slope.size(); s++) {
    nll -= dnorm(slope[s], priormean_slope, priorsd_slope, true);
  }
  
  if(iid) {
    // Likelihood of hyperparameter of polygon iid random effect.
    Type lambda = -log(prior_iideffect_sd_prob) / prior_iideffect_sd_max;
    Type pcdensityiid = lambda / 2 * pow(iideffect_tau, -3/2) * exp( - lambda * pow(iideffect_tau, -1/2));
    nll -= log(pcdensityiid);
    
    // Likelihood of random effect for polygons
    for(int p = 0; p < iideffect.size(); p++) {
      nll -= dnorm(iideffect[p], iideffect_mean, iideffect_sd, true);
    }
  }
  
  // Likelihood from the gaussian prior. 
  // log(prec) ~ loggamma
  // prec ~ gamma
  if(family == 0) {
    nll -= dgamma(tau_gaussian, prior_gamma_shape, prior_gamma_rate, true);
  }
  
  if(field) {
    // Likelihood of hyperparameters for field
    Type lambdatilde1 = -log(prior_rho_prob) * prior_rho_min;
    Type lambdatilde2 = -log(prior_sigma_prob) / prior_sigma_max;
    Type pcdensity = lambdatilde1 * lambdatilde2 * pow(rho, -2) * exp(-lambdatilde1 * pow(rho, -1) - lambdatilde2 * sigma);
    nll -= log(pcdensity);
    
    // Build spde matrix
    SparseMatrix<Type> Q = Q_spde(spde, kappa);
    
    // Likelihood of the random field.
    nll += SCALE(GMRF(Q), 1.0 / tau)(nodemean);
  }

  Type nll_priors = nll;
  
  // ------------------------------------------------------------------------ //
  // Likelihood from data
  // ------------------------------------------------------------------------ //
  
  vector<Type> pixel_linear_pred(n_pixels);
  pixel_linear_pred = intercept + x * slope;
  
  if(field) {
    // Calculate field for pixel data
    vector<Type> linear_pred_field(n_pixels);
    linear_pred_field = Apixel * nodemean;
    pixel_linear_pred += linear_pred_field.array();
  }
  
  // recalculate startendindices to be in the form start, n
  startendindex.col(1) = startendindex.col(1) - startendindex.col(0) + 1;
  
  Type polygon_response;
  Type normalised_polygon_response;
  Type normalisation_total;
  Type pred_polygoncases;
  Type pred_polygonrate;
  Type polygon_sd;
  vector<Type> pixel_pred;
  vector<Type> numerator_pixels;
  vector<Type> normalisation_pixels;
  vector<Type> reportnormalisation(n_polygons);
  vector<Type> reportprediction_cases(n_polygons);
  vector<Type> reportprediction_rate(n_polygons);
  vector<Type> reportnll(n_polygons);
  vector<Type> reportpolygonsd(n_polygons);
  
  // For each shape get pixel predictions within and aggregate to polygon level
  for (int polygon = 0; polygon < n_polygons; polygon++) {

    // Get pixel level predictions (rate)
    pixel_pred = pixel_linear_pred.segment(startendindex(polygon, 0), startendindex(polygon, 1)).array();
    if(iid) {
      pixel_pred += iideffect[polygon];
    }
    // Use correct link function
    if(link == 0) {
      pixel_pred = invlogit(pixel_pred);
    } else if(link == 1) {
      pixel_pred = exp(pixel_pred);
    } else if(link == 2){
      pixel_pred = pixel_pred;
    } else {
      error("Link function not implemented.");
    }
    
    // Aggregate to polygon prediction
    numerator_pixels = pixel_pred * aggregation_values.segment(startendindex(polygon, 0), startendindex(polygon, 1)).array();
    normalisation_pixels = aggregation_values.segment(startendindex(polygon, 0), startendindex(polygon, 1));
    normalisation_total = sum(normalisation_pixels);
    pred_polygoncases = sum(numerator_pixels);
    pred_polygonrate = pred_polygoncases/normalisation_total;
    
    reportnormalisation[polygon] = normalisation_total;
    reportprediction_cases[polygon] = pred_polygoncases;
    reportprediction_rate[polygon] = pred_polygonrate;
    
    // Use correct likelihood function
    if(family == 0) {
      // Scale the pixel sd to polygon level
      polygon_sd = gaussian_sd * sqrt((normalisation_pixels * normalisation_pixels).sum()) / normalisation_total;
      reportpolygonsd[polygon] = polygon_sd;
      // Calculate normal likelihood in rate space
      polygon_response = polygon_response_data(polygon);
      normalised_polygon_response = polygon_response/normalisation_total;
      nll -= dnorm(normalised_polygon_response, pred_polygonrate, polygon_sd, true);
      reportnll[polygon] = -dnorm(normalised_polygon_response, pred_polygonrate, polygon_sd, true);
    } else if(family == 1) {
      nll -= dbinom(polygon_response_data[polygon], response_sample_size[polygon], pred_polygonrate, true);
      reportnll[polygon] = -dbinom(polygon_response_data[polygon], response_sample_size[polygon], pred_polygonrate, true);
    } else if(family == 2) {
      nll -= dpois(polygon_response_data[polygon], pred_polygoncases, true);
      reportnll[polygon] = -dpois(polygon_response_data[polygon], pred_polygoncases, true);
    } else {
      error("Likelihood not implemented.");
    }
    
  }
  
  REPORT(reportprediction_cases);
  REPORT(reportprediction_rate);
  REPORT(reportnormalisation);
  REPORT(reportnll);
  REPORT(polygon_response_data);
  REPORT(nll_priors);
  REPORT(nll);
  if(family == 0) {
    REPORT(reportpolygonsd);
  }
  
  return nll;
}
