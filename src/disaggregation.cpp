//
// Author: Anita Nandi
// Date: 2019-02-14

// Data: Spatial field mesh and matrices, polygon data, covariate pixel data

// The model:
// response = inv.logit( raster covariates + spatial field 2016)
// polygon response = sum(pixel pop x response) / sum(pixel pop) + normal or gamma error

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
  
  // Priors for liklihood
  PARAMETER(polygon_sd);
  DATA_SCALAR(polygon_sd_mean);
  DATA_SCALAR(polygon_sd_sd);
  
  // iid random effect for each polygon
  PARAMETER_VECTOR(iideffect);
  
  Type priormean_iideffect = 0.0;
  DATA_SCALAR(priorsd_iideffect);
  
  // spde hyperparameters
  PARAMETER(log_kappa);
  PARAMETER(log_tau);
  
  // Priors on spde hyperparameters
  DATA_SCALAR(priormean_log_kappa);
  DATA_SCALAR(priorsd_log_kappa);
  DATA_SCALAR(priormean_log_tau);
  DATA_SCALAR(priorsd_log_tau);
  
  // Convert hyperparameters to natural scale
  Type tau = exp(log_tau);
  Type kappa = exp(log_kappa);
  
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
    for(int s = 0; s < iideffect.size(); s++){
      nll -= dnorm(iideffect[s], priormean_iideffect, priorsd_iideffect, true);
    } 
  }
  
  nll -= dnorm(polygon_sd, polygon_sd_mean, polygon_sd_sd, true);
  
  if(field) {
    // Likelihood of hyperparameters for field
    nll -= dnorm(log_kappa, priormean_log_kappa, priorsd_log_kappa, true);
    nll -= dnorm(log_tau, priormean_log_tau, priorsd_log_tau, true);
    
    // Build spde matrix
    SparseMatrix<Type> Q = Q_spde(spde, kappa);
    
    // Likelihood of the random field.
    nll += SCALE(GMRF(Q), 1.0 / tau)(nodemean);
  }

  Type nll1 = nll;
  
  // ------------------------------------------------------------------------ //
  // Likelihood from data
  // ------------------------------------------------------------------------ //
  
  vector<Type> pixel_linear_pred(n_pixels);
  pixel_linear_pred = intercept + x * slope;
  
  if(field) {
    // Calculate field for pixel data
    vector<Type> logit_prevalence_field;
    logit_prevalence_field = Apixel * nodemean;
    pixel_linear_pred += logit_prevalence_field.array();
  }
  
  // recalculate startendindices to be in the form start, n
  startendindex.col(1) = startendindex.col(1) - startendindex.col(0) + 1;
  
  Type polygon_response;
  Type normalised_polygon_response;
  vector<Type> pixel_pred;
  vector<Type> numerator_pixels;
  vector<Type> normalisation_pixels;
  vector<Type> reportprediction(n_polygons);
  vector<Type> reportnll(n_polygons);
  
  // For each shape get pixel predictions within and aggregate to polygon level
  for (int polygon = 0; polygon < n_polygons; polygon++) {

    // Get pixel level predictions
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

    // Use correct likelihood function
    if(family == 0) {
      // Calculate normal likelihood in rate space
      reportprediction[polygon] = sum(numerator_pixels) / sum(normalisation_pixels);
      polygon_response = polygon_response_data(polygon);
      normalised_polygon_response = polygon_response/sum(normalisation_pixels);
      nll -= dnorm(normalised_polygon_response, reportprediction[polygon], polygon_sd, true);
      reportnll[polygon] = -dnorm(normalised_polygon_response, reportprediction[polygon], polygon_sd, true);
    } else if(family == 1) {
      reportprediction[polygon] = sum(numerator_pixels) / sum(normalisation_pixels);
      nll -= dbinom(polygon_response_data[polygon], response_sample_size[polygon], reportprediction[polygon], true);
      reportnll[polygon] = -dbinom(polygon_response_data[polygon], response_sample_size[polygon], reportprediction[polygon], true);
    } else if(family == 2) {
      reportprediction[polygon] = sum(numerator_pixels);
      nll -= dpois(polygon_response_data[polygon], reportprediction[polygon], true);
      reportnll[polygon] = -dpois(polygon_response_data[polygon], reportprediction[polygon], true);
    } else {
      error("Likelihood not implemented.");
    }
    
  }
  
  REPORT(reportprediction);
  REPORT(reportnll);
  REPORT(polygon_response_data);
  if(iid) {
    REPORT(iideffect);
  }
  REPORT(nll1);
  REPORT(nll);
  
  return nll;
  }
