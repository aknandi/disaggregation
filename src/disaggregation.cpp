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
  
  // Pixel data.
  // All in long format so that if a pixel is in multiple polygons or multiple years, it will be represented by multiple rows.
  // Environmental/other covariates data matrix
  DATA_MATRIX(x);
  
  // ADMIN0 region id (pixel i is in ADMIN0 polygon j). Sorted by... shapefile id
  //DATA_IVECTOR(shapeadmin0);
  
  // two col matrix with start end indices for each shape case.
  DATA_IARRAY(startendindex);
  
  // Shape data. Cases and region id.
  DATA_VECTOR(polygon_response_data);
  
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
  
  /* // iid random effect for each country. Length is the number of countries in analysis.
  PARAMETER_VECTOR(admin0_slope);
  
  Type priormean_admin0_slope = 0.0;
  DATA_SCALAR(priorsd_admin0_slope); */
  
  // 2016 spde hyperparameters
  // tau defines strength of random field.
  // kappa defines distance within which points in field affect each other.
  PARAMETER(log_kappa);
  PARAMETER(log_tau);
  
  // Priors on spde hyperparameters
  //   kappa -- i.e. exp(priormean_log_kappa) -- set as approximately the width of the region being studied. This implies prior belief in a fairly flat field.
  //   tau -- exp(priormean_log_tau) -- set to close to zero. Betas on regression coefficients have priors of 0 so this is reasonable.
  
  DATA_SCALAR(priormean_log_kappa);
  DATA_SCALAR(priorsd_log_kappa);
  DATA_SCALAR(priormean_log_tau);
  DATA_SCALAR(priorsd_log_tau);
  
  // Convert hyperparameters to natural scale
  Type tau = exp(log_tau);
  Type kappa = exp(log_kappa);
  
  // Space-time random effect parameters
  // matrix logit_pr_offset [nrows = n_mesh, col=n_years].
  PARAMETER_VECTOR(nodemean);
  
  // get number of data points to loop over
  // y (cases) length
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
  
  /* for(int s = 0; s < admin0_slope.size(); s++){
  nll -= dnorm(admin0_slope[s], priormean_admin0_slope, priorsd_admin0_slope, true);
} */
  
  nll -= dnorm(polygon_sd, polygon_sd_mean, polygon_sd_sd, true);
  
  // Likelihood of hyperparameters for 2016 field
  nll -= dnorm(log_kappa, priormean_log_kappa, priorsd_log_kappa, true);
  nll -= dnorm(log_tau, priormean_log_tau, priorsd_log_tau, true);
  
  // Build spde matrix
  SparseMatrix<Type> Q = Q_spde(spde, kappa);
  
  // Likelihood of the random field.
  nll += SCALE(GMRF(Q), 1.0 / tau)(nodemean);
  
  Type nll1 = nll;
  
  // ------------------------------------------------------------------------ //
  // Calculate random field effects
  // ------------------------------------------------------------------------ //
  
  // Calculate field for pixel data
  vector<Type> logit_prevalence_field_2016;
  logit_prevalence_field_2016 = Apixel * nodemean;
  
  // ------------------------------------------------------------------------ //
  // Likelihood from data
  // ------------------------------------------------------------------------ //
  
  vector<Type> pixel_linear_pred(n_pixels);
  pixel_linear_pred = intercept + x * slope + logit_prevalence_field_2016.array();
  
  // recalculate startendindices to be in the form start, n
  startendindex.col(1) = startendindex.col(1) - startendindex.col(0) + 1;
  
  vector<Type> pixel_pred;
  
  vector<Type> reportprediction(n_polygons);
  vector<Type> reportnll(n_polygons);
  
  //For each shape use startendindex to find sum of pixel incidence rates
  for (int polygon = 0; polygon < n_polygons; polygon++) {
    // Sum pixel risks (raster + field
    
    // Get pixel level predictions
    pixel_pred = pixel_linear_pred.segment(startendindex(polygon, 0), startendindex(polygon, 1)).array(); // + admin0_slope[shapeadmin0[polygon] - 1];
    pixel_pred = invlogit(pixel_pred);
    
    reportprediction[polygon] = sum(pixel_pred);
    
    // Calculate likelihood from polygon prediction
    nll -= dnorm(polygon_response_data[polygon], reportprediction[polygon], polygon_sd, true);
    reportnll[polygon] = -dnorm(polygon_response_data[polygon], reportprediction[polygon], polygon_sd, true);
  }
  
  REPORT(reportprediction);
  REPORT(reportnll);
  REPORT(polygon_response_data);
  REPORT(nll1);
  REPORT(nll);
  
  return nll;
  }
