#' Function to fit the disaggregation model
#'
#' @param data disag.data object returned by prepare_data function that contains all the necessary objects for the model fitting
#' @param its number of iterations to run the optimisation for
#' 
#' @name fit_model
#'
#' @examples 
#'  \dontrun{
#'   fit_model(my_data)
#'  }
#'
#' @export

fit_model <- function(data, its = 10) {
  
  stopifnot(inherits(data, 'disag.data'))
  
  # Sort out mesh bits
  spde <- (INLA::inla.spde2.matern(data$mesh, alpha = 2)$param.inla)[c("M0", "M1", "M2")]	
  Apix <- INLA::inla.mesh.project(data$mesh, loc = data$coords)$A
  n_s <- nrow(spde$M0)
  
  cov_matrix <- as.matrix(data$covariate_data[, -c(1:2)])
  cov_matrix <- t(apply(cov_matrix, 1,as.numeric))

  parameters <- list(intercept = -5,
                     slope = rep(0, ncol(cov_matrix)),
                     polygon_sd_coverage = 0.1,
                     log_kappa = -3,
                     log_tau = -0.5,
                     nodemean = rep(0, n_s))
  
  input_data <- list(x = cov_matrix, 
                     Apixel = Apix,
                     spde = spde,
                     startendindex = data$startendindex,
                     polygon_coverage_data = data$polygon_data$response)
  
  
  obj <- TMB::MakeADFun(
    data = input_data, 
    parameters = parameters,
    random = c('nodemean'),
    DLL = "disaggregation")
  
  
  opt <- stats::nlminb(obj$par, obj$fn, obj$gr, control = list(iter.max = its, eval.max = 2*its, trace = 0))
  
  return(list(obj = obj,
              opt = opt))
  
  
}