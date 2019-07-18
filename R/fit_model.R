#' Function to fit the disaggregation model
#'
#' @param data disag.data object returned by prepare_data function that contains all the necessary objects for the model fitting
#' @param its number of iterations to run the optimisation for
#' 
#' @name fit_model
#'
#' @examples 
#' \dontrun{
#'  cds1 <- rbind(c(-120,-20), c(-100,5), c(-60, 0), c(-100,-60), c(-120,-20))
#'  cds2 <- rbind(c(30,0), c(50,60), c(70,0), c(70,-55), c(30,0))
#'  cds3 <- rbind(c(10,20), c(10,50), c(-60,40), c(-30,-30), c(30,10))
#'  polys <- raster::spPolygons(cds1, cds2, cds3)
#' 
#'  response_df <- data.frame(area_id = c('1', '2', '3'), response = c(4, 7, 2))
#' 
#'  spdf <- sp::SpatialPolygonsDataFrame(polys, response_df)
#' 
#'  # Create raster stack
#'  r <- raster::raster(ncol=36, nrow=18)
#'  r[] <- 1:raster::ncell(r)
#'  cov_rasters <- raster::stack(r, r)
#' 
#'  test_data <- prepare_data(polygon_shapefile = spdf, 
#'                            covariate_rasters = cov_rasters,
#'                            mesh.args = list(max.edge = c(50, 100)))
#'                          
#'  result <- fit_model(test_data, its = 2)
#'  }
#' 
#' @export

fit_model <- function(data, its = 10) {
  
  stopifnot(inherits(data, 'disag.data'))
  stopifnot(inherits(its, 'numeric'))
  
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
  
  sd_out <- TMB::sdreport(obj, getJointPrecision = TRUE)
  
  model_output <- list(obj = obj,
                       opt = opt,
                       sd_out = sd_out)
  
  class(model_output) <- c('fit.result', 'list')
  
  return(model_output)
  
  
}