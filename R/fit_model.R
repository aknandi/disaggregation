#' Function to fit the disaggregation model
#'
#' @param data disag.data object returned by prepare_data function that contains all the necessary objects for the model fitting
#' @param its number of iterations to run the optimisation for
#' 
#' @name fit_model
#'
#' @examples 
#' \dontrun{
#'  polygons <- list()
#'  for(i in 1:100) {
#'   row <- ceiling(i/10)
#'   col <- ifelse(i %% 10 != 0, i %% 10, 10)
#'   xmin = 2*(col - 1); xmax = 2*col; ymin = 2*(row - 1); ymax = 2*row
#'   polygons[[i]] <- rbind(c(xmin, ymax), c(xmax,ymax), c(xmax, ymin), c(xmin,ymin))
#'  }
#' 
#'  polys <- do.call(raster::spPolygons, polygons)
#'  response_df <- data.frame(area_id = 1:100, response = runif(100, min = 0, max = 10))
#'  spdf <- sp::SpatialPolygonsDataFrame(polys, response_df)
#' 
#'  r <- raster::raster(ncol=20, nrow=20)
#'  r[] <- sapply(1:raster::ncell(r), function(x) rnorm(1, ifelse(x %% 20 != 0, x %% 20, 20), 3))
#'  r2 <- raster::raster(ncol=20, nrow=20)
#'  r2[] <- sapply(1:raster::ncell(r), function(x) rnorm(1, ceiling(x/10), 3))
#'  cov_rasters <- raster::stack(r, r2)
#' 
#'  test_data <- prepare_data(polygon_shapefile = spdf, 
#'                            covariate_rasters = cov_rasters)
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