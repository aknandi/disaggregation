#' Parallel extraction of raster stack by shape file.
#' 
#' Parallelisation is performed across rasters, not shapes. 
#' So this function is only useful if you are extracting 
#' data from many raster layers.
#' As the overhead for parallel computation in windows is high
#' it only makes sense to parallelise in this way.
#'
#' 
#' @param raster A RasterBrick or RasterStack object.
#' @param shape A SpatialPolygons object.
#' @param fun The function used to aggregate the pixel data. If NULL, raw pixel data is returned.
#' @param id Name of column in shape object to be used to bind an ID column to output.
#' @param ... Other arguments to raster::extract.
#' 
#' @return A data.frame with columns of polygon id, cell id (if fun = NULL) and a column for each raster in the stack
#' 
#' @importFrom foreach %dopar%
#' @importFrom parallel stopCluster
#' @importFrom parallel makeCluster
#' @importFrom doParallel registerDoParallel
#'
#' @export
#' @examples 
#'  \dontrun{
#'   polygons <- list()
#'   for(i in 1:100) {
#'     row <- ceiling(i/10)
#'     col <- ifelse(i %% 10 != 0, i %% 10, 10)
#'     xmin = 2*(col - 1); xmax = 2*col; ymin = 2*(row - 1); ymax = 2*row
#'     polygons[[i]] <- rbind(c(xmin, ymax), c(xmax,ymax), c(xmax, ymin), c(xmin,ymin))
#'   }
#' 
#'   polys <- do.call(raster::spPolygons, polygons)
#'   response_df <- data.frame(area_id = 1:100, response = runif(100, min = 0, max = 10))
#'   spdf <- sp::SpatialPolygonsDataFrame(polys, response_df)
#' 
#'   r <- raster::raster(ncol=20, nrow=20)
#'   r <- raster::setExtent(r, raster::extent(spdf))
#'   r[] <- sapply(1:raster::ncell(r), function(x) rnorm(1, ifelse(x %% 20 != 0, x %% 20, 20), 3))
#'   r2 <- raster::raster(ncol=20, nrow=20)
#'   r2 <- raster::setExtent(r2, raster::extent(spdf))
#'   r2[] <- sapply(1:raster::ncell(r), function(x) rnorm(1, ceiling(x/10), 3))
#'   cov_rasters <- raster::stack(r, r2)
#' 
#'   cl <- parallel::makeCluster(2)
#'   doParallel::registerDoParallel(cl)
#'   result <- parallelExtract(cov_rasters, spdf, fun = NULL, id = 'area_id')
#'   parallel::stopCluster(cl)
#'   foreach::registerDoSEQ()
#'  }

parallelExtract <- function(raster, shape, fun = mean, id = 'OBJECTID',  ...){
  
  if (!requireNamespace("foreach", quietly = TRUE)) {
    stop("foreach needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  shape@data[, id] <- as.character(shape@data[, id])
  
  i <- NULL
  # Run extract in parallel.
  values <- foreach::foreach(i = seq_along(shape)) %dopar% {  
    raster::extract(raster, shape[i, ], fun = fun, na.rm = TRUE, cellnumbers = TRUE, ...)
  } 
  if(!is.null(fun)){
    # If a summary function was given, just bind everything together and add ID column
    df <- data.frame(do.call(rbind, values))
    if(inherits(shape, 'SpatialPolygonsDataFrame')){
      df <- cbind(ID = as.data.frame(shape)[, id], df)
    } else{
      df <- cbind(ID = names(shape), df)
      id <- 'id'
    }
    
    names(df) <- c(id, names(raster))
    
    return(df)
  } else {
    # If no summary was given we get a list of length n.shapes
    #   each entry in the list is a dataframe with n.covariates columns
    #   Want to make covariates columns, rbind shapes, and add shape and cell id columns.
    
    # list of vectors, one for each covariate
    values_id <- lapply(seq_along(values), function(x) data.frame(shape@data[, id][x], values[[x]][[1]]))
    
    
    df <- do.call(rbind, values_id)
    names(df) <- c(id, 'cellid', names(raster))
    
    return(df)
  }
  
}


#' Extract polygon id and response data into a data.frame from a SpatialPolygonsDataFrame
#' 
#' Returns a data.frame with a row for each polygon in the SpatialPolygonDataFrame and columns: area_id, response and N, containing the id of the
#' polygon, the values of the response for that polygon, and the sample size respectively. If the data is not survey data (the sample size does 
#' not exist), this column will contain NAs.
#' 
#' @param shape A SpatialPolygons object containing response data.
#' @param id_var Name of column in shape object with the polygon id. Default 'area_id'.
#' @param response_var Name of column in shape object with the response data. Default 'response'.
#' @param sample_size_var For survey data, name of column in SpatialPolygonDataFrame object (if it exists) with the sample size data. Default NULL.
#' 
#' @return A data.frame with a row for each polygon in the SpatialPolygonDataFrame and columns: area_id, response and N, containing the id of the
#' polygon, the values of the response for that polygon, and the sample size respectively. If the data is not survey data (the sample size does 
#' not exist), this column will contain NAs.
#' 
#' @export
#' @examples {
#'  polygons <- list()
#'  for(i in 1:100) {
#'    row <- ceiling(i/10)
#'    col <- ifelse(i %% 10 != 0, i %% 10, 10)
#'    xmin = 2*(col - 1); xmax = 2*col; ymin = 2*(row - 1); ymax = 2*row
#'    polygons[[i]] <- rbind(c(xmin, ymax), c(xmax,ymax), c(xmax, ymin), c(xmin,ymin))
#'  }
#' 
#'  polys <- do.call(raster::spPolygons, polygons)
#'  response_df <- data.frame(area_id = 1:100, response = runif(100, min = 0, max = 10))
#'  spdf <- sp::SpatialPolygonsDataFrame(polys, response_df)
#' 
#'  getPolygonData(spdf, id_var = 'area_id', response_var = 'response')
#' }
#' 
#' 

getPolygonData <- function(shape, id_var = 'area_id', response_var = 'response', sample_size_var = NULL) {
  
  if(is.null(sample_size_var)) {
    polygon_df <- shape@data[, c(id_var, response_var)]
    polygon_df$N <- rep(NA, nrow(polygon_df))
  } else {
    polygon_df <- shape@data[, c(id_var, response_var, sample_size_var)]
  }
  
  names(polygon_df) <- c('area_id', 'response', 'N')
  
  return(polygon_df)
}


#' Get a RasterStack of covariates from a folder containing .tif files
#' 
#' Looks in a specified folder for raster files. Returns a RasterStack of the rasters cropped to the extent specified by the shape parameter.
#' 
#' @param directory Filepath to the directory containing the rasters.
#' @param file_pattern Pattern the filenames must match. Default is all files ending in .tif .
#' @param shape An object with an extent that the rasters will be cropped to.
#' 
#' @return A RasterStack of the raster files in the directory
#' 
#' @export
#' @examples 
#' \dontrun{
#'   getCovariateRasters('/home/rasters', '.tif$', shape)
#'  }
#' 

getCovariateRasters <- function(directory, file_pattern = '.tif$', shape) {
  
  stopifnot(dir.exists(directory))
  
  covariate_files <- list.files(directory, pattern = file_pattern, full.names = TRUE)
  stopifnot(length(covariate_files) != 0)
  
  covariate_rasters <- lapply(covariate_files, function(x) raster::raster(x))
  covariate_stack <- raster::stack(covariate_rasters)
  
  covariate_stack <- raster::crop(covariate_stack, shape)
  
  return(covariate_stack)
}

# Extract coordinates from raster to use constructing the INLA mesh
# 
# @param cov_rasters RasterStack of the covariate rasters.
# @param selectIds numeric vector containing cell ids to retain. Default NULL retains all cell ids in the covariate rasters.
# 
# @return A matrix containing the coordinates used to make the mesh

extractCoordsForMesh <- function(cov_rasters, selectIds = NULL) {
  
  stopifnot(inherits(cov_rasters, c('RasterStack', 'RasterBrick')))
  if(!is.null(selectIds)) stopifnot(inherits(selectIds, 'numeric'))
  
  points_raster <- cov_rasters[[1]]
  points_raster[is.na(raster::values(points_raster))] <- -9999
  raster_pts <- raster::rasterToPoints(points_raster, spatial = TRUE)
  coords <- raster_pts@coords
  
  # If specified, only retain certain pixel ids
  if(!is.null(selectIds)) {
    coords <- coords[selectIds, ]
  }
  
  return(coords)
  
}



