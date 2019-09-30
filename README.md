Disaggregation
==============

[![Build Status](https://travis-ci.org/aknandi/disaggregation.svg)](https://travis-ci.org/aknandi/disaggregation)
[![codecov.io](https://codecov.io/github/aknandi/disaggregation/coverage.svg?branch=master)](https://codecov.io/github/aknandi/disaggregation?branch=master)

A package containing useful functions for disaggregation modelling

Installation
------------

```R
devtools::install_github('aknandi/disaggregation')
```

Overview
--------

## Data preparation

Function prepare_data takes in SpatialPolygonDataFrame (response) and RasterStack (covariates) to produce a data structure required for the disaggregation modelling. Calls functions to extract covariate data, polygon data, aggregation (population data), match points to polygons and build an INLA mesh for the spatial field (build_mesh)

```R
data_for_model <- prepare_data(polygon_shapefile = shps, 
                               covariate_rasters = covariate_stack, 
                               aggregation_raster = population_raster,
                               id_var = 'area_id',
                               response_var = 'inc_counts')
```

### Input data

* A RasterStack of covariate rasters to be used in the model (covariate_rasters)
* A SpatialPolygonsDataFrame (polygon_shapefile) containing at least two columns: one with the id for the polygons (id_var) and one with the response count data (response_var); for binomial data, i.e survey data, it can also contain a sample size column (sample_size_var).
* (Optional) Raster used to aggregate the pixel level predictions (aggregation_raster) to polygon level (usually population). If this is not supplied a uniform raster will be used

### Controlling the mesh

The argument mesh.args in prepare_data allows you to supply a list of INLA mesh parameter to control the mesh used for the spatial field

```R
data_for_model <- prepare_data(shps, covariate_stack, 
                               mesh.args = list(max.edge = c(0.2, 8), 
                                                cut = 0.05, 
                                                offset = c(1, 0.7)))

```

### Dealing with NAs

There is an na.action flag that is automatically off. If there are any NAs in your response or covariate data within the polygons the prepare_data method will fail. We advise you to sort out the NAs in your data yourself, however, if you want the function to automatically deal with NAs you can set na.action = TRUE. This will remove any polygons that have NAs as a response, set any aggregation pixels with NA to zero and set covariate NAs pixels to the median value for the that covariate.

```R
data_for_model <- prepare_data(shps, covariate_stack, 
                               aggregation_raster = population_raster,
                               na.action = TRUE)
```

## Model fitting

Function fit_model takes data structure returned by prepare_data and fits a TMB disaggregation model. Here you can specify priors, likelihood function, link function and whether to include a field or iid effect (default includes both)

```R
model_result <- fit_model(data_for_model, 
                          family = 'gaussian', 
                          link = 'identity')
```

## Model prediction

### Predict model

Function predict takes data structure returned by fit_model to predict model results and uncertainty.

```R
predictions <- predict(model_result)
```

## Plotting and summary functions

Plotting functions for input data, model results and predictions

```R
plot(data_for_model) # Plots polygon data, covariate rasters and INLa mesh
plot(model_result) # Plots of fixed effects parameters and in sample predictions
plot(predictions$mean_predictions) #  Plots the mean prediction, covariate contribution and the field and iid contribution (if they were used in the model)
plot(predictions$uncertainty_predictions) # Plots the upper and lower CI rasters
```

Summary functions for input data and model results

```R
summary(data_for_model)  
summary(model_result) 
```

