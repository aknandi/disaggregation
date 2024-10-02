Disaggregation
==============


[![CRANstatus](https://www.r-pkg.org/badges/version/disaggregation)](https://cran.r-project.org/package=disaggregation)
[![R-CMD-check](https://github.com/aknandi/disaggregation/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/aknandi/disaggregation/actions/workflows/R-CMD-check.yaml)


A package containing useful functions for disaggregation modelling.
An overview of the package is given in [our paper](https://www.jstatsoft.org/article/view/v106i11).

Installation
------------

```R
devtools::install_github('aknandi/disaggregation')
```

Overview
--------

## Data preparation

Function prepare_data takes in sf (response) and SpatRaster (covariates) to produce a data structure required for the disaggregation modelling. Calls functions to extract covariate data, polygon data, aggregation (population data), match points to polygons and build an INLA mesh for the spatial field (build_mesh)

```R
data_for_model <- prepare_data(polygon_shapefile = shps, 
                               covariate_rasters = covariate_stack, 
                               aggregation_raster = population_raster,
                               id_var = 'area_id',
                               response_var = 'inc_counts')
```

### Input data

* A SpatRaster of covariate rasters to be used in the model (covariate_rasters)
* A sf (polygon_shapefile) containing at least two columns: one with the id for the polygons (id_var) and one with the response count data (response_var); for binomial data, i.e survey data, it can also contain a sample size column (sample_size_var).
* (Optional) SpatRaster used to aggregate the pixel level predictions (aggregation_raster) to polygon level (usually population). If this is not supplied a uniform raster will be used

### Controlling the mesh

The argument mesh_args in prepare_data allows you to supply a list of INLA mesh parameters to control the mesh used for the spatial field

```R
data_for_model <- prepare_data(shps, covariate_stack, 
                               mesh_args = list(max.edge = c(0.2, 8), 
                                                cutoff = 0.05, 
                                                offset = c(1, 0.7)))

```

### Dealing with NAs

There is an na_action flag that is automatically off. If there are any NAs in your response or covariate data within the polygons the prepare_data method will fail. We advise you to sort out the NAs in your data yourself, however, if you want the function to automatically deal with NAs you can set na_action = TRUE. This will remove any polygons that have NAs as a response or a population of zero, set any aggregation pixels with NA to zero and set covariate NAs pixels to the median value for the that covariate.

```R
data_for_model <- prepare_data(shps, covariate_stack, 
                               aggregation_raster = population_raster,
                               na_action = TRUE)
```

## Model fitting

Function disag_model takes the data structure returned by prepare_data and fits a TMB disaggregation model. Here you can specify priors, likelihood function, link function and whether to include a field or iid effect (default includes both).

```R
model_result <- disag_model(data_for_model, 
                            family = 'gaussian', 
                            link = 'identity',
                            field = TRUE,
                            iid = FALSE)
```

### Specify priors

Priors can be specified for the regression parameters, field and iid effect as a single list.

```R
model_result <- disag_model(data_for_model, 
                            priors = list(priormean_intercept = -2.0,
                                          priorsd_intercept = 2.0,
                                          priormean_slope = 0.0,
                                          priorsd_slope = 0.3))
```

## Model prediction

### Predict model

Function predict takes data structure returned by disag_model to predict model results and uncertainty.

```R
prediction <- predict(model_result)
```

## Plotting and summary functions

Plotting functions for input data, model results and predictions

```R
plot(data_for_model) # Plots polygon data, covariate rasters and INLa mesh
plot(model_result) # Plots of fixed effects parameters and in sample predictions
plot(prediction) #  Plots the mean, upper CI and lower CI prediction rasters
```

Summary functions for input data and model results

```R
summary(data_for_model)  
summary(model_result) 
```
