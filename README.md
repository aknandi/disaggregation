Disaggregation
==============

[![Build Status](https://travis-ci.org/aknandi/disaggregation.svg)](https://travis-ci.org/aknandi/disaggregation)
[![codecov.io](https://codecov.io/github/aknandi/disaggregation/coverage.svg?branch=master)](https://codecov.io/github/aknandi/disaggregation?branch=master)

A package containing useful functions for disaggregation modelling


Overview
--------

### Extract

Functions to get covariate rasters, extract covariate data and setup polygon data to be used in model


### Matching

Function to match which pixels are contained within a given polygon


### Build mesh

Function to build INLA mesh

### Prepare data

Function takes in SpatialPolygonDataFrame (response) and RasterStack (covariates) to produce a data structure required for the disaggregation modelling. Calls functions from extract, matching and build_mesh



