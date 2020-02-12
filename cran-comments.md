## Update
This is a package update (version 0.1.3). The changes in this version are:

* Renamed fit_model function to disag_model. Deprecated fit_model, will be removed in the next version

* Renamed classes disag.data and fit.result to disag_data and disag_model

* Created a disag_predictions class which is returned by the predict function and contains the 
  mean and uncertainty predictions. This has replaced the predictions and uncertainty classes. 
  Plot, summary and print methods have been implemented for the disag_predictions class

* Extracted the function make_model_object to allow the user to make a TMB model object on its own, 
  so it can be used in different optimiser or a for MCMC
  
* Neatened up plot.disag_data function to produce 3 plots on the same canvas, with an optional which 
  argument for the user to choose which plots to display

* Made the summary and print function return different outputs. Print functions show minmial output, 
  summary function are more deatiled

## Test environments
* local Windows 10, R 3.6.1
* Ubuntu 16.04.6 LTS (on travis-ci, devel and release) 
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs. 

There were 3 NOTEs:

* checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Anita Nandi <anita.k.nandi@gmail.com>'

  Suggests or Enhances not in mainstream repositories:
    INLA
  Availability using Additional_repositories specification:
    INLA   yes   https://inla.r-inla-download.org/R/stable

  The Title field starts with the package name.
  
  I have not submitted a package before. These are not mis-spelled words. The package uses INLA, my understanding of this NOTE is that it is fine.

* checking installed package size ... NOTE
    installed size is 12.8Mb
    sub-directories of 1Mb or more:
      libs  12.5Mb

  Packages based on C++ can have large compiled libraries. This is as small as it can be, hope that is ok. I got a similar, but slightly different note when using R CMD check compared to devtools::check(). The gist was the same though.

* checking compilation flags used ... NOTE
  Compilation used the following non-portable flag(s):
    '-Wa,-mbig-obj'
    
  To compile large C++ source files on Windows a compilation flag is needed

## Downstream dependencies
There are currently no downstream dependencies for this package