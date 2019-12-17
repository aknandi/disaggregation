## Resubmission
This is a resubmission. In this version I have:

* Version 0.1.1 was accepted but subsequently failed a single build (r-devel-windows-ix86+x86_64-gcc8). These errors have now been fixed, details below. The version is now 0.1.2

  In tests, prepare_data function explicity only builds mesh if not on CRAN (if(identical(Sys.getenv("NOT_CRAN"), "true")))
  Affects: test-fit-model, test-plotting, test-predict-model
  
* Version 0.1.0 was accepted but subsequently failed some builds. These errors have now been fixed, details below. The version is now 0.1.1

  disaggregation.cpp:207:18: warning: explicitly assigning value of variable of type 'vector<Type>' to itself [-Wself-assign-overloaded]
  This line has been removed
  disaggregation.cpp:99:22: error: call of overloaded 'sqrt(int)' is ambiguous
    Type kappa = sqrt(8) / rho;
  Changed this to Type kappa = sqrt(8.0) / rho;

* Omitted the redundant and rather unspecific part "Useful Functions for" in the DESCRIPTION title.

* Added a useful reference for disaggregation modelling in the Description field of the DESCRIPTION file

  A useful reference for disaggregation modelling is Lucas et al. (2019) <doi:10.1101/548719>.
  
* Added missing Rd-tags: \value
  
  Affects: build_mesh, parallelExtract, getPolygonData, getCovariateRasters, getStartendindex, predict_model, predict_uncertatinty,
  as.disag.data, plotting functions and summary functions

* Remove Roxygen comments from functions that are not exported
  
  Affects: extractCoordsForMesh, plot_polygon_data, getCoords and getAmatrix

* Replaced the installed.packages() function with requireNamespace()
  
  Affects: prepare_data function and vignettes/disaggregation.Rmd

* Skipped tests on CRAN as they were taking too long

  The timings used to be:
  ** running tests for arch 'i386' ... [298s] OK
  ** running tests for arch 'x64' ... [305s] OK
  I have now reduced this by more than half:
  ** running tests for arch 'i386' ... [106s] OK
  ** running tests for arch 'x64' ... [101s] OK

* \dontrun{} is only used if the example really cannot be executed. 

  Examples that are the most interesting use INLA so we don't run these examples. 
  One example takes longer than 5 seconds and so \donttest{} is used instead.

## Test environments
* local Windows 10, R 3.6.1
* Ubuntu 16.04.6 LTS (on travis-ci, devel and release) 
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs. 

There were 3 NOTEs:

* checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Anita Nandi <anita.k.nandi@gmail.com>'

  New submission

  Possibly mis-spelled words in DESCRIPTION:
    Disaggregation (3:8)
    al (15:14)
    disaggregation (11:19, 14:82)
    et (15:11)

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