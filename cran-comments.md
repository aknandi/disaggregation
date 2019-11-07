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
    TMB (11:58)
    disaggregation (3:29, 11:19)

  Suggests or Enhances not in mainstream repositories:
    INLA
  Availability using Additional_repositories specification:
    INLA   yes   https://inla.r-inla-download.org/R/stable
  
  I have not submitted a package before. These are not mis-spelled words. The package uses INLA, my understanding of this NOTE is that it is fine.

* checking installed package size ... NOTE
    installed size is 12.8Mb
    sub-directories of 1Mb or more:
      libs  12.5Mb

  Packages based on C++ can have large compiled libraries

* checking compilation flags used ... NOTE
  Compilation used the following non-portable flag(s):
    '-Wa,-mbig-obj'
    
  To compile large C++ source files on Windows a compilation flag is needed

## Downstream dependencies
There are currently no downstream dependencies for this package