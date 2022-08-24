## Update
This is a package update (version 0.1.4). The changes in this version are:

* Change maintainer. Anita Nandi has emailed to confirm. Anita has moved industry and no longer has time to maintain this package.

* Fixed mistake in model definition. We were adjusting the jacobian for a change of variables incorrectly.

* Fixed predictions in models with no field

* Better documentation for priors.

* redocument to fix html5 issues.



## Test environments
Windows, R release
Ubuntu 20, R release
Ubuntu 20, r Oldrel
Ubuntu 20, R devel


## R CMD check results
There were no ERRORs or WARNINGs. 

There were 3 NOTEs:


  Suggests or Enhances not in mainstream repositories:
    INLA
  Availability using Additional_repositories specification:
    INLA   yes   https://inla.r-inla-download.org/R/stable

  The package uses INLA, my understanding of this NOTE is that it is fine.

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
