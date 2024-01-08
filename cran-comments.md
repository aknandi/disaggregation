## Update
This is a package update (version 0.2.0). The only real change in this version
is updating references to our Journal of Statistical Science paper that is in 
press.




## Test environments
Windows, R release
Ubuntu 20, R release
Ubuntu 20, r Oldrel
Ubuntu 20, R devel


## R CMD check results
There were no ERRORs or WARNINGs. 

There were 3 NOTES.


* checking CRAN incoming feasibility ... [14s] NOTE
Maintainer: 'Tim Lucas <timcdlucas@gmail.com>'

Possibly misspelled words in DESCRIPTION:
  Nandi (15:28)

Suggests or Enhances not in mainstream repositories:
  INLA
Availability using Additional_repositories specification:
  INLA   yes   https://inla.r-inla-download.org/R/stable

Found the following (possibly) invalid DOIs:
  DOI: 10.18637/jss.v106.i11
    From: DESCRIPTION
          inst/CITATION
    Status: 404
    Message: Not Found


Examples with CPU (user + system) or elapsed time > 10s
               user system elapsed
getPolygonData 9.89   0.17   10.08



Response: Anita Nandi's name is spelled correctly. The INLA availability
issue is the same as previous submissions. The doi is for our new Journal
of the Statistical Society paper and has been reserved but not registered yet.



* checking package dependencies ... NOTE
Package suggested but not available for checking: 'INLA'

Response: Same as above.


* checking examples ... [16s] NOTE
Examples with CPU (user + system) or elapsed time > 10s
               user system elapsed
getPolygonData 9.89   0.17   10.08



Response: As this is only just over the 10 second limit we hope it is ok. We 
have done our best to make the examples small throughout.


## Downstream dependencies
There are currently no downstream dependencies for this package
