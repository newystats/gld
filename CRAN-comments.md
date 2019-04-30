CRAN comments

## Test environment:
* local OS X install, R 3.5.1
* win-builder (release and devel)

## R CMD check results
There were no ERRORs or WARNINGs 

There was one NOTE
Found the following (possibly) invalid URLs:
  URL: http://nova.newcastle.edu.au/vital/access/manager/Repository/uon:13503
    From: man/fit.fkml.Rd
    Status: Error
    Message: libcurl error code 56:
      	SSL read: error:00000000:lib(0):func(0):reason(0), errno 10054
  URL: https://www.jstor.org/stable/1267789
    From: man/fit.fkml.Rd
    Status: 403
    Message: Forbidden
