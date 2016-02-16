## Resubmission
This is a resubmission. In this version I have:

* Updated inline citations to include the DOI when possible. In some cases the
  referenced paper has not yet been published (so no DOI exists), so instead of
  the DOI the ArXiv preprint identification number is given.

## Test environments
* local OS X install, R 3.2.3
* Using devtools::build_win()

## R CMD check results
There were no ERRORs or WARNINGs.

There were 2 NOTEs:

* checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Luca Weihs <lucaw@uw.edu>'
  
  Possibly mis-spelled words in DESCRIPTION:
    Bergsma (3:65, 11:37)
    Dassios (4:13, 11:49)
    Drton (14:29)
    Nandy (14:11)
    Weihs (14:18)
  
  The Title field should be in title case, current version then in title case:
  'Efficient Computation and Testing of the t* Statistic of Bergsma and Dassios'
  'Efficient Computation and Testing of the T* Statistic of Bergsma and Dassios'
  
  None of the above words are mis-spelled (they are all names) and t* is 
  intentionally lower-case as it is the name of a particular statistic.

## Downstream dependencies
There are currently no downstream dependencies for this package.