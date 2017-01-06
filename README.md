
### Updated R code from Gasparrini StatMed 2014

--------------------------------------------------------------------------------

An example illustrating the extension of DLNMs for modelling exposure-lag-response associations beyond time series analysis. The code completely reproduces the examples and simulation study described in the article:

Gasparrini A. Modeling exposure–lag–response associations with distributed lag non-linear models. *Statistics in Medicine*. 2014;**33**(5):881-899. [[freely available here](http://www.ag-myresearch.com/2014_gasparrini_statmed.html)]

This is complemented by the vignette *dlnmExtended* included in the R package `dlnm`, showing applications in alternative settings.

--------------------------------------------------------------------------------

The code:

  * *uminers.csv* stores the data from the Colorado Plateau uranium miners cohort, including individual information for 3,347 male subjects
  * the numbered files from *00.prep.R* to *08.simres.R* reproduce the results of the illustrative example and of the simulation study
  * *example.R* offers a simple example with of fitting some models and displaying/summarizing the results
  * *functions.R* creates some functions used in the other scripts
  * *simprep.R* and *simrun.R* are called from the other scripts for preparing and running the simulations
  
Download as a ZIP file using the green button *Clone or download* above