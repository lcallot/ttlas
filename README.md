# Replication material for: Sharp threshold detection based on sup-norm error rates in high-dimensional models.
### Laurent Callot, Mehmet Caner, Anders Bredahl Kock, and Juan Andres Riquelme.


---

Author: Laurent Callot (l.callot@vu.nl)

Date: 08/02/2015


This repository contains the replication material for the simulations and the empirical application in _Thresholded Lasso: Variable Selection in High Dimensions_. All the computations are carried using *R* and the *knitr* package for easy replication. The material is divided between 4 folders:

 - application: contains a *knitr* file, **application.Rnw**, to be compiled to replicate the results of the application. 
 - code: contains 3 files. **tlas.R** contains the code for the thresholded Lasso and the simulations. The other two files contain functions to print the result to a tex table. 
 - data: **BIS\_data.Rda** contains a 3 dimensional array with the data used in the application. 
 - mc: contains the files for the simulations. 5 files called *xp_XXX.R* contain the setup and call to *ttlas* for each of the 5 experiments. *xp_all.R* contains the global settings for the simulations and runs the experiments. It is the file that should be ran to compute the simulations, it takes around 2 hours to complete computations on 16 cores for 1000 replications. **mc.Rnw** should be compiled to generate the tables in the paper from the saved simulations statistics located in the mc/save/ folder. 