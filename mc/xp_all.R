library('glmnet')
library(parallel)
library(reshape2)
library(ggplot2)
source('code/tlas.R')


set.seed(42) #duh

# Number of cores for mclapply
ncores <- 8

# Search grids
taugrid <- seq(0.15,0.85,0.05)
Cgrid <- seq(0.1,5,0.1)

# Monte Carlo iterations
iter <-1000 

cat('\nTable Corr(X,Q)\n')
source('mc/xp_corr.R')

cat('\nTable psize\n')
#source('mc/xp_psize.R')

cat('\nTable 1\n')
#source('mc/xp_seo.R')

cat('\nTable no jump\n')
#source('mc/xp_nj.R')

cat('\nTable smpl\n')
#source('mc/xp_smpl.R')

cat('\nTable ones\n')
#source('mc/xp_ones.R')
