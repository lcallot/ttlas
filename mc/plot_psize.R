library('glmnet')
library(parallel)
library(reshape2)
library(ggplot2)
source('code/tlas.R')


set.seed(42) 

# Number of cores for mclapply
ncores <-16

# Search grids
taugrid <- seq(0.15,0.85,0.05)
Cgrid <- seq(0.1,5,0.1)

# Monte Carlo iterations
iter <-1000

# base param
pv1 <- c(0,1,1,1,1,1,rep(0,45))
pv2 <- c(1,-1,1,-1,1,rep(0,45))
# parameter list
lpar2 <- list(0.1*pv2,
			  0.2*pv2,
			  0.3*pv2,
			  0.4*pv2,
			  0.5*pv2,
			  0.75*pv2,
			  pv2,
			  1.5*pv2,
			  2*pv2)



lpar1 <- rep(list(pv1),length(lpar2))


vtau <- c(0.5)
vobs <- c(100,500)

mca <- mc_ttlas(lpar=lpar1,lpar2=lpar2,vobs=vobs,vtau=vtau,taugrid=taugrid,Cgrid=Cgrid,iter=iter,ncores = ncores,sig_eps=0.25)

mmc <- melt(mca)
mmc <- subset(mmc,stats%in%c('fnz','fz','ps'))
mmc$stats <- factor(mmc$stats)

# Names for the facets
levels(mmc$stats) <- c('False positive (nbr.)','False Negative (nbr.)','Perfect selection (pct.)')
mmc$nobs <- paste0('n=',mmc$nobs)

# Proper x-axis scale
tmp <- factor(mmc$par)
levels(tmp) <- c(0.1,0.2,0.3,0.4,0.5,0.75,1,1.5,2)
mmc$par <- as.numeric(as.vector(tmp))

# Estimator names
levels(mmc$Estimator) <- c('Lasso','Thresholded Lasso')

pplot <- ggplot(mmc,aes(x=par,y=value)) +
	geom_line(aes(colour=Estimator)) + 
	geom_point(aes(shape=Estimator,colour=Estimator)) + 
	scale_colour_manual(values=c('Lasso'='grey65','Thresholded Lasso'='black')) + 
	scale_shape_manual(values=c(44,19)) + 
	theme_bw() + labs(x=expression(paste('Scale of ',delta))) + 
	theme(legend.position="bottom") + 
	facet_grid(stats ~ nobs ,scale='free')
print(pplot)

ggsave(filename = 'mc/pplot.pdf',plot = pplot,width = 18,height=16,units='cm')
