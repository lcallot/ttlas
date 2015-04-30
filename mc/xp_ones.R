
lpar1 <- list(c(0,rep(2,1),rep(0,100)),
			  c(0,rep(2,5),rep(0,100)),
			  c(0,rep(2,10),rep(0,100)),
			  c(0,rep(2,25),rep(0,100))
				)
			  


vtau <- c(0.3,0.4,0.5)
vobs <- c(200)


mca <- mc_ttlas(lpar=lpar1,lpar2=NULL,vobs=vobs,vtau=vtau,taugrid=taugrid,Cgrid=Cgrid,iter=iter,ncores = ncores,sig_eps=0.25)

mmc <- melt(mca)
tmc <- acast(mmc,par+tau+nobs+Estimator~stats)


save(tmc,file = 'mc/save/xp_ones')
