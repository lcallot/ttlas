
lpar1 <- list(c(0,rep(2,5),rep(0,45)),
			  c(0,rep(2,5),rep(0,95)),
			  c(0,rep(2,5),rep(0,195)),
			  c(0,rep(2,5),rep(0,395))
			  )


lpar2 <- list(c(rep(0,50)),
			  c(rep(0,100)),
			  c(rep(0,200)),
			  c(rep(0,400))
			)

vtau <- c(0.5)
vobs <- c(200)


mca <- mc_ttlas(lpar=lpar1,lpar2=lpar2,vobs=vobs,vtau=vtau,taugrid=taugrid,Cgrid=Cgrid,iter=iter,ncores = ncores,sig_eps=0.25)

mmc <- melt(mca)
tmc <- acast(mmc,par+tau+nobs+Estimator~stats)


save(tmc,file = 'mc_code/save/xp_nj')
