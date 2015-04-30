
pv1 <- c(0,1,1,1,1,1,rep(0,95))
pv2 <- c(1,-1,1,-1,1,rep(0,95))

lpar1 <- list(0.5*pv1)
lpar2 <- list(0.5*pv2)


vtau <- c(0.5)
vobs <- c(100,200,1000)


mca <- mc_ttlas(lpar=lpar1,lpar2=lpar2,vobs=vobs,vtau=vtau,taugrid=taugrid,Cgrid=Cgrid,iter=iter,ncores = ncores,sig_eps=0.25)

mmc <- melt(mca)
tmc <- acast(mmc,par+tau+nobs+Estimator~stats)


save(tmc,file = 'mc/save/xp_psize5')

