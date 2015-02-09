
lpar1 <- list(c(0,2,2,2,2,2,rep(0,95)) )


lpar2 <- list(c(2,-2,2,-2,2,rep(0,95)))

vtau <- c(0.3,0.5)
vobs <- c(50,100,200,500,1000)



mca <- mc_ttlas(lpar=lpar1,lpar2=lpar2,vobs=vobs,vtau=vtau,taugrid=taugrid,Cgrid=Cgrid,iter=iter,ncores = ncores,sig_eps=0.25)

mmc <- melt(mca)
tmc <- acast(mmc,par+tau+nobs+Estimator~stats)


save(tmc,file = 'mc_code/save/xp_smpl')
