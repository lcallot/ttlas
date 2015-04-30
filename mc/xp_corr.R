
lpar1 <- list(c(0,rep(2,5),rep(0,45)))
lpar2 <- list(c(2,-2,2,-2,2,rep(0,45)))


vtau <- c(0.3,0.5)
vcorrXQ <- c(2,0,0.5,0.95)
vobs <- c(200)


mctmp<- NULL

for(corrXQ in vcorrXQ){
	mca <- mc_ttlas(lpar=lpar1,lpar2=lpar2,vobs=vobs,vtau=vtau,
					taugrid=taugrid,Cgrid=Cgrid,iter=iter,
					ncores = ncores,sig_eps=0.25,corrXQ=corrXQ)
					
	if(is.null(mctmp)) mctmp <- array(NA,dim=c(dim(mca),length(vcorrXQ)),
								  dimnames = c(dimnames(mca),list('corrXQ' = vcorrXQ)))
	
	mctmp[,,,,,as.character(corrXQ)] <- mca
	
}



mmc <- melt(mctmp)
tmc <- acast(mmc,par+corrXQ+tau+nobs+Estimator~stats)


save(tmc,file = 'mc/save/xp_corr')
