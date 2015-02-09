

# Simulation function for the threshdolded Lasso

mc_ttlas <- function( lpar,lpar2=NULL,vobs,vtau, taugrid , Cgrid , iter,ncores=1,sig_eps=1)
{	
	
	tt <- proc.time()
	
	# Names of the stats
	sn <- c('rmse','fnz','fz','ps','nnz','al1','alinf','tl1','C','lmbd')
	# Creating the storage array
	astat <- array(NA, dim=c(2,length(vobs),length(vtau),length(lpar),length(sn))
				   ,dimnames = list('Estimator' = c('las','tlas'),'nobs'=vobs,
				   					'tau'=vtau,'par'=NULL,'stats'=sn) )
	
	# Looping over the parameter specifications
	pcount <- 0 # specification counter
	for(par1 in lpar){
		pcount <- pcount + 1
		if(is.null(lpar2))par2 <- par1[-1]
		if(!is.null(lpar2))par2 <- lpar2[[pcount]]
		# Looping over the threshold locations:
		for(tau in vtau){
			# Looping over number of observations:
			for(nobs in vobs){
				# Estimate the model
				rt <- proc.time()
				#mctlas <- lapply(1:iter,ittlas, par1,par2,nobs,tau,y=y,xthd=xthd
				mctlas <- mclapply(1:iter,ittlas, par1,par2,nobs,tau,y=y,xthd=xthd
								   ,x=NULL,thdvar=thdvar,taugrid=taugrid,Cgrid=Cgrid,sig_eps=sig_eps
									#)
									,mc.cores=ncores)
				
				cat(paste0('\nn=',nobs,' par:',pcount,' XP time:',round(proc.time()-rt,2)[3]))
			
				# gather the results, compute the sstatistics, store in array
				vpar <- c(par1,par2)
				lstat <- mcstat(mctlas,vpar,tau)
				
				astat[,as.character(nobs),as.character(tau),pcount,'rmse'] <- rowMeans(lstat$rmse)
				astat[,as.character(nobs),as.character(tau),pcount,'fnz'] <- rowMeans(lstat$fnz)
				astat[,as.character(nobs),as.character(tau),pcount,'fz'] <- rowMeans(lstat$fz)
				astat[,as.character(nobs),as.character(tau),pcount,'ps'] <- 100*rowMeans(lstat$ps)
				astat[,as.character(nobs),as.character(tau),pcount,'nnz'] <- rowMeans(lstat$nnz)
				astat[,as.character(nobs),as.character(tau),pcount,'al1'] <- rowMeans(lstat$al1)
				astat[,as.character(nobs),as.character(tau),pcount,'alinf'] <- rowMeans(lstat$alinf)
				
				astat[1,as.character(nobs),as.character(tau),pcount,'tl1'] <- mean(lstat$tl1)
				astat[1,as.character(nobs),as.character(tau),pcount,'lmbd'] <- mean(lstat$lmbd)
				astat[2,as.character(nobs),as.character(tau),pcount,'C'] <- mean(lstat$C)
			}
		}
	}
	
	cat(paste0('\n\nTotal run time:',round((proc.time()-tt)/60,2)[3],' minutes. \n'))
	
	return(astat)
	
}

# A function taking a MC run, returning an array of stats
mcstat <- function(mctlas,vpar,tau)
{
	
	#init
	iter <- length(mctlas)
	lmbd <- Cest <- tl1  <- c()
	fnz <- nnz <- fz <- ps <- rmse <- al1 <- alinf <- matrix(0,2,iter)
	
	#loop over iterations
	for(i in 1:iter)
		{
		mc <- mctlas[[i]]
		# Estimation parameters
		gmin <- which.min(mc$Vobj)
		tl1 <- c(tl1,abs(mc$taugrid[gmin]-tau))
		lmbd <- c(lmbd,mc$lambda[gmin])
		Cest <- c(Cest,mc$Cmin[gmin])
	
		# Parameter storage
		npar  <- ncol(mc$X)
		alpha <- mc$alptau[1:npar,gmin]
		alpht <- mc$alpthd[1:npar,gmin]		
		
		# parameter l1 estimation error
		al1[1,i]   <- sum(abs(alpha-vpar))
		al1[2,i]   <- sum(abs(alpht-vpar))
		
		#l-infinity
		alinf[1,i] <- max(abs(alpha-vpar))
		alinf[2,i] <- max(abs(alpht-vpar))
		
		# root mean square prediction error
		rmse[1,i] <- (mean((mc$y-mc$X%*%alpha)^2))
		rmse[2,i] <- (mean((mc$y-mc$X%*%alpht)^2))
		
		# parameter selection
		# false non zero
		fnz[1,i] <- sum(((vpar[-1]==0)&(alpha[-1]!=0)))
		fnz[2,i] <- sum(((vpar[-1]==0)&(alpht[-1]!=0)))
		# false zero
		fz[1,i] <- sum(((vpar[-1]!=0)&(alpha[-1]==0)))
		fz[2,i] <- sum(((vpar[-1]!=0)&(alpht[-1]==0)))
		# perfect sel
		ps[1,i] <- sum((vpar[-1]==0)!=(alpha[-1]==0))==0
		ps[2,i] <- sum((vpar[-1]==0)!=(alpht[-1]==0))==0
		# total non zero
		nnz[1,i] <- sum(alpha[-1]!=0)
		nnz[2,i] <- sum(alpht[-1]!=0)
		
		}
	
	lstat <- list('rmse'=rmse,'fnz'=fnz,'fz'=fz,'ps'=ps,'nnz'=nnz,'alinf'=alinf,'al1'=al1,'tl1'=tl1,'C'=Cest,'lmbd'=lmbd)
	
	return(lstat)	
}


# generates the data
# calls the ttlas function
# called by the mc_ttlas function
ittlas <- function(i, par1,par2,nobs,tau , y , xthd , x=NULL , thdvar , taugrid , Cgrid , sig_eps)
{
	nvars <- length(par1) # nbr vars (incl. cste)
	
	thdvar <- matrix(rnorm(nobs),ncol=1) # generate the threshold variable
	X1 <- matrix(rnorm(nobs*(nvars-1)),ncol=nvars-1) # generating the regessors
	X2 <- X1
	X2[thdvar>tau,] <- 0 # threshold regressors
	
	
	y <- par1[1] + X1%*%par1[-1] + X2%*%par2 + sqrt(sig_eps)*rnorm(nobs) #generate the data
	xthd <- cbind(X1,X2) # regressor matrix	
	
	est <- ttlas(y , xthd = xthd, x=NULL , thdvar = thdvar, taugrid=taugrid, C=Cgrid, intercept = TRUE,
				 thd_intercept = FALSE, BIC = TRUE, standardize = FALSE)
	est$X <- cbind(1,X1,X2)
	est$y <- y
	return(est)
}

# core function to estimate the threshold model by Lasso and by thresholded Lasso. 
ttlas <- function(y , xthd , x , thdvar ,
				  taugrid , C=1 , nfolds=10 ,
				  intercept = TRUE , thd_intercept = TRUE , BIC = TRUE, standardize = FALSE)
{
	#init 
	n <- length(y)	
	# Storage
	lmbdcv <- c()
	alptau <- alpthd <- ols <- Cmin <- pcoef <- pstdE <- NULL
	Vobj <- c()
	
	# Looping over the taus
	for(cnt in 1:length(taugrid)){
		ta <- taugrid[cnt]
			
		# Construction of x tau
		if(thd_intercept) xtau <- cbind(1,xthd)		
		if(!thd_intercept) xtau <- cbind(0,xthd)		
		xtau <- xtau*matrix(thdvar<=ta,ncol=ncol(xtau),nrow=n)
		
		# Concatenation
		xall <- cbind(xthd,xtau)
		if(!is.null(x)) xall <- cbind(xall,x)	
		# Constructing the weights
		D <- apply(xall,2,function(x)sqrt(mean((x-mean(x))^2)))
		if(!is.null(x)) D[length(D)-c(ncol(x):1)+1] <- 0
		# print(D)
		# D <- NULL
		
		# cross validation for lambda selection
		if(!BIC)cgn	<- cv.glmnet(x = xall, y=y, penalty.factor = D, nfolds=nfolds, standardize = standardize, intercept = intercept)
		if(BIC){cgn	<- glmnet(x = xall, y=y, penalty.factor= D, intercept=intercept, standardize = standardize) 
		    yhat    <- predict(cgn, newx=xall, type = "response")
		    sigmah	<- colSums((matrix(y,nrow(yhat),ncol(yhat))-yhat)^2)/length(y)
		    bic     <- log(sigmah) + cgn$df *log(n)/n * log(log(ncol(xall))) 
		    cgn$lambda.min 	<- cgn$lambda[which.min(bic)] #lambda BIC
			}	
		    
		# Storing
		lmbdcv[cnt] <- cgn$lambda.min
		if(!BIC)alp <- matrix(coef(cgn,s='lambda.min'))
		if(BIC) alp <- cgn$beta[,which.min(bic)]
		if(intercept & BIC) alp<-c(cgn$a0[which.min(bic)],alp)
		alptau <- cbind(alptau,alp)	
		
		# Thresholding the estimate
		# Selecting threshold param by BIC
		rmse2 <- c()
		Cat <- Cbic <- c()
		for(iC in C){
			thd <- cgn$lambda.min *iC
			at  <- matrix(c(alp[1],alp[-1]*(abs(alp[-1])>D*thd)),ncol=1)
			Cat <- cbind(Cat,at) 
			rt  <- y - cbind(1,xall)%*%at
			Cbic<- c(Cbic,log(mean(rt^2)) + sum(at!=0) *log(n)/n * log(log(ncol(xall))))
		}
		# Selecting and storing
		alpthd <- cbind(alpthd,Cat[,which.min(Cbic)] )
		Cmin <- c(Cmin,C[which.min(Cbic)])
		#if(C[which.min(Cbic)]!=min(C)) browser()
		
		# Computing the value of the objective function:
		if(!BIC)yhat	 <- predict(cgn,newx = xall,s='lambda.min')
		RSS 	 <- sum((matrix(y,nrow(yhat),ncol(yhat))-yhat)^2)/(2*n)
		Vobj[cnt]<- RSS + cgn$lambda.min*sum(abs(alptau[-1,cnt]))
		
		# post Lasso
		plas <- summary(lm(y~cbind(1,xall[,Cat[,which.min(Cbic)][-1]!=0])))
		postcoef <- poststdE <- 0*Cat[,which.min(Cbic)]
		postcoef[Cat[,which.min(Cbic)]!=0] <- coefficients(plas)[,1]
		poststdE[Cat[,which.min(Cbic)]!=0] <- coefficients(plas)[,2]
		
		pcoef <- cbind(pcoef,postcoef)
		pstdE <- cbind(pstdE,poststdE)
		#postse
		
		# plain old OLS
		ols <- cbind(ols,coef(lm(y~xall)))
	}
	
	return(list('C'=C,'Vobj'=Vobj,'alptau'=alptau,'alpthd'=alpthd,
				'taugrid'=taugrid,'lambda'=lmbdcv,'ols'=ols,
				'Cmin'=Cmin,'pcoef'=pcoef,'pstdE'=pstdE))
}


