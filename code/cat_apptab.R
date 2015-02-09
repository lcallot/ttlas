
# function to print the tables for the application in the thresholded LAsso paper. Called from the knitr files after estimation of the model. 
prttab <- function(param,ndebt,tmin,lmbd,year,Cest,label,diag_dbt = TRUE, fe_line = NULL, smpl_line = NULL){
	
	coltype <- c('l','l',rep('r',2*length(ndebt)))
	
	cat('\\begin{table}\n')
	cat('\\begin{tabular}{',coltype, '}\n',sep=' ')
	cat('\\toprule\n')
	
	#column label 
	cat(' & Threshold:')
	for(n in ndebt){cat('& \\multicolumn{2}{c}{',n,'}')}
	cat('\\\\\n')
	cat('& & ')
	cat(paste(rep('L & T',length(ndebt)),collapse='&'))
	cat('\\\\\n')
	cat('\\midrule\n')
	
	# content
	np <- nrow(param)/2
	ptmp <- param[1:(np-ifelse(diag_dbt,1,0)),]	
	ptmp <- round(ptmp,3)
	ptmp[ptmp==0] <- ''
	rn <- rownames(ptmp)
	
	for(i in 1:nrow(ptmp)){
		if(i==1)cat('\\multirow{',nrow(ptmp) + ifelse(diag_dbt,-1+length(ndebt),0),'}{*}{$\\hat{\\beta}$} & ')
		if(i!=1)cat(' & ')
		cat(rn[i])
		cat(' & ')
		cat(paste(ptmp[i,],collapse = ' & '))
		cat('\\\\\n')	
	}
	
	if(diag_dbt){
		clmat <- kronecker(diag(length(ndebt)),matrix(rep(1,2),nrow=1))
		
		dd <- round(clmat%*%diag(param[np,]),3)
		dd[clmat==0] <- '-'
		dd[dd=='0'] <- ''
		for(i in 1:nrow(dd)){
			cat('&')
			cat(ndebt[i])
			cat('&')
			cat(paste(dd[i,],collapse = '&'))
			cat('\\\\\n')
		}
	}
	
	cat(paste('\n\\cmidrule(lr){1 - ',2*length(ndebt)+2,'}\n\n',sep=''))
	
	ptmp <- param[(np+1):(2*np-ifelse(diag_dbt,1,0)),]	
	ptmp <- round(ptmp,3)
	ptmp[ptmp==0] <- ''
	#rn <- rownames(ptmp)
	
	for(i in 1:nrow(ptmp)){
		if(i==1)cat('\\multirow{',nrow(ptmp) + ifelse(diag_dbt,-1+length(ndebt),0),'}{*}{$\\hat{\\delta}$} & ')
		if(i!=1)cat(' & ')
		cat(rn[i])
		cat(' & ')
		cat(paste(ptmp[i,],collapse = ' & '))
		cat('\\\\\n')	
	}
	
	if(diag_dbt){
		dd <- round(clmat%*%diag(param[2*np,]),3)
		dd[clmat==0] <- '-'
		dd[dd=='0'] <- ''
		for(i in 1:nrow(dd)){
			cat(paste(' & ',ndebt[i],sep=''))
			cat('&')
			cat(paste(dd[i,],collapse = '&'))
			cat('\\\\\n')
		}
	}
	
	cat('\\midrule\n')
	cat('& $\\widehat{\\tau}$ & ')
	cat(paste(round(rep(tmin,each=2),2),collapse = '&'))
	cat('\\\\\n')
	cat('& $\\widehat{\\lambda}$ & ')
	cat(paste(round(rep(lmbd,each=2),3),collapse = '&'))
	cat('\\\\\n')
	cat('& $\\widehat{C}$ & ')
	cat(paste(' - & ',Cest,collapse = '&'))
	cat('\\\\\n')
	if(!is.null(smpl_line))
	{cat(smpl_line)
	 cat('\\\\\n')}
	if(!is.null(fe_line))
	{cat(fe_line)
	 cat('\\\\\n')}
	
	cat('\\bottomrule\n')  		
	cat('\\end{tabular}')
	cat(paste('\\caption{ Estimated parameters, year: ',year,'. Empty cells are parameters set to zero, dashes indicate parameters not included in the model.}',sep=''))
	
	cat(paste('\\label{',label,'}',sep=''))
	
	cat('\\end{table}')
}



# function to print the tables for the application in the thresholded LAsso paper. Called from the knitr files after estimation of the model. 
prttab_nothdint <- function(param,ndebt,tmin,lmbd,year,Cest,label,diag_dbt = TRUE, fe_line = NULL, smpl_line = NULL){
	
	coltype <- c('l','l',rep('r',2*length(ndebt)))
	
	cat('\\begin{table}\n')
	cat('\\begin{tabular}{',coltype, '}\n',sep=' ')
	cat('\\toprule\n')
	
	#column label 
	cat(' & Threshold:')
	for(n in ndebt){cat('& \\multicolumn{2}{c}{',n,'}')}
	cat('\\\\\n')
	cat('& & ')
	cat(paste(rep('L & T',length(ndebt)),collapse='&'))
	cat('\\\\\n')
	cat('\\midrule\n')
	
	# content
	np <- (nrow(param)-1)/2
	ptmp <- param[1:(np+1-ifelse(diag_dbt,1,0)),]	
	ptmp <- round(ptmp,3)
	ptmp[ptmp==0] <- ''
	rn <- rownames(ptmp)
	
	for(i in 1:nrow(ptmp)){
		if(i==1)cat('\\multirow{',nrow(ptmp) + ifelse(diag_dbt,-1+length(ndebt),0),'}{*}{$\\hat{\\beta}$} & ')
		if(i!=1)cat(' & ')
		cat(rn[i])
		cat(' & ')
		cat(paste(ptmp[i,],collapse = ' & '))
		cat('\\\\\n')	
	}
		
	if(diag_dbt){
		clmat <- kronecker(diag(length(ndebt)),matrix(rep(1,2),nrow=1))
		
		dd <- round(clmat%*%diag(param[np,]),3)
		dd[clmat==0] <- '-'
		dd[dd=='0'] <- ''
		for(i in 1:nrow(dd)){
			cat('&')
			cat(ndebt[i])
			cat('&')
			cat(paste(dd[i,],collapse = '&'))
			cat('\\\\\n')
		}
	}
	
	cat(paste('\n\\cmidrule(lr){1 - ',2*length(ndebt)+2,'}\n\n',sep=''))
	
	ptmp <- param[(np+1):(2*np-ifelse(diag_dbt,1,0)),]	
	ptmp <- round(ptmp,3)
	ptmp[ptmp==0] <- ''
	#rn <- rownames(ptmp)
	
	for(i in 1:nrow(ptmp)){
		if(i==1)cat('\\multirow{',nrow(ptmp) + ifelse(diag_dbt,-1+length(ndebt),0),'}{*}{$\\hat{\\delta}$} & ')
		if(i!=1)cat(' & ')
		cat(rn[i+1])
		cat(' & ')
		cat(paste(ptmp[i,],collapse = ' & '))
		cat('\\\\\n')	
	}
	
	if(diag_dbt){
		dd <- round(clmat%*%diag(param[2*np,]),3)
		dd[clmat==0] <- '-'
		dd[dd=='0'] <- ''
		for(i in 1:nrow(dd)){
			cat(paste(' & ',ndebt[i],sep=''))
			cat('&')
			cat(paste(dd[i,],collapse = '&'))
			cat('\\\\\n')
		}
	}
		
	cat('\\midrule\n')
	cat('& $\\widehat{\\tau}$ & ')
	cat(paste(round(rep(tmin,each=2),2),collapse = '&'))
	cat('\\\\\n')
	cat('& $\\widehat{\\lambda}$ & ')
	cat(paste(round(rep(lmbd,each=2),3),collapse = '&'))
	cat('\\\\\n')
	cat('& $\\widehat{C}$ & ')
	cat(paste(' - & ',Cest,collapse = '&'))
	cat('\\\\\n')
	if(!is.null(smpl_line))
		{cat(smpl_line)
		cat('\\\\\n')}
	if(!is.null(fe_line))
		{cat(fe_line)
		cat('\\\\\n')}
	
	cat('\\bottomrule\n')  		
	cat('\\end{tabular}')
	cat(paste('\\caption{ Estimated parameters, year: ',year,'. Empty cells are parameters set to zero, dashes indicate parameters not included in the model.}',sep=''))
	
	cat(paste('\\label{',label,'}',sep=''))
	
	cat('\\end{table}')
}
