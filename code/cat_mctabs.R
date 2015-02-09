

# function to print the  Monte Carlo tables for the thresholded LAsso paper.
# Called from the knitr files using the saved mc statistics array.
# this version of the function is for a fixed tau
catmc_2col <- function(tmc,rfreq,clabs,statlabs=NULL,rlab1=NULL,rlab2=NULL,label,cdrop,rnd=2,tab_caption='',col_hide=NULL){
	
	

	# preparing the row names
	rn <- do.call(rbind,strsplit(rownames(tmc),'_'))
	rn <- rn[,-cdrop]
	
	if(is.null(rlab1)) rlab1 <- unique(rn[,1])
	if(is.null(rlab2)) rlab2 <- unique(rn[,2])
	
	coltype <- c('l','l','l',rep('r',ncol(tmc)))
	if(!is.null(col_hide))coltype[col_hide] <- 'H'
	
	cat('\\begin{table}\n')
	cat('\\centering\n')
	cat('\\begin{tabular}{',coltype, '}\n',sep=' ')
	cat('\\toprule\n')
	
	#column label 
#	cat(paste(clabs,collapse = ' & '))
	
	if(is.null(statlabs))statlabs <- colnames(tmc)

	cl <- c(head(statlabs,-2))
	cl <- paste('\\rotpar{',cl)
	cl <- paste(cl,'}')
	cl <- c(clabs,cl,tail(statlabs,2))
	cat(paste(cl,collapse = ' & '))

	#cat(paste(statlabs,collapse = ' & '))
	cat('\\\\\n')

	cat('\\midrule\n')
	
	c1 <- c2 <- 0
	
	for(i in 1:nrow(tmc)){
		
	# Grey Thresholded Lasso rows
	cat(ifelse((i%%2)==0,' \\rowcolor{black!10} ',' '))
	
	# booleans to check whether a multirow label should be printed.
	mr1 <- ((i%%rfreq[1])==0)
	mr2 <- ((i%%rfreq[2])==0)
	# Row label counters
	if(mr1)c1<-c1+1
	if(mr2)c2<-c2+1
	if((mr2)&(c2>length(rlab2))) c2 <- 1 # reset the nested counter if needed
	# printing multirow labels	
	cat(ifelse(mr1,paste0('\\cellcolor{white} \\multirow{-',rfreq[1],'}{*}{',rlab1[c1],'} & ')
			   ,'\\cellcolor{white}  & '))
	cat(ifelse(mr2,paste0('\\cellcolor{white} \\multirow{-',rfreq[2],'}{*}{',rlab2[c2],'} & ')
			   ,'\\cellcolor{white}  & '))
	# print estimator label
	cat(ifelse((i%%2)==0,' T & ',' L & '))
	
	# print content
	for(cl in 1:ncol(tmc)){
		if(cl>1)cat(' & ')
		if(is.na(tmc[i,cl]))cat(' - ')
		else {
			if(cl!=4) cat(formatC(tmc[i,cl],digits=2,format='f'))
			if(cl==4) cat(formatC(tmc[i,cl],format='d'))	
		}
	}

	cat('\\\\\n')
	
	if(mr1 & (i<nrow(tmc)) )	cat(paste('\\cmidrule(l){1 - ',ncol(tmc)+3,'}\n',sep=''))
		
	}



	cat('\\bottomrule\n')  		
	cat('\\end{tabular}')
	cat(paste('\\caption{',tab_caption,'}',sep=''))
	
	cat(paste('\\label{tab:',label,'}',sep=''))

	
	cat('\\end{table}')
}
