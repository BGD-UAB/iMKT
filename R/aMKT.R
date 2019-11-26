#' @title asymptotic MKT method
#' 
#' @description aMKT: MKT using asymptoticMKT method and estimation of negative selection fractions (d, b, f)
#'
#' @details The asymptotic MKT (aMKT) allows the estimation of the rate of adaptive evolution (alpha) and the diverse negative selection regimens. aMKT uses asymptotic MKT method (Messer and Petrov 2012 PNAS; Haller and Messer 2017 G3) to estimate alpha and the diverse negative selection fractions (d: strongly deleterious, b: weakly deleterious, f: neutral), based on the assumption that weakly deleterious mutations usually do not reach high allele frequencies and therefore, produce the underestimation of alpha at low DAF categories. The fraction of strongly deleterious mutations is estimated as the difference between neutral (0) and selected (i) polymorphic sites relative to the number of analyzed sites: d = 1 - (P0/m0 / Pi/mi). The fraction of weakly deleterious sites (b) corresponds to the relative proportion of selected polymorphic sites that cause the underestimation of alpha at low DAF categories. Finally, the fraction of neutral sites (f) is estimated as: f = 1 - d - b. iMKT() only fits an exponential model for the computation of alpha.
#' 
#' @param daf data frame containing DAF, Pi and P0 values
#' @param divergence data frame containing divergent and analyzed sites for selected (i) and neutral (0) classes
#' @param xlow lower limit for asymptotic alpha fit
#' @param xhigh higher limit for asymptotic alpha fit
#' @param seed seed value (optional). No seed by default
#' @param plot report plots of daf, alpha and negative selection fractions (optional). Default is FALSE
#'
#' @return iMKT method. List with asymptotic MK table and values, fractions of sites and graphs of DAF, asymptotic alpha model and negative selection fractions (optional).
#'
#' @examples
#' ## Without plot
#' iMKT(myDafData, myDivergenceData, xlow=0, xhigh=0.9)
#' ## With plot
#' iMKT(myDafData, myDivergenceData, xlow=0, xhigh=0.9, plot=TRUE)
#' 
#' @import utils
#' @import stats
#' @import ggplot2
#' @import nls2
#' @import MASS
#' @importFrom reshape2 melt
#' @importFrom ggthemes theme_foundation
#' @importFrom cowplot plot_grid
#' 
#' @keywords MKT
#' @export

aMKT = function(daf, divergence, xlow=0, xhigh=1, seed, plot=FALSE) {
	
	## Check data
	check = checkInput(daf, divergence, xlow, xhigh)
	
	if (check$data == FALSE) 
	{
		stop (check$print_errors) 
	}
	
	## Warning P0
	if (any(daf$P0 == 0)){
		warning('Input daf file contains P0 values = 0.\nThis can bias the function fitting and the estimation of alpha.') 
	}
	
	## Check seed
	if(missing(seed)) {
		seed = NULL
	} else {
		set.seed(seed)
	}

	## Create MKT table standard
	mktTableStandard = data.frame(
		Polymorphism = c(sum(daf[['P0']]), sum(daf[['Pi']])), 
		Divergence   = c(divergence[['D0']],divergence[['Di']]),
		row.names    = c('Neutral class','Selected class')
	)
	
	##Polymorphism and Divergence variables
	P0 = sum(daf[['P0']])
	Pi = sum(daf[['Pi']])
	D0 = divergence[['D0']]
	Di = divergence[['Di']]
	m0 = divergence[['m0']]
	mi = divergence[['mi']]

	## Run asymptotic MKT and retrieve alphas 
	## iMKT only fits exponential model in asymptotic alpha. it uses asymptoticMKExp()

	asymptoticMkTable = asymptoticMKExp(daf, divergence, xlow, xhigh)
	alphaAsymptotic   = as.numeric(asymptoticMkTable[['alphaAsymptotic']])
	alphaStandard     = as.numeric(asymptoticMkTable[['alphaOriginal']])
	alphaCiLow        = asymptoticMkTable[['ciLow']] 
	alphaCiHigh       = asymptoticMkTable[['ciHigh']] 
		
	## Estimate alpha for each DAF category
	daf[['alpha']] = 1-((as.numeric(D0)*as.numeric(daf[['Pi']]))/(as.numeric(Di)*as.numeric(daf[['P0']])))
	
	## Estimate the synonymous and non-synonymous ratio
	synonymousRatio    = P0/m0
	nonSynonymousRatio = Pi/mi 
	
	## Estimate the fraction of neutral sites incluiding weakly deleterious variants
	fb  = nonSynonymousRatio/synonymousRatio
	
	## Estimate the fraction of strongly deleleterious sites (d)
	d  = 1 - fb
	
	## Estimate the fraction of sligthly deleterious sites in each daf category (b)
	wd = daf[['Pi']] - (((1 - alphaAsymptotic) * divergence[['Di']] *  daf[['P0']])/divergence[['D0']])
	b  = (sum(wd)/sum(daf[['P0']]))*(m0/mi)

	## Re-estimate the truly number of neutral sites, removing the slightly deleterious 
	f = fb - b

	## Fraction of f, b and d sites

	# Estimate the relative proportion of non-synonymous and synonymous substitutions
	# daf[['N']] = daf[['Pi']]/sum(daf[['Pi']])
	# daf[['S']] = daf[['P0']]/sum(daf[['P0']])

	# wdS = 0
	# for (i in 1:nrow(daf)) {
	# 	row = daf[i,]
	# 	if (row$alpha < alphaCiLow) {
	# 		wdS = wdS + ((alphaAsymptotic-row$alpha)*row$N)
	# 	} else { break }
	# }
	# wdS = wdS/(alphaAsymptotic-min(daf$alpha,na.rm=T))
	# bS = wdS * f3

	fraction = data.frame(d = d, f = f, b = b)

	## Perform plots 
	if(plot == TRUE) {
		
		plotFraction = mutFractionsPlot(fraction)
		dafPoints = plotDaf(daf)
		plotAlpha = plotAsymptotic(asympProccesedDaf = daf,asymptoticDf = asymptoticMkTable,alpha='alpha',xl=xlow,xh=xhigh)

		## Render plots with labels
		plotsiMKT = plot_grid(dafPoints, plotAlpha, plotFraction, nrow=3,  labels=c('A','B','C'), rel_heights=c(2,2,1))
		
		## Store output
		## iMKT output
		
		asymptoticMkTable[2:8] = round(asymptoticMkTable[2:8],4)

		# Delete column necessary to plot in fraction
		output = list(asymptoticMkTable, fraction, plotsiMKT)
		names(output) = c('alphaCorrected', 'fractions', 'graphs')	

	}
	else if (plot == FALSE) {
		## iMKT output
		asymptoticMkTable[2:8] = round(asymptoticMkTable[2:8],4)

		# Delete column necessary to plot in fraction
		output = list(asymptoticMkTable, fraction)
		names(output) = c('alphaCorrected', 'fractions')
		
	}
	
	## Return output
	return(output) 
}


## Two-step nls2() model fit at a given level of precision (res)
fitMKmodel = function(alpha_trimmed, f_trimmed, res) {
	
	## First fitting using starting values (st)
	mod = tryCatch({
		
		## Starting values to fit the model  
		st = expand.grid(const_a=seq(-1,1,length.out=res + 1), const_b=seq(-1,1,length.out=res), const_c=seq(1,10,length.out=res + 1))
		
		## Fitting
		nls2::nls2(alpha_trimmed ~ const_a + const_b * exp(-const_c* f_trimmed), start=st, algorithm='brute-force', control=nls.control(maxiter=NROW(st)))
		
	}, error=function(cond) {}) ## Return condition of error when unable to fit
	
	## If mod fails...
	if (length(mod) == 0) { return(NULL) }
	
	## Second fitting, starting from previous fit (mod)
	mod2 = tryCatch({
		nls2::nls2(alpha_trimmed ~ const_a + const_b * exp(-const_c* f_trimmed), start = mod, control=nls.control(maxiter=200))
		
	}, error=function(cond) {}) ## Same error handling than the previous step
	
	## If mod2 fails...
	if (length(mod2) == 0) { return(NULL) }
	
	## Return mod2 if fitted
	return(mod2)
}

## Compute confidence intervals of alpha using predictNLS 
## Get a CI using Monte Carlo simulation based upon a fitted model.  
## Thanks to Andrej-Nikolai Spiess (http://www.dr-spiess.de) for this code.
predictNLS = function(object, newdata, level = 0.95, nsim = 10000) {
	  
	## get right-hand side of formula
	RHS  = as.list(object$call$formula)[[3]]
	EXPR = as.expression(RHS)
	
	## all variables in model
	VARS = all.vars(EXPR)
	
	## coefficients
	COEF = coef(object)
	
	## extract predictor variable    
	predNAME = setdiff(VARS, names(COEF))  
	
	## take fitted values, if 'newdata' is missing
	if (missing(newdata)) {
		newdata           = eval(object$data)[predNAME]
		colnames(newdata) = predNAME
	}
	  
	## check that 'newdata' has same name as predVAR
	if (names(newdata)[1] != predNAME) stop("newdata should have name '", predNAME, "'!")

	## get parameter coefficients
	COEF = coef(object)

	## get variance-covariance matrix
	VCOV = vcov(object)

	## augment variance-covariance matrix for 'mvrnorm' 
	## by adding a column/row for 'error in x'
	NCOL = ncol(VCOV)
	ADD1 = c(rep(0, NCOL))
	ADD1 = matrix(ADD1, ncol = 1)
	colnames(ADD1) = predNAME
	VCOV = cbind(VCOV, ADD1)
	ADD2 = c(rep(0, NCOL + 1))
	ADD2 = matrix(ADD2, nrow = 1)
	rownames(ADD2) = predNAME
	VCOV = rbind(VCOV, ADD2) 

	## iterate over all entries in 'newdata' as in usual 'predict.' functions
	NR = nrow(newdata)
	respVEC = numeric(NR)
	seVEC = numeric(NR)
	varPLACE = ncol(VCOV)   

	## define counter function
	counter = function(i) {
		if (i%%10 == 0) { cat(i) 
	} else { cat(".") }
		if (i%%50 == 0) { cat("\n") }
		flush.console()
	}
	  
	## create output matrix (df)
	outMAT = NULL 
	
	for (i in 1:NR) {
		
		## get predictor values and optional errors
		predVAL = newdata[i, 1]
		if (ncol(newdata) == 2) predERROR = newdata[i, 2] else predERROR = 0
		names(predVAL) = predNAME  
		names(predERROR) = predNAME  
		
		## create mean vector for 'mvrnorm'
		MU = c(COEF, predVAL)
		
		## create variance-covariance matrix for 'mvrnorm'
		## by putting error^2 in lower-right position of VCOV
		newVCOV = VCOV
		newVCOV[varPLACE, varPLACE] = predERROR^2
		
		## create MC simulation matrix
		simMAT = mvrnorm(n = nsim, mu = MU, Sigma = newVCOV, empirical = TRUE)
		
		## evaluate expression on rows of simMAT
		EVAL = try(eval(EXPR, envir = as.data.frame(simMAT)), silent = TRUE)
		if (inherits(EVAL, "try-error")) stop("There was an error evaluating the simulations!")
		
		## collect statistics
		PRED = data.frame(predVAL)
		colnames(PRED) = predNAME   
		FITTED = predict(object, newdata = data.frame(PRED))
		MEAN.sim = mean(EVAL, na.rm = TRUE)
		SD.sim = sd(EVAL, na.rm = TRUE)
		MEDIAN.sim = median(EVAL, na.rm = TRUE)
		MAD.sim = mad(EVAL, na.rm = TRUE)
		QUANT = quantile(EVAL, c((1 - level)/2, level + (1 - level)/2))
		RES = c(FITTED, MEAN.sim, SD.sim, MEDIAN.sim, MAD.sim, QUANT[1], QUANT[2])
		outMAT = rbind(outMAT, RES)
	}
	  
	colnames(outMAT) = c("fit", "mean", "sd", "median", "mad", names(QUANT[1]), names(QUANT[2]))
	rownames(outMAT) = NULL   
	return(outMAT)      
}
	
# Asymptotic execution
asymptoticMKExp = function(daf, divergence, xlow, xhigh, seed) {
	
	## Check data
	# check = checkInput(daf, divergence, xlow, xhigh)
	
	# if(check$data == FALSE) {
	# 	stop(check$print_errors) }
	
	if (any(daf$P0 == 0)){ ## Warning P0
		warning('Input daf file contains P0 values = 0.\nThis can bias the function fitting and the estimation of alpha.')}
	
	## Check seed
	if(missing(seed)) {
		seed = NULL
	} else {
		set.seed(seed)
	}
	
	## Parse the data from argument x
	f  = daf$daf #derived alelle frequencies
	p  = daf$Pi #non-synonymous polymorphism 
	p0 = daf$P0 #synonymous polymorphism
	
	## Parse the data from argument y
	m  = divergence$mi #number of non-synonymous analyzed positions   
	m0 = divergence$m0 ##number of synonymous analyzed positions
	d  = divergence$Di #non-synonymous divergence
	d0 = divergence$D0 #synonymous divergence
	
	## Compute alpha values and trim
	alpha         = 1 - (d0/d) * (p/p0)
	cutoff_f1     = xlow
	cutoff_f2     = xhigh
	trim          = ((f >= cutoff_f1) & (f <= cutoff_f2))
	f_trimmed     = f[trim]
	alpha_trimmed = alpha[trim]
	
	## Compute the original MK alpha
	alpha_nonasymp = 1 - (d0/d) * (sum(p[trim])/sum(p0[trim])) #using trimmed values
	
	## Two-step nls2() model fit at a given level of precision (res)
	mod1 = fitMKmodel(alpha_trimmed, f_trimmed, 10)
	
	## If mod1 did not work, try a deeper scan for a decent fit (res=20)
	if (length(mod1) == 0) {
		mod1 = fitMKmodel(alpha_trimmed, f_trimmed, 20)
	} 
	
	tryCatch({
		mod2 = lm(alpha_trimmed ~ f_trimmed)
	}, error=function(cond) {})
	

	tryCatch({
		ci_pred = predictNLS(mod1, newdata=data.frame(f_trimmed=1.0))
		alpha_1_low = ci_pred[6]
		alpha_1_high = ci_pred[7]

		## Preparation of ouput (alpha asym, a, b, c)
		alpha_1_est = predict(mod1, newdata=data.frame(f_trimmed=1.0))
		const_a = coef(mod1)['const_a']
		const_b = coef(mod1)['const_b']
		const_c = coef(mod1)['const_c']
		
		## Output table
		result_df = data.frame(model='exponential', a=const_a, b=const_b, c=const_c, alphaAsymptotic=alpha_1_est, ciLow=alpha_1_low, ciHigh=alpha_1_high, alphaOriginal=alpha_nonasymp, row.names=NULL)
		return(result_df)
	}, error=function(cond) {cat('Could not fit exponential model for the computation of asymptotic alpha.\n')})
}

mutFractionsPlot = function(fractionDf){

	fractionsMelt = melt(fractionDf, id.vars=NULL) 
	fractionsMelt[['test']] = rep(c('aMKT'),3)
	## Fractions graph
	pf = ggplot(fractionsMelt) + geom_bar(stat='identity', aes_string(x='test', y='value', fill='variable'), color='black') +
			coord_flip() + 
			themePublication() + 
			ylab(label='Fraction') + 
			xlab(label='Cut-off') +
			scale_fill_manual(values=c('#386cb0','#fdb462','#7fc97f','#ef3b2c','#662506','#a6cee3','#fb9a99','#984ea3','#ffff33'), breaks=c('f','d','b'), labels=c(expression(italic('f')),expression(italic('d')),expression(italic('b')))) +
			theme(axis.line=element_blank()) + scale_y_discrete(limit=seq(0,1,0.25), expand=c(0,0))

	# pf = ggplot(fractionDf) + 
	# 	geom_bar(stat='identity', aes_string(x='MKT', y='Fraction', fill='Type'), color='black') +
	# 	coord_flip() + 
	# 	ylab(label='Fraction') + 
	# 	xlab(label='Cut-off') +
	# 	scale_fill_manual(values=c('#386cb0','#fdb462','#7fc97f','#ef3b2c','#662506','#a6cee3','#fb9a99','#984ea3','#ffff33'), breaks=c('f','d','b'), labels=c(expression(italic('f')),expression(italic('d')),expression(italic('b')))) +
	# 	themePublication() + 
	# 	theme(axis.title.y=element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.line.x =element_blank()) + 
	# 	scale_y_discrete(limit=seq(0,1,0.25), expand=c(0,0))

	return(pf)
}

plotDaf = function(asympProccesedDaf){
	# Processing daf to input ggplot
	dafGraph = asympProccesedDaf[c('daf','Pi','P0')]
	dafGraph = melt(dafGraph,id.vars = 'daf')

	# Plotting values
	pD = ggplot(dafGraph) +
		geom_point(aes_string(x='daf', y='value', color='variable'), size=3) +
		scale_color_manual(values=c('#386cb0','#fdb462'), name='Type', breaks=c('Pi','P0'), labels=c('Non-synonymous','Synonymous')) +
		xlab('Derived Allele Frequency') + 
		ylab('Number of Sites') +
		themePublication() 
}

plotAsymptotic = function(asympProccesedDaf,asymptoticDf,alpha='alpha',xl,xh){

	# Defining asymptotic curve depending on daf values
	y1 = function(x) {
		asymptoticDf[['a']] + asymptoticDf[['b']]*exp(-asymptoticDf[['c']]*x)
	}
	
	# Creating values to shade weakly deleterious area
	xs       = seq(xl, xh, length.out = nrow(asympProccesedDaf)*5)
	ysmax    = rep(asymptoticDf$alphaAsymptotic, length(xs))
	ysmin    = y1(xs)
	shaderDf = data.frame(xs, ysmin, ysmax)
	
	
	pa       = ggplot(asympProccesedDaf, aes_string(x='daf', y=alpha))  +
		## Confidence intervals
		geom_rect(
			data = data.frame(
				xmin = -Inf, 
				xmax = Inf, 
				ymin = asymptoticDf$ciLow, 
				ymax = asymptoticDf$ciHigh
			),
			aes_string(xmin = 'xmin', xmax = 'xmax', ymin = 'ymin', ymax = 'ymax'), fill = 'gray30', alpha = 0.4, inherit.aes = F) +
		## Function curve
		stat_function(fun = y1, color = '#ef3b2c', size = 2) +
		## Asymptotic alpha
		geom_hline(yintercept = asymptoticDf$alphaAsymptotic, color = '#662506', linetype = 'dashed') +  
		## Alpha derived via classic MKT
		geom_hline(yintercept = asymptoticDf$alphaOriginal, color = '#386cb0', linetype = 'dashed') +
		geom_hline(yintercept = c(asympProccesedDaf$alpha[[1]],asymptoticDf$alphaAsymptotic), color = '#386cb0', linetype = 'dashed') +
		## Cut-offs
		geom_vline(xintercept = c(xl, xh), color = 'gray10', linetype='dotted') +    
		## Points
		geom_point(size=3, color='gray15') +
		## Shade the fraction of wd
		geom_ribbon(data=shaderDf, aes_string(x='xs', ymin='ysmin', ymax='ysmax'), fill='gray30', alpha=0.2, inherit.aes=F) + 
		## Customization
		xlab('Derived allele frequency') + ylab(expression(bold(paste('Adaptation (',alpha,')')))) +
		## Alphas labels
		annotate('text', x=xh-0.2, y=asymptoticDf$alphaAsymptotic-0.2, label=paste('alpha [asymptotic] ==', round(asymptoticDf$alphaAsymptotic, digits = 3)), color='#662506', size=5) +
		annotate('text', x=xh-0.2, y=0, label='a+b·e^(-cx)',color='#ef3b2c', size=5) +
		annotate('text', x = xh - 0.2, y = asymptoticDf$alphaOriginal - 0.1, label = paste0('alpha [standard] == ', round(asymptoticDf$alphaOriginal,digits = 3)), color = '#386cb0', size = 4) + 
		themePublication()

	return(pa)
}
