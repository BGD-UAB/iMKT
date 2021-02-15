#' @title imputedMKT correction method
#'
#' @description MKT calculation corrected using imputedMKT method (Mackay et al. 2012 Nature).
#'
#' @details In the standard McDonald and Kreitman test, the estimate of adaptive evolution (alpha) can be easily biased by the segregation of slightly deleterious non-synonymous substitutions. Specifically, slightly deleterious mutations contribute more to polymorphism than they do to Divergence, and thus, lead to an underestimation of alpha. Because adaptive mutations and weakly deleterious selection act in opposite Directions on the MKT, alpha and the fraction of substitutions that are slighlty deleterious, b, will be both underestimated when both selection regimes occur. To take adaptive and slighlty deleterious mutations mutually into account, Pi, the count off segregatning sites in class i, should be separated into the number of neutral variants and the number of weakly deleterious variants, Pi = Pineutral + Pi weak del. Alpha is then estimated as 1-(Pineutral/P0)(D0/Di). As weakly deleterious mutations tend to segregate at low frequencies, neutral and weakly deleterious fractions from Pi can be estimated based on any frequency cutoff established.
#'
#' @param daf data frame containing DAF, Pi and P0 values
#' @param Divergence data frame containing Divergent and analyzed sites for selected (i) and neutral (0) classes
#' @param listCutoffs list of cutoffs to use (optional). Default cutoffs are: 0, 0.05, 0.1
#' @param plot report plot (optional). Default is FALSE
#' 
#' @return MKT corrected by the imputedMKT method. List with alpha results, graph (optional), Divergence metrics, MKT tables and negative selection fractions
#'
#' @examples
#' ## Using default cutoffs
#' imputedMKT(myDafData, myDivergenceData)
#' ## Using custom cutoffs and rendering plot
#' imputedMKT(myDafData, myDivergenceData, c(0.05, 0.1, 0.15), plot=TRUE)
#'
#' @import utils
#' @import stats
#' @import ggplot2
#' @importFrom ggthemes theme_foundation
#' @importFrom cowplot plot_grid
#' @importFrom reshape2 melt 
#'
#' @keywords MKT
#' @export


####################################################
################# MKT-FWW function #################
####################################################

imputedMKT = function(daf, divergence, listCutoffs, plot=FALSE) {
	
	## Check data
	check = checkInput(daf, divergence, 0, 1)
	if(check$data == FALSE)
	{
		 stop(check$print_errors) 
	}
	
	## Declare output lists and data frames
	output     = list()
	mktTables  = list()
	divMetrics = list()
	divCutoff  = list()
	
	##Polymorphism and Divergence variables
	P0 = sum(daf[['P0']])
	Pi = sum(daf[['Pi']])
	D0 = divergence[['D0']]
	Di = divergence[['Di']]
	m0 = divergence[['m0']]
	mi = divergence[['mi']]
	
	## MKT tables
	mktTableStandard = data.frame(Polymorphism = c(sum(daf[['P0']]), sum(daf[['Pi']])), Divergence = c(D0,Di),row.names = c("Neutral class","Selected class"))
	
	## Divergence metrics
	Ka              = Di/mi
	Ks              = D0/m0
	omega           = Ka/Ks

	## Estimation of alpha
	# alpha = 1 - ((Pi/P0)*(D0/Di))
	# alphaC = 1 - ((PiNeutral/P0)*(mktTableStandard[1,2]/mktTableStandard[2,2]))
	alphaCorrected <- list()
	fractions <- list()
	for (c in listCutoffs) {
		
		## Estimating alpha with Pi/P0 ratio 
		PiMinus     = sum(daf[daf[['daf']] <= c,'Pi'])
		PiGreater   = sum(daf[daf[['daf']] > c,'Pi'])
		P0Minus     = sum(daf[daf[['daf']] <= c,'P0'])
		P0Greater   = sum(daf[daf[['daf']] > c,'P0'])
		
		ratioP0     = P0Minus/P0Greater
		deleterious = PiMinus - (PiGreater * ratioP0)
		PiNeutral = Pi - deleterious

		alphaC = 1 - (((Pi - deleterious)/P0)*(D0/Di))


		## Estimation of b: weakly deleterious
		b    = (deleterious/P0)*(m0/mi)
		
		## Estimation of f: neutral sites
		f = (m0*PiNeutral)/(as.numeric(mi)*as.numeric(P0))
		
		## Estimation of d, strongly deleterious sites
		d = 1 - (f+b)
		
		## Fisher exact test p-value from the MKT
		m      = matrix(c(P0,(Pi - deleterious),D0,Di), ncol=2)
		pvalue = fisher.test(round(m))$p.value
		
		## Omega A and Omega D
		omegaA = omega * alphaC
		omegaD = omega - omegaA
		
		## Store output  
		alphaCorrected[[paste0('cutoff=',c)]] = c(c, alphaC, pvalue)
		fractions[[paste0('cutoff=',c)]] = c(c, d, f, b)
		divCutoff[[paste0('cutoff=',c)]] = c(c, Ka,Ks,omegaA, omegaD)
		# mktTables[[paste0('cutoff=',c)]]  = mktTableCleaned
	}

	## Fractions
	# names(fraction) = 'Correction'
	
	## Store output 
	## Output format
	output[['mktTable']]                 = mktTableStandard 
	output[['alphaCorrected']]           = as.data.frame(do.call('rbind',alphaCorrected))
	colnames(output[['alphaCorrected']]) = c('cutoff', 'alphaCorrected', 'pvalue')
	## Divergence metricss
	divCutoff                            = as.data.frame(do.call('rbind',divCutoff))
	names(divCutoff)                     = c('cutoff', 'Ka','Ks','omegaA', 'omegaD')
	output[['divMetrics']]               = list('metricsByCutoff'=divCutoff)
	output[['fractions']]                = as.data.frame(do.call('rbind',fractions))
	names(output[['fractions']])                     = c('cutoff', 'd','f','b')
	
	# DivCutoff                     = data.frame('omegaA' = omegaA, 'omegaD' = omegaD)
	# output[['divMetrics']]        = list(DivTable, DivCutoff)
	# names(output[['divMetrics']]) = c("Global metrics", "Estimates by cutoff")
	# Results table
	# output = as.data.frame(do.call("rbind",output))
	# colnames(output) = c("cutoff", "alpha", "pvalue")
	
	# Divergence metrics
	# DivCutoff = as.data.frame(do.call("rbind",DivCutoff))
	# names(DivCutoff) = c("omegaA", "omegaD")

	## Render plot
	if (plot == TRUE) {

		## Cut-offs graph
		# plot = ggplot(output, aes(x=as.factor(cutoff), y=alpha, group=1)) +
		# 	geom_line(color="#386cb0") + 
		# 	geom_point(size=2.5, color="#386cb0")+
		# 	themePublication() +
		# 	xlab("Cut-off") + ylab(expression(bold(paste("Adaptation (",alpha,")")))) 
	
		## Re-format outputs
		# output = output[,c(2,3)]
		# names(output) = c("alpha.symbol","Fishers exact test P-value")
		# DivCutoff = DivCutoff[,c(2,3)]
		# colnames(DivCutoff) = c("omegaA.symbol", "omegaD.symbol")
		# DivMetrics = list(DivTable, DivCutoff)
		# names(DivMetrics) = c("Global metrics", "Estimates by cutoff")
	
		## Melt fractions data
		plotAlpha = ggplot(output[['alphaCorrected']], aes(x=as.factor(cutoff), y=alphaCorrected, group=1)) +
			geom_line(color="#386cb0") + 
			geom_point(size=2.5, color="#386cb0")+
			themePublication() +
			xlab("Cut-off") + ylab(expression(bold(paste("Adaptation (",alpha,")"))))
	
		## Fractions graph
		i = which.max(output$alphaCorrected$alphaCorrected)
		fractionsMelt = output[['fractions']][i,2:4]
		fractionsMelt = reshape2::melt(fractionsMelt, id.vars=NULL) 
		fractionsMelt[['test']] = rep(c('imputedMKT'),3)

		plotFraction = ggplot(fractionsMelt) + geom_bar(stat="identity", aes_string(x="test", y="value", fill="variable"), color="black") +
			coord_flip() + themePublication() + ylab(label="Fraction") + xlab(label="Cut-off") +
			scale_fill_manual(values=c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33"), breaks=c("f","d","b"), labels=c(expression(italic("f")),expression(italic("d")),expression(italic("b")))) +
			theme(axis.line=element_blank()) + scale_y_discrete(limit=seq(0,1,0.25), expand=c(0,0))
	
		plotEmkt = plot_grid(plotAlpha, plotFraction, nrow=2,  labels=c('A','B'), rel_heights=c(2,1))

		output[['graph']] = plotEmkt

		# plot = plot_grid(plot, plotfraction, nrow=2, labels=c("A","B"), rel_heights=c(2,1))
	
		## Return list output  
		# output[['Graph']] = plot
		# names(listOutput) = c("Results","Graph", "Divergence metrics", "MKT tables","Fractions")

		return(output)

	## If no plot to render
	} else if (plot==FALSE) {
			## Re-format outputs
			# output = output[,c(2,3)]
			# names(output) = c("alpha.symbol","Fishers exact test P-value")
			# DivCutoff = DivCutoff[,c(2,3)]
			# colnames(DivCutoff) = c("omegaA.symbol", "omegaD.symbol")
			# DivMetrics = list(DivTable, DivCutoff)
			# names(DivMetrics) = c("Global metrics", "Estimates by cutoff")
		
			# ## Melt fractions data
			# fractionsMelt = melt(fractions, id.vars=NULL) 
			# fractionsMelt$Fraction =  rep(c("d", "f", "b"),length(fractionsMelt$variable)/3)
			
			# ## Return list output  
			# listOutput = list(output, DivMetrics, mktTables, fractions)
			# names(listOutput) = c("Results", "Divergence metrics", "MKT tables","Fractions")

			return(output)
	}
	
}



# for (f in 2:nrow(daf)-1) {
# 	cleanedDafBellow = cleanedDaf[1, ]
# 	cleanedDafAbove  = cleanedDaf[2:nrow(cleanedDaf), ] ## over cleanedDaf 
	
# 	## Create MKT table 
# 	mktTable  = data.frame(`cleanedDaf below cutoff` = c(sum(cleanedDafBellow$P0), sum(cleanedDafBellow$Pi)), `cleanedDaf above cutoff`=c(sum(cleanedDafAbove$P0), sum(cleanedDafAbove$Pi)),row.names = c("Neutral class","Selected class"))
	
# 	cleanedPi = sum(cleanedDaf$Pi)
	
# 	## Estimate fractions
# 	fNeutral  = mktTable[1,1]/sum(cleanedDaf$P0)
# 	PiNeutralBelowCutoff = cleanedPi * fNeutral

# 	# If we haven't deleterious in Pi then PiNeutral = Pi_in_current_cutoff
# 	if(PiNeutralBelowCutoff > cleanedDafBellow[['Pi']]){
# 		PiNeutral = PiNeutral + cleanedDafBellow[['Pi']]
# 	}
# 	else{
# 		PiNeutral = PiNeutral + PiNeutralBelowCutoff
# 	}
	
# 	## Deleting current frequency 
# 	cleanedDaf = cleanedDaf[2:nrow(cleanedDaf),]

# } 
