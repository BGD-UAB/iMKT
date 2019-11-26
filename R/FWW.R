#' @title FWW correction method
#'
#' @description MKT calculation corrected using FWW method (Fay et al. 2001 Genetics).
#'
#' @details In the standard McDonald and Kreitman test, the estimate of adaptive evolution (alpha) can be easily biased by the segregation of slightly deleterious non-synonymous substitutions. Specifically, slightly deleterious mutations contribute more to polymorphism than they do to divergence, and thus, lead to an underestimation of alpha. Because they tend to segregate at lower frequencies than do neutral mutations, they can be partially controled by removing low frequency polymorphisms from the analysis. This is known as the FWW method.
#'
#' @param daf data frame containing DAF, Pi and P0 values
#' @param divergence data frame containing divergent and analyzed sites for selected (i) and neutral (0) classes
#' @param listCutoffs list of cutoffs to use (optional). Default cutoffs are: 0, 0.05, 0.1
#' @param plot report plot (optional). Default is FALSE
#' 
#' @return MKT corrected by the FWW method. List with alpha results, graph (optional), divergence metrics, MKT tables and negative selection fractions
#'
#' @examples
#' ## Using default cutoffs
#' FWW(myDafData, myDivergenceData)
#' ## Using custom cutoffs and rendering plot
#' FWW(myDafData, myDivergenceData, c(0.05, 0.1, 0.15), plot=TRUE)
#' 
#' @import utils
#' @import stats
#' @import ggplot2
#' @importFrom ggthemes theme_foundation
#' @importFrom cowplot plot_grid
#'
#' @keywords MKT
#' @export

FWW = function(daf, divergence, listCutoffs, plot=FALSE) {
	
	## Check data
	check = checkInput(daf, divergence, 0, 1)
	if(check$data == FALSE) {
		 stop(check$print_errors) }

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

	## Divergence metrics
	Ka              = Di/mi
	Ks              = D0/m0
	omega           = Ka/Ks
	divTable        = data.frame('Ka' = Ka, 'Ks' = Ks, 'omega' = omega)
	
	## Iterate along cutoffs
	# Opening necessary variables to iter
	P0         = sum(daf$P0) 
	Pi         = sum(daf$Pi) 
	cleanedDaf = daf
	PiNeutral  = 0

	mktTableStandard = data.frame(Polymorphism = c(sum(P0), sum(Pi)), Divergence=c(D0,Di),row.names = c("Neutral class","Selected class"))

	## Iterate along cutoffs
	alphaCorrected <- list()
	for (c in listCutoffs) {
		
		dafRemove = daf[daf[['daf']] > c,] 
	
		## Create MKT table 
		mktTableCleaned = data.frame(Polymorphism=c(sum(dafRemove[['P0']]), sum(dafRemove[['Pi']])), Divergence=c(D0,Di), row.names=c("Neutral class","Selected class"))
	
		## Estimate of alpha
		alphaC = 1-(mktTableCleaned[2,1]/mktTableCleaned[1,1])*(mktTableCleaned[1,2]/mktTableCleaned[2,2])
		
		## Fisher exact test p-value from the MKT
		pvalue = fisher.test(mktTableCleaned)$p.value
		
		## Omega A and Omega D
		omegaA = omega * alphaC
		omegaD = omega - omegaA
		
		## Store output  
		alphaCorrected[[paste0('cutoff=',c)]] = c(c, alphaC, pvalue)
		divCutoff[[paste0('cutoff=',c)]] = c(c, Ka,Ks,omegaA, omegaD)
		# mktTables[[paste0('cutoff=',c)]]  = mktTableCleaned
	}
	
	## Output format
	output[['mktTable']]                 = mktTableStandard 
	output[['alphaCorrected']]           = as.data.frame(do.call('rbind',alphaCorrected))
	colnames(output[['alphaCorrected']]) = c('cutoff', 'alphaCorrected', 'pvalue')
	## Divergence metricss
	divCutoff                            = as.data.frame(do.call('rbind',divCutoff))
	names(divCutoff)                     = c('cutoff', 'Ka','Ks','omegaA', 'omegaD')
	output[['divMetrics']]               = list('metricsByCutoff'=divCutoff)

	## Perform plot
	if(plot == TRUE) {
		## Cut-offs graphs
		plot = ggplot(output[['alphaCorrected']], aes(x=as.factor(cutoff), y=alphaCorrected, group=1)) +
			geom_line(color="#386cb0") + 
			geom_point(size=2.5, color="#386cb0")+
			themePublication() +
			xlab("Cut-off") + ylab(expression(bold(paste("Adaptation (",alpha,")"))))
	
		## Re-format outputs
		output[['graph']] = plot
		
		return(output)

	## If no plot to perform  
	} else if (plot == FALSE) {
		return(output)
	}
}
