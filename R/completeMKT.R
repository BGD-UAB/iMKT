#' @title Complete MK methodologies
#'
#' @description MKT calculation using all methodologies included in the package: standardMKT, FWW, eMKT, imputedMKT, aMKT.
#' 
#' @details Perform all MKT derived methodologies at once using the same input data and parameters.
#'
#' @param daf data frame containing DAF, Pi and P0 values
#' @param divergence data frame containing divergent and analyzed sites for selected (i) and neutral (0) classes
#' @param listCutoffs list of cutoffs to use for FWW/eMKT/imputedMKT (optional). Default cutoffs are: 0, 0.05, 0.1
#' @param xlow lower limit for asymptotic alpha fit
#' @param xhigh higher limit for asymptotic alpha fit
#' @param seed seed value (optional). No seed by default
#'
#' @return List with the diverse MKT results: standardMKT, FWW, eMKT, imputedMKT, aMKT
#'
#' @examples 
#' completimputedMKT(myDafData, myDivergenceData, xlow=0, xhigh=0.9)
#' 
#' @import utils
#' @import stats
#'
#' @keywords MKT
#' @export

completimputedMKT = function(daf, divergence, listCutoffs=c(0, 0.05, 0.1), xlow=0, xhigh=1, seed) {
	
	## Check data
	check = checkInput(daf, divergence, xlow, xhigh)
	if(check$data == FALSE){
		 stop(check$print_errors) 
	}
	
	## Check seed
	if(missing(seed)) 
	{
		seed = NULL
	} else 
	{
		set.seed(seed)
	}

	## Create output list
	fullResults = list()
	
	## Standard MKT
	fullResults[['standardMKT']] = standardMKT(daf,divergence)
	
	## FWW MKT
	fullResults[['FWW']] = FWW(daf, divergence, listCutoffs=listCutoffs)
	
	## eMKT
	fullResults[['eMKT']] = eMKT(daf, divergence, listCutoffs=listCutoffs)
	
	## imputedMKT
	fullResults[['imputedMKT']] = imputedMKT(daf, divergence, listCutoffs=listCutoffs)
	
	## Asymptotic MKT (aMKT)
	fullResults[['aMKT']] = aMKT(daf, divergence, xlow, xhigh)
	
	## Output
	return(fullResults)
}

