#' @title iMKT using PopHuman data
#'
#' @description Perform any MKT method using a subset of PopHuman data defined by custom genes and populations lists
#'
#' @details Execute any MKT method (standardMKT, FWW, imputedMKT, eMKT, aMKT) using a subset of PopHuman data defined by custom genes and populations lists. It uses the dataframe PopHumanData, which can be already loaded in the workspace (using loadPopHuman()) or is directly loaded when executing this function. It also allows deciding whether to analyze genes groupped by recombination bins or not, using recombination rate values corresponding to the sex average estimates from Bhérer et al. 2017 Nature Commun. 
#'
#' @param genes list of genes to analyze
#' @param pops list of populations to analyze
#' @param cutoff list of cutoffs to perform FWW, eMKT and/or imputedMKT
#' @param recomb group genes according to recombination values (TRUE/FALSE)
#' @param bins number of recombination bins to compute (mandatory if recomb=TRUE)
#' @param test which test to perform. Options include: standardMKT (default), imputedMKT, eMKT, FWW, aMKT
#' @param xlow lower limit for asymptotic alpha fit (default=0)
#' @param xhigh higher limit for asymptotic alpha fit (default=1)
#' @param plot report plot (optional). Default is FALSE
#' 
#' @return List of lists with the default test output for each selected population (and recombination bin when defined)
#'
#' @examples
#' ## List of genes (gene symbols, not Ensembl IDs)
#' mygenes <- c("SCYL3","C1orf112","FGR","CFH","STPG1","NIPAL3",
#'              "AK2","KDM1A","TTC22","ST7L","SELE","DNAJC11")
#' ## Perform analyses
#' PopHumanAnalysis(genes=mygenes, pops=c("CEU","YRI"), recomb=FALSE, test="standardMKT")
#' PopHumanAnalysis(genes=mygenes, pops="CEU", recomb=TRUE, bins=3, test="imputedMKT")
#' 
#' @import utils
#' @import stats
#'
#' @keywords PopData
#' @export
			
PopHumanAnalysis <- function(genes=c("gene1","gene2","..."), pops=c("pop1","pop2","..."), cutoff=0.05, recomb=TRUE/FALSE, bins=0, test=c("standardMKT","imputedMKT","eMKT","FWW","aMKT"), xlow=0, xhigh=1, plot=FALSE) { 
	
	## Get PopHuman data
	if (exists("PopHumanData") == TRUE) {
	data <- get("PopHumanData")
	} else {
	loadPopHuman()
	data <- get("PopHumanData") }
	
	## Check input variables
	## Numer of arguments
	if (nargs() < 3 && nargs()) {
	stop("You must specify 3 arguments at least: genes, pops, recomb (T/F).\nIf test = asymptoticMKT or test = aMKT, you must specify xlow and xhigh values.") }
	
	## Argument genes
	if (length(genes) == 0 || all(genes == "") || !is.character(genes)) {
	stop("You must specify at least one gene.") }
	if (!all(genes %in% data$symbol) == TRUE) {
	difGenes <- setdiff(genes, data$symbol)
	difGenes <- paste(difGenes, collapse=", ")
	stopMssg <- paste0("MKT data is not available for the requested gene(s).\nRemember to use gene symbols.\nThe genes that caused the error are: ", difGenes, ".")
	stop(stopMssg) }
	
	## Argument pops
	if (length(pops) == 0 || all(pops == "") || !is.character(pops)) {
	stop("You must specify at least one population.") }
	if (!all(pops %in% data$pop) == TRUE) {
	correctPops <- c("ACB","ASW","BEB","CDX","CEU","CHB","CHS","CLM","ESN","FIN","GBR","GIH","GWD","IBS","ITU","JPT","KHV","LWK","MSL","MXL","PEL","PJL","PUR","STU","TSI","YRI")
	difPops <- setdiff(pops, correctPops)
	difPops <- paste(difPops, collapse=", ")
	stopMssg <- paste0("MKT data is not available for the sequested populations(s).\nSelect among the following populations:\nACB, ASW, BEB, CDX, CEU, CHB, CHS, CLM, ESN, FIN, GBR, GIH, GWD, IBS, ITU, JPT, KHV, LWK, MSL, MXL, PEL, PJL, PUR, STU, TSI, YRI!.\nThe populations that caused the error are: ", difPops, ".")
	stop(stopMssg) }
	
	## Argument recomb
	if (recomb != TRUE && recomb != FALSE) {
	stop("Parameter recomb must be TRUE or FALSE.") }
	
	## Argument bins
	if (recomb == TRUE) {
	if (!is.numeric(bins) || bins == 0	|| bins == 1) {
		stop("If recomb = TRUE, you must specify the number of bins to use (> 1).") }
	if (bins > round(length(genes)/2)) {
		stop("Parameter bins > (genes/2). At least 2 genes for each bin are required.") }
	}
	if (recomb == FALSE && bins != 0) {
	warning("Parameter bins not used! (recomb=F selected)") }
	
	## Argument test and xlow + xhigh (when necessary)
	if(missing(test)) {
	test <- "standardMKT"
	}
	else if (test != "standardMKT" && test != "imputedMKT" && test != "eMKT" && test != "FWW" && test != "aMKT") {
	stop("Parameter test must be one of the following: standardMKT, imputedMKT, eMKT, FWW, aMKT")
	}
	if (length(test) > 1) {
	stop("Select only one of the following tests to perform: standardMKT, imputedMKT, eMKT, FWW, aMKT") }
	if ((test == "standardMKT" || test == "imputedMKT" || test == "eMKT" || test == "FWW") && (xlow != 0 || xhigh != 1)) {
	warningMssgTest <- paste0("Parameters xlow and xhigh not used! (test = ",test," selected)")
	warning(warningMssgTest) }
	
	## Arguments xlow, xhigh features (numeric, bounds...) checked in checkInput()
	
	## Perform subset
	subsetGenes <- data[(data$symbol %in% genes & data$pop %in% pops), ]
	subsetGenes$symbol <- as.factor(subsetGenes$symbol)
	subsetGenes <- droplevels(subsetGenes)
	
	## If recomb analysis is selected
	if (recomb == TRUE) {
	
	## Declare output list (each element 1 pop)
	outputList <- list()
	
	for (k in levels(subsetGenes$pop)) {
		print(paste0("Population = ", k))
		
		## Declare bins output list (each element 1 bin)
		outputListBins <- list()
		
		x <- subsetGenes[subsetGenes$pop == k, ]
		x <- x[order(x$recomb), ]
		
		## create bins
		## NOTE: replaced the previous manual binsize/modulo loop, which
		## silently dropped an entire bin's worth of genes whenever
		## nrow(x) was not exactly divisible by 'bins' (same confirmed
		## bug as in PopFlyAnalysis.R). cut() assigns every gene to
		## exactly one of 'bins' groups, with no gene left out.
		x$Group <- cut(seq_len(nrow(x)), breaks = bins, labels = FALSE)
		dat <- x
		dat$Group <- as.factor(dat$Group)
		
		## Iterate through each recomb bin
		for (j in levels(dat$Group)) {
		print(paste0("Recombination bin = ", j))
		x1 <- dat[dat$Group == j, ]
		
		## Recomb stats from bin j
		numGenes <- nrow(x1)
		minRecomb <- min(x1$recomb, na.rm=T)
		medianRecomb <- median(x1$recomb, na.rm=T)
		meanRecomb <- mean(x1$recomb, na.rm=T)
		maxRecomb <- max(x1$recomb, na.rm=T)
		recStats <- cbind(numGenes,minRecomb,medianRecomb,meanRecomb,maxRecomb)
		recStats <- as.data.frame(recStats)
		recStats <- list("Recombination bin Summary"=recStats)
		
		## Set counters to 0
		Pi <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
		P0 <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
		f <- seq(0.025,0.975,0.05)
		mi <- 0; m0 <- 0
		Di <- 0; D0 <- 0
		
		## Group genes
		x1 <- droplevels(x1)
		for (l in levels(x1$symbol)) {
			x2 <- x1[x1$symbol == l, ]
			
			## DAF
			## DAF
			x2$DAF0f <- as.character(x2$DAF0f); x2$DAF4f <- as.character(x2$DAF4f)
			daf0f <- Reduce(`+`, lapply(x2$DAF0f, function(z) as.numeric(unlist(strsplit(z, split=";")))))
			daf4f <- Reduce(`+`, lapply(x2$DAF4f, function(z) as.numeric(unlist(strsplit(z, split=";")))))
			Pi <- Pi + daf0f; P0 <- P0 + daf4f
			
			## Divergence
			mi <- mi + sum(x2$mi); m0 <- m0 + sum(x2$m0)
			Di <- Di + sum(x2$di); D0 <- D0 + sum(x2$d0)
		}
		
		## Proper formats
		daf <- cbind(f, Pi, P0); daf <- as.data.frame(daf)
		names(daf) <- c("daf","Pi","P0")
		div <- cbind(mi, Di, m0, D0); div <- as.data.frame(div)
		names(div) <- c("mi","Di","m0","D0")
		
		## Check data inside each test!
		
		## Transform daf20 to daf10 (faster fitting) for asymptoticMKT and aMKT
		if (nrow(daf) == 20) {
			daf1 <- daf
			daf1$daf10 <- sort(rep(seq(0.05,0.95,0.1),2)) ## Add column with the daf10 frequencies
			daf1 <- daf1[c("daf10","Pi","P0")] ## Keep new frequencies, Pi and P0
			daf1 <- aggregate(. ~ daf10, data=daf1, FUN=sum)	## Sum Pi and P0 two by two based on daf
			colnames(daf1)<-c("daf","Pi","P0") ## Final daf columns name
		}
		
		## Perform test
		if(test == "standardMKT") {
			output <- standardMKT(daf, div) 
			output <- c(output, recStats) } ## Add recomb summary for bin j
		else if(test == "imputedMKT" && plot == FALSE) {
			output <- imputedMKT(daf, div,listCutoffs=cutoff) 
			output <- c(output, recStats) }
		else if(test == "imputedMKT" && plot == TRUE) {
			output <- imputedMKT(daf, div, listCutoffs=cutoff ,plot=TRUE) 
			output <- c(output, recStats) }
		else if(test == "eMKT" && plot == FALSE) {
			output <- eMKT(daf, div, listCutoffs=cutoff)
			output <- c(output, recStats) }
		else if(test == "eMKT" && plot == TRUE) {
			output <- eMKT(daf, div, listCutoffs=cutoff, plot=TRUE)
			output <- c(output, recStats) }
		else if(test == "FWW" && plot == FALSE) {
			output <- FWW(daf, div, listCutoffs=cutoff)					 
			output <- c(output, recStats) }
		else if(test == "FWW" && plot == TRUE) {
			output <- FWW(daf, div, listCutoffs=cutoff, plot=TRUE)					 
				output <- c(output, recStats) }
		else if(test == "aMKT" && plot == FALSE) {
			output <- aMKT(daf1, div, xlow, xhigh)
			output <- c(output, recStats) }
		else if(test == "aMKT" && plot == TRUE) {
			output <- aMKT(daf1, div, xlow, xhigh, plot=TRUE)
			output <- c(output, recStats) }
		
		## Fill list with each bin
		outputListBins[[paste("Recombination bin = ",j)]] <- output
		}
		
		## Fill list with each pop
		outputList[[paste("Population = ",k)]] <- outputListBins
	}
	
	## Warning if some genes are lost. Bins must be equally sized.
	if (nrow(dat) != length(genes)) {
		missingGenes <- round(length(genes) - nrow(dat))
		genesNames <- as.vector(tail(x, missingGenes)$Name)
		genesNames <- paste(genesNames, collapse=", ")
		warningMssg <- paste0("The ",missingGenes," gene(s) with highest recombination rate estimates (", genesNames, ") was/were excluded from the analysis in order to get equally sized bins.\n")
		warning(warningMssg) }
	
	## Return output
	cat("\n")
	return(outputList)
	}
	
	## If NO recombination analysis selected
	else if (recomb == FALSE) {
	
	## Declare output list (each element 1 pop)
	outputList <- list()
	
	for (i in levels(subsetGenes$pop)) {
		print(paste0("Population = ", i))
		x <- subsetGenes[subsetGenes$pop == i, ]
		
		## Set counters to 0
		Pi <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
		P0 <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
		f <- seq(0.025,0.975,0.05)
		mi <- 0; m0 <- 0
		Di <- 0; D0 <- 0
		
		## Group genes
		for (j in levels(x$symbol)) {
		x1 <- x[x$symbol == j, ]
		
		## DAF
		x1$DAF0f <- as.character(x1$DAF0f); x1$DAF4f <- as.character(x1$DAF4f)
		daf0f <- Reduce(`+`, lapply(x1$DAF0f, function(z) as.numeric(unlist(strsplit(z, split=";")))))
		daf4f <- Reduce(`+`, lapply(x1$DAF4f, function(z) as.numeric(unlist(strsplit(z, split=";")))))
		Pi <- Pi + daf0f; P0 <- P0 + daf4f
		
		## Divergence
		mi <- mi + sum(x1$mi); m0 <- m0 + sum(x1$m0)
		Di <- Di + sum(x1$di); D0 <- D0 + sum(x1$d0)
		}
		
		## Proper formats
		daf <- cbind(f, Pi, P0); daf <- as.data.frame(daf)
		names(daf) <- c("daf","Pi","P0")
		div <- cbind(mi, Di, m0, D0); div <- as.data.frame(div)
		names(div) <- c("mi","Di","m0","D0")
		
		## Check data inside each test!
		
		## Transform daf20 to daf10 (faster fitting) for asymptoticMKT and aMKT
		if (nrow(daf) == 20) {
		daf1 <- daf
		daf1$daf10 <- sort(rep(seq(0.05,0.95,0.1),2)) ## Add column with the daf10 frequencies
		daf1 <- daf1[c("daf10","Pi","P0")] ## Keep new frequencies, Pi and P0
		daf1 <- aggregate(. ~ daf10, data=daf1, FUN=sum)	## Sum Pi and P0 two by two based on daf
		colnames(daf1)<-c("daf","Pi","P0") ## Final daf columns name
		}
		
		## Perform test
		if(test == "standardMKT") {
		output <- standardMKT(daf, div) }
		else if(test == "imputedMKT" && plot == FALSE) {
		output <- imputedMKT(daf, div,listCutoffs=cutoff) }
		else if(test == "imputedMKT" && plot == TRUE) {
		output <- imputedMKT(daf, div,listCutoffs=cutoff, plot=TRUE) }
		else if(test == "eMKT" && plot == FALSE) {
		output <- eMKT(daf, div, listCutoffs=cutoff) }
		else if(test == "eMKT" && plot == TRUE) {
		output <- eMKT(daf, div, listCutoffs=cutoff, plot=TRUE) }
		else if(test == "FWW" && plot == FALSE) {
		output <- FWW(daf, div, listCutoffs=cutoff) }
		else if(test == "FWW" && plot == TRUE) {
		output <- FWW(daf, div, listCutoffs=cutoff, plot=TRUE) }
		else if(test == "aMKT" && plot == FALSE) {
		output <- aMKT(daf1, div, xlow, xhigh) }
		else if(test == "aMKT" && plot == TRUE) {
		output <- aMKT(daf1, div, xlow, xhigh, plot=TRUE) }
		
		## Fill list with each pop
		outputList[[paste("Population = ",i)]] <- output
	}
	
	## Return output
	cat("\n")
	return(outputList)
	}
}

