## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
	echo = TRUE,
	fig.align = "center",
	fig.height = 7,
	fig.width = 7,
	collapse = TRUE,
	comment = "#>"
)

## ----load package and see sample data, echo=T, message=F, warning=F-----------
## Install devtools if needed
# install.packages("devtools")
# devtools::install_github("BGD-UAB/iMKT")
library(iMKT)

## Sample daf data
head(myDafData)

## Sample divergence data
myDivergenceData

## ----Standard MKT, echo=TRUE--------------------------------------------------
standardMKT(daf = myDafData, divergence = myDivergenceData)

## ----FWW, echo=TRUE-----------------------------------------------------------
FWW(daf = myDafData, divergence = myDivergenceData, listCutoffs = c(0, 0.05, 0.1))

## ----FWW plot, echo=TRUE, fig.width=6, fig.height=4---------------------------
FWW(daf = myDafData, divergence = myDivergenceData, listCutoffs=c(0.05, 0.15, 0.25, 0.35), plot = TRUE)

## ----imputedMKT, echo=TRUE----------------------------------------------------
imputedMKT(daf = myDafData, divergence = myDivergenceData, listCutoffs = c(0.05))

## ----echo=TRUE, fig.width=6, fig.height=6-------------------------------------
imputedMKT(daf = myDafData, divergence = myDivergenceData, listCutoffs = c(0.05, 0.15,0.25,0.35), plot = TRUE)

## ----eMKT, echo=TRUE----------------------------------------------------------
eMKT(daf = myDafData, divergence = myDivergenceData, listCutoffs = c(0.05))

## ----eMKT plot, echo=TRUE, fig.width=6, fig.height=6--------------------------
eMKT(daf = myDafData, divergence = myDivergenceData, listCutoffs = c(0.05, 0.15,0.25,0.35), plot=TRUE)

## ----Asymptotic MKT, echo=TRUE------------------------------------------------
asymptoticMKT(daf = myDafData, divergence = myDivergenceData, xlow = 0, xhigh = 0.9)

## ----aMKT, echo=TRUE, fig.width=6, fig.height=9-------------------------------
aMKT(daf=myDafData, divergence=myDivergenceData, xlow=0, xhigh=0.9, plot=TRUE)

## ----summary, echo=FALSE------------------------------------------------------
results <- data.frame("Standard"=0.2364, "FWW_0.05"=0.5410, "imputedMKT_0.05"=0.5410, "eMKT_0.05"=0.4249, "aMKT"=0.6572)
knitr::kable(results, align="c")
rm(results)

