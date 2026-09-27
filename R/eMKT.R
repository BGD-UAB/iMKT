#' @title eMKT correction method
#'
#' @description MKT calculation corrected using the eMKT method (Mackay et al. 2012 Nature).
#'
#' @details In the standard McDonald and Kreitman test, the estimate of adaptive evolution (alpha) can be easily biased by the segregation of slightly deleterious non-synonymous substitutions. Specifically, slightly deleterious mutations contribute more to polymorphism than they do to Divergence, and thus, lead to an underestimation of alpha. Because adaptive mutations and weakly deleterious selection act in opposite directions on the MKT, alpha and the fraction of substitutions that are slightly deleterious, b, will be both underestimated when both selection regimes occur. To take adaptive and slightly deleterious mutations mutually into account, Pi, the count of segregating sites in class i, should be separated into the number of neutral variants and the number of weakly deleterious variants, Pi = Pineutral + Pi weak del. Alpha is then estimated as 1-(Pineutral/P0)(D0/Di). As weakly deleterious mutations tend to segregate at low frequencies, the neutral and weakly deleterious fractions from Pi can be estimated based on any frequency cutoff established.
#'
#' @param daf data frame containing DAF, Pi and P0 values
#' @param divergence data frame containing divergent and analyzed sites for selected (i) and neutral (0) classes
#' @param listCutoffs list of cutoffs to use (optional). Default cutoffs are: 0, 0.05, 0.1
#' @param plot report plot (optional). Default is FALSE
#'
#' @return MKT corrected by the eMKT method. List with alpha results, graph (optional), Divergence metrics, MKT tables and negative selection fractions
#'
#' @examples
#' ## Using default cutoffs
#' eMKT(myDafData, myDivergenceData)
#' ## Using custom cutoffs and rendering plot
#' eMKT(myDafData, myDivergenceData, c(0.05, 0.1, 0.15), plot=TRUE)
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

eMKT <- function (daf, divergence, listCutoffs=c(0, 0.05, 0.1), plot = FALSE)
{
  output = list()
  mktTables = list()
  divMetrics = list()
  divCutoff = list()
  P0 = sum(daf[["P0"]])
  Pi = sum(daf[["Pi"]])
  D0 = divergence[["D0"]]
  Di = divergence[["Di"]]
  m0 = divergence[["m0"]]
  mi = divergence[["mi"]]
  mktTableStandard = data.frame(Polymorphism = c(sum(daf[["P0"]]),
                                                  sum(daf[["Pi"]])), Divergence = c(D0, Di), row.names = c("Neutral class",
                                                                                                            "Selected class"))
  Ka = Di/mi
  Ks = D0/m0
  omega = Ka/Ks
  alphaCorrected <- list()
  fractions <- list()
  for (c in listCutoffs) {
    PiMinus = sum(daf[daf[["daf"]] <= c, "Pi"])
    PiGreater = sum(daf[daf[["daf"]] > c, "Pi"])
    P0Minus = sum(daf[daf[["daf"]] <= c, "P0"])
    P0Greater = sum(daf[daf[["daf"]] > c, "P0"])
    ratioP0 = P0Minus/P0Greater
    deleterious = PiMinus - (PiGreater * ratioP0)
    PiNeutral = Pi - deleterious
    f_neutral <- P0Minus/sum(daf$P0)
    Pi_neutral_below_cutoff <- Pi * f_neutral
    Pi_wd <- PiMinus - Pi_neutral_below_cutoff
    Pi_neutral <- round(Pi_neutral_below_cutoff + PiGreater)
    alphaC <- 1 - ((Pi_neutral/P0) * (D0/Di))
    b = (Pi_wd/P0) * (m0/mi)
    f = (m0 * Pi_neutral)/(as.numeric(mi) * as.numeric(P0))
    d = 1 - (f + b)
    m = matrix(c(P0, (Pi - Pi_wd), D0, Di), ncol = 2)
    pvalue = fisher.test(round(m))$p.value
    omegaA = omega * alphaC
    omegaD = omega - omegaA
    alphaCorrected[[paste0("cutoff=", c)]] = c(c, alphaC,
                                                pvalue)
    fractions[[paste0("cutoff=", c)]] = c(c, d, f, b)
    divCutoff[[paste0("cutoff=", c)]] = c(c, Ka, Ks, omegaA,
                                           omegaD)
  }
  output[["mktTable"]] = mktTableStandard
  output[["alphaCorrected"]] = as.data.frame(do.call("rbind",
                                                       alphaCorrected))
  colnames(output[["alphaCorrected"]]) = c("cutoff", "alphaCorrected",
                                            "pvalue")
  divCutoff = as.data.frame(do.call("rbind", divCutoff))
  names(divCutoff) = c("cutoff", "Ka", "Ks", "omegaA", "omegaD")
  output[["divMetrics"]] = list(metricsByCutoff = divCutoff)
  output[["fractions"]] = as.data.frame(do.call("rbind", fractions))
  names(output[["fractions"]]) = c("cutoff", "d", "f", "b")
  if (plot == TRUE) {
    plotAlpha = ggplot(output[["alphaCorrected"]], aes(x = as.factor(cutoff),
                                                        y = alphaCorrected, group = 1)) + geom_line(color = "#386cb0") +
      geom_point(size = 2.5, color = "#386cb0") + themePublication() +
      xlab("Cut-off") + ylab(expression(bold(paste("Adaptation (",
                                                     alpha, ")"))))
    i = which.max(output$alphaCorrected$alphaCorrected)
    fractionsMelt = output[["fractions"]][i, 2:4]
    fractionsMelt = reshape2::melt(fractionsMelt, id.vars = NULL)
    fractionsMelt[["test"]] = rep(c("eMKT"), 3)
    plotFraction = ggplot(fractionsMelt) + geom_bar(stat = "identity",
                                                     aes(x = test, y = value, fill = variable),
                                                     color = "black") + coord_flip() + themePublication() +
      ylab(label = "Fraction") + xlab(label = "Cut-off") +
      scale_fill_manual(values = c("#386cb0", "#fdb462",
                                    "#7fc97f", "#ef3b2c", "#662506", "#a6cee3", "#fb9a99",
                                    "#984ea3", "#ffff33"), breaks = c("f", "d", "b"),
                         labels = c(expression(italic("f")), expression(italic("d")),
                                    expression(italic("b")))) + theme(axis.line = element_blank()) +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = c(0,
                                                            0))
    plotEmkt = plot_grid(plotAlpha, plotFraction, nrow = 2,
                          labels = c("A", "B"), rel_heights = c(2, 1))
    output[["graph"]] = plotEmkt
    return(output)
  }
  else if (plot == FALSE) {
    return(output)
  }
}
