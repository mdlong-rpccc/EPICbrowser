#' If input file name is found in the working directory, the file is loaded directly
#'
#' @details This approach uses the ChAMP package first described by Morris TJ et al:
#' Morris TJ, Butcher LM, Teschendorff AE, Chakravarthy AR, Wojdacz TK and Beck S (2014). “ChAMP: 450k Chip Analysis Methylation
#' Pipeline.” _Bioinformatics_, *30*(3), pp. 428-430. doi: 10.1093/bioinformatics/btt684
#' (URL:http://doi.org/10.1093/bioinformatics/btt684).
#'
#' @param limma.comparison data.frame(1) A data frame returned from the limma.MVP function
#' @param deltaB numeric(1) A numeric value equal to the fraction
#'
#' @return A table of
#'
#' @export


beta.thresh <- function(limma.comparison, deltaB){
    up <- length(limma.comparison$probeID[limma.comparison$adj.P.Val < 0.05 & limma.comparison$deltaBeta > deltaB])
    dn <- length(limma.comparison$probeID[limma.comparison$adj.P.Val < 0.05 & limma.comparison$deltaBeta < -deltaB])
    up.perc <- up / length(limma.comparison$probeID[limma.comparison$adj.P.Val < 0.05 & limma.comparison$deltaBeta > 0]) * 100
    dn.perc <- dn / length(limma.comparison$probeID[limma.comparison$adj.P.Val < 0.05 & limma.comparison$deltaBeta < 0]) * 100
    mvp.table <- data.frame(Count = c(up, dn), PercentOfTotal = c(up.perc, dn.perc))
    row.names(mvp.table) <- c("UP", "DN")
    print(mvp.table)
}
