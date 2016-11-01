#'Calling Methylation Variable Positions (MVPs)
#' Implements the limma package to calculate the p-value for differential methylation using a linear model
#' Can only do one comparison at a time - use the 'compute.group' argument to make comparisons from the sample sheet variables (Sample_Group)
#' Compare Groups: i.e. C1 (C42_shCTL-ETOH) and C4 (LNCaP_shCTL_ETOH)
#' If input file name is found in the working directory, the file is loaded directly
#'
#' @details This approach uses the ChAMP package first described by Morris TJ et al:
#' Morris TJ, Butcher LM, Teschendorff AE, Chakravarthy AR, Wojdacz TK and Beck S (2014). “ChAMP: 450k Chip Analysis Methylation
#' Pipeline.” _Bioinformatics_, *30*(3), pp. 428-430. doi: 10.1093/bioinformatics/btt684
#' (URL:http://doi.org/10.1093/bioinformatics/btt684).
#'
#' @param file_name character(1) The name of the output file returned from the champ.MVP function.
#' @param group1 character(1) The name of the first comparison group.
#' @param group2 character(1) The name of the second comparison group.
#'
#' @return A loaded object containing calculated MVP information from given comparison.
#'
#' @importFrom ChAMP champ.MVP
#'
#' @export


limma.MVP <- function(file_name, group1, group2, ...){

    # Check for normalized methylation information
    if(!exists("batchNorm")){
        stop("Load ComBat normalized methylation information")
    }

    # Run champ.MVP function to determine MVPs
    if(!file.exists(file_name)){
        limma.comparison <- champ.MVP(beta.norm = batchNorm$beta,
                                      arraytype="EPIC",
                                      compare.group = c(group1, group2),
                                      adjust.method = "fdr",
                                      adjPVal = 0.05,
                                      ...
                                      )
    }

    # Save MVP data frame objects
    if(!file.exists(file_name)){
        save(limma.comparison, file = file_name)
    }

    # Load the previously saved MVP data frame objects
    load(file_name, envir = .GlobalEnv)

}
