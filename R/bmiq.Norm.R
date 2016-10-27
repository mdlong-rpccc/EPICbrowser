#' Performs beta-mixture quantile (BMIQ) normalization using the champ.norm() function
#' BMIQ is an intra-sample normalization procedure, correcting for type-2 probe values (methylation arrays have two types of probes that yield different distributions).
#' 3-step procedure: (i) fitting of a 3-state beta mixture model, (ii) transformation of state-membership probabilities of type2 probes
#' into quantiles of the type1 distribution, and (iii) a conformal transformation for the hemi-methylated probes.#' Normalize Imported IDAT files using 'champ.norm'
#'
#' @details This approach uses the ChAMP package first described by Morris TJ et al:
#' Morris TJ, Butcher LM, Teschendorff AE, Chakravarthy AR, Wojdacz TK and Beck S (2014). “ChAMP: 450k Chip Analysis Methylation
#' Pipeline.” _Bioinformatics_, *30*(3), pp. 428-430. doi: 10.1093/bioinformatics/btt684
#' (URL:http://doi.org/10.1093/bioinformatics/btt684).
#'
#' @param myNorm_name character(1) The name of the output file returned from the champ.norm function.
#'
#' @return A loaded object containing BMIQ normalized beta values.
#'
#' @importFrom ChAMP champ.norm
#'
#' @export


bmiq.Norm <- function(myNorm_name = "myNorm_BMIQ_normalized_EPIC_betaVals"){

    if(!exists("myLoad", envir = .GlobalEnv) | !class(myLoad) == "list"){
        stop("Run idat.Load() to import IDAT files")
    }

    if(!file.exists(paste(myNorm_name, ".RData", sep = ""))){
        myNorm <- champ.norm(arraytype="EPIC")
        save(myNorm, file = paste(myNorm_name, ".RData", sep = ""))
    }

    # If champ.norm() previously utilized to normalize data and subsequently
    # saved you can load from here.
    load(paste(myNorm_name, ".RData", sep = ""), envir = .GlobalEnv)

    # Convert myNorm to data.frame
    myNorm.table <- as.data.frame(myNorm)

    # Plot hierarchical clustering (with euclidean distance) of normalized beta-values
    if(!file.exists(paste(myNorm_name, "_clustering.pdf", sep = ""))){
    pdf(file = paste(myNorm_name, "_clustering.pdf", sep = ""))
    plot(hclust(dist(t(myNorm.table))))
    dev.off()
    }
}
