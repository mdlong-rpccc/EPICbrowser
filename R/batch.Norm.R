#' Implement the ComBat normalization method developed for microarrays via the sva package.
#' This corrects for batch effects related to the slide (Sentrix_ID)
#' Uses an empirical Bayes method to correct for technical variation (after logit transformation, required as beta values are 0-1, and followed by reverse logit transformation)
#' Will abort if all samples are from a single slide
#'
#' @details This approach uses the ChAMP package first described by Morris TJ et al:
#' Morris TJ, Butcher LM, Teschendorff AE, Chakravarthy AR, Wojdacz TK and Beck S (2014). “ChAMP: 450k Chip Analysis Methylation
#' Pipeline.” _Bioinformatics_, *30*(3), pp. 428-430. doi: 10.1093/bioinformatics/btt684
#' (URL:http://doi.org/10.1093/bioinformatics/btt684).
#'
#' @param myNorm_name character(1) The name of the output file returned from the champ.norm function.
#'
#' @return A loaded object containing ComBat normalized beta values.
#'
#' @importFrom ChAMP champ.runCombat
#' @importFrom pheatmap pheatmap
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette pdf
#'
#' @export


batch.Norm <- function(batchNorm_name = "batchNorm_ComBat_normalized_EPIC_betaVals"){


    if(!file.exists(paste(batchNorm_name, ".RData", sep = ""))){
        # Create stop if myLoad isn't available
        if(!exists("myLoad", envir = .GlobalEnv) | !class(myLoad) == "list"){
            stop("Run idat.Load() to import IDAT files")
        }
        # Create stop if myNorm isn't available
        if(!exists("myNorm", envir = .GlobalEnv) | !class(myNorm) == "list"){
            stop("Run bmiq.Norm() to run BMIQ")
        }
        # Run/Save batchNorm using champ.runCombat
        batchNorm <- champ.runCombat()
        save(batchNorm, file = paste(batchNorm_name, ".RData", sep = ""))
    }

    # If champ.runCombat() previously utilized to ComBat normalize data and subsequently
    # saved you can load from here.
    load(paste(batchNorm_name, ".RData", sep = ""), envir = .GlobalEnv)

    # Convert myNorm to data.frame
    batchNorm.table <- as.data.frame(batchNorm)

    # Write out BMIQ normalized, ComBat adjusted beta values
    write.csv(batchNorm, "ChAMP-normalized_ComBat-adjusted_EPIC-beta-values_A.csv")

    # Read in BMIQ normalized beta values for each sample
    batchNorm.table <- read.csv("ChAMP-normalized_ComBat-adjusted_EPIC-beta-values_A.csv", header = TRUE, row.names = 1)

    # Plot hierarchical clustering (with euclidean distance) of normalized beta-values
    if(!file.exists(paste(batchNorm_name, "_clustering.pdf", sep = ""))){
    pdf(file = paste(batchNorm_name, "_clustering.pdf", sep = ""))
    plot(hclust(dist(t(batchNorm.table))))
    dev.off()
    }

    # Plot Heatmap of euclidean distance
    if(!file.exists(paste(batchNorm_name, "_distanceMatrix.pdf", sep = ""))){
    colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
    labels <- rownames(myLoad$pd)
    pdf(file = paste(batchNorm_name, "_distanceMatrix.pdf", sep = ""))
    pheatmap(as.matrix(dist(t(batchNorm.table))),
             labels_col = labels,
             labels_row = labels,
             col=colors)
    dev.off()
    }
}

