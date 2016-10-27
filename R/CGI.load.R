#' Implements the 'CGI.load' function to import UCSC CpG island (hg19) annotations
#'
#' @details The CGI.load() function is used to import and organize the UCSC CpG island annotation files (hg19/grch37)
#'          further analyses. The function looks for the original file, 'cpgIslandExt.txt' (by default, set by the CGI_file argument),
#'          in the default directory (set by the file_path argument). It filters information for only CGIs
#'          on chromosomes 1-22, X, Y & MT. The file is sorted by location and resaved as "refGene_hg19_info.txt" (by default,
#'          set the refSeq.post argument). If this file has been saved, the function will directly load it instead.
#'
#' @param file_path character(1) The name of the directory containing the UCSC CGI annotation file.
#' @param CGI_file.pre character(1) The name of the original UCSC CGI annotation file.
#' @param CGI_file.post character(1) The name of the saved UCSC CGI file.
#'
#' @return A data frame containing UCSC CpG island (hg19).
#'
#'
#' @importFrom biomaRt getBM
#'
#' @export


CGI.load <- function(file_path = getwd(), CGI_file.pre = "cpgIslandExt", CGI_file.post = "UCSC_hg19_cpgIsland"){

    # Create and/or Set working directory
    if(!dir.exists(paste(file_path, "/grch37", sep = ""))){
        dir.create(paste(file_path, "/grch37", sep = ""))
    }
    setwd(paste(file_path, "/grch37", sep = ""))

    if(!file.exists(paste(CGI_file.post, ".txt", sep = ""))){
        # Read in CGI file (hg19)
        cpg_islands <- read.table(paste(CGI_file.pre, ".txt", sep = ""), header = FALSE, sep = "\t")
        cpg_islands[,1] <- NULL
        # Add column head information - http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/cpgIslandExt.sql
        colnames(cpg_islands) <- c("chr", "start", "end", "ID", "length", "CpgNum", "GcNum", "perCpG", "perGc", "obsExp")
        # Remove 'chr' string from chromosome information
        cpg_islands$chr <- gsub("chr", "", cpg_islands$chr)
        # Filter for only full chromosomes and remove unmappable/haplotype chromosomes
        chromosomes <- as.factor(c(1:22, "X", "Y", "MT"))
        cpg_islands <- subset(cpg_islands, cpg_islands$chr %in% chromosomes)
        cpg_islands <- droplevels(cpg_islands)
        # Save organized CGI file
        write.table(cpg_islands, paste(CGI_file.post, ".txt", sep = ""), row.names = FALSE)
    }

    # Import saved/organized UCSC CGI annotation file
    read.table(paste(CGI_file.post, ".txt", sep = ""), header = TRUE)
}

