#' Implements the 'refSeq.load' function to import refSeq (hg19) gene annotations
#'
#' @details The refseq.load() function is used to import and organize the refSeq annotation file (hg19) for further
#'          processing. Specifically, TSS location is determined based on strand information and TxStart/TxEnd.
#'          The function looks for the original file, 'refGene_nochr.txt' (by default, set by the refSeq.pre argument),
#'          in the current directory (by default, set bu the file_path argument). It filters information for only genes
#'          on chromosomes 1-22, X, Y & MT, then identifies TSS locations as described above. The file is sorted by
#'          location and resaved as "refGene_hg19_info.txt" (by default, set the refSeq.post argument). If this file
#'          has been saved, the function will directly load it instead.
#'
#' @param file_path character(1) The name of the directory containing the refSeq annotation file.
#' @param refSeq.pre character(1) The name of the file containing the original refSeq annotation file.
#' @param refSeq.post character(1) The name of the file containing the organized/saved refSeq
#'        annotation file with TSS information.
#'
#' @return A data frame containing refSeq (hg19) annotations with corrected TSS information
#'
#' @export



refSeq.load <- function(file_path = getwd(), refSeq.pre = "refGene_nochr.txt", refSeq.post = "refGene_hg19_info.txt"){

    # Set directory
    setwd(file_path)

    # Create a stop if refseq file is not present in directory
    if(!file.exists(refSeq.pre) & !file.exists(refSeq.post)){
        stop("refSeq information is not present in current directory. Use file_path argument to set directory
             and the refSeq.pre and/or refSeq.post arguments to set the file names.")
    }


    if(!file.exists(refSeq.post)){
        # Read in reference gene information (refseq, hg19)
        refGene <- read.table(refSeq.pre, header = FALSE)
        # Add column head information - https://genome.ucsc.edu/FAQ/FAQformat.html#format9
        colnames(refGene) <- c("num", "refseq_transcript_ID", "chr", "strand", "txStart", "txEnd", "CdsStart",
                               "CdsEnd", "NumExons", "ExonStarts", "ExonEnds", "Score", "geneID", "CdsStartStat",
                               "CdsEndStat", "ExonFrames")
        # Filter for only full chromosomes and remove unmappable/haplotype chromosomes
        chromosomes <- as.factor(c(1:22, "X", "Y", "MT"))
        refGene <- subset(refGene, refGene$chr %in% chromosomes)
        refGene <- droplevels(refGene)
        # Indicate RefSeq TSS sites
        refGene.TSS <- numeric()
        for(i in 1:length(refGene$refseq_transcript_ID)){
            if(refGene$strand[i] == "+"){
                refGene.TSS[i] = refGene$txStart[i]
            } else if (refGene$strand[i] == "-"){
                refGene.TSS[i] = refGene$txEnd[i]
            }
        }
        refGene$TSS <- refGene.TSS
        # Order gene information by location
        refGene <- refGene[order(refGene$chr, refGene$TSS) ,]
        # Save new refGene file
        write.table(refGene, "refGene_hg19_info.txt")
    }

    # Otherwise, load organized refGene (hg19) information file
    read.table("refGene_hg19_info.txt", header = TRUE)
    }

