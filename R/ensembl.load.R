#' Implements the 'ensembl.load' function to import ensembl (hg19) gene annotations
#'
#' @details The ensembl.load() function is used to import and organize the ensembl annotation files (hg19/grch37) using
#'          biomaRt for further processing. Specifically, TSS location is determined based on strand information and TxStart/TxEnd.
#'          The function looks for the original file, 'refGene_nochr.txt' (by default, set by the refSeq.pre argument),
#'          in the current directory (by default, set bu the file_path argument). It filters information for only genes
#'          on chromosomes 1-22, X, Y & MT, then identifies TSS locations as described above. The file is sorted by
#'          location and resaved as "refGene_hg19_info.txt" (by default, set the refSeq.post argument). If this file
#'          has been saved, the function will directly load it instead.
#'
#' @param file_path character(1) The name of the directory containing the refSeq annotation file.
#' @param file_name character(1) The name of the saved ensembl annotation file.
#' @param ensembl_attributes character(1) The name of the ensembl attributes imported through biomaRt.
#'
#' @return A data frame containing ensembl (hg19) annotations with TSS information
#'
#'
#' @importFrom biomaRt getBM
#'
#' @export


ensembl.load <- function(file_path = getwd(), file_name = "grch37_geneInfo.txt",
                         ensembl_attributes = c("ensembl_gene_id", "ensembl_transcript_id",
                                                "chromosome_name", "start_position",
                                                "end_position", "strand",
                                                "transcript_start", "transcript_end",
                                                "transcription_start_site", "transcript_length",
                                                "external_gene_name", "percentage_gc_content")){

    # Create and/or Set working directory
    if(!dir.exists(paste(file_path, "/grch37", sep = ""))){
        dir.create(paste(file_path, "/grch37", sep = ""))
    }
    setwd(paste(file_path, "/grch37", sep = ""))

    # IMPORTING ENSEMBL TSS site/GENE/EXON LENGTH INFORMATION FROM ENSEMBL USING BIOMART.
    # ONCE COMPLETED, FILES ARE SAVED AND THIS IS SKIPPED #
    # getBM function is the main query function in biomaRt having 4 main arguments: attributes, filters, values, mart
    if(!file.exists(file_name)){
    grch37_data <- getBM(attributes = c(ensembl_attributes),
                         mart = grch37)

    # Convert grch37 PATCH information to correct chromosome number
    PATCH <- read.csv("grch37_PATCH_conversion.csv")
    for(i in 1:length(grch37_data$chromosome_name)){
        if(grch37_data$chromosome_name[i] %in% PATCH$patch){
            grch37_data$chromosome_name[i] <- as.character(PATCH$chromosome[PATCH$patch == grch37_data$chromosome_name[i]])
        }
    }
    # Order gene information by location
    grch37_data <- grch37_data[order(grch37_data$chromosome_name, grch37_data$start_position) ,]
    # Change strand information to '+' and '-'
    grch37_data$strand <- gsub("^1", "+", grch37_data$strand)
    grch37_data$strand <- gsub("-1$", "-", grch37_data$strand)
    # Filter for only genes with a TSS
    grch37_data <- subset(grch37_data, grepl("[0-9]", grch37_data$transcription_start_site))
    # Filter out contigs and haplotype chromosomes
    grch37_data <- subset(grch37_data, grepl("^[^H|G]", grch37_data$chromosome_name))
    # Filter for only full chromosomes and remove unmappable/haplotype chromosomes
    chromosomes <- as.factor(c(1:22, "X", "Y", "MT"))
    grch37_data <- subset(grch37_data, grch37_data$chromosome_name %in% chromosomes)
    grch37_data <- droplevels(grch37_data)
    # Save grch37 gene information as .txt file
    write.table(grch37_data, "grch37_geneInfo.txt")
    }


    # Import ALL gene/data length information as downloaded through Ensembl Biomart -
    # http://useast.ensembl.org/biomart/martview/5142e6a944dfd9765baff37d45d84f17
    read.table("grch37_geneInfo.txt", header=TRUE)


}
