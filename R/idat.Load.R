#' Import IDAT files using 'champ.load'
#'
#' @details This approach uses the ChAMP package first described by Morris TJ et al:
#' Morris TJ, Butcher LM, Teschendorff AE, Chakravarthy AR, Wojdacz TK and Beck S (2014). “ChAMP: 450k Chip Analysis Methylation
#' Pipeline.” _Bioinformatics_, *30*(3), pp. 428-430. doi: 10.1093/bioinformatics/btt684
#' (URL:http://doi.org/10.1093/bioinformatics/btt684).
#'
#' @param IDAT_path character(1) The path to the IDAT data files.
#' @param myLoad_name character(1) The name of the output file returned from the champ.load function.
#'
#' @return A loaded object (list) containing imported IDAT information for further processing.
#'
#' @importFrom ChAMP champ.load
#'
#' @export


idat.Load <- function(IDAT_path, myLoad_name = "currentStudyloadedData"){

    if(!file.exists(paste(myLoad_name, ".RData", sep = ""))){
        # Load Data from idat files using minfi, and filter for failed probes using
        # minfi method (p < 0.01 in at least 1 sample)
        myLoad <- champ.load(directory = IDAT_path, arraytype = "EPIC")
        # The champ.load() function uses the most memory. Save this for future analyses:
        save(myLoad, file = paste(myLoad_name, ".RData", sep = ""))
    }

    # If champ.load() previously utilized to import data and subsequently
    # saved you can load from here.
    load(file = paste(myLoad_name, ".RData", sep = ""), envir = .GlobalEnv)
}
