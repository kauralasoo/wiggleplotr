#' wiggleplotr
#' 
#' wiggleplotr package provides tools to visualise transcript annotations (\code{\link[wiggleplotr]{plotTranscripts}}) and plot 
#' sequencing read coverage over annotated transcripts (\code{\link[wiggleplotr]{plotCoverage}}). 
#' 
#' You can also use covenient wrapper functions 
#' (\code{\link[wiggleplotr]{plotTranscriptsFromEnsembldb}}), (\code{\link[wiggleplotr]{plotCoverageFromEnsembldb}}), 
#' (\code{\link[wiggleplotr]{plotTranscriptsFromUCSC}}) and (\code{\link[wiggleplotr]{plotCoverageFromUCSC}}).
#' 
#' To learn more about wiggleplotr, start with the vignette: 
#' \code{browseVignettes(package = "wiggleplotr")}
#'
#' @name wiggleplotr
#' @docType package
#' @import ggplot2
#' @importFrom dplyr "%>%"
#' @importFrom dplyr "row_number"
utils::globalVariables(c("strand","gene_name","transcript_id", "tx_id"))
