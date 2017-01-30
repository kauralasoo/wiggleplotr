#' Exons from 9 protein coding transcripts of NCOA7
#'
#' A dataset containing start and end coordinates of exons from nine 
#' protein coding transcripts of NCOA7.
#'
#' @format A GRangesList object with 9 elements:
#' \describe{
#'   \item{element}{Exon start and end coordinates for a single transcript (GRanges object)}
#'   ...
#' }
#' @source \url{http://www.ensembl.org/}
"ncoa7_exons"

#' Coding sequences from 9 protein coding transcripts of NCOA7
#'
#' A dataset containing start and end coordinates of coding sequences (CDS) from nine 
#' protein coding transcripts of NCOA7.
#'
#' @format A GRangesList object with 9 elements:
#' \describe{
#'   \item{element}{CDS start and end coordinates for a single transcript (GRanges object)}
#'   ...
#' }
#' @source \url{http://www.ensembl.org/}
"ncoa7_cdss"

#' Gene metadata for NCOA7
#'
#' A a list of transcripts for NCOA7.
#'
#' @format A data.frame object with 4 columns:
#' \describe{
#'   \item{transcript_id}{Ensembl transcript id.}
#'   \item{gene_id}{Ensembl gene id.}
#'   \item{gene_name}{Human readable gene name.}
#'   \item{strand}{Strand of the transcript (either +1 or -1).}
#'   ...
#' }
#' @source \url{http://www.ensembl.org/}
"ncoa7_metadata"