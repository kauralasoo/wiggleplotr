extractTranscriptAnnotationsFromEnsembldb <- function(ensembldb, gene_names, transcript_ids){

  #Fetch gene metadata
  gene_filter = AnnotationFilter::GeneNameFilter(gene_names)
  gene_metadata = ensembldb::transcripts(ensembldb, filter = gene_filter) %>%
    as.data.frame() %>%
    dplyr::transmute(transcript_id = tx_id, gene_name, strand)
  
  #Fetch gene exons and cdss
  exons = ensembldb::exonsBy(ensembldb, filter = gene_filter)
  cdss = ensembldb::cdsBy(ensembldb, filter = gene_filter)
  
  if(!is.null(transcript_ids)){
    gene_metadata = dplyr::filter(gene_metadata, transcript_id %in% transcript_ids)
    exons = exons[intersect(names(exons), transcript_ids)]
    cdss = cdss[intersect(names(cdss), transcript_ids)]
  }
  return(list(exons = exons, cdss = cdss, transcript_annotations = gene_metadata))
}

extractTranscriptAnnotationsFromUCSC <- function(orgdb, txdb, gene_names, transcript_ids = NULL){
  gene_meta = AnnotationDbi::select(orgdb, keys = gene_names, columns = c("SYMBOL", "UCSCKG"), keytype = "SYMBOL")
  colnames(gene_meta) = c("gene_name", "transcript_id")
  
  if(!is.null(transcript_ids)){
    gene_meta = dplyr::filter(gene_meta, transcript_id %in% transcript_ids)
  }
  
  #Extract exons and cdss
  tx_list = stats::setNames(as.list(gene_meta$transcript_id), gene_meta$transcript_id)
  exons_list = purrr::map(tx_list, ~GenomicFeatures::exons(txdb, filter = list(tx_name = .)))
  cds_list = purrr::map(tx_list, ~GenomicFeatures::cds(txdb, filter = list(tx_name = .)))
  
  #Add strand to gene gene metada
  gene_meta = dplyr::mutate(gene_meta, strand = extractStrandsFromGrangesList(exons_list))
  
  return(list(exons = exons_list, cdss = cds_list, transcript_annotations = gene_meta))
  
}


#' Plot transcripts directly from ensembldb object.
#' 
#' A wrapper around the plotTranscripts function. See the documentation for (\code{\link[wiggleplotr]{plotTranscripts}}) 
#' for more information.
#'
#' @param ensembldb ensembldb object.
#' @param gene_names List of gene names to be plotted.
#' @param transcript_ids Optional list of transcript ids to be plotted.
#' @param ... Additional parameters to be passed to plotTranscripts
#'
#' @return ggplot2 object
#' @export
#'
#' @examples
#' require("EnsDb.Hsapiens.v86")
#' plotTranscriptsFromEnsembldb(EnsDb.Hsapiens.v86, "NCOA7", transcript_ids = c("ENST00000438495", "ENST00000392477"))
plotTranscriptsFromEnsembldb <- function(ensembldb, gene_names, transcript_ids = NULL, ...){
  tx_annot = extractTranscriptAnnotationsFromEnsembldb(ensembldb, gene_names, transcript_ids)
  plotTranscripts(exons = tx_annot$exons, 
                  cdss = tx_annot$cdss, 
                  transcript_annotations = tx_annot$transcript_annotations, ...)
}

#' Plot read coverage directly from ensembldb object.
#' 
#' A wrapper around the plotCoverage function. See the documentation for (\code{\link[wiggleplotr]{plotCoverage}}) 
#' for more information.
#'
#' @param ensembldb ensembldb object.
#' @param gene_names List of gene names to be plotted.
#' @param transcript_ids Optional list of transcript ids to be plotted.
#' @param ... Additional parameters to be passed to plotCoverage.
#'
#' @return ggplot2 object
#' @export
#'
#' @examples
#' require("EnsDb.Hsapiens.v86")
#' require("dplyr")
#' require("GenomicRanges")
#' sample_data = dplyr::data_frame(sample_id = c("aipt_A", "aipt_C", "bima_A", "bima_C"), 
#'  condition = factor(c("Naive", "LPS", "Naive", "LPS"), levels = c("Naive", "LPS")), 
#'  scaling_factor = 1) %>%
#'  dplyr::mutate(bigWig = system.file("extdata",  paste0(sample_id, ".str2.bw"), package = "wiggleplotr"))
#'  
#' track_data = dplyr::mutate(sample_data, track_id = condition, colour_group = condition)
#' \dontrun{
#' plotCoverageFromEnsembldb(EnsDb.Hsapiens.v86, "NCOA7", transcript_ids = c("ENST00000438495", "ENST00000392477"), 
#' track_data, heights = c(2,1), fill_palette = getGenotypePalette())
#' }
plotCoverageFromEnsembldb <- function(ensembldb, gene_names, transcript_ids = NULL, ...){
  tx_annot = extractTranscriptAnnotationsFromEnsembldb(ensembldb, gene_names, transcript_ids)
  plotCoverage(exons = tx_annot$exons, 
                  cdss = tx_annot$cdss, 
                  transcript_annotations = tx_annot$transcript_annotations, ...)
}


#' Plot transcripts directly from UCSC OrgDb and TxDb objects.
#'
#' A wrapper around the plotTranscripts function. See the documentation for (\code{\link[wiggleplotr]{plotTranscripts}}) 
#' for more information. Note that this function is much slower than (\code{\link[wiggleplotr]{plotTranscripts}}) or 
#' (\code{\link[wiggleplotr]{plotTranscriptsFromEnsembldb}}) functions,  because indivudally extracting exon 
#' coordinates from txdb objects is quite inefficient.
#' @param orgdb UCSC OrgDb object.
#' @param txdb UCSC TxDb obejct.
#' @param gene_names List of gene genaes to be plot.
#' @param transcript_ids Optional list of transcript ids to be plot. (default = NULL)
#' @param ... Additional parameters to be passed to plotTranscripts
#'
#' @return Transcript plot.
#' @export
#'
#' @examples
#' #Load OrgDb and TxDb objects with UCSC gene annotations
#' require("org.Hs.eg.db")
#' require("TxDb.Hsapiens.UCSC.hg38.knownGene")
#' orgdb = org.Hs.eg.db
#' txdb = TxDb.Hsapiens.UCSC.hg38.knownGene
#' 
#' plotTranscriptsFromUCSC(orgdb, txdb, "NCOA7", transcript_ids = c("ENST00000438495.6", "ENST00000392477.6"))
plotTranscriptsFromUCSC <- function(orgdb, txdb, gene_names, transcript_ids = NULL, ...){
  tx_annot = extractTranscriptAnnotationsFromUCSC(orgdb, txdb, gene_names, transcript_ids)
  plotTranscripts(exons = tx_annot$exons, 
                  cdss = tx_annot$cdss, 
                  transcript_annotations = tx_annot$transcript_annotations, ...)
}


#' Plot read coverage directly from UCSC OrgDb and TxDb objects.
#' 
#' A wrapper around the plotCoverage function. See the documentation for (\code{\link[wiggleplotr]{plotCoverage}}) 
#' for more information.
#'
#' @param orgdb UCSC OrgDb object.
#' @param txdb UCSC TxDb obejct.
#' @param gene_names List of gene names to be plotted.
#' @param transcript_ids Optional list of transcript ids to be plotted.
#' @param ... Additional parameters to be passed to plotCoverage.
#'
#' @return ggplot2 object
#' @export
#'
#' @examples
#' require("dplyr")
#' require("GenomicRanges")
#' require("org.Hs.eg.db")
#' require("TxDb.Hsapiens.UCSC.hg38.knownGene")
#' 
#' orgdb = org.Hs.eg.db
#' txdb = TxDb.Hsapiens.UCSC.hg38.knownGene
#' 
#' sample_data = dplyr::data_frame(sample_id = c("aipt_A", "aipt_C", "bima_A", "bima_C"), 
#'  condition = factor(c("Naive", "LPS", "Naive", "LPS"), levels = c("Naive", "LPS")), 
#'  scaling_factor = 1) %>%
#'  dplyr::mutate(bigWig = system.file("extdata",  paste0(sample_id, ".str2.bw"), package = "wiggleplotr"))
#'  
#' track_data = dplyr::mutate(sample_data, track_id = condition, colour_group = condition)
#' \dontrun{
#' #Note: This example does not work, becasue UCSC and Ensembl use different chromosome names
#' plotCoverageFromUCSC(orgdb, txdb, "NCOA7", transcript_ids = c("ENST00000438495.6", "ENST00000392477.6"), 
#' track_data, heights = c(2,1), fill_palette = getGenotypePalette())
#' }
plotCoverageFromUCSC <- function(orgdb, txdb, gene_names, transcript_ids = NULL, ...){
  tx_annot = extractTranscriptAnnotationsFromUCSC(orgdb, txdb, gene_names, transcript_ids)
  plotCoverage(exons = tx_annot$exons, 
                  cdss = tx_annot$cdss, 
                  transcript_annotations = tx_annot$transcript_annotations, ...)
}
