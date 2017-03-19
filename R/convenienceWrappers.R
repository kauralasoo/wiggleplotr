extractTranscriptAnnotationsFromEnsembldb <- function(ensembldb, gene_ids, transcript_ids){
  require(ensembldb)
  
  #Fetch gene metadata
  gene_filter = ensembldb::GenenameFilter(gene_ids)
  gene_metadata = ensembldb::transcripts(ensembldb, filter = gene_filter) %>%
    as.data.frame() %>%
    dplyr::transmute(transcript_id = tx_id, gene_name, strand)
  
  #Fetch gene exons and cdss
  exons = exonsBy(ensembldb, filter = gene_filter)
  cdss = cdsBy(ensembldb, filter = gene_filter)
  
  if(!is.null(transcript_ids)){
    gene_metadata = dplyr::filter(gene_metadata, transcript_id %in% transcript_ids)
    exons = exons[intersect(names(exons), transcript_ids)]
    cdss = cdss[intersect(names(cdss), transcript_ids)]
  }
  return(list(exons = exons, cdss = cdss, transcript_annotations = gene_metadata))
}


plotTranscriptsFromEnsembldb <- function(ensembldb, gene_ids, transcript_ids = NULL, ...){
  tx_annot = extractTranscriptAnnotationsFromEnsembldb(ensembldb, gene_ids, transcript_ids)
  plotTranscripts(exons = tx_annot$exons, 
                  cdss = tx_annot$cdss, 
                  transcript_annotations = tx_annot$gene_metadata, ...)
}

plotCoverageFromEnsembldb <- function(ensembldb, gene_ids, transcript_ids = NULL, ...){
  tx_annot = extractTranscriptAnnotationsFromEnsembldb(ensembldb, gene_ids, transcript_ids)
  plotCoverage(exons = tx_annot$exons, 
                  cdss = tx_annot$cdss, 
                  transcript_annotations = tx_annot$gene_metadata, ...)
}

