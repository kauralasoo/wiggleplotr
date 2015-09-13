#Helper function to make wiggle plots

readCoverageFromBigWig <- function(bigwig_path, gene_range, flanking_length){
  #Read coverage over a region from a bigWig file
  sel = BigWigSelection(gene_range)
  coverage_ranges = rtracklayer::import.bw(bigwig_path, selection = sel)
  seqlevels(coverage_ranges) = IRanges::as.vector(seqnames(gene_range))
  coverage_rle = coverage(coverage_ranges, weight = score(coverage_ranges))[[1]]
  coverage_rle = coverage_rle[(start(gene_range)-flanking_length[1]):(end(gene_range)+flanking_length[2])] #Keep the region of interest
}

joinExons <- function(exons) {
  #Join a list of exons into one GRanges object
  
  #Test that all transcripts are on the same chromosome
  chrs = unlist(lapply(exons, function(x) IRanges::as.vector(seqnames(x)[1])))
  if (!all(chrs == chrs[1])){
    stop("Some transcripts are on different chromosomes.")
  }
  
  #Join all exons together
  transcript_ids = names(exons)
  joint_exons = c()
  for(tx_id in transcript_ids){
    tx = exons[[tx_id]]
    if(length(joint_exons) == 0){
      joint_exons = tx
    }
    else{
      joint_exons = c(joint_exons, tx)
    }
  }
  joint_exons = reduce(joint_exons)
  return(joint_exons)
}

plotGenesCoverage <- function(transcript_list, tx_annot, exons, cdss, sample_data, heights = c(0.75, 0.25)){
  gene_names = names(transcript_list)
  results = list()
  for (gene in gene_names){
    print(gene)
    transcripts = transcript_list[[gene]]
    gene_exons = exons[intersect(transcripts, names(exons))]
    gene_cdss = cdss[intersect(transcripts, names(cdss))]
    coverage_plot = wiggleplotr(gene_exons, gene_cdss, sample_data, tx_annot, plot_fraction = 0.2, heights = heights)
    results[[gene]] = coverage_plot
  }
  return(results)
}

saveCoveragePlots <- function(plot_list, path, width, height){
  #Save a list of plots into the folder specified by path
  gene_names = names(plot_list)
  for (gene in gene_names){
    file_name = file.path(path, paste(gene, ".pdf", sep = ""))
    ggsave(file_name, plot_list[[gene]], width = width, height = height)
  }
}

prepareTranscriptStructureForPlotting <- function(exon_ranges, cds_ranges, transcript_annotations){
  #Combine exon_ranges and cds_ranges into a single data.frame that also contains transcript rank
  
  #Convert exon ranges into data.frame and add transcript rank
  exons_df = plyr::ldply(exon_ranges, data.frame)
  colnames(exons_df)[colnames(exons_df) == ".id"] = "transcript_id"
  exons_df = dplyr::mutate(exons_df, transcript_rank = as.numeric(factor(exons_df$transcript_id)), type = "isoforms")
  transcript_rank = unique(exons_df[,c("transcript_id", "transcript_rank", "type")])
  
  #Convert CDS ranges into a data.frame
  cds_df = plyr::ldply(cds_ranges, data.frame)
  colnames(cds_df)[colnames(cds_df) == ".id"] = "transcript_id"
  cds_df = dplyr::left_join(cds_df, transcript_rank, by = "transcript_id") #Add matching transcript rank
  
  #Join exons and cdss together
  exons_df = dplyr::mutate(exons_df, feature_type = "exon")
  cds_df = dplyr::mutate(cds_df, feature_type = "cds")
  transcript_struct = rbind(exons_df, cds_df)
  
  #Add additional metadata
  transcript_struct = dplyr::left_join(transcript_struct, transcript_annotations, by = "transcript_id") %>% #Add gene name
    #Construct a label for each transcript
    dplyr::mutate(transcript_label = ifelse(strand == 1, 
                    paste(paste(gene_name, transcript_id, sep = ":")," >",sep =""), 
                    paste("< ",paste(gene_name, transcript_id, sep = ":"),sep =""))) 
  
  return(transcript_struct)
}

intronsFromJointExonRanges <- function(joint_exon_ranges, flanking_length){
  #Construct intron ranges from joint exon ranges
  introns = gaps(joint_exon_ranges, 
                     start = min(start(joint_exon_ranges)) - flanking_length[1], 
                     end = max(end(joint_exon_ranges)) + flanking_length[2])
  return(introns)
}