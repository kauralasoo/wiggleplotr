readCoverageFromBigWig <- function(bigwig_path, gene_range, flanking = 50){
  #Read coverage over a region from a bigWig file
  sel = BigWigSelection(gene_range)
  coverage_ranges = rtracklayer::import.bw(bigwig_path, selection = sel)
  seqlevels(coverage_ranges) = IRanges::as.vector(seqnames(gene_range))
  coverage_rle = coverage(coverage_ranges, weight = score(coverage_ranges))[[1]]
  coverage_rle = coverage_rle[(start(gene_range)-flanking):(end(gene_range)+flanking)] #Keep the region of interest
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

wiggleplotr <- function(exons, cdss, sample_data, transcript_annotations, rescale_introns = TRUE,
                        new_intron_length = 50, plot_fraction = 0.1, heights = c(0.75, 0.25)){
  #transcript_annotations is a data.frame with 3 columns: transcript_id, gene_id, gene_name.
  
  #Join exons together
  joint_exons = joinExons(exons)
  
  #Read coverage tracks from BigWig file
  sample_list = as.list(sample_data$bigWig)
  names(sample_list) = sample_data$sample_id
  gene_range = reduce(c(joint_exons,gaps(joint_exons, start = NA, end = NA)))
  seqlevels(gene_range) = IRanges::as.vector(seqnames(gene_range))[1]
  coverage_list = lapply(sample_list, readCoverageFromBigWig, gene_range, flanking = new_intron_length)
  
  #Shorten introns and translate exons into the new introns
  if(rescale_introns){
    #Recale transcript annotations
    tx_annotations = rescaleIntrons(exons, cdss, joint_exons, new_intron_length = new_intron_length)
    #Rescale read coverage vectors
    coverage_list = lapply(coverage_list, shrinkIntronsCoverage, 
                                    tx_annotations$old_introns, tx_annotations$new_introns)  
  }
  else{
    #TODO: this behaviour not yet implemented
    tx_annotations = list(exon_ranges = lapply(exons, ranges), cds_ranges = lapply(cdss, ranges))
    coverage_list = coverage_list #TODO:fix this
  }
  
  #Take a subsample of points that's easier to plot
  n_total = nrow(coverage_list[[1]])
  points = sample(n_total, floor(n_total*plot_fraction))
  exon_starts = unique(unlist(lapply(tx_annotations$exon_ranges, start)))
  exon_ends = unique(unlist(lapply(tx_annotations$exon_ranges, end)))
  points = unique(sort(c(points, exon_starts, exon_ends, 
                         exon_starts -3, exon_starts +3, 
                         exon_ends + 3, exon_ends -3)))
  points = points[points >= 0]
  coverage_list = lapply(coverage_list, function(x) {x[points,]} )
  
  #Covert to data frame and plot
  coverage_df = plyr::ldply(coverage_list, data.frame)
  colnames(coverage_df)[colnames(coverage_df) == ".id"] = "sample_id"
  coverage_df = dplyr::left_join(coverage_df, sample_data, by = "sample_id") %>%
    dplyr::mutate(coverage = coverage/library_size) #Normalize by library size

  #Make plots
  limits = c(0,n_total)
  #Construct transcript structure data.frame from ranges lists
  transcript_struct = prepareTranscriptStructureForPlotting(tx_annotations$exon_ranges, 
                                                            tx_annotations$cds_ranges, transcript_annotations)
  tx_structure = plotTranscriptStructure(transcript_struct, limits)
  
  coverage_plot = plotCoverage(coverage_df, limits)
  plot = gridExtra::arrangeGrob(coverage_plot, tx_structure, heights = heights, ncol = 1, nrow = 2)
  return(plot)
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
    dplyr::mutate(transcript_label = ifelse(strand == "+", 
                    paste(paste(gene_name, transcript_id, sep = ":")," >",sep =""), 
                    paste("< ",paste(gene_name, transcript_id, sep = ":"),sep =""))) 
  
  return(transcript_struct)
}

plotTranscripts <- function(exons, cdss, annotations, rescale_introns = TRUE, new_intron_length = 50){
  #Function to quickly plot transcript structure without any read coverage
  # annotations: data.frame with 4 columns: transcript_id, gene_id, gene_name, strand
  
  #Join exons together
  joint_exons = joinExons(exons)
  
  #Rescale introns
  if (rescale_introns){
    tx_annotations = rescaleIntrons(exons, cdss, joint_exons, new_intron_length = new_intron_length)
  } else {
    tx_annotations = list(exon_ranges = lapply(exons, ranges), cds_ranges = lapply(cdss, ranges))
  }
  
  plot = prepareTranscriptStructureForPlotting(tx_annotations$exon_ranges, 
                                               tx_annotations$cds_ranges, annotations) %>%
    plotTranscriptStructure()
  return(plot)
}
