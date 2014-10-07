readCoverageFromBigWig <- function(bigwig_path, joint_exons, flanking = 50){
  #Read coverage over a region from a bigWig file
  chromosome = IRanges::as.vector(seqnames(joint_exons))[1]
  seqlevels(joint_exons) = chromosome
  sel = BigWigSelection(joint_exons)
  coverage_ranges = rtracklayer::import.bw(bigwig_path, selection = sel)
  seqlevels(coverage_ranges) = chromosome
  coverage_rle = coverage(coverage_ranges, weight = score(coverage_ranges))[[1]]
  coverage_rle = coverage_rle[(min(start(joint_exons))-flanking):(max(end(joint_exons))+flanking)] #Keep the region of interest
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

wiggleplotr <- function(exons, cdss, sample_data, transcript_annotations, new_intron_length = 50, plot_fraction = 0.1){
  #Plot read coverage over exons
  transcript_ids = names(exons)
  exon_ranges = lapply(exons, ranges)
  
  #Join exons together into single GRanges object
  joint_exons = joinExons(exons)
  joint_ranges = ranges(joint_exons)
  
  #Shorten introns and translate exons into the new exons
  old_introns = gaps(joint_ranges, start = min(start(joint_ranges)) - new_intron_length, end = max(end(joint_ranges)) + new_intron_length)
  new_introns = shortenIntrons(old_introns,new_intron_length)
  new_exon_ranges = lapply(exon_ranges, translateExonCoordinates, old_introns, new_introns)
  
  #Convert exons list to data frame
  exons_df = plyr::ldply(new_exon_ranges, data.frame)
  colnames(exons_df)[colnames(exons_df) == ".id"] = "transcript_id"
  exons_df = dplyr::mutate(exons_df, transcript_rank = as.numeric(factor(exons_df$transcript_id)), type = "isoforms")
  transcript_rank = unique(exons_df[,c("transcript_id", "transcript_rank", "type")])
  
  #Translate cds coordinates
  cds_ranges = lapply(cdss, ranges)
  new_cds_ranges = lapply(cds_ranges, translateExonCoordinates, old_introns, new_introns)
  cds_df = plyr::ldply(new_cds_ranges, data.frame)
  colnames(cds_df)[colnames(cds_df) == ".id"] = "transcript_id"
  cds_df = plyr::join(cds_df, transcript_rank, by = "transcript_id") #Add matching transcript rank
  
  #Join exons and cdss together
  exons_df = dplyr::mutate(exons_df, feature_type = "exon")
  cds_df = dplyr::mutate(cds_df, feature_type = "cds")
  transcript_struct = rbind(exons_df, cds_df)
  transcript_struct = dplyr::left_join(transcript_struct, transcript_annotations, by = "transcript_id") #Add gene name
  
  #Read coverage tracks from BigWig file
  sample_list = as.list(sample_data$bigWig)
  names(sample_list) = sample_data$sample_id
  coverage_list = lapply(sample_list, readCoverageFromBigWig, joint_exons, new_intron_length)
  shrunken_coverage_list = lapply(coverage_list, shrinkIntronsCoverage, old_introns, new_introns)
  
  #Take a subsample of points that's easier to plot
  n_total = nrow(shrunken_coverage_list[[1]])
  points = sample(n_total, floor(n_total*plot_fraction))
  points = unique(sort(c(points, start(new_exon_ranges), end(new_exon_ranges), 
                         start(new_exon_ranges) -3, start(new_exon_ranges) +3, 
                         end(new_exon_ranges) + 3, end(new_exon_ranges) -3)))
  points = points[points >= 0]
  shrunken_coverage_list = lapply(shrunken_coverage_list, function(x) {x[points,]} )
  
  #Covert to data frame and plot
  coverage_df = plyr::ldply(shrunken_coverage_list, data.frame)
  colnames(coverage_df)[colnames(coverage_df) == ".id"] = "sample_id"
  coverage_df = plyr::join(coverage_df, sample_data, by = "sample_id")
  coverage_df = dplyr::mutate(coverage_df, coverage = coverage/library_size)

  #Make plots
  limits = c(0,n_total)
  tx_structure = plotTranscriptStructure(transcript_struct, limits)
  coverage_plot = plotCoverage(coverage_df, limits)
  plot = gridExtra::arrangeGrob(coverage_plot, tx_structure, heights = c(0.75, 0.25), ncol = 1, nrow = 2)
  return(plot)
}

plotGenesCoverage <- function(transcript_list, tx_annot, exons, cdss, sample_data){
  gene_names = names(transcript_list)
  results = list()
  for (gene in gene_names){
    print(gene)
    transcripts = transcript_list[[gene]]
    gene_exons = exons[intersect(transcripts, names(exons))]
    gene_cdss = cdss[intersect(transcripts, names(cdss))]
    coverage_plot = wiggleplotr(gene_exons, gene_cdss, sample_data, tx_annot, new_intron_length = 50, plot_fraction = 0.2)
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



