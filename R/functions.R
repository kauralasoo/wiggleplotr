#Helper function to make wiggle plots

readCoverageFromBigWig <- function(bigwig_path, gene_range){
  #Read coverage over a region from a bigWig file
  sel = rtracklayer::BigWigSelection(gene_range)
  coverage_ranges = rtracklayer::import.bw(bigwig_path, selection = sel)
  GenomeInfoDb::seqlevels(coverage_ranges) = IRanges::as.vector(GenomicRanges::seqnames(gene_range))
  coverage_rle = GenomicRanges::coverage(coverage_ranges, weight = GenomicRanges::score(coverage_ranges))[[1]]
  coverage_rle = coverage_rle[(GenomicRanges::start(gene_range)):(GenomicRanges::end(gene_range))] #Keep the region of interest
}

joinExons <- function(exons) {
  #Join a list of exons into one GRanges object
  
  #Test that all transcripts are on the same chromosome
  chrs = unlist(lapply(exons, function(x) IRanges::as.vector(GenomicRanges::seqnames(x)[1])))
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
  joint_exons = GenomicRanges::reduce(joint_exons)
  return(joint_exons)
}

saveCoveragePlots <- function(plot_list, path, ...){
  #Save a list of plots into the folder specified by path
  gene_names = names(plot_list)
  for (gene in gene_names){
    file_name = file.path(path, paste(gene, ".pdf", sep = ""))
    ggplot2::ggsave(file_name, plot_list[[gene]], ...)
  }
}

prepareTranscriptStructureForPlotting <- function(exon_ranges, cds_ranges, transcript_annotations, label_type){
  #Combine exon_ranges and cds_ranges into a single data.frame that also contains transcript rank
  
  #Convert exon ranges into data.frame and add transcript rank
  exons_df = plyr::ldply(exon_ranges, data.frame)
  colnames(exons_df)[colnames(exons_df) == ".id"] = "transcript_id"
  exons_df = dplyr::mutate(exons_df, transcript_rank = as.numeric(factor(exons_df$transcript_id)), type = "annotations")
  transcript_rank = unique(exons_df[,c("transcript_id", "transcript_rank", "type")])
  
  #Convert CDS ranges into a data.frame
  cds_df = plyr::ldply(cds_ranges, data.frame)
  colnames(cds_df)[colnames(cds_df) == ".id"] = "transcript_id"
  cds_df = dplyr::left_join(cds_df, transcript_rank, by = "transcript_id") #Add matching transcript rank
  
  #Join exons and cdss together
  exons_df = dplyr::mutate(exons_df, feature_type = "exon")
  cds_df = dplyr::mutate(cds_df, feature_type = "cds")
  transcript_struct = rbind(exons_df, cds_df)
  
  #Prepare transcript annotations for plotting:
  #Keep only required columns
  transcript_annotations = dplyr::select(transcript_annotations, transcript_id, gene_id, gene_name, strand) %>%
    dplyr::mutate(strand = ifelse(strand == "+" | strand == 1, 1, -1)) #Change strand indicator to number is specified by character

  #Add additional metadata
  transcript_struct = dplyr::left_join(transcript_struct, transcript_annotations, by = "transcript_id") #Add gene name
  #Construct a label for each transcript
  if(label_type == "transcript"){
    transcript_struct = dplyr::mutate(transcript_struct, transcript_label = ifelse(strand == 1, 
                        paste(paste(gene_name, transcript_id, sep = ":")," >",sep =""), 
                        paste("< ",paste(gene_name, transcript_id, sep = ":"),sep =""))) 
  } else if(label_type == "peak"){
    transcript_struct = dplyr::mutate(transcript_struct, transcript_label = gene_name)
  }
  
  return(transcript_struct)
}

intronsFromJointExonRanges <- function(joint_exon_ranges, flanking_length){
  #Construct intron ranges from joint exon ranges
  introns = IRanges::gaps(joint_exon_ranges, 
                     start = min(IRanges::start(joint_exon_ranges)) - flanking_length[1], 
                     end = max(IRanges::end(joint_exon_ranges)) + flanking_length[2])
  return(introns)
}

#' Find the start and end coordinates of the whole gene form joint exons. 
constructGeneRange <- function(joint_exon_ranges, flanking_length){
  gene_range = GenomicRanges::reduce(c(joint_exon_ranges, GenomicRanges::gaps(joint_exon_ranges, start = NA, end = NA)))
  GenomeInfoDb::seqlevels(gene_range) = IRanges::as.vector(GenomicRanges::seqnames(gene_range))[1]
  GenomicRanges::start(gene_range) = GenomicRanges::start(gene_range) - flanking_length[1]
  GenomicRanges::end(gene_range) = GenomicRanges::end(gene_range) + flanking_length[2]
  return(gene_range)
}

#' Paste two factors together and preserved their joint order.
#'
#' @param factor1 First factor
#' @param factor2 Second factor
pasteFactors <- function(factor1, factor2){
  #Extract levels
  levels1 = levels(factor1)
  levels2 = levels(factor2)
  
  #Construct joint levels
  new1 = rep(levels1, length(levels2))
  new2 = rep(levels2, each = length(levels1))
  new_levels = paste(new1, new2, sep = "_")
  
  new_factor = factor(paste(factor1, factor2, sep = "_"), levels = new_levels)
  return(new_factor)
}

#' Calculate mean coverage within each track_id and colour_group
meanCoverage <- function(coverage_df){
  coverage_df = group_by(coverage_df, track_id, colour_group, bins) %>% 
    dplyr::summarise(coverage = mean(coverage)) %>%
    dplyr::mutate(sample_id = pasteFactors(track_id, colour_group)) %>% #Construct a new sample id for mean vector
    dplyr::ungroup()
  return(coverage_df)
}

#' Choose a subsample of points to make plotting faster
#' 
#' Makes sure that intron-exon boundaries are well samples.
subsamplePoints <- function(tx_annotations, plot_fraction){
  #Define the start and end coorinates of the region
  region_start = min(IRanges::start(tx_annotations$new_introns))
  region_end = max(IRanges::end(tx_annotations$new_introns))
  region_length = region_end - region_start

  #Take a subsample of points that's easier to plot
  points = sample(region_length, floor(region_length*plot_fraction))
  #Subtract the start coordinate of the region
  exon_starts = unique(unlist(lapply(tx_annotations$exon_ranges, start))) - (region_start -1)
  exon_ends = unique(unlist(lapply(tx_annotations$exon_ranges, end))) - (region_start - 1)
  points = unique(sort(c(points, exon_starts, exon_ends, 
                         exon_starts -3, exon_starts +3, 
                         exon_ends + 3, exon_ends -3)))
  points = points[points >= 0]
  return(points)
}

#' Returns a three-colour colour palette suitable for plotting coverage stratified by genotype
getGenotypePalette <- function(){
  c("#d7191c","#fdae61","#1a9641")
}
