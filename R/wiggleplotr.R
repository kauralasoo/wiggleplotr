### External functions to make the actual plots

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

wiggleplotr <- function(exons, cdss, sample_data, transcript_annotations, rescale_introns = TRUE,
                        new_intron_length = 50, plot_fraction = 0.1, heights = c(0.75, 0.25),
                        flanking_length = c(50,50)){
  #transcript_annotations is a data.frame with 3 columns: transcript_id, gene_id, gene_name.
  
  #Join exons together
  joint_exons = joinExons(exons)
  
  #Read coverage tracks from BigWig file
  sample_list = as.list(sample_data$bigWig)
  names(sample_list) = sample_data$sample_id
  gene_range = reduce(c(joint_exons,gaps(joint_exons, start = NA, end = NA)))
  seqlevels(gene_range) = IRanges::as.vector(seqnames(gene_range))[1]
  coverage_list = lapply(sample_list, readCoverageFromBigWig, gene_range, flanking_length = flanking_length)
  
  #Shorten introns and translate exons into the new introns
  if(rescale_introns){
    #Recale transcript annotations
    tx_annotations = rescaleIntrons(exons, cdss, joint_exons, 
                    new_intron_length = new_intron_length, flanking_length = flanking_length)
  }
  else{ #Do not rescale transcript annotationn
    #Need to calculate joint intron coordinates for transcript annotations
    old_introns = intronsFromJointExonRanges(ranges(joint_exons), flanking_length = flanking_length)
    tx_annotations = list(exon_ranges = lapply(exons, ranges), cds_ranges = lapply(cdss, ranges),
                          old_introns = old_introns, new_introns = old_introns)
  }
  
  #Shrink intron coverage and convert coverage vectors into data frames
  coverage_list = lapply(coverage_list, shrinkIntronsCoverage, 
                         tx_annotations$old_introns, tx_annotations$new_introns)
  
  #Define the start and end coorinates of the region
  region_start = min(start(tx_annotations$new_introns))
  region_end = max(end(tx_annotations$new_introns))
  region_length = region_end - region_start
  limits = c(region_start, region_end)
  
  #Take a subsample of points that's easier to plot
  points = sample(region_length, floor(region_length*plot_fraction))
  #Subtract the start coordinate of the region
  exon_starts = unique(unlist(lapply(tx_annotations$exon_ranges, start))) - (region_start -1)
  exon_ends = unique(unlist(lapply(tx_annotations$exon_ranges, end))) - (region_start - 1)
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
  #Construct transcript structure data.frame from ranges lists
  transcript_struct = prepareTranscriptStructureForPlotting(tx_annotations$exon_ranges, 
                                                            tx_annotations$cds_ranges, transcript_annotations)
  tx_structure = plotTranscriptStructure(transcript_struct, limits)
  
  coverage_plot = plotCoverage(coverage_df, limits)
  plot = gridExtra::arrangeGrob(coverage_plot, tx_structure, heights = heights, ncol = 1, nrow = 2)
  return(plot)
}