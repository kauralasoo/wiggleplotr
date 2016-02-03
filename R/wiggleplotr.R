#' Function to quickly plot transcript structure without any read coverage
#'
#' @param exons 
#' @param cdss 
#' @param annotations 
#' @param rescale_introns 
#' @param new_intron_length 
#' @param flanking_length 
#'
#' @return ggplot2 object
#' @export
plotTranscripts <- function(exons, cdss, annotations, rescale_introns = TRUE, new_intron_length = 50, flanking_length = c(50,50)){
  #
  # annotations: data.frame with 4 columns: transcript_id, gene_id, gene_name, strand
  
  #Join exons together
  joint_exons = joinExons(exons)
  
  #Rescale introns
  if (rescale_introns){
    tx_annotations = rescaleIntrons(exons, cdss, joint_exons, new_intron_length = new_intron_length, flanking_length)
  } else {
    tx_annotations = list(exon_ranges = lapply(exons, ranges), cds_ranges = lapply(cdss, ranges))
  }
  
  plot = prepareTranscriptStructureForPlotting(tx_annotations$exon_ranges, 
                                               tx_annotations$cds_ranges, annotations) %>%
    plotTranscriptStructure()
  return(plot)
}

#' Plot read coverage across genomic regions
#' 
#' Also supports rescaling introns to constant length.
#' 
#' @param exons List of GRanges objects, each object containing exons for one transcript. 
#' The list must have names that correspond to transcript_id column in transript_annotations data.frame.
#' @param cdss list of GRanges objects, each object containing the coding regions (CDS) of a single transcript. 
#' The list must have names that correspond to transcript_id column in transript_annotations data.frame. 
#' If specified then coding and non-coding regions will be plotted in a different colour.
#' @param track_data data.frame with the metadata for the bigWig read coverage files. Must contain the following columns:
#' \itemize{
#'  \item sample_id - unique id for each sample.
#'  \item track_id - if multiple samples (bigWig files) have the same track_id they will be overlayed on the same 
#' plot, track_id is also used as the facet label on the right.
#'  \item bigWig - path to the bigWig file.
#'  \item scaling_factor - normalisation factor for each sample, useful if different samples sequenced to different 
#' depth and bigWig files not normalised for that.
#'  \item colour_group - additional column to group samples into, is used as the colour of the coverage track.
#' }
#' @param transcript_annotations Data.frame with three columns: transcript_id, gene_id, gene_name.
#' @param rescale_introns Specifies if the introns should be scaled to fixed length or not. (default: TRUE)
#' @param new_intron_length length (bp) of introns after scaling. (default: 50)
#' @param flanking_length Lengths of the flanking regions upstream and downstream of the gene. (default: c(50,50))
#' @param plot_fraction Size of the random sub-sample of points used to plot coverage (between 0 and 1). 
#' Smaller values make plotting significantly faster. (default: 0.1)
#' @param heights  Specifies the proportion of the height that is dedicated to coverage plots (first value) 
#' relative to transcript annotations (second value). (default: c(0.75,0.25))
#' @param alpha Transparency (alpha) value for the read coverage tracks. 
#' Useful to set to something < 1 when overlaying multiple tracks (see track_id). (default: 1)
#' @param fill_palette Colour palette used for the coverage tracks.
#' @param mean_only Plot only mean coverage within each combination of track_id and colour_group values. 
#' Useful for example for plotting mean coverage stratified by genotype (which is specified in the colour_group column) (default: TRUE).
#'
#' @return Object form gridExtra::arrangeGrob() function. Can we viewed using plot command or saved using ggsave.
#' @export
wiggleplotr <- function(exons, cdss, track_data, transcript_annotations, rescale_introns = TRUE,
                        new_intron_length = 50, flanking_length = c(50,50),
                        plot_fraction = 0.1, heights = c(0.75, 0.25), alpha = 1,
                        fill_palette = c("#a1dab4","#41b6c4","#225ea8"), mean_only = TRUE){
  #Join exons together
  joint_exons = joinExons(exons)
  
  #Read coverage tracks from BigWig file
  sample_list = as.list(track_data$bigWig)
  names(sample_list) = track_data$sample_id
  gene_range = reduce(c(joint_exons,gaps(joint_exons, start = NA, end = NA)))
  seqlevels(gene_range) = IRanges::as.vector(seqnames(gene_range))[1]
  coverage_list = lapply(sample_list, readCoverageFromBigWig, 
                         gene_range, flanking_length = flanking_length)
  
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
  #return(coverage_list)
  
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
  coverage_df = plyr::ldply(coverage_list, data.frame, .id = "sample_id") %>% 
    dplyr::mutate(sample_id = as.character(sample_id)) #Convert factor to character
  coverage_df = dplyr::left_join(coverage_df, track_data, by = "sample_id") %>%
    dplyr::mutate(coverage = coverage/scaling_factor) #Normalize by library size
  
  #Calculate mean coverage within each track and colour group
  if(mean_only){
    coverage_df = meanCoverage(coverage_df)
  }
  
  #Make plots
  #Construct transcript structure data.frame from ranges lists
  transcript_struct = prepareTranscriptStructureForPlotting(tx_annotations$exon_ranges, 
                                                            tx_annotations$cds_ranges, transcript_annotations)
  tx_structure = plotTranscriptStructure(transcript_struct, limits)
  
  coverage_plot = plotCoverage(coverage_df, limits, alpha, fill_palette)
  plot = gridExtra::arrangeGrob(coverage_plot, tx_structure, heights = heights, ncol = 1, nrow = 2)
  return(plot)
  }
  