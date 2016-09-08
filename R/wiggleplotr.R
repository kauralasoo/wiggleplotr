#' Function to quickly plot transcript structure without any read coverage
#'
#' @param exons 
#' @param cdss 
#' @param annotations data.frame with 4 columns: transcript_id, gene_id, gene_name, strand
#' @param rescale_introns 
#' @param new_intron_length 
#' @param flanking_length 
#' @param connect_exons 
#' @param label_type Specifies the format for annotation labels. If set to "transcript" then plots both gene_name and transcript_id, 
#' if set to "peak" then plots only gene_name form transcript_annotations data.frame (default: "transcript"). 
#' @param region_coords 
#'
#' @return ggplot2 object
#' @export
plotTranscripts <- function(exons, cdss, annotations, rescale_introns = TRUE, new_intron_length = 50, 
                            flanking_length = c(50,50), connect_exons = TRUE, label_type = "transcript", 
                            region_coords = NULL){
  #Join exons together
  joint_exons = joinExons(exons)
  
  #Extract chromosome name
  chromosome_name = as.vector(GenomicRanges::seqnames(joint_exons)[1])
  
  #If region_coords is specificed, then ignore the flanking_length attrbute and compute
  # flanking_length form region_coords
  if(!is.null(region_coords)){
    gene_range = constructGeneRange(joint_exons, c(0,0))
    flanking_length = c(GenomicRanges::start(gene_range) - region_coords[1],
                        region_coords[2] - GenomicRanges::end(gene_range))
  }
  
  #Rescale introns
  if (rescale_introns){
    tx_annotations = rescaleIntrons(exons, cdss, joint_exons, new_intron_length = new_intron_length, flanking_length)
    xlabel = "Distance from region start (bp)"
  } else {
    tx_annotations = list(exon_ranges = lapply(exons, GenomicRanges::ranges), cds_ranges = lapply(cdss, GenomicRanges::ranges))
    xlabel = paste("Chromosome", chromosome_name, "position (bp)")
  }
  structure = prepareTranscriptStructureForPlotting(tx_annotations$exon_ranges, 
                                               tx_annotations$cds_ranges, annotations, label_type)
  plot = plotTranscriptStructure(structure, connect_exons = connect_exons, xlabel = xlabel)
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
#' @param transcript_annotations Data.frame with four columns: transcript_id, gene_id, gene_name, strand.
#' @param rescale_introns Specifies if the introns should be scaled to fixed length or not. (default: TRUE)
#' @param new_intron_length length (bp) of introns after scaling. (default: 50)
#' @param flanking_length Lengths of the flanking regions upstream and downstream of the gene. (default: c(50,50))
#' @param plot_fraction Size of the random sub-sample of points used to plot coverage (between 0 and 1). 
#' Smaller values make plotting significantly faster. (default: 0.1)
#' @param heights  Specifies the proportion of the height that is dedicated to coverage plots (first value) 
#' relative to transcript annotations (second value). (default: c(0.75,0.25))
#' @param alpha Transparency (alpha) value for the read coverage tracks. 
#' Useful to set to something < 1 when overlaying multiple tracks (see track_id). (default: 1)
#' @param fill_palette Vector of fill colours used for the coverage tracks. Length must be equal to the number of 
#' unique values in track_data$colour_group column.
#' @param mean_only Plot only mean coverage within each combination of track_id and colour_group values. 
#' Useful for example for plotting mean coverage stratified by genotype (which is specified in the colour_group column) (default: TRUE).
#' @param connect_exons Print lines that connect exons together. Set to FALSE when plotting peaks (default: TRUE).
#' @param label_type Specifies the format for annotation labels. If set to "transcript" then plots both gene_name and transcript_id, 
#' if set to "peak" then plots only gene_name form transcript_annotations data.frame (default: "transcript"). 
#' @param return_subplots_list Instead of a joint plot return a list of subplots that can be joined together manually. 
#' @param region_coords Start and end coordinates of the region, overrides flanking_length parameter. 
#'
#' @return Either object from cow_plot::plot_grid() function or a list of subplots (if return_subplots_list == TRUE)
#' @export
plotCoverage <- function(exons, cdss, track_data, transcript_annotations, rescale_introns = TRUE,
                        new_intron_length = 50, flanking_length = c(50,50),
                        plot_fraction = 0.1, heights = c(0.75, 0.25), alpha = 1,
                        fill_palette = c("#a1dab4","#41b6c4","#225ea8"), mean_only = TRUE, 
                        connect_exons = TRUE, label_type = "transcript", return_subplots_list = FALSE,
                        region_coords = NULL){
  
  #Make some assertions about the input data
  #Check track_data
  assertthat::assert_that(assertthat::has_name(track_data, "sample_id"))
  assertthat::assert_that(assertthat::has_name(track_data, "track_id"))
  assertthat::assert_that(assertthat::has_name(track_data, "bigWig"))
  assertthat::assert_that(assertthat::has_name(track_data, "scaling_factor"))
  assertthat::assert_that(assertthat::has_name(track_data, "colour_group"))
  
  #Make sure that bigWig column is not a factor
  if(is.factor(track_data$bigWig)){
    warning("bigWig column in track_data data.frame is a factor, coverting to a character vector.")
    track_data = dplyr::mutate(track_data, bigWig = as.character(bigWig))
  }
  
  #Check transcript annotation
  assertthat::assert_that(assertthat::has_name(transcript_annotations, "transcript_id"))
  assertthat::assert_that(assertthat::has_name(transcript_annotations, "gene_id"))
  assertthat::assert_that(assertthat::has_name(transcript_annotations, "gene_name"))
  assertthat::assert_that(assertthat::has_name(transcript_annotations, "strand"))
  
  
  #Find the start and end cooridinates of the whole region spanning the gene
  joint_exons = joinExons(exons)
  
  #If region_coords is specificed, then ignore the flanking_length attrbute and compute
  # flanking_length form region_coords
  if(!is.null(region_coords)){
    gene_range = constructGeneRange(joint_exons, c(0,0))
    flanking_length = c(GenomicRanges::start(gene_range) - region_coords[1],
                        region_coords[2] - GenomicRanges::end(gene_range))
    gene_range = constructGeneRange(joint_exons, flanking_length)
  } else{
    gene_range = constructGeneRange(joint_exons, flanking_length)
  }

  #Extract chromosome name
  chromosome_name = as.vector(GenomicRanges::seqnames(gene_range)[1])

  #Read coverage tracks from BigWig file
  sample_list = as.character(track_data$bigWig)
  names(sample_list) = track_data$sample_id
  coverage_list = lapply(sample_list, readCoverageFromBigWig, gene_range)

  #Shorten introns and translate exons into the new introns
  if(rescale_introns){
    #Recale transcript annotations
    tx_annotations = rescaleIntrons(exons, cdss, joint_exons, 
                                    new_intron_length = new_intron_length, flanking_length = flanking_length)
    #Make a label for gene structure plot
    xlabel = "Distance from region start (bp)"
  }
  else{ #Do not rescale transcript annotationn
    #Need to calculate joint intron coordinates for transcript annotations
    old_introns = intronsFromJointExonRanges(GenomicRanges::ranges(joint_exons), flanking_length = flanking_length)
    tx_annotations = list(exon_ranges = lapply(exons, GenomicRanges::ranges), cds_ranges = lapply(cdss, GenomicRanges::ranges),
                          old_introns = old_introns, new_introns = old_introns)
    
    #Make a label for gene structure plot
    xlabel = paste("Chromosome", chromosome_name, "position (bp)")
  }
  #Shrink intron coverage and convert coverage vectors into data frames
  coverage_list = lapply(coverage_list, shrinkIntronsCoverage, tx_annotations$old_introns, tx_annotations$new_introns)

  #Take a subsample of points that is easier to plot
  points = subsamplePoints(tx_annotations, plot_fraction)
  coverage_list = lapply(coverage_list, function(x) {x[points,]} )

  #Covert to data frame and plot
  coverage_df = plyr::ldply(coverage_list, data.frame, .id = "sample_id") %>% 
    dplyr::mutate(sample_id = as.character(sample_id)) #Convert factor to character
  coverage_df = dplyr::left_join(coverage_df, track_data, by = "sample_id") %>%
    dplyr::mutate(coverage = coverage/scaling_factor) #Normalize by library size

  #Calculate mean coverage within each track and colour group
  if(mean_only){  coverage_df = meanCoverage(coverage_df) }
  
  #Make plots
  #Construct transcript structure data.frame from ranges lists
  limits = c( min(IRanges::start(tx_annotations$new_introns)), max(IRanges::end(tx_annotations$new_introns)))
  transcript_struct = prepareTranscriptStructureForPlotting(tx_annotations$exon_ranges, 
                       tx_annotations$cds_ranges, transcript_annotations, label_type)
  tx_structure = plotTranscriptStructure(transcript_struct, limits, connect_exons, xlabel)
  
  coverage_plot = makeCoveragePlot(coverage_df, limits, alpha, fill_palette)
  
  #Choose between returning plot list or a joint plot using plot_grid
  if(return_subplots_list){
    plot_list = list(coverage_plot = coverage_plot, tx_structure = tx_structure)
    return(plot_list)
  } else {
    plot = cowplot::plot_grid(coverage_plot, tx_structure, align = "v", rel_heights = heights, ncol = 1)
    return(plot)
  }
}
  