#' Quickly plot transcript structure without read coverage tracks
#'
#' @param exons list of GRanges objects, each object containing exons for one transcript.
#' The list must have names that correspond to transcript_id column in transript_annotations data.frame.
#' @param cdss list of GRanges objects, each object containing the coding regions (CDS) of a single transcript. 
#' The list must have names that correspond to transcript_id column in transript_annotations data.frame. 
#' If cdss is not specified then exons list will be used for both arguments. (default: NULL)
#' @param transcript_annotations Data frame with at least three columns: transcript_id, gene_name, strand.
#' Used to construct transcript labels. (default: NULL)
#' @param rescale_introns Specifies if the introns should be scaled to fixed length or not. (default: TRUE)
#' @param new_intron_length length (bp) of introns after scaling. (default: 50)
#' @param flanking_length Lengths of the flanking regions upstream and downstream of the gene. (default: c(50,50))
#' @param connect_exons Print lines that connect exons together. Set to FALSE when plotting peaks (default: TRUE).
#' @param transcript_label If TRUE then transcript labels are printed above each transcript. (default: TRUE). 
#' @param region_coords Start and end coordinates of the region to plot, overrides flanking_length parameter.
#'
#' @return ggplot2 object
#' @examples
#' plotTranscripts(ncoa7_exons, ncoa7_cdss, ncoa7_metadata, rescale_introns = FALSE)
#' 
#' @export
plotTranscripts <- function(exons, cdss = NULL, transcript_annotations = NULL, 
                            rescale_introns = TRUE, new_intron_length = 50, 
                            flanking_length = c(50,50), connect_exons = TRUE, 
                            transcript_label = TRUE, region_coords = NULL){
  
  #IF cdss is not specified then use exons instead on cdss
  if(is.null(cdss) || length(cdss) == 0){
    cdss = exons
  }
  
  #Check exons and cdss
  assertthat::assert_that(is.list(exons)|| is(exons, "GRangesList")) #Check that exons and cdss objects are lists
  assertthat::assert_that(is.list(cdss) || is(exons, "GRangesList"))
  
  #Join exons together
  joint_exons = joinExons(exons)
  
  #Extract chromosome name
  chromosome_name = as.vector(GenomicRanges::seqnames(joint_exons)[1])
  
  #If region_coords is specificed, then ignore the flanking_length attrbute and compute
  # flanking_length form region_coords
  if(!is.null(region_coords)){
    gene_range = constructGeneRange(joint_exons, c(0,0))
    min_start = min(GenomicRanges::start(gene_range))
    max_end = max(GenomicRanges::end(gene_range))
    flanking_length = c(min_start - region_coords[1], region_coords[2] - max_end)
  }
  #Make sure that flanking_length is a vector of two elements
  assertthat::assert_that(length(flanking_length) == 2) 

  #Rescale introns
  if (rescale_introns){
    tx_annotations = rescaleIntrons(exons, cdss, joint_exons, new_intron_length = new_intron_length, flanking_length)
    xlabel = "Distance from region start (bp)"
  } else {
    old_introns = intronsFromJointExonRanges(GenomicRanges::ranges(joint_exons), flanking_length = flanking_length)
    tx_annotations = list(exon_ranges = lapply(exons, GenomicRanges::ranges), cds_ranges = lapply(cdss, GenomicRanges::ranges),
                          old_introns = old_introns, new_introns = old_introns)
    
    xlabel = paste("Chromosome", chromosome_name, "position (bp)")
  }
  
  #If transcript annotations are not supplied then construct them manually from the GRanges list
  if(is.null(transcript_annotations)){
    plotting_annotations = dplyr::tibble(transcript_id = names(exons),
                                             strand = extractStrandsFromGrangesList(exons)) %>%
      prepareTranscriptAnnotations()
  } else{
    plotting_annotations = prepareTranscriptAnnotations(transcript_annotations)
  }
  
  #Plot transcript structures
  limits = c( min(IRanges::start(tx_annotations$new_introns)), max(IRanges::end(tx_annotations$new_introns)))
  structure = prepareTranscriptStructureForPlotting(tx_annotations$exon_ranges, 
                                               tx_annotations$cds_ranges, plotting_annotations)
  plot = plotTranscriptStructure(structure, limits, connect_exons = connect_exons, xlabel = xlabel, 
                                 transcript_label = transcript_label)
  return(plot)
}

#' Plot read coverage across a genomic region
#' 
#' Also supports rescaling introns to constant length. Extracts read coverage from bigWig files with 
#' extractCoverageData and plots it with plotCoverageData. Custom visualisations can be created by 
#' modifying the plotCoverageData function. Does not work on Windows, because rtracklayer cannot 
#' read BigWig files on Windows.
#' 
#' @param exons list of GRanges objects, each object containing exons for one transcript. 
#' The list must have names that correspond to transcript_id column in transript_annotations data.frame.
#' @param cdss list of GRanges objects, each object containing the coding regions (CDS) of a single transcript. 
#' The list must have names that correspond to transcript_id column in transript_annotations data.frame. 
#' If cdss is not specified then exons list will be used for both arguments. (default: NULL).
#' @param transcript_annotations Data frame with at least three columns: transcript_id, gene_name, strand. 
#' Used to construct transcript labels. (default: NULL)
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
#' @param transcript_label If TRUE then transcript labels are printed above each transcript. (default: TRUE). 
#' @param return_subplots_list Instead of a joint plot return a list of subplots that can be joined together manually. 
#' @param region_coords Start and end coordinates of the region to plot, overrides flanking_length parameter.
#' @param coverage_type Specifies if the read coverage is represented by either 'line', 'area' or 'both'. 
#' The 'both' option tends to give better results for wide regions. (default: area).
#' @param show_legend display legend for the colour_group next to the read coverage plot (default: FALSE).
#'
#' @return Either object from cow_plot::plot_grid() function or a list of subplots (if return_subplots_list == TRUE)
#' @examples
#' require("dplyr")
#' require("GenomicRanges")
#' sample_data = dplyr::data_frame(sample_id = c("aipt_A", "aipt_C", "bima_A", "bima_C"), 
#'     condition = factor(c("Naive", "LPS", "Naive", "LPS"), levels = c("Naive", "LPS")), 
#'     scaling_factor = 1) %>%
#'     dplyr::mutate(bigWig = system.file("extdata",  paste0(sample_id, ".str2.bw"), package = "wiggleplotr"))
#' 
#' track_data = dplyr::mutate(sample_data, track_id = condition, colour_group = condition)
#' 
#' selected_transcripts = c("ENST00000438495", "ENST00000392477") #Plot only two transcripts of the gens
#' \dontrun{
#' plotCoverage(ncoa7_exons[selected_transcripts], ncoa7_cdss[selected_transcripts], 
#'    ncoa7_metadata, track_data, 
#'    heights = c(2,1), fill_palette = getGenotypePalette())
#' }
#' 
#' @export
plotCoverage <- function(exons, cdss = NULL, transcript_annotations = NULL, track_data, rescale_introns = TRUE,
                        new_intron_length = 50, flanking_length = c(50,50),
                        plot_fraction = 0.1, heights = c(0.75, 0.25), alpha = 1,
                        fill_palette = c("#a1dab4","#41b6c4","#225ea8"), mean_only = TRUE, 
                        connect_exons = TRUE, transcript_label = TRUE, return_subplots_list = FALSE,
                        region_coords = NULL, coverage_type = "area", show_legend = FALSE){
  
  #Extract coverage data from bigWig files (and rescale introns)
  coverage_data_list = extractCoverageData(exons, cdss, transcript_annotations, track_data, rescale_introns,
                                           new_intron_length, flanking_length, plot_fraction, mean_only, region_coords)
  
  plot = plotCoverageData(coverage_data_list, heights, alpha,fill_palette,
                          connect_exons, transcript_label, return_subplots_list, coverage_type, show_legend)
  return(plot)
  
}

#' Extract read coverage data from the bigWig files
#' 
#' Does not work on Windows, because rtracklayer cannot read BigWig files on Windows.
#' 
#' @param exons list of GRanges objects, each object containing exons for one transcript. 
#' The list must have names that correspond to transcript_id column in transcript_annotations data.frame.
#' @param cdss list of GRanges objects, each object containing the coding regions (CDS) of a single transcript. 
#' The list must have names that correspond to transcript_id column in trancsript_annotations data.frame. 
#' If cdss is not specified then exons list will be used for both arguments. (default: NULL).
#' @param transcript_annotations Data frame with at least three columns: transcript_id, gene_name, strand. 
#' Used to construct transcript labels. (default: NULL)
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
#' @param rescale_introns Specifies if the introns should be scaled to fixed length or not. (default: TRUE)
#' @param new_intron_length length (bp) of introns after scaling. (default: 50)
#' @param flanking_length Lengths of the flanking regions upstream and downstream of the gene. (default: c(50,50))
#' @param plot_fraction Size of the random sub-sample of points used to plot coverage (between 0 and 1). 
#' Smaller values make plotting significantly faster. (default: 0.1)
#' @param mean_only Plot only mean coverage within each combination of track_id and colour_group values. 
#' Useful for example for plotting mean coverage stratified by genotype (which is specified in the colour_group column) (default: TRUE).
#' @param region_coords Start and end coordinates of the region to plot, overrides flanking_length parameter.
#' The 'both' option tends to give better results for wide regions. (default: area). 
#'
#' @return List containing all of the necessary data for the plotCoverageData function ()
#' @examples
#' require("dplyr")
#' require("GenomicRanges")
#' sample_data = dplyr::data_frame(sample_id = c("aipt_A", "aipt_C", "bima_A", "bima_C"), 
#'     condition = factor(c("Naive", "LPS", "Naive", "LPS"), levels = c("Naive", "LPS")), 
#'     scaling_factor = 1) %>%
#'     dplyr::mutate(bigWig = system.file("extdata",  paste0(sample_id, ".str2.bw"), package = "wiggleplotr"))
#' 
#' track_data = dplyr::mutate(sample_data, track_id = condition, colour_group = condition)
#' 
#' selected_transcripts = c("ENST00000438495", "ENST00000392477") #Plot only two transcripts of the gens
#' \dontrun{
#' extractCoverageData(ncoa7_exons[selected_transcripts], ncoa7_cdss[selected_transcripts], ncoa7_metadata, track_data)
#' }
#' 
#' @export
extractCoverageData <- function(exons, cdss = NULL, transcript_annotations = NULL, track_data, rescale_introns = TRUE,
                                new_intron_length = 50, flanking_length = c(50,50),
                                plot_fraction = 0.1, mean_only = TRUE, region_coords = NULL){
  
  #IF cdss is not specified then use exons instead on cdss
  if(is.null(cdss)){
    cdss = exons
  }
  
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
    track_data = dplyr::mutate_(track_data, .dots = stats::setNames(list(~as.character(bigWig)), c("bigWig")))
  }
  
  #Check transcript annotation
  #If transcript annotations are not supplied then construct them manually from the GRanges list
  if(is.null(transcript_annotations)){
    plotting_annotations = dplyr::data_frame(transcript_id = names(exons),
                                             strand = extractStrandsFromGrangesList(exons)) %>%
      prepareTranscriptAnnotations()
  } else{
    assertthat::assert_that(assertthat::has_name(transcript_annotations, "transcript_id"))
    assertthat::assert_that(assertthat::has_name(transcript_annotations, "gene_name"))
    assertthat::assert_that(assertthat::has_name(transcript_annotations, "strand"))
    plotting_annotations = prepareTranscriptAnnotations(transcript_annotations)
  }
  
  #Check exons and cdss
  assertthat::assert_that(is.list(exons) || is(exons, "GRangesList")) #Check that exons and cdss objects are lists
  assertthat::assert_that(is.list(cdss) || is(exons, "GRangesList"))
  #TODO: Check that the names of the exons and cdss list match that of the transcript_annotations data.frame
  
  #Find the start and end coordinates of the whole region spanning the gene
  joint_exons = joinExons(exons)
  
  #If region_coords is specified, then ignore the flanking_length attribute and compute
  # flanking_length form region_coords
  if(!is.null(region_coords)){
    gene_range = constructGeneRange(joint_exons, c(0,0))
    min_start = min(GenomicRanges::start(gene_range))
    max_end = max(GenomicRanges::end(gene_range))
    flanking_length = c(min_start - region_coords[1], region_coords[2] - max_end)
    
    gene_range = constructGeneRange(joint_exons, flanking_length)
  } else{
    gene_range = constructGeneRange(joint_exons, flanking_length)
  }
  assertthat::assert_that(length(flanking_length) == 2) #flanking_length is a vector of two elements
  
  #Extract chromosome name
  chromosome_name = as.vector(GenomicRanges::seqnames(gene_range)[1])
  
  #Read coverage tracks from BigWig file
  sample_list = as.list(track_data$bigWig)
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
  
  #Convert to data frame and plot
  coverage_df = purrr::map_df(coverage_list, identity, .id = "sample_id") %>% 
    as.data.frame() %>%
    dplyr::mutate_(.dots = stats::setNames(list(~as.character(sample_id)), c("sample_id")) ) #Convert factor to character
  coverage_df = dplyr::left_join(coverage_df, track_data, by = "sample_id") %>%
    dplyr::mutate_(.dots = stats::setNames(list(~coverage/scaling_factor), c("coverage")) ) #Normalize by library size
  
  #Calculate mean coverage within each track and colour group
  if(mean_only){  coverage_df = meanCoverage(coverage_df) }
  
  #Extract plot limits
  limits = c( min(IRanges::start(tx_annotations$new_introns)), max(IRanges::end(tx_annotations$new_introns)))
  
  #Make a result list
  result = list(exons = exons, cdss = cdss,
                tx_annotations = tx_annotations,
                plotting_annotations = plotting_annotations,
                xlabel = xlabel,
                coverage_df = coverage_df,
                limits = limits)
}

#' Plot read coverage across a genomic region
#' 
#' Does not work on Windows, because rtracklayer cannot read BigWig files on Windows.
#' 
#' @param coverage_data_list List of required from the extractCoverageData function:
#' \itemize{
#'  \item exons - list of GRanges objects, each object containing exons for one transcript.
#'  \item cdss - list of GRanges objects, each object containing the coding regions (CDS) of a single transcript. 
#'  \item plotting_annotations - Transcript labels for plotting.
#'  \item tx_annotations - Transcript coordinates for plotting.
#'  \item xlabel - Label of the x-axis.
#'  \item coverage_df - Read coverage data frame.
#'  \item limits - x-axis limits.
#' }
#' @param heights  Specifies the proportion of the height that is dedicated to coverage plots (first value) 
#' relative to transcript annotations (second value). (default: c(0.75,0.25))
#' @param alpha Transparency (alpha) value for the read coverage tracks. 
#' Useful to set to something < 1 when overlaying multiple tracks (see track_id). (default: 1)
#' @param fill_palette Vector of fill colours used for the coverage tracks. Length must be equal to the number of 
#' unique values in track_data$colour_group column.
#' @param connect_exons Print lines that connect exons together. Set to FALSE when plotting peaks (default: TRUE).
#' @param transcript_label If TRUE then transcript labels are printed above each transcript. (default: TRUE). 
#' @param return_subplots_list Instead of a joint plot return a list of subplots that can be joined together manually. 
#' @param coverage_type Specifies if the read coverage is represented by either 'line', 'area' or 'both'. 
#' The 'both' option tends to give better results for wide regions. (default: area). 
#' @param show_legend display legend for the colour_group next to the read coverage plot (default: FALSE).
#'
#' @return Either object from cow_plot::plot_grid() function or a list of subplots (if return_subplots_list == TRUE)
#' @examples
#' require("dplyr")
#' require("GenomicRanges")
#' sample_data = dplyr::data_frame(sample_id = c("aipt_A", "aipt_C", "bima_A", "bima_C"), 
#'     condition = factor(c("Naive", "LPS", "Naive", "LPS"), levels = c("Naive", "LPS")), 
#'     scaling_factor = 1) %>%
#'     dplyr::mutate(bigWig = system.file("extdata",  paste0(sample_id, ".str2.bw"), package = "wiggleplotr"))
#' 
#' track_data = dplyr::mutate(sample_data, track_id = condition, colour_group = condition)
#' 
#' selected_transcripts = c("ENST00000438495", "ENST00000392477") #Plot only two transcripts of the gens
#' \dontrun{
#' cov_data = extractCoverageData(ncoa7_exons[selected_transcripts], ncoa7_cdss[selected_transcripts], ncoa7_metadata, track_data)
#' plotCoverageData(cov_data, heights = c(2,1), fill_palette = getGenotypePalette())
#' }
#' 
#' @export
plotCoverageData <- function(coverage_data_list, heights = c(0.75, 0.25), alpha = 1,
                             fill_palette = c("#a1dab4","#41b6c4","#225ea8"),
                             connect_exons = TRUE, transcript_label = TRUE, 
                             return_subplots_list = FALSE, coverage_type = "area", 
                             show_legend = FALSE){
  
  #Make plots
  #Construct transcript structure data.frame from ranges lists
  transcript_struct = prepareTranscriptStructureForPlotting(coverage_data_list$tx_annotations$exon_ranges, 
                                                            coverage_data_list$tx_annotations$cds_ranges, coverage_data_list$plotting_annotations)
  tx_structure = plotTranscriptStructure(transcript_struct, coverage_data_list$limits, connect_exons, coverage_data_list$xlabel, transcript_label)
  
  coverage_plot = makeCoveragePlot(coverage_data_list$coverage_df, coverage_data_list$limits, alpha, fill_palette, coverage_type, show_legend)
  
  #Choose between returning plot list or a joint plot using plot_grid
  if(return_subplots_list){
    plot_list = list(coverage_plot = coverage_plot, tx_structure = tx_structure)
    return(plot_list)
  } else {
    plot = cowplot::plot_grid(coverage_plot, tx_structure, align = "v", 
                              axis = "lr", rel_heights = heights, ncol = 1)
    return(plot)
  }
}
