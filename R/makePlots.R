#' Make a Transcript structure plot 
#' 
#' The transcript structure plots is compatible with wiggpleplotr read coverage and manhattan plots. 
#' Can be appended to those using the cowplot::plot_grid() function. 
#'
#' @param exons_df the dataframe. Result of generateTxStructurePlotData function
#' @param limits The cevtor of limits. Result of generateTxStructurePlotData function (default: NA)
#' @param connect_exons Print lines that connect exons together. Set to FALSE when plotting peaks (default: TRUE).
#' @param xlabel Label to be shown in X-axis (default: "Distance from gene start (bp)"). 
#' @param transcript_label If TRUE then transcript labels are printed above each transcript. (default: TRUE). 
#'
#' @return gglot2 object
#' @examples
#' exon_plot_data <- wiggleplotr::generateTxStructurePlotData(..)
#' exon_plot <- wiggleplotr::plotTranscriptStructure(exons_df = exon_plot_data$transcript_struct_df, limits = exon_plot_data$limits)
#'
#' @export
plotTranscriptStructure <- function(exons_df, limits = NA, connect_exons = TRUE,  
                                    xlabel = "Distance from gene start (bp)", transcript_label = TRUE){
  
  #Extract the position for plotting transcript name
  transcript_annot = dplyr::group_by_(exons_df, ~transcript_id) %>% 
    dplyr::filter_(~feature_type == "exon") %>%
    dplyr::arrange_('transcript_id', 'start') %>%
    dplyr::filter(row_number() == 1)

  #Create a plot of transcript structure
  plot = ggplot(exons_df) + geom_blank()
  if(connect_exons){ #Print line connecting exons
    plot = plot + geom_line(aes_(x = ~start, y = ~transcript_rank, group = ~transcript_rank, color = ~feature_type))
  }
  plot = plot + 
    geom_rect(aes_(xmin = ~start, 
                   xmax = ~end, 
                   ymax = ~transcript_rank + 0.25, 
                   ymin = ~transcript_rank - 0.25, 
                   fill = ~feature_type)) + 
    theme_light() +
    theme(plot.margin=unit(c(0,1,1,1),"line"), 
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position="none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.text.y = element_text(colour = "grey10"),
          strip.background = element_rect(fill = "grey85")) +
    xlab(xlabel) +
    facet_grid(type~.) +
    scale_y_continuous(expand = c(0.2,0.15)) +
    scale_fill_manual(values = c("#2c7bb6","#abd9e9")) + 
    scale_colour_manual(values = c("#2c7bb6","#abd9e9"))
  if(all(!is.na(limits))){
    plot = plot + scale_x_continuous(expand = c(0,0)) +
      coord_cartesian(xlim = limits)
  }
  if(transcript_label){
    plot = plot + geom_text(aes_(x = ~start, 
                                 y = ~transcript_rank + 0.30, 
                                 label = ~transcript_label), 
                            data = transcript_annot, hjust = 0, vjust = 0, size = 4)

  }
  return(plot)
}

#' Make a Coverage plot of the region
#' 
#' The coverage plot is compatible with wiggpleplotr Manhattan and transcript strucutre plots. 
#' Can be appended to those using the cowplot::plot_grid() function. 
#'
#' @param coverage_df the dataframe. Result of generateCoveragePlotData function
#' @param limits The cevtor of limits. Result of generateCoveragePlotData function
#' @param alpha Transparency (alpha) value for the read coverage tracks. 
#' Useful to set to something < 1 when overlaying multiple tracks (see track_id). (default: 1)
#' @param fill_palette Vector of fill colours used for the coverage tracks. Length must be equal to the number of 
#' unique values in track_data$colour_group column.
#' @param coverage_type Specifies if the read coverage is represented by either 'line', 'area' or 'both'. 
#' The 'both' option tends to give better results for wide regions. (default: area). 
#' @param show_genotype_legend If TRUE, shows the legend of genotypes on the right of the plot. (default: false)
#'
#' @return gglot2 object
#' @examples
#' coverage_plot_data = wiggleplotr::generateCoveragePlotData(...)
#'   coverage_plot = wiggleplotr::makeCoveragePlot(coverage_df = coverage_plot_data$coverage_df, 
#'                                                 limits = coverage_plot_data$limits, 
#'                                                 alpha = 1, 
#'                                                 fill_palette = wiggleplotr::getGenotypePalette(), 
#'                                                 coverage_type = "line", 
#'                                                 show_genotype_legend = TRUE)
#' @export
makeCoveragePlot <- function(coverage_df, limits,  alpha = 1, fill_palette = c("#a1dab4","#41b6c4","#225ea8"), coverage_type = "area", show_genotype_legend = FALSE){
  #Plot coverage over a region
  coverage_plot = ggplot(coverage_df, aes_(~bins, ~coverage, group = ~sample_id, alpha = ~alpha)) + 
    geom_blank() +
    theme_light()
  #Choose between plotting a line and plotting area
  if(coverage_type == "line"){
    coverage_plot = coverage_plot + 
      geom_line(aes_(colour = ~colour_group), alpha = alpha, position = "identity") 
  } else if (coverage_type == "area"){
    coverage_plot = coverage_plot + 
      geom_area(aes_(fill = ~colour_group), alpha = alpha, position = "identity")
  } else if (coverage_type == "both"){
    coverage_plot = coverage_plot + 
      geom_area(aes_(fill = ~colour_group), alpha = alpha, position = "identity") +
      geom_line(aes_(colour = ~colour_group), alpha = alpha, position = "identity") 
  } else{
    stop("Coverage type not supported.")
  }
  
  coverage_df <- coverage_df %>% dplyr::filter(!is.na(coverage))
  coverage_max_y = max(coverage_df$coverage) + max(coverage_df$coverage) * 0.05
  coverage_plot = coverage_plot +
    facet_grid(track_id~.) +
    dataTrackTheme() + 
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0), limits = c(0, coverage_max_y)) +
    coord_cartesian(xlim = limits) +
    scale_color_manual(values = fill_palette) +
    scale_fill_manual(values = fill_palette) +
    ylab("FPM")
  
  # Add genotype legend if show_genotype_legend
  if(show_genotype_legend) {
    coverage_plot = coverage_plot + 
      labs(colour = NULL) +
      theme(legend.justification = c(1, 1), 
            legend.position = "right", 
            legend.direction = "vertical", 
            legend.title.align = 0, legend.box = "vertical",
            legend.key = element_rect(colour = "transparent", fill = "white"),
            legend.margin=margin(t=0, r=0, b=0, l=0))  
  }
  
  return(coverage_plot)
}

#' Make a Manahattan plot of p-values
#' 
#' The Manhattan plots is compatible with wiggpleplotr read coverage and transcript strucutre plots. 
#' Can be appended to those using the cowplot::plot_grid() function. 
#'
#' @param pvalues_df Data frame of association p-values (required columns: track_id, p_nominal, pos)
#' @param region_coords Start and end coordinates of the region to plot. 
#' @param color_R2 Color the points according to R2 from the lead variant. Require R2 column in the pvalues_df data frame.
#' @param data_track If TRUE, then remove all information from x-axis. 
#' Makes it easy to append to read coverage or transcript strcture plots using cowplot::plot_grid(). 
#'
#' @return gglot2 object
#' @examples
#' data = dplyr::data_frame(track_id = "GWAS", pos = sample(c(1:1000), 200), p_nominal = runif(200, min = 0.0000001, 1))
#' makeManhattanPlot(data, c(1,1000), data_track = FALSE)
#' @export
makeManhattanPlot <- function(pvalues_df, region_coords, color_R2 = FALSE, data_track = TRUE){
  
  #Make assertions
  assertthat::assert_that(assertthat::has_name(pvalues_df, "track_id"))
  assertthat::assert_that(assertthat::has_name(pvalues_df, "p_nominal"))
  assertthat::assert_that(assertthat::has_name(pvalues_df, "pos"))
  
  #If R2 is specified
  if(color_R2){
    assertthat::assert_that(assertthat::has_name(pvalues_df, "R2"))
    plot_base = ggplot(pvalues_df, aes_(x = ~pos, y = ~-log(p_nominal, 10), colour = ~R2)) + geom_blank()
  } else{
    #Else do not colour
    plot_base = ggplot(pvalues_df, aes_(x = ~pos, y = ~-log(p_nominal, 10))) + geom_blank()
  }
  
  #Make the rest of the plot
  plot = plot_base + 
    facet_grid(track_id ~ .) +
    geom_point() + 
    theme_light() + 
    ylab(expression(paste("-",log[10], " p-value"))) +
    scale_x_continuous(expand = c(0,0)) +
    coord_cartesian(xlim = region_coords)
  
  #Apply data track theme so that plots can later be pasted together with cowplot
  if(data_track){
    plot = plot + dataTrackTheme()
  }
  return(plot)
}
