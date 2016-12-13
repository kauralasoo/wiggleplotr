plotTranscriptStructure <- function(exons_df, limits = NA, connect_exons = TRUE, 
                                    xlabel = "Distance from gene start (bp)"){
  
  #Extract the position for plotting transcript name
  transcript_annot = dplyr::group_by(exons_df, transcript_id) %>% 
    dplyr::filter(feature_type == "exon") %>%
    dplyr::arrange(transcript_id, start) %>%
    dplyr::filter(row_number() == 1)

  #Create a plot of transcript structure
  plot = ggplot(exons_df) + geom_blank()
  if(connect_exons){ #Print line connecting exons
    plot = plot + geom_line(aes(x = start, y = transcript_rank, group = transcript_rank, color = feature_type))
  }
  plot = plot + geom_rect(aes(xmin = start, xmax = end, ymax = transcript_rank + 0.25, ymin = transcript_rank - 0.25, fill = feature_type)) + 
    geom_text(aes(x = start, y = transcript_rank + 0.30, label = transcript_label), data = transcript_annot, hjust = 0, vjust = 0, size  =4) +
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
    plot = plot + scale_x_continuous(limits = limits, expand = c(0,0))
  }
  return(plot)
}

makeCoveragePlot <- function(coverage_df, limits, alpha, fill_palette, coverage_type){
  #Plot coverage over a region
  coverage_plot = ggplot(coverage_df, aes(bins, coverage, group = sample_id, alpha = alpha)) + 
    geom_blank() +
    theme_light()
  #Choose between plotting a line and plotting area
  if(coverage_type == "line"){
    coverage_plot = coverage_plot + 
      geom_line(aes(colour = colour_group), alpha = alpha, position = "identity") 
  } else if (coverage_type == "area"){
    coverage_plot = coverage_plot + 
      geom_area(aes(fill = colour_group), alpha = alpha, position = "identity")
  } else if (coverage_type == "both"){
    coverage_plot = coverage_plot + 
      geom_area(aes(fill = colour_group), alpha = alpha, position = "identity") +
      geom_line(aes(colour = colour_group), alpha = alpha, position = "identity") 
  } else{
    error("Coverage type not supported.")
  }
  coverage_plot = coverage_plot +
    facet_grid(track_id~.) +
    dataTrackTheme() + 
    scale_x_continuous(limits = limits, expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    scale_color_manual(values = fill_palette) +
    scale_fill_manual(values = fill_palette) +
    ylab("FPM")
  return(coverage_plot)
}

makeManhattanPlot <- function(pvalues_df, limits, color_R2 = FALSE, data_track = TRUE){
  
  #Make assertions
  assertthat::assert_that(assertthat::has_name(pvalues_df, "track_id"))
  assertthat::assert_that(assertthat::has_name(pvalues_df, "p_nominal"))
  assertthat::assert_that(assertthat::has_name(pvalues_df, "pos"))
  
  #If R2 is specified
  if(color_R2){
    assertthat::assert_that(assertthat::has_name(pvalues_df, "R2"))
    plot_base = ggplot(pvalues_df, aes(x = pos, y = -log(p_nominal, 10), colour = R2)) + geom_blank()
  } else{
    #Else do not colour
    plot_base = ggplot(pvalues_df, aes(x = pos, y = -log(p_nominal, 10))) + geom_blank()
  }
  
  #Make the rest of the plot
  plot = plot_base + 
    facet_grid(track_id ~ .) +
    geom_point() + 
    theme_light() + 
    ylab(expression(paste("-",log[10], " p-value"))) +
    scale_x_continuous(limits = limits, expand = c(0,0))
  
  #Apply data track theme so that plots can later be pasted together with cowplot
  if(data_track){
    plot = plot + dataTrackTheme()
  }
  return(plot)
}
