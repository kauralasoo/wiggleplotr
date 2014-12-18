plotTranscriptStructure <- function(exons_df, limits = NA){
  
  #Extract the position for plotting transcript name
  transcript_annot = dplyr::group_by(exons_df, transcript_id) %>% 
    dplyr::filter(feature_type == "exon") %>%
    dplyr::arrange(transcript_id, start) %>%
    dplyr::filter(row_number() == 1)

  #Create a plot of transcript structure
  plot = ggplot(exons_df) + 
    geom_line(aes(x = start, y = transcript_rank, group = transcript_rank, color = feature_type)) +
    geom_rect(aes(xmin = start, xmax = end, ymax = transcript_rank + 0.2, ymin = transcript_rank - 0.2, fill = feature_type)) + 
    geom_text(aes(x = start, y = transcript_rank + 0.25, label = transcript_label), data = transcript_annot, hjust = 0, vjust = 0, size  =4) +
    theme(plot.margin=unit(c(0,1,1,1),"line"), 
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position="none",
          panel.background = element_blank()) +
    xlab("Distance from gene start (bp)") +
    facet_grid(type~.) +
    scale_y_continuous(expand = c(0,0.5)) +
    scale_fill_manual(values = c("#2c7bb6","#abd9e9")) + 
    scale_colour_manual(values = c("#2c7bb6","#abd9e9"))
  if(all(!is.na(limits))){
    plot = plot + scale_x_continuous(limits = limits, expand = c(0,0))
  }
  return(plot)
}

plotCoverage <- function(coverage_df, limits){
  #Plot coverage over a region
  coverage_plot = ggplot(coverage_df, aes(bins, coverage, group = sample_id, fill = tissue)) + 
    geom_area(alpha = 0.5, position = "identity") + 
    facet_grid(track_id~.) +
    theme(axis.text.x = element_blank(), 
          axis.title.x = element_blank(), 
          axis.ticks.x = element_blank(),
          plot.margin=unit(c(1,1,0,1),"line"),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position="none",
          panel.background = element_blank()) +
    scale_x_continuous(limits = limits, expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    scale_fill_manual(values = c("#1b9e77","#7570b3"))
  return(coverage_plot)
}
