plotTranscriptStructure <- function(exons_df, limits){
  #Create a plot of transcript structure
  plot = ggplot(exons_df, aes(xmin = start, xmax = end,
                              ymax = transcript_rank + 0.3, ymin = transcript_rank - 0.3,
                              x = start, y = transcript_rank, group = transcript_rank)) + 
    geom_line() +
    geom_rect(fill = "white", color = 1) + 
    theme(plot.margin=unit(c(0,1,1,1),"line"), axis.title.y = element_blank()) +
    scale_x_continuous(limits = limits) +
    xlab("Distance from gene start (bp)") +
    facet_grid(type~.)
  return(plot)
}

plotCoverage <- function(coverage_df, limits){
  #Plot coverage over a region
  coverage_plot = ggplot(coverage_df, aes(bins, coverage, group = sample_id, fill = treatment)) + 
    geom_area(alpha = 0.5, position = "identity") + 
    facet_grid(track_id~.) +
    theme(axis.text.x = element_blank(), 
          axis.title.x = element_blank(), 
          axis.ticks.x = element_blank(),
          plot.margin=unit(c(1,1,0,1),"line"),
          axis.title.y = element_blank(),
          legend.position="top") +
    scale_x_continuous(limits = limits)    
}
