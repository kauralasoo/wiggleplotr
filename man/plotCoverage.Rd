% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wiggleplotr.R
\name{plotCoverage}
\alias{plotCoverage}
\title{Plot read coverage across a genomic region}
\usage{
plotCoverage(
  exons,
  cdss = NULL,
  transcript_annotations = NULL,
  track_data,
  rescale_introns = TRUE,
  new_intron_length = 50,
  flanking_length = c(50, 50),
  plot_fraction = 0.1,
  heights = c(0.75, 0.25),
  alpha = 1,
  fill_palette = c("#a1dab4", "#41b6c4", "#225ea8"),
  mean_only = TRUE,
  connect_exons = TRUE,
  transcript_label = TRUE,
  return_subplots_list = FALSE,
  region_coords = NULL,
  coverage_type = "area",
  show_legend = FALSE
)
}
\arguments{
\item{exons}{list of GRanges objects, each object containing exons for one transcript. 
The list must have names that correspond to transcript_id column in transript_annotations data.frame.}

\item{cdss}{list of GRanges objects, each object containing the coding regions (CDS) of a single transcript. 
The list must have names that correspond to transcript_id column in transript_annotations data.frame. 
If cdss is not specified then exons list will be used for both arguments. (default: NULL).}

\item{transcript_annotations}{Data frame with at least three columns: transcript_id, gene_name, strand. 
Used to construct transcript labels. (default: NULL)}

\item{track_data}{data.frame with the metadata for the bigWig read coverage files. Must contain the following columns:
\itemize{
 \item sample_id - unique id for each sample.
 \item track_id - if multiple samples (bigWig files) have the same track_id they will be overlayed on the same 
plot, track_id is also used as the facet label on the right.
 \item bigWig - path to the bigWig file.
 \item scaling_factor - normalisation factor for each sample, useful if different samples sequenced to different 
depth and bigWig files not normalised for that.
 \item colour_group - additional column to group samples into, is used as the colour of the coverage track.
}}

\item{rescale_introns}{Specifies if the introns should be scaled to fixed length or not. (default: TRUE)}

\item{new_intron_length}{length (bp) of introns after scaling. (default: 50)}

\item{flanking_length}{Lengths of the flanking regions upstream and downstream of the gene. (default: c(50,50))}

\item{plot_fraction}{Size of the random sub-sample of points used to plot coverage (between 0 and 1). 
Smaller values make plotting significantly faster. (default: 0.1)}

\item{heights}{Specifies the proportion of the height that is dedicated to coverage plots (first value) 
relative to transcript annotations (second value). (default: c(0.75,0.25))}

\item{alpha}{Transparency (alpha) value for the read coverage tracks. 
Useful to set to something < 1 when overlaying multiple tracks (see track_id). (default: 1)}

\item{fill_palette}{Vector of fill colours used for the coverage tracks. Length must be equal to the number of 
unique values in track_data$colour_group column.}

\item{mean_only}{Plot only mean coverage within each combination of track_id and colour_group values. 
Useful for example for plotting mean coverage stratified by genotype (which is specified in the colour_group column) (default: TRUE).}

\item{connect_exons}{Print lines that connect exons together. Set to FALSE when plotting peaks (default: TRUE).}

\item{transcript_label}{If TRUE then transcript labels are printed above each transcript. (default: TRUE).}

\item{return_subplots_list}{Instead of a joint plot return a list of subplots that can be joined together manually.}

\item{region_coords}{Start and end coordinates of the region to plot, overrides flanking_length parameter.}

\item{coverage_type}{Specifies if the read coverage is represented by either 'line', 'area' or 'both'. 
The 'both' option tends to give better results for wide regions. (default: area).}

\item{show_legend}{display legend for the colour_group next to the read coverage plot (default: FALSE).}
}
\value{
Either object from cow_plot::plot_grid() function or a list of subplots (if return_subplots_list == TRUE)
}
\description{
Also supports rescaling introns to constant length. Extracts read coverage from bigWig files with 
extractCoverageData and plots it with plotCoverageData. Custom visualisations can be created by 
modifying the plotCoverageData function. Does not work on Windows, because rtracklayer cannot 
read BigWig files on Windows.
}
\examples{
require("dplyr")
require("GenomicRanges")
sample_data = dplyr::data_frame(sample_id = c("aipt_A", "aipt_C", "bima_A", "bima_C"), 
    condition = factor(c("Naive", "LPS", "Naive", "LPS"), levels = c("Naive", "LPS")), 
    scaling_factor = 1) \%>\%
    dplyr::mutate(bigWig = system.file("extdata",  paste0(sample_id, ".str2.bw"), package = "wiggleplotr"))

track_data = dplyr::mutate(sample_data, track_id = condition, colour_group = condition)

selected_transcripts = c("ENST00000438495", "ENST00000392477") #Plot only two transcripts of the gens
\dontrun{
plotCoverage(ncoa7_exons[selected_transcripts], ncoa7_cdss[selected_transcripts], 
   ncoa7_metadata, track_data, 
   heights = c(2,1), fill_palette = getGenotypePalette())
}

}
