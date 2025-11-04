convertGRangesToCoverageRle <- function(granges_object, gene_range){
  coverage_rle = GenomicRanges::coverage(granges_object, weight = GenomicRanges::score(granges_object))[[1]]
  
  #Keep coverage in the relevant region
  filter_start = GenomicRanges::start(gene_range)
  filter_end = min(GenomicRanges::end(gene_range), length(coverage_rle))
  coverage_rle = coverage_rle[filter_start:filter_end]
  
  #Convert to from Rle to standard vector and pad end with zeros
  region_width = GenomicRanges::end(gene_range) - GenomicRanges::start(gene_range) + 1
  coverage_vector = S4Vectors::as.vector.Rle(coverage_rle, mode = "double")
  length_diff = region_width - length(coverage_vector)
  
  if(length_diff > 0){
    padding = rep(0, length_diff)
    coverage_vector = c(coverage_vector, padding)
  }
  return(coverage_vector)
}

coverageParquetToGRangesList <- function(coverage_ranges_df, gene_range){
  #Keep ranges from a single region
  coverage_df = dplyr::filter(coverage_ranges_df, 
                              seqnames == as.character(GenomicRanges::seqnames(gene_range)), 
                              end > GenomicRanges::start(gene_range), 
                              start <= GenomicRanges::end(gene_range))
  
  #Group coverage by sample ids and make GRanges
  grouped_df = dplyr::group_by(coverage_df, sample_id)
  sample_ids = dplyr::group_keys(grouped_df)$sample_id
  df_list = dplyr::group_split(grouped_df)
  names(df_list) = sample_ids
  granges_list = purrr::map(df_list, GenomicRanges::GRanges)
  return(granges_list)
}