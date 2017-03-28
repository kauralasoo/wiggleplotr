shortenIntrons <- function(introns, intron_length){
  #Shorten introns from a fixed length to a variable length
  
  #Calculate neccesary parameters
  exons = IRanges::gaps(introns)
  n_introns = length(introns)
  n_exons = length(exons)
  
  #Calculate cumulative with of introns
  intron_cum_width = seq(intron_length,(n_introns-1)*intron_length,intron_length)
  #Calculate new exon starts ignoring introns
  new_intron_starts = c(1,IRanges::start(introns)[2:n_introns] - (IRanges::end(introns)[1:n_introns-1] - intron_cum_width))
  #Add exon widths to the introns
  new_intron_starts = new_intron_starts + c(0,cumsum(IRanges::width(exons)) - IRanges::width(exons))
  
  new_introns = IRanges::IRanges(start = new_intron_starts, width = rep(intron_length, n_introns))
  return(new_introns)
}

shrinkIntronsCoverage <- function(coverage, old_introns, new_introns){
  
  #Covert coverage vector from Rle to normal vector
  coverage = S4Vectors::as.vector.Rle(coverage, mode = "double")

  #Calculate full annotations
  old_annot = sort(c(old_introns, IRanges::gaps(old_introns)))
  new_annot = sort(c(new_introns, IRanges::gaps(new_introns)))
  
  #If new and old annotations are identical then return coverage as data frame
  if(all(IRanges::width(old_annot) == IRanges::width(new_annot))){
    bins = seq(min(IRanges::start(new_annot)), max(IRanges::end(new_annot)))
   
    #Make sure that coverage vector and bins vector have equal length
    assertthat::assert_that(assertthat::are_equal(length(bins), length(coverage)))
    new_coverage = dplyr::data_frame(bins = bins, coverage = coverage)
    return(new_coverage)
    
  } else{ #Otherwise shrink intron converage
    
    #Calculate the width of each annotation bin
    bin_width = ceiling(IRanges::width(old_annot)/IRanges::width(new_annot))
    
    #Build summarisation groups
    s_coord = IRanges::start(new_annot)
    e_coord = IRanges::end(new_annot)
    w_old = IRanges::width(old_annot)
    
    bins = c()
    
    for (i in seq_along(new_annot)){
      bin_id = rep(c(s_coord[i]:e_coord[i]),each = bin_width[i])[1:w_old[i]]
      bins = c(bins, bin_id)
    }
    
    #Calculate mean coverage in bins
    df = data.frame(coverage, bins)
    new_coverage = dplyr::summarize(dplyr::group_by(df, bins), coverage = mean(coverage))
    return(new_coverage)
  }
}

translateExonCoordinates <- function(exons, old_introns, new_introns){
  #Tranlate exon coordinates by shortening introns
  old_exon_starts = IRanges::start(exons)
  old_intron_ends = IRanges::end(old_introns)
  new_intron_ends = IRanges::end(new_introns)
  
  #Translate old exon coordinates to new exon coordinates
  new_exon_starts = rep(0,length(old_exon_starts))
  for (i in seq_along(old_exon_starts)){
    #Find the nearest upstream intron for the current gene
    nearest_intron_number = max(which(old_exon_starts[i] > old_intron_ends))
    new_exon_starts[i] = old_exon_starts[i] - old_intron_ends[nearest_intron_number] + new_intron_ends[nearest_intron_number]
  }
  
  #Create new exon coordinates
  new_exons = IRanges::IRanges(start = new_exon_starts, width = IRanges::width(exons))
  return(new_exons)
}

rescaleIntrons <- function(exons, cdss, joint_exons, new_intron_length, flanking_length){
  
  #Convert exons and cds objects to ranges
  exon_ranges = lapply(exons, GenomicRanges::ranges)
  cds_ranges = lapply(cdss, GenomicRanges::ranges)
  
  #Shorten introns and translate exons into the new exons
  old_introns = intronsFromJointExonRanges(GenomicRanges::ranges(joint_exons), flanking_length = flanking_length)
  new_introns = shortenIntrons(old_introns,new_intron_length)
  new_exon_ranges = lapply(exon_ranges, translateExonCoordinates, old_introns, new_introns)
  new_cds_ranges = lapply(cds_ranges, translateExonCoordinates, old_introns, new_introns)
  
  return(list(exon_ranges = new_exon_ranges, cds_ranges = new_cds_ranges, 
              old_introns = old_introns, new_introns = new_introns))
}