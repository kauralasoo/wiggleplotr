shortenIntrons <- function(introns, intron_length){
  #Shorten introns from a fixed length to a variable length
  
  #Calculate neccesary parameters
  exons = gaps(introns)
  n_introns = length(introns)
  n_exons = length(exons)
  
  #Calculate cumulative with of introns
  intron_cum_width = seq(intron_length,(n_introns-1)*intron_length,intron_length)
  #Calculate new exon starts ignoring introns
  new_intron_starts = c(1,start(introns)[2:n_introns] - (end(introns)[1:n_introns-1] - intron_cum_width))
  #Add exon widths to the introns
  new_intron_starts = new_intron_starts + c(0,cumsum(width(exons)) - width(exons))
  
  new_introns = IRanges(start = new_intron_starts, width = rep(intron_length, n_introns))
  return(new_introns)
}

shrinkIntronsCoverage <- function(coverage, old_introns, new_introns){
  
  #Covert coverage vector from Rle to normal vector
  coverage = IRanges::as.vector(coverage)  
  
  #Calculate full annotations
  old_annot = sort(c(old_introns, gaps(old_introns)))
  new_annot = sort(c(new_introns, gaps(new_introns)))
  
  #Calculate the width of each annotation bin
  bin_width = ceiling(width(old_annot)/width(new_annot))
  
  #Build summarisation groups
  s_coord = start(new_annot)
  e_coord = end(new_annot)
  w_old = width(old_annot)
  
  bins = c()
  for (i in 1:length(new_annot)){
    bin_id = rep(c(s_coord[i]:e_coord[i]),each = bin_width[i])[1:w_old[i]]
    bins = c(bins, bin_id)
  }
  
  #Calculate mean coverage in bins
  df = data.frame(coverage, bins)
  new_coverage = dplyr::summarize(dplyr::group_by(df, bins), coverage = mean(coverage))
  return(new_coverage)
}

translateExonCoordinates <- function(exons, old_introns, new_introns){
  #Tranlate exon coordinates by shortening introns
  old_exon_starts = start(exons)
  old_intron_ends = end(old_introns)
  new_intron_ends = end(new_introns)
  
  #Translate old exon coordinates to new exon coordinates
  new_exon_starts = rep(0,length(old_exon_starts))
  for (i in 1:length(old_exon_starts)){
    #Find the nearest upstream intron for the current gene
    nearest_intron_number = max(which(old_exon_starts[i] > old_intron_ends))
    new_exon_starts[i] = old_exon_starts[i] - old_intron_ends[nearest_intron_number] + new_intron_ends[nearest_intron_number]
  }
  
  #Create new exon coordinates
  new_exons = IRanges(start = new_exon_starts, width = width(exons))
  return(new_exons)
}

rescaleIntrons <- function(exons, cdss, joint_exons, new_intron_length){

  #Join exons together into single GRanges object
  joint_ranges = ranges(joint_exons)
  
  #Convert exons and cds objects to ranges
  exon_ranges = lapply(exons, ranges)
  cds_ranges = lapply(cdss, ranges)
  
  #Shorten introns and translate exons into the new exons
  old_introns = gaps(joint_ranges, 
                start = min(start(joint_ranges)) - new_intron_length, 
                end = max(end(joint_ranges)) + new_intron_length)
  new_introns = shortenIntrons(old_introns,new_intron_length)
  new_exon_ranges = lapply(exon_ranges, translateExonCoordinates, old_introns, new_introns)
  new_cds_ranges = lapply(cds_ranges, translateExonCoordinates, old_introns, new_introns)
  
  return(list(exon_ranges = new_exon_ranges, cds_ranges = new_cds_ranges, 
              old_introns = old_introns, new_introns = new_introns))
}