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

readCoverageFromBigWig <- function(bigwig_path, joint_exons, flanking = 50){
  #Read coverage over a region from a bigWig file
  chromosome = IRanges::as.vector(seqnames(joint_exons))[1]
  seqlevels(joint_exons) = chromosome
  sel = BigWigSelection(joint_exons)
  coverage_ranges = rtracklayer::import.bw(bigwig_path, selection = sel)
  seqlevels(coverage_ranges) = chromosome
  coverage_rle = coverage(coverage_ranges, weight = score(coverage_ranges))[[1]]
  coverage_rle = coverage_rle[(min(start(joint_exons))-flanking):(max(end(joint_exons))+flanking)] #Keep the region of interest
}

wiggleplotr <- function(exons, sample_data, new_intron_length = 50, plot_fraction = 0.1){
  #Plot read coverage over exons
  transcript_ids = names(exons)
  exon_ranges = lapply(exons, ranges)
  
  #Test that all transcripts are on the same chromosome
  chrs = unlist(lapply(exons, function(x) IRanges::as.vector(seqnames(x)[1])))
  if (!all(chrs == chrs[1])){
    stop("Some transcripts are on different chromosomes.")
  }
  
  #Join all exons together
  joint_exons = c()
  for(tx_id in transcript_ids){
    tx = exons[[tx_id]]
    if(length(joint_exons) == 0){
      joint_exons = tx
    }
    else{
      joint_exons = c(joint_exons, tx)
    }
  }
  joint_exons = reduce(joint_exons)
  joint_ranges = ranges(joint_exons)
  
  #Shorten introns and translate exons into the new exons
  old_introns = gaps(joint_ranges, start = min(start(joint_ranges)) - new_intron_length, end = max(end(joint_ranges)) + new_intron_length)
  new_introns = shortenIntrons(old_introns,new_intron_length)
  new_exon_ranges = lapply(exon_ranges, translateExonCoordinates, old_introns, new_introns)
  
  #Convert exons list to data frame
  exons_df = plyr::ldply(new_exon_ranges, data.frame, .id = "transcript_id")
  exons_df = dplyr::mutate(exons_df, transcript_rank = as.numeric(exons_df$transcript_id), type = "isoforms")
  
  #Read coverage tracks from BigWig file
  sample_list = as.list(sample_data$bigWig)
  names(sample_list) = sample_data$sample_id
  coverage_list = lapply(sample_list, readCoverageFromBigWig, joint_exons, new_intron_length)
  shrunken_coverage_list = lapply(coverage_list, shrinkIntronsCoverage, old_introns, new_introns)
  
  #Take a subsample of points that's easier to plot
  n_total = nrow(shrunken_coverage_list[[1]])
  points = sample(n_total, floor(n_total*plot_fraction))
  points = unique(sort(c(points, start(new_exon_ranges), end(new_exon_ranges), 
                         start(new_exon_ranges) -3, start(new_exon_ranges) +3, 
                         end(new_exon_ranges) + 3, end(new_exon_ranges) -3)))
  points = points[points >= 0]
  shrunken_coverage_list = lapply(shrunken_coverage_list, function(x) {x[points,]} )
  
  #Covert to data frame and plot
  coverage_df = plyr::ldply(shrunken_coverage_list, data.frame, .id = "sample_id")
  coverage_df = plyr::join(coverage_df, sample_data, by = "sample_id")
  coverage_df = dplyr::mutate(coverage_df, coverage = coverage/library_size)

  
  #Make plots
  limits = c(0,n_total)
  tx_structure = plotTranscriptStructure(exons_df, limits)
  coverage_plot = plotCoverage(coverage_df, limits)
  plot = arrangeGrob(coverage_plot, tx_structure, heights = c(0.80, 0.20), ncol = 1, nrow = 2)
  return(plot)
}




