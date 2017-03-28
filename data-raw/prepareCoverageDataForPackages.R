library("rtracklayer")
library("GenomicFeatures")
library("dplyr")

#Import transcript annotations and metadata
txdb = loadDb("../../annotations/GRCh38/genes/Ensembl_79/TranscriptDb_GRCh38_79.db")
tx_metadata = readRDS("../../annotations/GRCh38/genes/Ensembl_79/Homo_sapiens.GRCh38.79.transcript_data.rds") %>%
  dplyr::rename(transcript_id = ensembl_transcript_id,
                gene_id = ensembl_gene_id,
                gene_name = external_gene_name)
exons = exonsBy(txdb, by = "tx", use.names = TRUE)
cdss = cdsBy(txdb, by = "tx", use.names = TRUE)


#Extract subset of the data needed for the vignette
ncoa7_metadata = dplyr::filter(tx_metadata, gene_name == "NCOA7") %>% tbl_df() %>%
  dplyr::filter(transcript_biotype == "protein_coding") %>%
  dplyr::select(transcript_id, gene_id, gene_name, strand)
ncoa7_exons = exons[transcript_metadata$transcript_id]
ncoa7_cdss = cdss[transcript_metadata$transcript_id]
devtools::use_data(ncoa7_metadata, ncoa7_exons, ncoa7_cdss, overwrite = TRUE)


#Make small versions of some example bigwig files
region = range(exons["ENST00000368357"])
selection = rtracklayer::BigWigSelection(region[[1]])
bw1 = rtracklayer::import.bw("data-raw/bigWig/aipt_A.str2.bw", selection = selection)
bw2 = rtracklayer::import.bw("data-raw/bigWig/aipt_C.str2.bw", selection = selection)
bw3 = rtracklayer::import.bw("data-raw/bigWig/bima_A.str2.bw", selection = selection)
bw4 = rtracklayer::import.bw("data-raw/bigWig/bima_C.str2.bw", selection = selection)

#Export small bigWig files
rtracklayer::export.bw(bw1, "inst/extdata/aipt_A.str2.bw")
rtracklayer::export.bw(bw2, "inst/extdata/aipt_C.str2.bw")
rtracklayer::export.bw(bw3, "inst/extdata/bima_A.str2.bw")
rtracklayer::export.bw(bw4, "inst/extdata/bima_C.str2.bw")

