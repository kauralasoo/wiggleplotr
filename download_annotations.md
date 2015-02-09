
Using biomart to download genomic annotations in R
===================================================

This is document provides instructions on how to download Ensembl gene annotations from [biomart](http://www.ensembl.org/biomart/) using R. 

First, we need to load all neccessary R packages.


```r
library("biomaRt")
library("dplyr")
```

## Downloading transcript metadata
First, let's define which mart and dataset we want to use. 

```r
ensembl_mart = useMart("ENSEMBL_MART_ENSEMBL", host = "dec2014.archive.ensembl.org")
ensembl_dataset = useDataset("hsapiens_gene_ensembl",mart=ensembl_mart)
ensembl_dataset
```

```
## Object of class 'Mart':
##  Using the ENSEMBL_MART_ENSEMBL BioMart database
##  Using the hsapiens_gene_ensembl dataset
```
The `host` helps to make sure that we get the annotations from a specific ensembl version. For example, Ensembl 78 correseponds to `host="dec2014.archive.ensembl.org"`. You can use the Ensembl Archives [website](http://www.ensembl.org/info/website/archives/index.html) to check which host name corresponds to desired Ensembl version. More information using specific ensembl versions with biomaRt can be found in the [biomaRt vignette].

We can see all available attributes with the `listAttributes` command. 

```r
attributes = listAttributes(ensembl_dataset)
head(attributes)
```

```
##                    name           description
## 1       ensembl_gene_id       Ensembl Gene ID
## 2 ensembl_transcript_id Ensembl Transcript ID
## 3    ensembl_peptide_id    Ensembl Protein ID
## 4       ensembl_exon_id       Ensembl Exon ID
## 5           description           Description
## 6       chromosome_name       Chromosome Name
```

Now, let's select gene id, gene name, transcript_id and strand from the biomart and download the corresponding columns.

```r
selected_attributes = c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name", "strand")
data = getBM(attributes = selected_attributes, mart = ensembl_dataset)
head(data)
```

```
##   ensembl_transcript_id ensembl_gene_id external_gene_name strand
## 1       ENST00000508957 ENSG00000197468      RP11-747H12.1      1
## 2       ENST00000435337 ENSG00000231049            OR52B5P      1
## 3       ENST00000618935 ENSG00000276385       RP5-859I17.3     -1
## 4       ENST00000614589 ENSG00000275151       RP11-42L13.2      1
## 5       ENST00000432676 ENSG00000228913                UBD     -1
## 6       ENST00000449391 ENSG00000231948         HS1BP3-IT1     -1
```

Finally, we need to rename the columns

```r
data = dplyr::rename(data, transcript_id = ensembl_transcript_id, gene_id = ensembl_gene_id, gene_name = external_gene_name)
head(data)
```

```
##     transcript_id         gene_id     gene_name strand
## 1 ENST00000508957 ENSG00000197468 RP11-747H12.1      1
## 2 ENST00000435337 ENSG00000231049       OR52B5P      1
## 3 ENST00000618935 ENSG00000276385  RP5-859I17.3     -1
## 4 ENST00000614589 ENSG00000275151  RP11-42L13.2      1
## 5 ENST00000432676 ENSG00000228913           UBD     -1
## 6 ENST00000449391 ENSG00000231948    HS1BP3-IT1     -1
```

We can now save the metadata into a file to avoid downloading it every time we need to use it.

```r
saveRDS(data, "transcript_metadata.rds")
```

Next time that we need to access the data we can load it directly from disk.

```r
transcript_metadata = readRDS("transcript_metadata.rds")
head(transcript_metadata)
```

```
##     transcript_id         gene_id     gene_name strand
## 1 ENST00000508957 ENSG00000197468 RP11-747H12.1      1
## 2 ENST00000435337 ENSG00000231049       OR52B5P      1
## 3 ENST00000618935 ENSG00000276385  RP5-859I17.3     -1
## 4 ENST00000614589 ENSG00000275151  RP11-42L13.2      1
## 5 ENST00000432676 ENSG00000228913           UBD     -1
## 6 ENST00000449391 ENSG00000231948    HS1BP3-IT1     -1
```


## Downloading transcript and exon coordinates
We use the [GenomicFeatures] packages to download the transcript and exon coordinates directly from biomaRt.

```r
library("GenomicFeatures")
```




## References
1. [biomaRt vignette]
2. [GenomicFeatures]

[biomaRt vignette]:http://www.bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/biomaRt.pdf
[GenomicFeatures]:http://www.bioconductor.org/packages/release/bioc/html/GenomicFeatures.html
