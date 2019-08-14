
# _wiggleplotr_
_wiggleplotr_ is a tool to visualise RNA-seq read overage overlapping gene annotations. A key feature of _wiggleplotr_ is that it is able rescale all introns of a gene to fixed length, making it easier to see differences in read coverage between neighbouring exons that can otherwise be too far away. Since _wiggleplotr_ takes standard BigWig files as input, it can also be used to visualise read overage from other sequencing-based assays such as ATAC-seq and ChIP-seq. 

<img src="PTK2B.png" width="450">

## Installation
This repostitory contains the development version of _wiggleplotr_. The latest stable version can be installed directly from [Bioconductor](https://bioconductor.org/packages/wiggleplotr/):
```r
## try http:// if https:// URLs are not supported
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("wiggleplotr")
```

Alternatively, you can still install the development version using devtools: 
```r
library("devtools")
devtools::install_github("kauralasoo/wiggleplotr")
```
However, the stable Bioconductor version is likely to be the best option for most people.

## Getting started
See the [vignette](https://htmlpreview.github.io/?https://github.com/kauralasoo/wiggleplotr/blob/master/vignettes/wiggleplotr.html) for instructions on how to get started.

## Citation
If you use wiggleplotr for research, please cite the Bioconductor package directly: [Alasoo K (2019). wiggleplotr: Make read coverage plots from BigWig files. R package version 1.8.0.](https://doi.org/doi:10.18129/B9.bioc.wiggleplotr)
