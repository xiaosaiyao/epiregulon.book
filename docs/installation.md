# Installation


## Installing from github
All the epiregulon components are available on [github](https://github.com/xiaosaiyao?tab=repositories). 
There are three separate epiregulon packages. The core epiregulon package supports input in the form of `SingleCellExperiment` objects. If the users would like to start from `ArchR` projects, they may choose to use `epiregulon.archr` package, which allows for the seamless integration with `ArchR` package through accepting its output to be used in the downstream workflow. 


``` r
# install devtools
if(!require(devtools)) install.packages("devtools")

# install basic epiregulon package
devtools::install_github(repo='xiaosaiyao/epiregulon')

# install extended version of epiregulon
devtools::install_github(repo='xiaosaiyao/epiregulon.archr')
```
 
Moreover, we provide a suite of tools for the enrichment analysis, visualization, and network analysis which can be run on the `epireglon` or `epiregulon.archr` output.


``` r
# install extended version of epiregulon
devtools::install_github(repo='xiaosaiyao/epiregulon.extra')
```

The data package scMultiome is for storing TF ChIP-seq data and pre-processed datasets.


``` r
# install extended version of epiregulon
devtools::install_github(repo='xiaosaiyao/scMultiome')
```

## Installing from Bioconductor

Epiregulon is now available through [Bioconductor](https://www.bioconductor.org/).


``` r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
 
BiocManager::install("epiregulon")
BiocManager::install("epiregulon.extra")
BiocManager::install("scMultiome")
```

