# Data preparation

Epiregulon operates on the `SingleCellExperiment` class. We assume that gene expression, chromatin accessibility and dimension reduction have been obtained by users' favorite packages prior to the use of Epiregulon. This chapter provides instructions on how to convert gene expression matrix and peak matrix into `SingleCellExperiment` objects from other formats including ArchR project, Seurat objects, AnnData and 10x genomics output. The first section provides a quick primer on the components of a `SingleCellExperiment` object necessary to run Epiregulon. It is thus recommended for all users to go through it. 

## Constructing a SingleCellExperiment object

Let's construct a GeneExpressionMatrix from scratch. First we will create the count matrix. 

``` r
library(SingleCellExperiment)
counts <- matrix(rpois(100000, lambda = 2), ncol=1000, nrow=100)
GeneExpressionMatrix <- SingleCellExperiment(list(counts=counts))
rownames(GeneExpressionMatrix) <- paste("Gene",1:100, sep="-")
colnames(GeneExpressionMatrix) <- paste("Cell",1:1000, sep="-")
```
 

It is important (and efficient) to convert the count matrix to `dgCMatrix` format at the beginning of the workflow.


``` r
library(Matrix)
```

```
## 
## Attaching package: 'Matrix'
```

```
## The following object is masked from 'package:S4Vectors':
## 
##     expand
```

``` r
counts(GeneExpressionMatrix) <- as(counts(GeneExpressionMatrix), "dgCMatrix")
```
  
 
Next we will add the cell information to `colData`

``` r
colData(GeneExpressionMatrix) <- DataFrame(Cluster = paste("cluster", sample(1:3,1000, TRUE)))
```


For the purpose of Epiregulon, it is important to provide the start and end position of the genes so that we can link genes to the peak regions.


``` r
seqnames <- paste0("chr", sample(1:2, 100, TRUE))
start <- sample(1:100000, 100, TRUE)
end <- start + sample(100:500, 100, TRUE)
strand <- sample(c("+", "-", "*"), 100, TRUE)

gr <- GRanges(
    seqnames = seqnames,
    ranges = IRanges(start = start, end = end),
    strand = strand
)
```

Next we provide additional information about the genes. We can add them into the `mcols` as `DataFrame`. The `GRanges` become the `rowRanges` of the GeneExpressionMatrix.


``` r
mcols(gr) <- DataFrame(name = paste("Gene", 1:100, sep="-"),
                       ID = paste0("ID", 1:100))

rowRanges(GeneExpressionMatrix) <- gr
rowRanges(GeneExpressionMatrix)
```

```
## GRanges object with 100 ranges and 2 metadata columns:
##         seqnames       ranges strand |        name          ID
##            <Rle>    <IRanges>  <Rle> | <character> <character>
##     [1]     chr1 99727-100177      * |      Gene-1         ID1
##     [2]     chr2    9362-9700      + |      Gene-2         ID2
##     [3]     chr2  15287-15674      + |      Gene-3         ID3
##     [4]     chr2    7816-8139      * |      Gene-4         ID4
##     [5]     chr1  54227-54387      * |      Gene-5         ID5
##     ...      ...          ...    ... .         ...         ...
##    [96]     chr1  36777-37191      + |     Gene-96        ID96
##    [97]     chr1  63321-63588      - |     Gene-97        ID97
##    [98]     chr2  26599-27032      + |     Gene-98        ID98
##    [99]     chr1  82008-82259      * |     Gene-99        ID99
##   [100]     chr2  11651-12071      * |    Gene-100       ID100
##   -------
##   seqinfo: 2 sequences from an unspecified genome; no seqlengths
```

The additional information about the genes will also appear in `rowData`

``` r
rowData(GeneExpressionMatrix)
```

```
## DataFrame with 100 rows and 2 columns
##            name          ID
##     <character> <character>
## 1        Gene-1         ID1
## 2        Gene-2         ID2
## 3        Gene-3         ID3
## 4        Gene-4         ID4
## 5        Gene-5         ID5
## ...         ...         ...
## 96      Gene-96        ID96
## 97      Gene-97        ID97
## 98      Gene-98        ID98
## 99      Gene-99        ID99
## 100    Gene-100       ID100
```

Finally we add some reduced dimension data which is needed for clustering

``` r
reducedDim(GeneExpressionMatrix, "PCA") <- matrix(data=rnorm(20000), nrow=1000, ncol=20)
GeneExpressionMatrix
```

```
## class: SingleCellExperiment 
## dim: 100 1000 
## metadata(0):
## assays(1): counts
## rownames: NULL
## rowData names(2): name ID
## colnames: NULL
## colData names(1): Cluster
## reducedDimNames(1): PCA
## mainExpName: NULL
## altExpNames(0):
```

Repeat the process to create a `PeakMatrix`


``` r
# add counts
counts <- matrix(rpois(1000000, lambda = 1), ncol=1000, nrow=1000)
PeakMatrix <- SingleCellExperiment(list(counts=counts))
rownames(PeakMatrix) <- paste("Peak",1:1000, sep="-")
colnames(PeakMatrix) <- paste("Cell",1:1000, sep="-")

# convert count matrix to dgCMatrix
counts(PeakMatrix) <- as(counts(PeakMatrix), "dgCMatrix")

# add rowRanges
seqnames <- paste0("chr", sample(1:2, 1000, TRUE))
start <- sample(1:100000, 1000, TRUE)
end <- start + sample(100:500, 1000, TRUE)
strand <- sample(c("+", "-", "*"), 1000, TRUE)

gr <- GRanges(
    seqnames = seqnames,
    ranges = IRanges(start = start, end = end),
    strand = strand
)

# add rowData
mcols(gr) <- DataFrame(name = paste("Peak", 1:1000, sep="-"))

rowRanges(PeakMatrix) <- gr

# add reduced dimensionality matrix

reducedDim(PeakMatrix, "ATAC_LSI") <- matrix(data=rnorm(20000), nrow=1000, ncol=20)
PeakMatrix
```

```
## class: SingleCellExperiment 
## dim: 1000 1000 
## metadata(0):
## assays(1): counts
## rownames: NULL
## rowData names(1): name
## colnames(1000): Cell-1 Cell-2 ... Cell-999 Cell-1000
## colData names(0):
## reducedDimNames(1): ATAC_LSI
## mainExpName: NULL
## altExpNames(0):
```

For more information on SingleCellExperiment, please refer to the [documentation]("https://bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html") on bioconductor.

For a real example of a properly formatted `SingleCellExperiment` object, check out datasets from the [scMultiome](https://www.bioconductor.org/packages/release/data/experiment/html/scMultiome.html) package. Both the GeneExpressionMatrix and PeakMatrix are combined into a `MultiAssayExperiment` object.

``` r
mae <- scMultiome::reprogramSeq()
```

```
## see ?scMultiome and browseVignettes('scMultiome') for documentation
```

```
## loading from cache
```

``` r
GeneExpressMatrix <- mae[["GeneExpressMatrix "]]
GeneExpressMatrix
```

```
## NULL
```

``` r
PeakMatrix <- mae[["PeakMatrix"]]
PeakMatrix
```

```
## class: SingleCellExperiment 
## dim: 126602 3903 
## metadata(1): .internal
## assays(1): counts
## rownames: NULL
## rowData names(1): idx
## colnames(3903): reprogram#TTAGGAACAAGGTACG-1
##   reprogram#GAGCGGTCAACCTGGT-1 ... reprogram#GGTTACTAGACACCGC-1
##   reprogram#CGCTATGAGTGAACAG-1
## colData names(23): BlacklistRatio DoubletEnrichment ... ReadsInPeaks
##   FRIP
## reducedDimNames(0):
## mainExpName: NULL
## altExpNames(0):
```
## Starting from an ArchR object

Epiregulon is designed to work seamlessly with ArchR. GeneExpressionMatrix and PeakMatrix can be easily exported from ArchR using ArchR's build-in function.

Download a test ArchR project

``` r
library(ArchR)
archr.proj <- getTestProject()
```

Check available matrices in ArchR project

``` r
getAvailableMatrices(ArchRProj = archr.proj)
```

```
## [1] "GeneIntegrationMatrix" "GeneScoreMatrix"       "MotifMatrix"          
## [4] "PeakMatrix"            "TileMatrix"
```

Export the GeneExpressionMatrix and PeakMatrix from the ArchR project.


``` r
GeneExpressionMatrix <- getMatrixFromProject(ArchRProj = archr.proj, useMatrix = "GeneIntegrationMatrix")
PeakMatrix <- getMatrixFromProject(ArchRProj = archr.proj, useMatrix = "PeakMatrix")
GeneExpressionMatrix
```

```
## class: SummarizedExperiment 
## dim: 2051 127 
## metadata(0):
## assays(1): GeneIntegrationMatrix
## rownames: NULL
## rowData names(6): seqnames start ... name idx
## colnames(127): PBSmall#B.43 PBSmall#T.24 ... PBSmall#B.37 PBSmall#B.4
## colData names(20): BlacklistRatio nDiFrags ... predictedGroup_Un
##   predictedScore_Un
```

``` r
PeakMatrix
```

```
## class: RangedSummarizedExperiment 
## dim: 2142 127 
## metadata(0):
## assays(1): PeakMatrix
## rownames: NULL
## rowData names(1): idx
## colnames(127): PBSmall#B.43 PBSmall#T.24 ... PBSmall#B.37 PBSmall#B.4
## colData names(20): BlacklistRatio nDiFrags ... predictedGroup_Un
##   predictedScore_Un
```

The GeneExpressionMatrix and PeakMatrix exported from ArchR project are in the form of `SummarizedExperiment` and `RangedSummarizedExperiment` respectively. We provide a helper function to convert both to `SingleCellExperiment` class. Furthermore, the genomic location of the genes are transferred from `rowData` in the `RangedSummarizedExperiment` to the `rowRanges` of the `SingleCellExperiment`.

Please note that GeneExpressionMatrix extracted from ArchR project contain normalized counts (not logged) and for clarity, we rename the assay as "normalizedCounts"

``` r
library(epiregulon.archr)
```

```
## Loading required package: epiregulon
```

```
## 
## Attaching package: 'epiregulon.archr'
```

```
## The following objects are masked from 'package:epiregulon':
## 
##     addMotifScore, addTFMotifInfo, calculateP2G, getTFMotifInfo
```

``` r
GeneExpressionMatrix <- ArchRMatrix2SCE(rse = GeneExpressionMatrix, rename = "normalizedCounts")
GeneExpressionMatrix
```

```
## class: SingleCellExperiment 
## dim: 2051 127 
## metadata(0):
## assays(1): normalizedCounts
## rownames: NULL
## rowData names(2): name idx
## colnames(127): PBSmall#B.43 PBSmall#T.24 ... PBSmall#B.37 PBSmall#B.4
## colData names(20): BlacklistRatio nDiFrags ... predictedGroup_Un
##   predictedScore_Un
## reducedDimNames(0):
## mainExpName: NULL
## altExpNames(0):
```

We can also transform the counts to logcounts. This will add a new assay and we name it as "logcounts".

``` r
GeneExpressionMatrix <- ArchRMatrix2SCE(rse = GeneExpressionMatrix, transform = TRUE, transform_method = "log", log_name = "logcounts")
GeneExpressionMatrix
```

```
## class: SingleCellExperiment 
## dim: 2051 127 
## metadata(0):
## assays(2): counts logcounts
## rownames: NULL
## rowData names(2): name idx
## colnames(127): PBSmall#B.43 PBSmall#T.24 ... PBSmall#B.37 PBSmall#B.4
## colData names(20): BlacklistRatio nDiFrags ... predictedGroup_Un
##   predictedScore_Un
## reducedDimNames(0):
## mainExpName: NULL
## altExpNames(0):
```

Check that rowRanges have been transferred 

``` r
rowRanges(GeneExpressionMatrix)
```

```
## GRanges object with 2051 ranges and 2 metadata columns:
##          seqnames              ranges strand |        name       idx
##             <Rle>           <IRanges>  <Rle> | <character> <integer>
##      [1]    chr11       126987-139152      - |   LINC01001         1
##      [2]    chr11       193080-194573      + |     SCGB1C1         2
##      [3]    chr11       196761-200258      + |        ODF3         3
##      [4]    chr11       202924-207422      - |       BET1L         4
##      [5]    chr11       208530-215110      + |       RIC8A         5
##      ...      ...                 ...    ... .         ...       ...
##   [2047]     chr5 180551357-180552304      - |       OR2V1       827
##   [2048]     chr5 180581943-180582890      + |       OR2V2       828
##   [2049]     chr5 180620924-180632177      - |       TRIM7       829
##   [2050]     chr5 180650263-180662808      + |      TRIM41       830
##   [2051]     chr5 180683386-180688119      - |      TRIM52       831
##   -------
##   seqinfo: 2 sequences from an unspecified genome; no seqlengths
```

We next convert the PeakMatrix to  a `SingleCellExperiment` object. The counts exported from ArchR are raw counts and thus we rename the assay as "counts"

``` r
PeakMatrix <- ArchRMatrix2SCE(PeakMatrix, rename = "counts")
PeakMatrix
```

```
## class: SingleCellExperiment 
## dim: 2142 127 
## metadata(0):
## assays(1): counts
## rownames: NULL
## rowData names(1): idx
## colnames(127): PBSmall#B.43 PBSmall#T.24 ... PBSmall#B.37 PBSmall#B.4
## colData names(20): BlacklistRatio nDiFrags ... predictedGroup_Un
##   predictedScore_Un
## reducedDimNames(0):
## mainExpName: NULL
## altExpNames(0):
```

Neither SummarizedExperiment object or RangedSummarizedExperiment object contains reduced dimension, so we must extract the reduced dimensionality matrix from the ArchR project and add it to GeneExpressionMatrix and/or PeakMatrix


``` r
reducedDim(PeakMatrix, "IterativeLSI") <- getReducedDims(ArchRProj = archr.proj, reducedDims = "IterativeLSI")
```

Refer to ArchR [manual](https://www.archrproject.com/bookdown/index.html) for full documentation.

## Starting from a Seurat/Signac object

We download an example multimodel dataset from [SeuratData](https://github.com/satijalab/seurat-data). 

``` r
library(SeuratData)

InstallData("pbmcMultiome")
```

This is a PBMC dataset consisting of both RNAseq and ATACseq data and we load each modality as a separate Seurat object

``` r
library(Seurat)
```

```
## Loading required package: SeuratObject
```

```
## Loading required package: sp
```

```
## 
## Attaching package: 'sp'
```

```
## The following object is masked from 'package:IRanges':
## 
##     %over%
```

```
## 
## Attaching package: 'SeuratObject'
```

```
## The following object is masked from 'package:SummarizedExperiment':
## 
##     Assays
```

```
## The following object is masked from 'package:GenomicRanges':
## 
##     intersect
```

```
## The following object is masked from 'package:GenomeInfoDb':
## 
##     intersect
```

```
## The following object is masked from 'package:IRanges':
## 
##     intersect
```

```
## The following object is masked from 'package:S4Vectors':
## 
##     intersect
```

```
## The following object is masked from 'package:BiocGenerics':
## 
##     intersect
```

```
## The following objects are masked from 'package:base':
## 
##     intersect, t
```

```
## 
## Attaching package: 'Seurat'
```

```
## The following object is masked from 'package:SummarizedExperiment':
## 
##     Assays
```

``` r
library(Signac)
library(SeuratData)
```

```
## ── Installed datasets ──────────────────────────────── SeuratData v0.2.2.9001 ──
```

```
## ✔ pbmcMultiome 0.1.4
```

```
## ────────────────────────────────────── Key ─────────────────────────────────────
```

```
## ✔ Dataset loaded successfully
## ❯ Dataset built with a newer version of Seurat than installed
## ❓ Unknown version of Seurat installed
```

``` r
pbmc.rna <- LoadData("pbmcMultiome", "pbmc.rna")
```

```
## Validating object structure
```

```
## Updating object slots
```

```
## Ensuring keys are in the proper structure
## Ensuring keys are in the proper structure
```

```
## Ensuring feature names don't have underscores or pipes
```

```
## Updating slots in RNA
```

```
## Validating object structure for Assay 'RNA'
```

```
## Object representation is consistent with the most current Seurat version
```

```
## Warning: Assay RNA changing from Assay to Assay5
```

``` r
pbmc.atac <- LoadData("pbmcMultiome", "pbmc.atac")
```

```
## Validating object structure
```

```
## Updating object slots
```

```
## Ensuring keys are in the proper structure
## Ensuring keys are in the proper structure
```

```
## Ensuring feature names don't have underscores or pipes
```

```
## Updating slots in ATAC
```

```
## Validating object structure for ChromatinAssay 'ATAC'
```

```
## Object representation is consistent with the most current Seurat version
```
We first convert the RNA Seurat object to GeneExpressionMatrix `SingleCellExperiment` class. 

``` r
GeneExpressionMatrix <- as.SingleCellExperiment(pbmc.rna, assay="RNA")
```

```
## Warning: Layer 'scale.data' is empty
```

Because these genes are missing genomic positions, we must first annotate them with a genomic database, for example ensembl

``` r
library(AnnotationHub)
ah <- AnnotationHub()
edb <- ah[["AH98047"]]  #"EnsDb.Hsapiens.v105"
```

```
## loading from cache
```

```
## require("ensembldb")
```

``` r
gr <- genes(edb, columns = c("gene_id", "gene_name"))
```

We retain only the genes that have genomic positions. The genomic positions are necessary to link peak positions to the nearby target genes.


``` r
common_genes <- na.omit(intersect(gr$gene_name, rownames(GeneExpressionMatrix)))
GeneExpressionMatrix <- GeneExpressionMatrix[common_genes,]
rowRanges(GeneExpressionMatrix) <- gr[match(common_genes, gr$gene_name)]
rowRanges(GeneExpressionMatrix)
```

```
## GRanges object with 23644 ranges and 2 metadata columns:
##                   seqnames            ranges strand |         gene_id
##                      <Rle>         <IRanges>  <Rle> |     <character>
##   ENSG00000243485        1       29554-31109      + | ENSG00000243485
##   ENSG00000237613        1       34554-36081      - | ENSG00000237613
##   ENSG00000186092        1       65419-71585      + | ENSG00000186092
##   ENSG00000284733        1     450740-451678      - | ENSG00000284733
##   ENSG00000284662        1     685716-686654      - | ENSG00000284662
##               ...      ...               ...    ... .             ...
##   ENSG00000228296        Y 25063083-25099892      - | ENSG00000228296
##   ENSG00000223641        Y 25182277-25213389      - | ENSG00000223641
##   ENSG00000228786        Y 25378300-25394719      - | ENSG00000228786
##   ENSG00000172288        Y 25622162-25624902      + | ENSG00000172288
##   ENSG00000231141        Y 25728490-25733388      + | ENSG00000231141
##                      gene_name
##                    <character>
##   ENSG00000243485  MIR1302-2HG
##   ENSG00000237613      FAM138A
##   ENSG00000186092        OR4F5
##   ENSG00000284733       OR4F29
##   ENSG00000284662       OR4F16
##               ...          ...
##   ENSG00000228296       TTTY4C
##   ENSG00000223641      TTTY17C
##   ENSG00000228786 LINC00266-4P
##   ENSG00000172288         CDY1
##   ENSG00000231141        TTTY3
##   -------
##   seqinfo: 456 sequences (1 circular) from GRCh38 genome
```

We then convert ATAC matrix to PeakMatrix. After conversion to `SingleCellExperiment`, the peak positions appear as rownames and must be converted  to `GRanges`.


``` r
PeakMatrix <- as.SingleCellExperiment(pbmc.atac, assay="ATAC")

peak_position <- strsplit(rownames(PeakMatrix), split = "-")
gr <- GRanges(
    seqnames = sapply(peak_position,"[",1),
    ranges = IRanges(start = as.numeric(sapply(peak_position,"[",2)), 
                     end = as.numeric(sapply(peak_position,"[",3)))
)
rowRanges(PeakMatrix) <- gr
PeakMatrix
```

```
## class: SingleCellExperiment 
## dim: 108377 11909 
## metadata(0):
## assays(2): counts logcounts
## rownames: NULL
## rowData names(0):
## colnames(11909): AAACAGCCAAGGAATC-1 AAACAGCCAATCCCTT-1 ...
##   TTTGTTGGTTGGTTAG-1 TTTGTTGGTTTGCAGA-1
## colData names(5): orig.ident nCount_ATAC nFeature_ATAC
##   seurat_annotations ident
## reducedDimNames(0):
## mainExpName: ATAC
## altExpNames(0):
```

For more information, refer to Signac [tutorial](https://stuartlab.org/signac/articles/overview)


## Starting from AnnData

We download an example PBMC anndata dataset from [scglue](https://scglue.readthedocs.io/en/latest/data.html).

We first import the GeneExpressionMatrix.

``` r
library(zellkonverter)
```

```
## Registered S3 method overwritten by 'zellkonverter':
##   method                                             from      
##   py_to_r.pandas.core.arrays.categorical.Categorical reticulate
```

``` r
library(GenomicRanges)
library(SingleCellExperiment)
url <- "http://download.gao-lab.org/GLUE/dataset/10x-Multiome-Pbmc10k-RNA.h5ad"
destfile <- tempfile(fileext = ".h5ad")

# Download the file
download.file(url, destfile, mode = "wb")
GeneExpressionMatrix <- readH5AD(destfile)
```

The field "chrom" "chromStart"  "chromEnd" correspond to genomic positions of the genes.

``` r
rowData(GeneExpressionMatrix)
```

```
## DataFrame with 29095 rows and 26 columns
##                   gene_ids   feature_types   genome      chrom chromStart
##                <character>        <factor> <factor>   <factor>  <numeric>
## AL627309.1 ENSG00000238009 Gene Expression   GRCh38       chr1      89294
## AL627309.5 ENSG00000241860 Gene Expression   GRCh38       chr1     141473
## AL627309.4 ENSG00000241599 Gene Expression   GRCh38       chr1     160445
## AP006222.2 ENSG00000286448 Gene Expression   GRCh38       chr1     266854
## AL669831.2 ENSG00000229905 Gene Expression   GRCh38       chr1     760910
## ...                    ...             ...      ...        ...        ...
## AC004556.3 ENSG00000276345 Gene Expression   GRCh38 KI270721.1       2584
## AC233755.2 ENSG00000277856 Gene Expression   GRCh38 KI270726.1      26240
## AC233755.1 ENSG00000275063 Gene Expression   GRCh38 KI270726.1      41443
## AC007325.1 ENSG00000276017 Gene Expression   GRCh38 KI270734.1      72410
## AC007325.4 ENSG00000278817 Gene Expression   GRCh38 KI270734.1     131493
##             chromEnd            name    score   strand thickStart thickEnd
##            <numeric>     <character> <factor> <factor>   <factor> <factor>
## AL627309.1    133723 ENSG00000238009        .        -          .        .
## AL627309.5    173862 ENSG00000241860        .        -          .        .
## AL627309.4    161525 ENSG00000241599        .        +          .        .
## AP006222.2    268655 ENSG00000286448        .        +          .        .
## AL669831.2    761989 ENSG00000229905        .        +          .        .
## ...              ...             ...      ...      ...        ...      ...
## AC004556.3     11802 ENSG00000276345        .        +          .        .
## AC233755.2     26534 ENSG00000277856        .        +          .        .
## AC233755.1     41876 ENSG00000275063        .        +          .        .
## AC007325.1     74814 ENSG00000276017        .        +          .        .
## AC007325.4    137392 ENSG00000278817        .        +          .        .
##             itemRgb blockCount blockSizes blockStarts      gene_type  gene_name
##            <factor>   <factor>   <factor>    <factor>       <factor>   <factor>
## AL627309.1        .          .          .           .         lncRNA AL627309.1
## AL627309.5        .          .          .           .         lncRNA AL627309.5
## AL627309.4        .          .          .           .         lncRNA AL627309.4
## AP006222.2        .          .          .           .         lncRNA AP006222.2
## AL669831.2        .          .          .           .         lncRNA AL669831.2
## ...             ...        ...        ...         ...            ...        ...
## AC004556.3        .          .          .           . protein_coding AC004556.3
## AC233755.2        .          .          .           . protein_coding AC233755.2
## AC233755.1        .          .          .           . protein_coding AC233755.1
## AC007325.1        .          .          .           . protein_coding AC007325.1
## AC007325.4        .          .          .           . protein_coding AC007325.4
##             hgnc_id          havana_gene               tag  n_counts
##            <factor>             <factor>          <factor> <numeric>
## AL627309.1       NA OTTHUMG00000001096.2 overlapping_locus        70
## AL627309.5       NA OTTHUMG00000002480.4 ncRNA_host              442
## AL627309.4       NA OTTHUMG00000002525.1 NA                       44
## AP006222.2       NA OTTHUMG00000194680.1 NA                        1
## AL669831.2       NA OTTHUMG00000002408.1 NA                       10
## ...             ...                  ...               ...       ...
## AC004556.3       NA                   NA                NA       320
## AC233755.2       NA                   NA                NA         1
## AC233755.1       NA                   NA                NA         1
## AC007325.1       NA                   NA                NA         3
## AC007325.4       NA                   NA                NA        43
##            highly_variable highly_variable_rank       means   variances
##                  <logical>            <numeric>   <numeric>   <numeric>
## AL627309.1           FALSE                  NaN 0.007268196 0.008462225
## AL627309.5           FALSE                  NaN 0.045893469 0.050853072
## AL627309.4           FALSE                  NaN 0.004568581 0.004755865
## AP006222.2           FALSE                  NaN 0.000103831 0.000103831
## AL669831.2           FALSE                  NaN 0.001038314 0.001037343
## ...                    ...                  ...         ...         ...
## AC004556.3           FALSE                  NaN 0.033226041 0.035448356
## AC233755.2           FALSE                  NaN 0.000103831 0.000103831
## AC233755.1           FALSE                  NaN 0.000103831 0.000103831
## AC007325.1           FALSE                  NaN 0.000311494 0.000311429
## AC007325.4           FALSE                  NaN 0.004464749 0.004445277
##            variances_norm
##                 <numeric>
## AL627309.1       0.971895
## AL627309.5       0.888672
## AL627309.4       0.891707
## AP006222.2       0.999904
## AL669831.2       0.933916
## ...                   ...
## AC004556.3       0.856870
## AC233755.2       0.999904
## AC233755.1       0.999904
## AC007325.1       0.975766
## AC007325.4       0.854162
```

They must be renamed to "seqnames" "start" and "end" before conversion to `GRanges`

``` r
index_to_rename <- match(c("chrom", "chromStart",  "chromEnd"), colnames(rowData(GeneExpressionMatrix)))
colnames(rowData(GeneExpressionMatrix))[index_to_rename] <- c("seqnames", "start",   "end")
rowRanges(GeneExpressionMatrix) <- GRanges(rowData(GeneExpressionMatrix))
rowRanges(GeneExpressionMatrix)
```

```
## GRanges object with 29095 ranges and 22 metadata columns:
##                seqnames        ranges strand |        gene_ids   feature_types
##                   <Rle>     <IRanges>  <Rle> |     <character>        <factor>
##   AL627309.1       chr1  89294-133723      - | ENSG00000238009 Gene Expression
##   AL627309.5       chr1 141473-173862      - | ENSG00000241860 Gene Expression
##   AL627309.4       chr1 160445-161525      + | ENSG00000241599 Gene Expression
##   AP006222.2       chr1 266854-268655      + | ENSG00000286448 Gene Expression
##   AL669831.2       chr1 760910-761989      + | ENSG00000229905 Gene Expression
##          ...        ...           ...    ... .             ...             ...
##   AC004556.3 KI270721.1    2584-11802      + | ENSG00000276345 Gene Expression
##   AC233755.2 KI270726.1   26240-26534      + | ENSG00000277856 Gene Expression
##   AC233755.1 KI270726.1   41443-41876      + | ENSG00000275063 Gene Expression
##   AC007325.1 KI270734.1   72410-74814      + | ENSG00000276017 Gene Expression
##   AC007325.4 KI270734.1 131493-137392      + | ENSG00000278817 Gene Expression
##                genome            name    score thickStart thickEnd  itemRgb
##              <factor>     <character> <factor>   <factor> <factor> <factor>
##   AL627309.1   GRCh38 ENSG00000238009        .          .        .        .
##   AL627309.5   GRCh38 ENSG00000241860        .          .        .        .
##   AL627309.4   GRCh38 ENSG00000241599        .          .        .        .
##   AP006222.2   GRCh38 ENSG00000286448        .          .        .        .
##   AL669831.2   GRCh38 ENSG00000229905        .          .        .        .
##          ...      ...             ...      ...        ...      ...      ...
##   AC004556.3   GRCh38 ENSG00000276345        .          .        .        .
##   AC233755.2   GRCh38 ENSG00000277856        .          .        .        .
##   AC233755.1   GRCh38 ENSG00000275063        .          .        .        .
##   AC007325.1   GRCh38 ENSG00000276017        .          .        .        .
##   AC007325.4   GRCh38 ENSG00000278817        .          .        .        .
##              blockCount blockSizes blockStarts      gene_type  gene_name
##                <factor>   <factor>    <factor>       <factor>   <factor>
##   AL627309.1          .          .           .         lncRNA AL627309.1
##   AL627309.5          .          .           .         lncRNA AL627309.5
##   AL627309.4          .          .           .         lncRNA AL627309.4
##   AP006222.2          .          .           .         lncRNA AP006222.2
##   AL669831.2          .          .           .         lncRNA AL669831.2
##          ...        ...        ...         ...            ...        ...
##   AC004556.3          .          .           . protein_coding AC004556.3
##   AC233755.2          .          .           . protein_coding AC233755.2
##   AC233755.1          .          .           . protein_coding AC233755.1
##   AC007325.1          .          .           . protein_coding AC007325.1
##   AC007325.4          .          .           . protein_coding AC007325.4
##               hgnc_id          havana_gene               tag  n_counts
##              <factor>             <factor>          <factor> <numeric>
##   AL627309.1       NA OTTHUMG00000001096.2 overlapping_locus        70
##   AL627309.5       NA OTTHUMG00000002480.4 ncRNA_host              442
##   AL627309.4       NA OTTHUMG00000002525.1 NA                       44
##   AP006222.2       NA OTTHUMG00000194680.1 NA                        1
##   AL669831.2       NA OTTHUMG00000002408.1 NA                       10
##          ...      ...                  ...               ...       ...
##   AC004556.3       NA                   NA                NA       320
##   AC233755.2       NA                   NA                NA         1
##   AC233755.1       NA                   NA                NA         1
##   AC007325.1       NA                   NA                NA         3
##   AC007325.4       NA                   NA                NA        43
##              highly_variable highly_variable_rank       means   variances
##                    <logical>            <numeric>   <numeric>   <numeric>
##   AL627309.1           FALSE                  NaN 0.007268196 0.008462225
##   AL627309.5           FALSE                  NaN 0.045893469 0.050853072
##   AL627309.4           FALSE                  NaN 0.004568581 0.004755865
##   AP006222.2           FALSE                  NaN 0.000103831 0.000103831
##   AL669831.2           FALSE                  NaN 0.001038314 0.001037343
##          ...             ...                  ...         ...         ...
##   AC004556.3           FALSE                  NaN 0.033226041 0.035448356
##   AC233755.2           FALSE                  NaN 0.000103831 0.000103831
##   AC233755.1           FALSE                  NaN 0.000103831 0.000103831
##   AC007325.1           FALSE                  NaN 0.000311494 0.000311429
##   AC007325.4           FALSE                  NaN 0.004464749 0.004445277
##              variances_norm
##                   <numeric>
##   AL627309.1       0.971895
##   AL627309.5       0.888672
##   AL627309.4       0.891707
##   AP006222.2       0.999904
##   AL669831.2       0.933916
##          ...            ...
##   AC004556.3       0.856870
##   AC233755.2       0.999904
##   AC233755.1       0.999904
##   AC007325.1       0.975766
##   AC007325.4       0.854162
##   -------
##   seqinfo: 34 sequences from an unspecified genome; no seqlengths
```
We next import the PeakMatrix. The rownames are already in the format of seqnames:start-end, so rowRanges can be extracted directly from rowData.


``` r
library(zellkonverter)
url <- "http://download.gao-lab.org/GLUE/dataset/10x-Multiome-Pbmc10k-ATAC.h5ad"
destfile <- tempfile(fileext = ".h5ad")

# Download the file
download.file(url, destfile, mode = "wb")

PeakMatrix <- readH5AD(destfile)
rowRanges(PeakMatrix) <- GRanges(rowData(PeakMatrix))
rowRanges(PeakMatrix)
```

```
## GRanges object with 107194 ranges and 3 metadata columns:
##                            seqnames        ranges strand | feature_types
##                               <Rle>     <IRanges>  <Rle> |      <factor>
##       chr1:816881-817647       chr1 816881-817647      * |         Peaks
##       chr1:819912-823500       chr1 819912-823500      * |         Peaks
##       chr1:825827-825889       chr1 825827-825889      * |         Peaks
##       chr1:826612-827979       chr1 826612-827979      * |         Peaks
##       chr1:841243-843059       chr1 841243-843059      * |         Peaks
##                      ...        ...           ...    ... .           ...
##   KI270713.1:20444-22615 KI270713.1   20444-22615      * |         Peaks
##   KI270713.1:27118-28927 KI270713.1   27118-28927      * |         Peaks
##   KI270713.1:29485-30706 KI270713.1   29485-30706      * |         Peaks
##   KI270713.1:31511-32072 KI270713.1   31511-32072      * |         Peaks
##   KI270713.1:37129-37638 KI270713.1   37129-37638      * |         Peaks
##                            genome  n_counts
##                          <factor> <numeric>
##       chr1:816881-817647   GRCh38      1025
##       chr1:819912-823500   GRCh38      1384
##       chr1:825827-825889   GRCh38        20
##       chr1:826612-827979   GRCh38      4555
##       chr1:841243-843059   GRCh38       555
##                      ...      ...       ...
##   KI270713.1:20444-22615   GRCh38     12640
##   KI270713.1:27118-28927   GRCh38       533
##   KI270713.1:29485-30706   GRCh38       622
##   KI270713.1:31511-32072   GRCh38       207
##   KI270713.1:37129-37638   GRCh38       388
##   -------
##   seqinfo: 32 sequences from an unspecified genome; no seqlengths
```

This example SingleCellExperiment dataset is missing the dimensionality reduction information. Users should have been this precalculated. Refer to the section on how to import from [csv files](### Reading directly from CSV),  or refer to [OSCA](https://bioconductor.org/books/3.20/OSCA.basic/dimensionality-reduction.html) book on how to perform dimensionality reduction. 

## Starting from 10x data formats 
### Reading from .h5 file 

We first download the publicly available PBMC data from 10x Genomics


``` r
url <- "https://cf.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5"
destfile <- tempfile(fileext = ".h5")  

# Download the file
download.file(url, destfile, mode = "wb")
```
This .h5 file contains both genes and peak regions as features, so we read in all the features into a single `SingleCellExperiment`, and we can specify "Gene Expression" as the main experiment and "Peaks" as the alternative experiment. We then split this single object into 2 separate  `SingleCellExperiment` objects for consistency with other sections.

``` r
library(DropletUtils)
sce10x <- read10xCounts(destfile)

sce10x <- splitAltExps(sce10x, f=rowData(sce10x)$Type, ref="Gene Expression")

peakMatrix <- altExp(sce10x)
rowRanges(peakMatrix) <- GRanges(rownames(peakMatrix))

GeneExpressionMatrix <- removeAltExps(sce10x)
```

### Reading directly from a 10x directory

We can also create a `SingleCellExperiment` object from a directory containing "matrix.mtx.gz", "barcodes.tsv.gz" and "features.tsv.gz"

``` r
example(write10xCounts)
```

```
## 
## wrt10C> # Mocking up some count data.
## wrt10C> library(Matrix)
## 
## wrt10C> my.counts <- matrix(rpois(1000, lambda=5), ncol=10, nrow=100)
## 
## wrt10C> my.counts <- as(my.counts, "CsparseMatrix")
## 
## wrt10C> cell.ids <- paste0("BARCODE-", seq_len(ncol(my.counts)))
## 
## wrt10C> ngenes <- nrow(my.counts)
## 
## wrt10C> gene.ids <- paste0("ENSG0000", seq_len(ngenes))
## 
## wrt10C> gene.symb <- paste0("GENE", seq_len(ngenes))
## 
## wrt10C> # Writing this to file:
## wrt10C> tmpdir <- tempfile()
## 
## wrt10C> write10xCounts(tmpdir, my.counts, gene.id=gene.ids, 
## wrt10C+     gene.symbol=gene.symb, barcodes=cell.ids)
## 
## wrt10C> list.files(tmpdir)
## [1] "barcodes.tsv" "genes.tsv"    "matrix.mtx"  
## 
## wrt10C> # Creating a version 3 HDF5 file:
## wrt10C> tmph5 <- tempfile(fileext=".h5")
## 
## wrt10C> write10xCounts(tmph5, my.counts, gene.id=gene.ids, 
## wrt10C+     gene.symbol=gene.symb, barcodes=cell.ids, version='3')
```

``` r
list.files(tmpdir)
```

```
## [1] "barcodes.tsv" "genes.tsv"    "matrix.mtx"
```

``` r
sce10x <- read10xCounts(tmpdir)
```

Follow the section on [Constructing a SingleCellExperiment object](## Starting from a SingleCellExperiment object) on how to construct SingleCellExperiment, add colData, rowData and reduced dimensions

## Reading directly from CSV

If the count matrices are in the form of csv files, read in the count matrix as a sparse matrix and convert it to a `dgCMatrix`. 


``` r
library(SparseArray)
counts <- readSparseCSV("counts.csv")
counts <- as(counts, "dgCMatrix")
```


Then read in the genomic positions as data.frame and convert it to a `GRanges`.

``` r
peak_position <- read.csv("peaks.csv")
peak_position <- data.frame(peak_position)
gr <- GenomicRanges::makeGRangesFromDataFrame(peak_position)
```

Refer to the section on [Constructing a SingleCellExperiment object](## Constructing a SingleCellExperiment object) on how to construct SingleCellExperiment, add colData, rowData and reduced dimensions.

## Important points

Finally, we would like to emphasize again the importance of converting all count matrices to dgCMatrix for speed and compatibility. Both gene expression and peak matrices require `rowRanges` to indicate genomic positions of target gene position and peak position.
