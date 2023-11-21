
           ___      .______        ______  __    __  .______      
          /   \     |   _  \      /      ||  |  |  | |   _  \     
         /  ^  \    |  |_)  |    |  ,----'|  |__|  | |  |_)  |    
        /  /_\  \   |      /     |  |     |   __   | |      /     
       /  _____  \  |  |\  \\___ |  `----.|  |  |  | |  |\  \\___.
      /__/     \__\ | _| `._____| \______||__|  |__| | _| `._____|
    
Logging With ArchR!

Start Time : 2023-11-13 04:19:48.849849

------- ArchR Info

ArchRThreads = 1

------- System Info

Computer OS = unix

------- Session Info

R version 4.3.0 (2023-04-21)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 18.04.6 LTS

Matrix products: default
BLAS:   /usr/local/lib/R/lib/libRblas.so 
LAPACK: /usr/local/lib/R/lib/libRlapack.so;  LAPACK version 3.11.0

locale:
 [1] LC_CTYPE=en_US.UTF-8 LC_NUMERIC=C         LC_TIME=C           
 [4] LC_COLLATE=C         LC_MONETARY=C        LC_MESSAGES=C       
 [7] LC_PAPER=C           LC_NAME=C            LC_ADDRESS=C        
[10] LC_TELEPHONE=C       LC_MEASUREMENT=C     LC_IDENTIFICATION=C 

time zone: Etc/UTC
tzcode source: system (glibc)

attached base packages:
[1] grid      stats4    stats     graphics  grDevices utils     datasets 
[8] methods   base     

other attached packages:
 [1] rhdf5_2.45.1                RcppArmadillo_0.12.6.4.0   
 [3] Rcpp_1.0.11                 Matrix_1.6-1.1             
 [5] sparseMatrixStats_1.13.4    data.table_1.14.8          
 [7] stringr_1.5.0               plyr_1.8.9                 
 [9] magrittr_2.0.3              ggplot2_3.4.4              
[11] gtable_0.3.4                gtools_3.9.4               
[13] gridExtra_2.3               devtools_2.4.5             
[15] usethis_2.1.6               ArchR_1.0.3                
[17] epiregulon_1.0.34           SingleCellExperiment_1.23.0
[19] SummarizedExperiment_1.31.1 Biobase_2.61.0             
[21] GenomicRanges_1.53.3        GenomeInfoDb_1.37.7        
[23] IRanges_2.35.3              S4Vectors_0.39.3           
[25] BiocGenerics_0.47.1         MatrixGenerics_1.13.2      
[27] matrixStats_1.0.0          

loaded via a namespace (and not attached):
  [1] rstudioapi_0.15.0                  jsonlite_1.8.7                    
  [3] MultiAssayExperiment_1.27.5        ggbeeswarm_0.7.2                  
  [5] rmarkdown_2.25                     fs_1.6.2                          
  [7] BiocIO_1.11.0                      zlibbioc_1.47.0                   
  [9] vctrs_0.6.4                        memoise_2.0.1                     
 [11] Rsamtools_2.17.0                   DelayedMatrixStats_1.23.9         
 [13] RCurl_1.98-1.12                    htmltools_0.5.5                   
 [15] S4Arrays_1.1.6                     Rhdf5lib_1.23.2                   
 [17] SparseArray_1.1.12                 sass_0.4.6                        
 [19] bslib_0.5.1                        htmlwidgets_1.6.2                 
 [21] cachem_1.0.8                       GenomicAlignments_1.37.0          
 [23] igraph_1.5.1                       mime_0.12                         
 [25] lifecycle_1.0.3                    pkgconfig_2.0.3                   
 [27] rsvd_1.0.5                         R6_2.5.1                          
 [29] fastmap_1.1.1                      GenomeInfoDbData_1.2.11           
 [31] shiny_1.7.5.1                      digest_0.6.31                     
 [33] colorspace_2.1-0                   AnnotationDbi_1.63.2              
 [35] ps_1.7.5                           irlba_2.3.5.1                     
 [37] pkgload_1.3.2                      RSQLite_2.3.1                     
 [39] beachmat_2.17.17                   fansi_1.0.4                       
 [41] httr_1.4.6                         abind_1.4-5                       
 [43] compiler_4.3.0                     remotes_2.4.2                     
 [45] withr_2.5.0                        bit64_4.0.5                       
 [47] BiocParallel_1.35.4                DBI_1.1.3                         
 [49] pkgbuild_1.4.0                     HDF5Array_1.29.3                  
 [51] sessioninfo_1.2.2                  DelayedArray_0.27.10              
 [53] rjson_0.2.21                       BSgenome.Hsapiens.UCSC.hg19_1.4.3 
 [55] tools_4.3.0                        vipor_0.4.5                       
 [57] beeswarm_0.4.0                     httpuv_1.6.11                     
 [59] glue_1.6.2                         restfulr_0.0.15                   
 [61] callr_3.7.3                        rhdf5filters_1.13.5               
 [63] promises_1.2.0.1                   generics_0.1.3                    
 [65] BSgenome_1.69.1                    BiocSingular_1.17.1               
 [67] ScaledMatrix_1.9.1                 utf8_1.2.3                        
 [69] XVector_0.41.2                     pillar_1.9.0                      
 [71] GSVA_1.49.8                        later_1.3.1                       
 [73] dplyr_1.1.3                        lattice_0.22-5                    
 [75] rtracklayer_1.61.2                 bit_4.0.5                         
 [77] annotate_1.79.0                    tidyselect_1.2.0                  
 [79] BSgenome.Hsapiens.UCSC.hg38_1.4.5  Biostrings_2.69.2                 
 [81] miniUI_0.1.1.1                     knitr_1.44                        
 [83] bookdown_0.36                      xfun_0.39                         
 [85] stringi_1.7.12                     yaml_2.3.7                        
 [87] evaluate_0.21                      codetools_0.2-19                  
 [89] BSgenome.Mmusculus.UCSC.mm10_1.4.3 tibble_3.2.1                      
 [91] graph_1.79.4                       cli_3.6.1                         
 [93] xtable_1.8-4                       munsell_0.5.0                     
 [95] processx_3.8.1                     jquerylib_0.1.4                   
 [97] png_0.1-8                          XML_3.99-0.14                     
 [99] parallel_4.3.0                     ellipsis_0.3.2                    
[101] blob_1.2.4                         prettyunits_1.1.1                 
[103] profvis_0.3.8                      urlchecker_1.0.1                  
[105] bitops_1.0-7                       GSEABase_1.63.0                   
[107] scales_1.2.1                       purrr_1.0.2                       
[109] crayon_1.5.2                       rlang_1.1.1                       
[111] KEGGREST_1.41.4                   


------- Log Info


2023-11-13 04:19:48.98504 : getMatrixFromProject Input-Parameters, Class = list

getMatrixFromProject Input-Parameters$ArchRProj: length = 1

getMatrixFromProject Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromProject Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromProject Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromProject Input-Parameters$verbose: length = 1
[1] TRUE


getMatrixFromProject Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromProject Input-Parameters$threads: length = 1
[1] 1


getMatrixFromProject Input-Parameters$logFile: length = 1
[1] "x"


2023-11-13 04:19:52.941444 : Reading GeneIntegrationMatrix : SRR13927737(1 of 4), 0.068 mins elapsed.

2023-11-13 04:19:52.980292 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927737 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927737.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 4439
[1] "SRR13927737#GTCACGGGTCAACAGG-1" "SRR13927737#TGCCTGTTCTTACTCA-1"
[3] "SRR13927737#TTGCGGGTCAGGTCTA-1" "SRR13927737#GGTTGCGCAAAGCTTC-1"
[5] "SRR13927737#CTTAATCTCGCGCTGA-1" "SRR13927737#TGTACGATCGGAAAGC-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-13 04:19:54.161611 : featureDF SRR13927737, Class = DFrame
DataFrame with 23525 rows and 6 columns
      seqnames     start       end    strand        name       idx
         <Rle> <integer> <integer> <integer> <character> <integer>
1         chr1     69091     70008         1       OR4F5         1
2         chr1    817371    819837         1      FAM87B         2
3         chr1    825138    859446         1   LINC01128         3
4         chr1    827522    826206         2   LINC00115         4
5         chr1    876903    868071         2      FAM41C         5
...        ...       ...       ...       ...         ...       ...
23521     chrX 155264589 155258241         2      RAB39B       923
23522     chrX 155334657 155276211         2       CLIC2       924
23523     chrX 155380787 155381134         1      H2AFB1       925
23524     chrX 155382115 155383230         1        F8A2       926
23525     chrX 155669944 155490115         2       TMLHE       927

2023-11-13 04:19:54.232423 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927737.arrow, 0.021 mins elapsed.

2023-11-13 04:20:32.833171 : mat SRR13927737, Class = dgCMatrix
mat SRR13927737: nRows = 23525, nCols = 4439
mat SRR13927737: NonZeroEntries = 50277774, EntryRange = [ 0.000148130054904683 , 376.345870564644 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927737#GTCACGGGTCAACAGG-1 SRR13927737#TGCCTGTTCTTACTCA-1
[1,]                      .                             .         
[2,]                      .                             .         
[3,]                      0.1057576                     0.09630696
[4,]                      .                             0.01806508
[5,]                      .                             .         
     SRR13927737#TTGCGGGTCAGGTCTA-1 SRR13927737#GGTTGCGCAAAGCTTC-1
[1,]                     .                              .         
[2,]                     .                              .         
[3,]                     0.09113499                     0.08396973
[4,]                     0.02096304                     0.01727603
[5,]                     .                              .         
     SRR13927737#CTTAATCTCGCGCTGA-1
[1,]                     .         
[2,]                     .         
[3,]                     0.08327649
[4,]                     0.01368917
[5,]                     .         


2023-11-13 04:20:33.31441 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927737.arrow, 0.672 mins elapsed.

2023-11-13 04:20:34.006922 : se SRR13927737, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 4439 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(4439): SRR13927737#GTCACGGGTCAACAGG-1
  SRR13927737#TGCCTGTTCTTACTCA-1 ... SRR13927737#TTGCTATAGATGTTCC-1
  SRR13927737#GCTCACTAGTTAGCGG-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-13 04:20:34.015947 : Completed GeneIntegrationMatrix : SRR13927737(1 of 4), 0.753 mins elapsed.
2023-11-13 04:20:34.016556 : Reading GeneIntegrationMatrix : SRR13927738(2 of 4), 0.753 mins elapsed.

2023-11-13 04:20:34.040802 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927738 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927738.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3736
[1] "SRR13927738#TAGCATGAGGCTAAAT-1" "SRR13927738#TGTTAGGTCGGTCAGC-1"
[3] "SRR13927738#CTCTCAGCATGTCCCT-1" "SRR13927738#ACTATTCGTCCAACCG-1"
[5] "SRR13927738#TGGTCCTTCTACTGCC-1" "SRR13927738#ATTGTCTAGACTAGGC-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-13 04:20:35.00245 : featureDF SRR13927738, Class = DFrame
DataFrame with 23525 rows and 6 columns
      seqnames     start       end    strand        name       idx
         <Rle> <integer> <integer> <integer> <character> <integer>
1         chr1     69091     70008         1       OR4F5         1
2         chr1    817371    819837         1      FAM87B         2
3         chr1    825138    859446         1   LINC01128         3
4         chr1    827522    826206         2   LINC00115         4
5         chr1    876903    868071         2      FAM41C         5
...        ...       ...       ...       ...         ...       ...
23521     chrX 155264589 155258241         2      RAB39B       923
23522     chrX 155334657 155276211         2       CLIC2       924
23523     chrX 155380787 155381134         1      H2AFB1       925
23524     chrX 155382115 155383230         1        F8A2       926
23525     chrX 155669944 155490115         2       TMLHE       927

2023-11-13 04:20:35.034617 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927738.arrow, 0.017 mins elapsed.

2023-11-13 04:20:56.774514 : mat SRR13927738, Class = dgCMatrix
mat SRR13927738: nRows = 23525, nCols = 3736
mat SRR13927738: NonZeroEntries = 41898897, EntryRange = [ 9.73975135001421e-05 , 758.763031303812 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927738#TAGCATGAGGCTAAAT-1 SRR13927738#TGTTAGGTCGGTCAGC-1
[1,]                    .                                .        
[2,]                    .                                .        
[3,]                    0.002346352                      0.0084575
[4,]                    0.005836413                      .        
[5,]                    .                                .        
     SRR13927738#CTCTCAGCATGTCCCT-1 SRR13927738#ACTATTCGTCCAACCG-1
[1,]                     .                             .          
[2,]                     .                             .          
[3,]                     .                             0.004081077
[4,]                     0.02001583                    .          
[5,]                     .                             .          
     SRR13927738#TGGTCCTTCTACTGCC-1
[1,]                    .          
[2,]                    .          
[3,]                    0.007640877
[4,]                    0.015048796
[5,]                    .          


2023-11-13 04:20:57.107709 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927738.arrow, 0.384 mins elapsed.

2023-11-13 04:20:57.409441 : se SRR13927738, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 3736 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(3736): SRR13927738#TAGCATGAGGCTAAAT-1
  SRR13927738#TGTTAGGTCGGTCAGC-1 ... SRR13927738#GCTCCTATCCGTGCGA-1
  SRR13927738#GGCACGTAGTCGATAA-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-13 04:20:57.414642 : Completed GeneIntegrationMatrix : SRR13927738(2 of 4), 1.143 mins elapsed.
2023-11-13 04:20:57.415082 : Reading GeneIntegrationMatrix : SRR13927735(3 of 4), 1.143 mins elapsed.

2023-11-13 04:20:57.448454 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927735 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927735.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3813
[1] "SRR13927735#TTATGTCTCCAGGTAT-1" "SRR13927735#TATTGCTCATCAGAAA-1"
[3] "SRR13927735#TTCGATTGTAGGGTTG-1" "SRR13927735#CATTCATTCGGATGTT-1"
[5] "SRR13927735#ACGTTAGGTCAACTGT-1" "SRR13927735#AAATGCCCAGCAATGG-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-13 04:20:58.083169 : featureDF SRR13927735, Class = DFrame
DataFrame with 23525 rows and 6 columns
      seqnames     start       end    strand        name       idx
         <Rle> <integer> <integer> <integer> <character> <integer>
1         chr1     69091     70008         1       OR4F5         1
2         chr1    817371    819837         1      FAM87B         2
3         chr1    825138    859446         1   LINC01128         3
4         chr1    827522    826206         2   LINC00115         4
5         chr1    876903    868071         2      FAM41C         5
...        ...       ...       ...       ...         ...       ...
23521     chrX 155264589 155258241         2      RAB39B       923
23522     chrX 155334657 155276211         2       CLIC2       924
23523     chrX 155380787 155381134         1      H2AFB1       925
23524     chrX 155382115 155383230         1        F8A2       926
23525     chrX 155669944 155490115         2       TMLHE       927

2023-11-13 04:20:58.09939 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927735.arrow, 0.011 mins elapsed.

2023-11-13 04:21:20.779682 : mat SRR13927735, Class = dgCMatrix
mat SRR13927735: nRows = 23525, nCols = 3813
mat SRR13927735: NonZeroEntries = 45535634, EntryRange = [ 7.20715653308254e-05 , 673.990433360078 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927735#TTATGTCTCCAGGTAT-1 SRR13927735#TATTGCTCATCAGAAA-1
[1,]                    .                              .          
[2,]                    .                              .          
[3,]                    0.070892732                    0.057009268
[4,]                    0.008800039                    0.002978408
[5,]                    .                              .          
     SRR13927735#TTCGATTGTAGGGTTG-1 SRR13927735#CATTCATTCGGATGTT-1
[1,]                     .                               .        
[2,]                     .                               .        
[3,]                     0.05523011                      0.0248301
[4,]                     .                               .        
[5,]                     .                               .        
     SRR13927735#ACGTTAGGTCAACTGT-1
[1,]                     .         
[2,]                     .         
[3,]                     0.05384571
[4,]                     0.01228688
[5,]                     .         


2023-11-13 04:21:21.146522 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927735.arrow, 0.395 mins elapsed.

2023-11-13 04:21:21.447566 : se SRR13927735, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 3813 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(3813): SRR13927735#TTATGTCTCCAGGTAT-1
  SRR13927735#TATTGCTCATCAGAAA-1 ... SRR13927735#GATCATGGTGTTGTTG-1
  SRR13927735#GTCCATCCATTTCACT-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-13 04:21:21.453427 : Completed GeneIntegrationMatrix : SRR13927735(3 of 4), 1.543 mins elapsed.
2023-11-13 04:21:21.453822 : Reading GeneIntegrationMatrix : SRR13927736(4 of 4), 1.543 mins elapsed.

2023-11-13 04:21:21.476689 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927736 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927736.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3534
[1] "SRR13927736#GCTGTTCTCCTTGACC-1" "SRR13927736#AGTTACGAGTAGGTCG-1"
[3] "SRR13927736#AGTTTGGAGCATTCCA-1" "SRR13927736#GCACGGTTCGGCAATT-1"
[5] "SRR13927736#TCAGCTCTCTCTATTG-1" "SRR13927736#CCTTAATTCTACTTTG-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-13 04:21:22.373833 : featureDF SRR13927736, Class = DFrame
DataFrame with 23525 rows and 6 columns
      seqnames     start       end    strand        name       idx
         <Rle> <integer> <integer> <integer> <character> <integer>
1         chr1     69091     70008         1       OR4F5         1
2         chr1    817371    819837         1      FAM87B         2
3         chr1    825138    859446         1   LINC01128         3
4         chr1    827522    826206         2   LINC00115         4
5         chr1    876903    868071         2      FAM41C         5
...        ...       ...       ...       ...         ...       ...
23521     chrX 155264589 155258241         2      RAB39B       923
23522     chrX 155334657 155276211         2       CLIC2       924
23523     chrX 155380787 155381134         1      H2AFB1       925
23524     chrX 155382115 155383230         1        F8A2       926
23525     chrX 155669944 155490115         2       TMLHE       927

2023-11-13 04:21:22.390387 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927736.arrow, 0.015 mins elapsed.

2023-11-13 04:21:39.816405 : mat SRR13927736, Class = dgCMatrix
mat SRR13927736: nRows = 23525, nCols = 3534
mat SRR13927736: NonZeroEntries = 39003616, EntryRange = [ 8.77627031078293e-05 , 770.387276519725 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927736#GCTGTTCTCCTTGACC-1 SRR13927736#AGTTACGAGTAGGTCG-1
[1,]                     .                             .          
[2,]                     .                             .          
[3,]                     0.07845381                    0.028344539
[4,]                     0.01395289                    0.008499189
[5,]                     .                             .          
     SRR13927736#AGTTTGGAGCATTCCA-1 SRR13927736#GCACGGTTCGGCAATT-1
[1,]                     .                              .         
[2,]                     .                              .         
[3,]                     0.06535130                     0.06576222
[4,]                     0.01810177                     0.01315728
[5,]                     .                              .         
     SRR13927736#TCAGCTCTCTCTATTG-1
[1,]                    .          
[2,]                    .          
[3,]                    0.016229526
[4,]                    0.003862592
[5,]                    .          


2023-11-13 04:21:40.115661 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927736.arrow, 0.311 mins elapsed.

2023-11-13 04:21:40.363046 : se SRR13927736, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 3534 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(3534): SRR13927736#GCTGTTCTCCTTGACC-1
  SRR13927736#AGTTACGAGTAGGTCG-1 ... SRR13927736#TAAACCGCAAGTCTGT-1
  SRR13927736#TGCATTTGTTCTACCC-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-13 04:21:40.36806 : Completed GeneIntegrationMatrix : SRR13927736(4 of 4), 1.859 mins elapsed.
2023-11-13 04:21:40.36897 : Organizing colData, 1.859 mins elapsed.
2023-11-13 04:21:40.51463 : Organizing rowData, 1.861 mins elapsed.
2023-11-13 04:21:40.517638 : Organizing rowRanges, 1.861 mins elapsed.
2023-11-13 04:21:40.523145 : Organizing Assays (1 of 1), 1.861 mins elapsed.
2023-11-13 04:21:47.747395 : Constructing SummarizedExperiment, 1.982 mins elapsed.
2023-11-13 04:21:48.979943 : Finished Matrix Creation, 2.002 mins elapsed.

2023-11-13 04:21:48.982364 : getMatrixFromProject Input-Parameters, Class = list

getMatrixFromProject Input-Parameters$ArchRProj: length = 1

getMatrixFromProject Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromProject Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromProject Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromProject Input-Parameters$verbose: length = 1
[1] TRUE


getMatrixFromProject Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromProject Input-Parameters$threads: length = 1
[1] 1


getMatrixFromProject Input-Parameters$logFile: length = 1
[1] "x"


2023-11-13 04:21:50.641209 : Reading PeakMatrix : SRR13927737(1 of 4), 0.028 mins elapsed.

2023-11-13 04:21:50.673947 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927737 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927737.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 4439
[1] "SRR13927737#GTCACGGGTCAACAGG-1" "SRR13927737#TGCCTGTTCTTACTCA-1"
[3] "SRR13927737#TTGCGGGTCAGGTCTA-1" "SRR13927737#GGTTGCGCAAAGCTTC-1"
[5] "SRR13927737#CTTAATCTCGCGCTGA-1" "SRR13927737#TGTACGATCGGAAAGC-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-13 04:21:51.454308 : featureDF SRR13927737, Class = DFrame
DataFrame with 210721 rows and 4 columns
       seqnames       idx     start       end
          <Rle> <integer> <integer> <integer>
1          chr1         1    804677    805177
2          chr1         2    811866    812366
3          chr1         3    817090    817590
4          chr1         4    818831    819331
5          chr1         5    826543    827043
...         ...       ...       ...       ...
210717     chrX      5777 155898525 155899025
210718     chrX      5778 155899804 155900304
210719     chrX      5779 155956608 155957108
210720     chrX      5780 155966827 155967327
210721     chrX      5781 155971625 155972125

2023-11-13 04:21:51.522882 : Getting PeakMatrix from ArrowFile : SRR13927737.arrow, 0.014 mins elapsed.

2023-11-13 04:22:03.478211 : mat SRR13927737, Class = dgCMatrix
mat SRR13927737: nRows = 210721, nCols = 4439
mat SRR13927737: NonZeroEntries = 29963541, EntryRange = [ 1 , 4 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927737#GTCACGGGTCAACAGG-1 SRR13927737#TGCCTGTTCTTACTCA-1
[1,]                              .                              .
[2,]                              .                              2
[3,]                              .                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927737#TTGCGGGTCAGGTCTA-1 SRR13927737#GGTTGCGCAAAGCTTC-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              2                              .
[4,]                              .                              .
[5,]                              4                              .
     SRR13927737#CTTAATCTCGCGCTGA-1
[1,]                              .
[2,]                              .
[3,]                              .
[4,]                              .
[5,]                              .


2023-11-13 04:22:03.707798 : Organizing SE PeakMatrix from ArrowFile : SRR13927737.arrow, 0.217 mins elapsed.

2023-11-13 04:22:04.009673 : se SRR13927737, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 4439 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(4439): SRR13927737#GTCACGGGTCAACAGG-1
  SRR13927737#TGCCTGTTCTTACTCA-1 ... SRR13927737#TTGCTATAGATGTTCC-1
  SRR13927737#GCTCACTAGTTAGCGG-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-13 04:22:04.019036 : Completed PeakMatrix : SRR13927737(1 of 4), 0.251 mins elapsed.
2023-11-13 04:22:04.019432 : Reading PeakMatrix : SRR13927738(2 of 4), 0.251 mins elapsed.

2023-11-13 04:22:04.032003 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927738 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927738.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3736
[1] "SRR13927738#TAGCATGAGGCTAAAT-1" "SRR13927738#TGTTAGGTCGGTCAGC-1"
[3] "SRR13927738#CTCTCAGCATGTCCCT-1" "SRR13927738#ACTATTCGTCCAACCG-1"
[5] "SRR13927738#TGGTCCTTCTACTGCC-1" "SRR13927738#ATTGTCTAGACTAGGC-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-13 04:22:04.640627 : featureDF SRR13927738, Class = DFrame
DataFrame with 210721 rows and 4 columns
       seqnames       idx     start       end
          <Rle> <integer> <integer> <integer>
1          chr1         1    804677    805177
2          chr1         2    811866    812366
3          chr1         3    817090    817590
4          chr1         4    818831    819331
5          chr1         5    826543    827043
...         ...       ...       ...       ...
210717     chrX      5777 155898525 155899025
210718     chrX      5778 155899804 155900304
210719     chrX      5779 155956608 155957108
210720     chrX      5780 155966827 155967327
210721     chrX      5781 155971625 155972125

2023-11-13 04:22:04.691145 : Getting PeakMatrix from ArrowFile : SRR13927738.arrow, 0.011 mins elapsed.

2023-11-13 04:22:16.380662 : mat SRR13927738, Class = dgCMatrix
mat SRR13927738: nRows = 210721, nCols = 3736
mat SRR13927738: NonZeroEntries = 30165921, EntryRange = [ 1 , 4 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927738#TAGCATGAGGCTAAAT-1 SRR13927738#TGTTAGGTCGGTCAGC-1
[1,]                              2                              .
[2,]                              .                              .
[3,]                              .                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927738#CTCTCAGCATGTCCCT-1 SRR13927738#ACTATTCGTCCAACCG-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              .                              .
[4,]                              .                              .
[5,]                              2                              .
     SRR13927738#TGGTCCTTCTACTGCC-1
[1,]                              .
[2,]                              .
[3,]                              2
[4,]                              .
[5,]                              .


2023-11-13 04:22:16.616275 : Organizing SE PeakMatrix from ArrowFile : SRR13927738.arrow, 0.21 mins elapsed.

2023-11-13 04:22:16.929423 : se SRR13927738, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 3736 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(3736): SRR13927738#TAGCATGAGGCTAAAT-1
  SRR13927738#TGTTAGGTCGGTCAGC-1 ... SRR13927738#GCTCCTATCCGTGCGA-1
  SRR13927738#GGCACGTAGTCGATAA-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-13 04:22:16.936591 : Completed PeakMatrix : SRR13927738(2 of 4), 0.466 mins elapsed.
2023-11-13 04:22:16.93702 : Reading PeakMatrix : SRR13927735(3 of 4), 0.466 mins elapsed.

2023-11-13 04:22:16.96289 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927735 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927735.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3813
[1] "SRR13927735#TTATGTCTCCAGGTAT-1" "SRR13927735#TATTGCTCATCAGAAA-1"
[3] "SRR13927735#TTCGATTGTAGGGTTG-1" "SRR13927735#CATTCATTCGGATGTT-1"
[5] "SRR13927735#ACGTTAGGTCAACTGT-1" "SRR13927735#AAATGCCCAGCAATGG-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-13 04:22:17.760486 : featureDF SRR13927735, Class = DFrame
DataFrame with 210721 rows and 4 columns
       seqnames       idx     start       end
          <Rle> <integer> <integer> <integer>
1          chr1         1    804677    805177
2          chr1         2    811866    812366
3          chr1         3    817090    817590
4          chr1         4    818831    819331
5          chr1         5    826543    827043
...         ...       ...       ...       ...
210717     chrX      5777 155898525 155899025
210718     chrX      5778 155899804 155900304
210719     chrX      5779 155956608 155957108
210720     chrX      5780 155966827 155967327
210721     chrX      5781 155971625 155972125

2023-11-13 04:22:17.860524 : Getting PeakMatrix from ArrowFile : SRR13927735.arrow, 0.015 mins elapsed.

2023-11-13 04:22:28.420969 : mat SRR13927735, Class = dgCMatrix
mat SRR13927735: nRows = 210721, nCols = 3813
mat SRR13927735: NonZeroEntries = 26038084, EntryRange = [ 1 , 4 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927735#TTATGTCTCCAGGTAT-1 SRR13927735#TATTGCTCATCAGAAA-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              1                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927735#TTCGATTGTAGGGTTG-1 SRR13927735#CATTCATTCGGATGTT-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              .                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927735#ACGTTAGGTCAACTGT-1
[1,]                              2
[2,]                              .
[3,]                              .
[4,]                              .
[5,]                              .


2023-11-13 04:22:28.618948 : Organizing SE PeakMatrix from ArrowFile : SRR13927735.arrow, 0.194 mins elapsed.

2023-11-13 04:22:28.900783 : se SRR13927735, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 3813 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(3813): SRR13927735#TTATGTCTCCAGGTAT-1
  SRR13927735#TATTGCTCATCAGAAA-1 ... SRR13927735#GATCATGGTGTTGTTG-1
  SRR13927735#GTCCATCCATTTCACT-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-13 04:22:28.908024 : Completed PeakMatrix : SRR13927735(3 of 4), 0.665 mins elapsed.
2023-11-13 04:22:28.908472 : Reading PeakMatrix : SRR13927736(4 of 4), 0.665 mins elapsed.

2023-11-13 04:22:28.921191 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927736 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927736.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3534
[1] "SRR13927736#GCTGTTCTCCTTGACC-1" "SRR13927736#AGTTACGAGTAGGTCG-1"
[3] "SRR13927736#AGTTTGGAGCATTCCA-1" "SRR13927736#GCACGGTTCGGCAATT-1"
[5] "SRR13927736#TCAGCTCTCTCTATTG-1" "SRR13927736#CCTTAATTCTACTTTG-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-13 04:22:29.713809 : featureDF SRR13927736, Class = DFrame
DataFrame with 210721 rows and 4 columns
       seqnames       idx     start       end
          <Rle> <integer> <integer> <integer>
1          chr1         1    804677    805177
2          chr1         2    811866    812366
3          chr1         3    817090    817590
4          chr1         4    818831    819331
5          chr1         5    826543    827043
...         ...       ...       ...       ...
210717     chrX      5777 155898525 155899025
210718     chrX      5778 155899804 155900304
210719     chrX      5779 155956608 155957108
210720     chrX      5780 155966827 155967327
210721     chrX      5781 155971625 155972125

2023-11-13 04:22:29.764508 : Getting PeakMatrix from ArrowFile : SRR13927736.arrow, 0.014 mins elapsed.

2023-11-13 04:22:41.114779 : mat SRR13927736, Class = dgCMatrix
mat SRR13927736: nRows = 210721, nCols = 3534
mat SRR13927736: NonZeroEntries = 25583439, EntryRange = [ 1 , 4 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927736#GCTGTTCTCCTTGACC-1 SRR13927736#AGTTACGAGTAGGTCG-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              .                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927736#AGTTTGGAGCATTCCA-1 SRR13927736#GCACGGTTCGGCAATT-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              1                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927736#TCAGCTCTCTCTATTG-1
[1,]                              .
[2,]                              .
[3,]                              .
[4,]                              .
[5,]                              .


2023-11-13 04:22:41.310856 : Organizing SE PeakMatrix from ArrowFile : SRR13927736.arrow, 0.206 mins elapsed.

2023-11-13 04:22:41.597687 : se SRR13927736, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 3534 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(3534): SRR13927736#GCTGTTCTCCTTGACC-1
  SRR13927736#AGTTACGAGTAGGTCG-1 ... SRR13927736#TAAACCGCAAGTCTGT-1
  SRR13927736#TGCATTTGTTCTACCC-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-13 04:22:41.604352 : Completed PeakMatrix : SRR13927736(4 of 4), 0.877 mins elapsed.
2023-11-13 04:22:41.605227 : Organizing colData, 0.877 mins elapsed.
2023-11-13 04:22:41.742789 : Organizing rowData, 0.879 mins elapsed.
2023-11-13 04:22:41.751229 : Organizing rowRanges, 0.879 mins elapsed.
2023-11-13 04:22:41.762007 : Organizing Assays (1 of 1), 0.88 mins elapsed.
2023-11-13 04:22:46.144614 : Constructing SummarizedExperiment, 0.953 mins elapsed.
2023-11-13 04:22:59.664671 : Finished Matrix Creation, 1.178 mins elapsed.

2023-11-13 04:29:03.848744 : getMatrixFromProject Input-Parameters, Class = list

getMatrixFromProject Input-Parameters$ArchRProj: length = 1

getMatrixFromProject Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromProject Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromProject Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromProject Input-Parameters$verbose: length = 1
[1] TRUE


getMatrixFromProject Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromProject Input-Parameters$threads: length = 1
[1] 1


getMatrixFromProject Input-Parameters$logFile: length = 1
[1] "x"


2023-11-13 04:29:05.110051 : Reading GeneIntegrationMatrix : SRR13927737(1 of 4), 0.021 mins elapsed.

2023-11-13 04:29:05.140921 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927737 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927737.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 4439
[1] "SRR13927737#GTCACGGGTCAACAGG-1" "SRR13927737#TGCCTGTTCTTACTCA-1"
[3] "SRR13927737#TTGCGGGTCAGGTCTA-1" "SRR13927737#GGTTGCGCAAAGCTTC-1"
[5] "SRR13927737#CTTAATCTCGCGCTGA-1" "SRR13927737#TGTACGATCGGAAAGC-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-13 04:29:05.615549 : featureDF SRR13927737, Class = DFrame
DataFrame with 23525 rows and 6 columns
      seqnames     start       end    strand        name       idx
         <Rle> <integer> <integer> <integer> <character> <integer>
1         chr1     69091     70008         1       OR4F5         1
2         chr1    817371    819837         1      FAM87B         2
3         chr1    825138    859446         1   LINC01128         3
4         chr1    827522    826206         2   LINC00115         4
5         chr1    876903    868071         2      FAM41C         5
...        ...       ...       ...       ...         ...       ...
23521     chrX 155264589 155258241         2      RAB39B       923
23522     chrX 155334657 155276211         2       CLIC2       924
23523     chrX 155380787 155381134         1      H2AFB1       925
23524     chrX 155382115 155383230         1        F8A2       926
23525     chrX 155669944 155490115         2       TMLHE       927

2023-11-13 04:29:05.656595 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927737.arrow, 0.009 mins elapsed.

2023-11-13 04:29:41.440584 : mat SRR13927737, Class = dgCMatrix
mat SRR13927737: nRows = 23525, nCols = 4439
mat SRR13927737: NonZeroEntries = 50277774, EntryRange = [ 0.000148130054904683 , 376.345870564644 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927737#GTCACGGGTCAACAGG-1 SRR13927737#TGCCTGTTCTTACTCA-1
[1,]                      .                             .         
[2,]                      .                             .         
[3,]                      0.1057576                     0.09630696
[4,]                      .                             0.01806508
[5,]                      .                             .         
     SRR13927737#TTGCGGGTCAGGTCTA-1 SRR13927737#GGTTGCGCAAAGCTTC-1
[1,]                     .                              .         
[2,]                     .                              .         
[3,]                     0.09113499                     0.08396973
[4,]                     0.02096304                     0.01727603
[5,]                     .                              .         
     SRR13927737#CTTAATCTCGCGCTGA-1
[1,]                     .         
[2,]                     .         
[3,]                     0.08327649
[4,]                     0.01368917
[5,]                     .         


2023-11-13 04:29:41.905041 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927737.arrow, 0.613 mins elapsed.

2023-11-13 04:29:42.306553 : se SRR13927737, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 4439 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(4439): SRR13927737#GTCACGGGTCAACAGG-1
  SRR13927737#TGCCTGTTCTTACTCA-1 ... SRR13927737#TTGCTATAGATGTTCC-1
  SRR13927737#GCTCACTAGTTAGCGG-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-13 04:29:42.312463 : Completed GeneIntegrationMatrix : SRR13927737(1 of 4), 0.641 mins elapsed.
2023-11-13 04:29:42.312849 : Reading GeneIntegrationMatrix : SRR13927738(2 of 4), 0.641 mins elapsed.

2023-11-13 04:29:42.324679 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927738 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927738.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3736
[1] "SRR13927738#TAGCATGAGGCTAAAT-1" "SRR13927738#TGTTAGGTCGGTCAGC-1"
[3] "SRR13927738#CTCTCAGCATGTCCCT-1" "SRR13927738#ACTATTCGTCCAACCG-1"
[5] "SRR13927738#TGGTCCTTCTACTGCC-1" "SRR13927738#ATTGTCTAGACTAGGC-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-13 04:29:42.600272 : featureDF SRR13927738, Class = DFrame
DataFrame with 23525 rows and 6 columns
      seqnames     start       end    strand        name       idx
         <Rle> <integer> <integer> <integer> <character> <integer>
1         chr1     69091     70008         1       OR4F5         1
2         chr1    817371    819837         1      FAM87B         2
3         chr1    825138    859446         1   LINC01128         3
4         chr1    827522    826206         2   LINC00115         4
5         chr1    876903    868071         2      FAM41C         5
...        ...       ...       ...       ...         ...       ...
23521     chrX 155264589 155258241         2      RAB39B       923
23522     chrX 155334657 155276211         2       CLIC2       924
23523     chrX 155380787 155381134         1      H2AFB1       925
23524     chrX 155382115 155383230         1        F8A2       926
23525     chrX 155669944 155490115         2       TMLHE       927

2023-11-13 04:29:42.616799 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927738.arrow, 0.005 mins elapsed.

2023-11-13 04:29:58.449424 : mat SRR13927738, Class = dgCMatrix
mat SRR13927738: nRows = 23525, nCols = 3736
mat SRR13927738: NonZeroEntries = 41898897, EntryRange = [ 9.73975135001421e-05 , 758.763031303812 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927738#TAGCATGAGGCTAAAT-1 SRR13927738#TGTTAGGTCGGTCAGC-1
[1,]                    .                                .        
[2,]                    .                                .        
[3,]                    0.002346352                      0.0084575
[4,]                    0.005836413                      .        
[5,]                    .                                .        
     SRR13927738#CTCTCAGCATGTCCCT-1 SRR13927738#ACTATTCGTCCAACCG-1
[1,]                     .                             .          
[2,]                     .                             .          
[3,]                     .                             0.004081077
[4,]                     0.02001583                    .          
[5,]                     .                             .          
     SRR13927738#TGGTCCTTCTACTGCC-1
[1,]                    .          
[2,]                    .          
[3,]                    0.007640877
[4,]                    0.015048796
[5,]                    .          


2023-11-13 04:29:58.774061 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927738.arrow, 0.274 mins elapsed.

2023-11-13 04:29:59.085448 : se SRR13927738, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 3736 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(3736): SRR13927738#TAGCATGAGGCTAAAT-1
  SRR13927738#TGTTAGGTCGGTCAGC-1 ... SRR13927738#GCTCCTATCCGTGCGA-1
  SRR13927738#GGCACGTAGTCGATAA-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-13 04:29:59.090318 : Completed GeneIntegrationMatrix : SRR13927738(2 of 4), 0.921 mins elapsed.
2023-11-13 04:29:59.09071 : Reading GeneIntegrationMatrix : SRR13927735(3 of 4), 0.921 mins elapsed.

2023-11-13 04:29:59.10399 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927735 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927735.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3813
[1] "SRR13927735#TTATGTCTCCAGGTAT-1" "SRR13927735#TATTGCTCATCAGAAA-1"
[3] "SRR13927735#TTCGATTGTAGGGTTG-1" "SRR13927735#CATTCATTCGGATGTT-1"
[5] "SRR13927735#ACGTTAGGTCAACTGT-1" "SRR13927735#AAATGCCCAGCAATGG-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-13 04:29:59.319517 : featureDF SRR13927735, Class = DFrame
DataFrame with 23525 rows and 6 columns
      seqnames     start       end    strand        name       idx
         <Rle> <integer> <integer> <integer> <character> <integer>
1         chr1     69091     70008         1       OR4F5         1
2         chr1    817371    819837         1      FAM87B         2
3         chr1    825138    859446         1   LINC01128         3
4         chr1    827522    826206         2   LINC00115         4
5         chr1    876903    868071         2      FAM41C         5
...        ...       ...       ...       ...         ...       ...
23521     chrX 155264589 155258241         2      RAB39B       923
23522     chrX 155334657 155276211         2       CLIC2       924
23523     chrX 155380787 155381134         1      H2AFB1       925
23524     chrX 155382115 155383230         1        F8A2       926
23525     chrX 155669944 155490115         2       TMLHE       927

2023-11-13 04:29:59.332984 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927735.arrow, 0.004 mins elapsed.

2023-11-13 04:30:15.757064 : mat SRR13927735, Class = dgCMatrix
mat SRR13927735: nRows = 23525, nCols = 3813
mat SRR13927735: NonZeroEntries = 45535634, EntryRange = [ 7.20715653308254e-05 , 673.990433360078 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927735#TTATGTCTCCAGGTAT-1 SRR13927735#TATTGCTCATCAGAAA-1
[1,]                    .                              .          
[2,]                    .                              .          
[3,]                    0.070892732                    0.057009268
[4,]                    0.008800039                    0.002978408
[5,]                    .                              .          
     SRR13927735#TTCGATTGTAGGGTTG-1 SRR13927735#CATTCATTCGGATGTT-1
[1,]                     .                               .        
[2,]                     .                               .        
[3,]                     0.05523011                      0.0248301
[4,]                     .                               .        
[5,]                     .                               .        
     SRR13927735#ACGTTAGGTCAACTGT-1
[1,]                     .         
[2,]                     .         
[3,]                     0.05384571
[4,]                     0.01228688
[5,]                     .         


2023-11-13 04:30:16.147089 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927735.arrow, 0.284 mins elapsed.

2023-11-13 04:30:16.440343 : se SRR13927735, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 3813 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(3813): SRR13927735#TTATGTCTCCAGGTAT-1
  SRR13927735#TATTGCTCATCAGAAA-1 ... SRR13927735#GATCATGGTGTTGTTG-1
  SRR13927735#GTCCATCCATTTCACT-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-13 04:30:16.446122 : Completed GeneIntegrationMatrix : SRR13927735(3 of 4), 1.21 mins elapsed.
2023-11-13 04:30:16.44652 : Reading GeneIntegrationMatrix : SRR13927736(4 of 4), 1.21 mins elapsed.

2023-11-13 04:30:16.458937 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927736 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927736.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3534
[1] "SRR13927736#GCTGTTCTCCTTGACC-1" "SRR13927736#AGTTACGAGTAGGTCG-1"
[3] "SRR13927736#AGTTTGGAGCATTCCA-1" "SRR13927736#GCACGGTTCGGCAATT-1"
[5] "SRR13927736#TCAGCTCTCTCTATTG-1" "SRR13927736#CCTTAATTCTACTTTG-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-13 04:30:16.775187 : featureDF SRR13927736, Class = DFrame
DataFrame with 23525 rows and 6 columns
      seqnames     start       end    strand        name       idx
         <Rle> <integer> <integer> <integer> <character> <integer>
1         chr1     69091     70008         1       OR4F5         1
2         chr1    817371    819837         1      FAM87B         2
3         chr1    825138    859446         1   LINC01128         3
4         chr1    827522    826206         2   LINC00115         4
5         chr1    876903    868071         2      FAM41C         5
...        ...       ...       ...       ...         ...       ...
23521     chrX 155264589 155258241         2      RAB39B       923
23522     chrX 155334657 155276211         2       CLIC2       924
23523     chrX 155380787 155381134         1      H2AFB1       925
23524     chrX 155382115 155383230         1        F8A2       926
23525     chrX 155669944 155490115         2       TMLHE       927

2023-11-13 04:30:16.788753 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927736.arrow, 0.005 mins elapsed.

2023-11-13 04:30:29.648088 : mat SRR13927736, Class = dgCMatrix
mat SRR13927736: nRows = 23525, nCols = 3534
mat SRR13927736: NonZeroEntries = 39003616, EntryRange = [ 8.77627031078293e-05 , 770.387276519725 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927736#GCTGTTCTCCTTGACC-1 SRR13927736#AGTTACGAGTAGGTCG-1
[1,]                     .                             .          
[2,]                     .                             .          
[3,]                     0.07845381                    0.028344539
[4,]                     0.01395289                    0.008499189
[5,]                     .                             .          
     SRR13927736#AGTTTGGAGCATTCCA-1 SRR13927736#GCACGGTTCGGCAATT-1
[1,]                     .                              .         
[2,]                     .                              .         
[3,]                     0.06535130                     0.06576222
[4,]                     0.01810177                     0.01315728
[5,]                     .                              .         
     SRR13927736#TCAGCTCTCTCTATTG-1
[1,]                    .          
[2,]                    .          
[3,]                    0.016229526
[4,]                    0.003862592
[5,]                    .          


2023-11-13 04:30:29.942569 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927736.arrow, 0.225 mins elapsed.

2023-11-13 04:30:30.209999 : se SRR13927736, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 3534 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(3534): SRR13927736#GCTGTTCTCCTTGACC-1
  SRR13927736#AGTTACGAGTAGGTCG-1 ... SRR13927736#TAAACCGCAAGTCTGT-1
  SRR13927736#TGCATTTGTTCTACCC-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-13 04:30:30.214782 : Completed GeneIntegrationMatrix : SRR13927736(4 of 4), 1.439 mins elapsed.
2023-11-13 04:30:30.215711 : Organizing colData, 1.439 mins elapsed.
2023-11-13 04:30:30.363652 : Organizing rowData, 1.442 mins elapsed.
2023-11-13 04:30:30.366734 : Organizing rowRanges, 1.442 mins elapsed.
2023-11-13 04:30:30.372027 : Organizing Assays (1 of 1), 1.442 mins elapsed.
2023-11-13 04:30:35.033452 : Constructing SummarizedExperiment, 1.52 mins elapsed.
2023-11-13 04:30:36.294945 : Finished Matrix Creation, 1.541 mins elapsed.

2023-11-13 04:30:36.297351 : getMatrixFromProject Input-Parameters, Class = list

getMatrixFromProject Input-Parameters$ArchRProj: length = 1

getMatrixFromProject Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromProject Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromProject Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromProject Input-Parameters$verbose: length = 1
[1] TRUE


getMatrixFromProject Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromProject Input-Parameters$threads: length = 1
[1] 1


getMatrixFromProject Input-Parameters$logFile: length = 1
[1] "x"


2023-11-13 04:30:37.409393 : Reading PeakMatrix : SRR13927737(1 of 4), 0.019 mins elapsed.

2023-11-13 04:30:37.433427 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927737 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927737.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 4439
[1] "SRR13927737#GTCACGGGTCAACAGG-1" "SRR13927737#TGCCTGTTCTTACTCA-1"
[3] "SRR13927737#TTGCGGGTCAGGTCTA-1" "SRR13927737#GGTTGCGCAAAGCTTC-1"
[5] "SRR13927737#CTTAATCTCGCGCTGA-1" "SRR13927737#TGTACGATCGGAAAGC-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-13 04:30:38.067569 : featureDF SRR13927737, Class = DFrame
DataFrame with 210721 rows and 4 columns
       seqnames       idx     start       end
          <Rle> <integer> <integer> <integer>
1          chr1         1    804677    805177
2          chr1         2    811866    812366
3          chr1         3    817090    817590
4          chr1         4    818831    819331
5          chr1         5    826543    827043
...         ...       ...       ...       ...
210717     chrX      5777 155898525 155899025
210718     chrX      5778 155899804 155900304
210719     chrX      5779 155956608 155957108
210720     chrX      5780 155966827 155967327
210721     chrX      5781 155971625 155972125

2023-11-13 04:30:38.140647 : Getting PeakMatrix from ArrowFile : SRR13927737.arrow, 0.012 mins elapsed.

2023-11-13 04:30:48.753088 : mat SRR13927737, Class = dgCMatrix
mat SRR13927737: nRows = 210721, nCols = 4439
mat SRR13927737: NonZeroEntries = 29963541, EntryRange = [ 1 , 4 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927737#GTCACGGGTCAACAGG-1 SRR13927737#TGCCTGTTCTTACTCA-1
[1,]                              .                              .
[2,]                              .                              2
[3,]                              .                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927737#TTGCGGGTCAGGTCTA-1 SRR13927737#GGTTGCGCAAAGCTTC-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              2                              .
[4,]                              .                              .
[5,]                              4                              .
     SRR13927737#CTTAATCTCGCGCTGA-1
[1,]                              .
[2,]                              .
[3,]                              .
[4,]                              .
[5,]                              .


2023-11-13 04:30:48.978143 : Organizing SE PeakMatrix from ArrowFile : SRR13927737.arrow, 0.192 mins elapsed.

2023-11-13 04:30:49.274369 : se SRR13927737, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 4439 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(4439): SRR13927737#GTCACGGGTCAACAGG-1
  SRR13927737#TGCCTGTTCTTACTCA-1 ... SRR13927737#TTGCTATAGATGTTCC-1
  SRR13927737#GCTCACTAGTTAGCGG-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-13 04:30:49.28407 : Completed PeakMatrix : SRR13927737(1 of 4), 0.216 mins elapsed.
2023-11-13 04:30:49.284514 : Reading PeakMatrix : SRR13927738(2 of 4), 0.216 mins elapsed.

2023-11-13 04:30:49.297028 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927738 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927738.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3736
[1] "SRR13927738#TAGCATGAGGCTAAAT-1" "SRR13927738#TGTTAGGTCGGTCAGC-1"
[3] "SRR13927738#CTCTCAGCATGTCCCT-1" "SRR13927738#ACTATTCGTCCAACCG-1"
[5] "SRR13927738#TGGTCCTTCTACTGCC-1" "SRR13927738#ATTGTCTAGACTAGGC-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-13 04:30:49.940258 : featureDF SRR13927738, Class = DFrame
DataFrame with 210721 rows and 4 columns
       seqnames       idx     start       end
          <Rle> <integer> <integer> <integer>
1          chr1         1    804677    805177
2          chr1         2    811866    812366
3          chr1         3    817090    817590
4          chr1         4    818831    819331
5          chr1         5    826543    827043
...         ...       ...       ...       ...
210717     chrX      5777 155898525 155899025
210718     chrX      5778 155899804 155900304
210719     chrX      5779 155956608 155957108
210720     chrX      5780 155966827 155967327
210721     chrX      5781 155971625 155972125

2023-11-13 04:30:49.992424 : Getting PeakMatrix from ArrowFile : SRR13927738.arrow, 0.012 mins elapsed.

2023-11-13 04:30:59.998836 : mat SRR13927738, Class = dgCMatrix
mat SRR13927738: nRows = 210721, nCols = 3736
mat SRR13927738: NonZeroEntries = 30165921, EntryRange = [ 1 , 4 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927738#TAGCATGAGGCTAAAT-1 SRR13927738#TGTTAGGTCGGTCAGC-1
[1,]                              2                              .
[2,]                              .                              .
[3,]                              .                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927738#CTCTCAGCATGTCCCT-1 SRR13927738#ACTATTCGTCCAACCG-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              .                              .
[4,]                              .                              .
[5,]                              2                              .
     SRR13927738#TGGTCCTTCTACTGCC-1
[1,]                              .
[2,]                              .
[3,]                              2
[4,]                              .
[5,]                              .


2023-11-13 04:31:00.225322 : Organizing SE PeakMatrix from ArrowFile : SRR13927738.arrow, 0.182 mins elapsed.

2023-11-13 04:31:00.515602 : se SRR13927738, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 3736 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(3736): SRR13927738#TAGCATGAGGCTAAAT-1
  SRR13927738#TGTTAGGTCGGTCAGC-1 ... SRR13927738#GCTCCTATCCGTGCGA-1
  SRR13927738#GGCACGTAGTCGATAA-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-13 04:31:00.522265 : Completed PeakMatrix : SRR13927738(2 of 4), 0.404 mins elapsed.
2023-11-13 04:31:00.522656 : Reading PeakMatrix : SRR13927735(3 of 4), 0.404 mins elapsed.

2023-11-13 04:31:00.534678 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927735 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927735.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3813
[1] "SRR13927735#TTATGTCTCCAGGTAT-1" "SRR13927735#TATTGCTCATCAGAAA-1"
[3] "SRR13927735#TTCGATTGTAGGGTTG-1" "SRR13927735#CATTCATTCGGATGTT-1"
[5] "SRR13927735#ACGTTAGGTCAACTGT-1" "SRR13927735#AAATGCCCAGCAATGG-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-13 04:31:00.888191 : featureDF SRR13927735, Class = DFrame
DataFrame with 210721 rows and 4 columns
       seqnames       idx     start       end
          <Rle> <integer> <integer> <integer>
1          chr1         1    804677    805177
2          chr1         2    811866    812366
3          chr1         3    817090    817590
4          chr1         4    818831    819331
5          chr1         5    826543    827043
...         ...       ...       ...       ...
210717     chrX      5777 155898525 155899025
210718     chrX      5778 155899804 155900304
210719     chrX      5779 155956608 155957108
210720     chrX      5780 155966827 155967327
210721     chrX      5781 155971625 155972125

2023-11-13 04:31:00.927737 : Getting PeakMatrix from ArrowFile : SRR13927735.arrow, 0.007 mins elapsed.

2023-11-13 04:31:09.614441 : mat SRR13927735, Class = dgCMatrix
mat SRR13927735: nRows = 210721, nCols = 3813
mat SRR13927735: NonZeroEntries = 26038084, EntryRange = [ 1 , 4 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927735#TTATGTCTCCAGGTAT-1 SRR13927735#TATTGCTCATCAGAAA-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              1                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927735#TTCGATTGTAGGGTTG-1 SRR13927735#CATTCATTCGGATGTT-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              .                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927735#ACGTTAGGTCAACTGT-1
[1,]                              2
[2,]                              .
[3,]                              .
[4,]                              .
[5,]                              .


2023-11-13 04:31:09.810177 : Organizing SE PeakMatrix from ArrowFile : SRR13927735.arrow, 0.155 mins elapsed.

2023-11-13 04:31:10.098029 : se SRR13927735, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 3813 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(3813): SRR13927735#TTATGTCTCCAGGTAT-1
  SRR13927735#TATTGCTCATCAGAAA-1 ... SRR13927735#GATCATGGTGTTGTTG-1
  SRR13927735#GTCCATCCATTTCACT-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-13 04:31:10.104911 : Completed PeakMatrix : SRR13927735(3 of 4), 0.563 mins elapsed.
2023-11-13 04:31:10.105299 : Reading PeakMatrix : SRR13927736(4 of 4), 0.563 mins elapsed.

2023-11-13 04:31:10.11905 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927736 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927736.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3534
[1] "SRR13927736#GCTGTTCTCCTTGACC-1" "SRR13927736#AGTTACGAGTAGGTCG-1"
[3] "SRR13927736#AGTTTGGAGCATTCCA-1" "SRR13927736#GCACGGTTCGGCAATT-1"
[5] "SRR13927736#TCAGCTCTCTCTATTG-1" "SRR13927736#CCTTAATTCTACTTTG-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-13 04:31:10.763555 : featureDF SRR13927736, Class = DFrame
DataFrame with 210721 rows and 4 columns
       seqnames       idx     start       end
          <Rle> <integer> <integer> <integer>
1          chr1         1    804677    805177
2          chr1         2    811866    812366
3          chr1         3    817090    817590
4          chr1         4    818831    819331
5          chr1         5    826543    827043
...         ...       ...       ...       ...
210717     chrX      5777 155898525 155899025
210718     chrX      5778 155899804 155900304
210719     chrX      5779 155956608 155957108
210720     chrX      5780 155966827 155967327
210721     chrX      5781 155971625 155972125

2023-11-13 04:31:10.815126 : Getting PeakMatrix from ArrowFile : SRR13927736.arrow, 0.012 mins elapsed.

2023-11-13 04:31:19.311929 : mat SRR13927736, Class = dgCMatrix
mat SRR13927736: nRows = 210721, nCols = 3534
mat SRR13927736: NonZeroEntries = 25583439, EntryRange = [ 1 , 4 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927736#GCTGTTCTCCTTGACC-1 SRR13927736#AGTTACGAGTAGGTCG-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              .                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927736#AGTTTGGAGCATTCCA-1 SRR13927736#GCACGGTTCGGCAATT-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              1                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927736#TCAGCTCTCTCTATTG-1
[1,]                              .
[2,]                              .
[3,]                              .
[4,]                              .
[5,]                              .


2023-11-13 04:31:19.505013 : Organizing SE PeakMatrix from ArrowFile : SRR13927736.arrow, 0.156 mins elapsed.

2023-11-13 04:31:19.785132 : se SRR13927736, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 3534 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(3534): SRR13927736#GCTGTTCTCCTTGACC-1
  SRR13927736#AGTTACGAGTAGGTCG-1 ... SRR13927736#TAAACCGCAAGTCTGT-1
  SRR13927736#TGCATTTGTTCTACCC-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-13 04:31:19.79189 : Completed PeakMatrix : SRR13927736(4 of 4), 0.725 mins elapsed.
2023-11-13 04:31:19.792737 : Organizing colData, 0.725 mins elapsed.
2023-11-13 04:31:19.92461 : Organizing rowData, 0.727 mins elapsed.
2023-11-13 04:31:19.932677 : Organizing rowRanges, 0.727 mins elapsed.
2023-11-13 04:31:19.943039 : Organizing Assays (1 of 1), 0.727 mins elapsed.
2023-11-13 04:31:22.6196 : Constructing SummarizedExperiment, 0.772 mins elapsed.
2023-11-13 04:31:36.475508 : Finished Matrix Creation, 1.003 mins elapsed.

2023-11-13 04:53:33.052626 : getMatrixFromProject Input-Parameters, Class = list

getMatrixFromProject Input-Parameters$ArchRProj: length = 1

getMatrixFromProject Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromProject Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromProject Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromProject Input-Parameters$verbose: length = 1
[1] TRUE


getMatrixFromProject Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromProject Input-Parameters$threads: length = 1
[1] 1


getMatrixFromProject Input-Parameters$logFile: length = 1
[1] "x"


2023-11-13 04:53:47.138391 : Reading GeneIntegrationMatrix : SRR13927737(1 of 4), 0.235 mins elapsed.

2023-11-13 04:53:47.403928 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927737 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927737.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 4439
[1] "SRR13927737#GTCACGGGTCAACAGG-1" "SRR13927737#TGCCTGTTCTTACTCA-1" "SRR13927737#TTGCGGGTCAGGTCTA-1" "SRR13927737#GGTTGCGCAAAGCTTC-1"
[5] "SRR13927737#CTTAATCTCGCGCTGA-1" "SRR13927737#TGTACGATCGGAAAGC-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-13 04:53:48.30107 : featureDF SRR13927737, Class = DFrame
DataFrame with 23525 rows and 6 columns
      seqnames     start       end    strand        name       idx
         <Rle> <integer> <integer> <integer> <character> <integer>
1         chr1     69091     70008         1       OR4F5         1
2         chr1    817371    819837         1      FAM87B         2
3         chr1    825138    859446         1   LINC01128         3
4         chr1    827522    826206         2   LINC00115         4
5         chr1    876903    868071         2      FAM41C         5
...        ...       ...       ...       ...         ...       ...
23521     chrX 155264589 155258241         2      RAB39B       923
23522     chrX 155334657 155276211         2       CLIC2       924
23523     chrX 155380787 155381134         1      H2AFB1       925
23524     chrX 155382115 155383230         1        F8A2       926
23525     chrX 155669944 155490115         2       TMLHE       927

2023-11-13 04:53:48.364951 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927737.arrow, 0.016 mins elapsed.

2023-11-13 04:54:33.620575 : mat SRR13927737, Class = dgCMatrix
mat SRR13927737: nRows = 23525, nCols = 4439
mat SRR13927737: NonZeroEntries = 50277774, EntryRange = [ 0.000148130054904683 , 376.345870564644 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927737#GTCACGGGTCAACAGG-1 SRR13927737#TGCCTGTTCTTACTCA-1 SRR13927737#TTGCGGGTCAGGTCTA-1 SRR13927737#GGTTGCGCAAAGCTTC-1
[1,]                      .                             .                              .                              .         
[2,]                      .                             .                              .                              .         
[3,]                      0.1057576                     0.09630696                     0.09113499                     0.08396973
[4,]                      .                             0.01806508                     0.02096304                     0.01727603
[5,]                      .                             .                              .                              .         
     SRR13927737#CTTAATCTCGCGCTGA-1
[1,]                     .         
[2,]                     .         
[3,]                     0.08327649
[4,]                     0.01368917
[5,]                     .         


2023-11-13 04:54:34.229988 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927737.arrow, 0.78 mins elapsed.

2023-11-13 04:54:34.805047 : se SRR13927737, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 4439 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(4439): SRR13927737#GTCACGGGTCAACAGG-1 SRR13927737#TGCCTGTTCTTACTCA-1 ... SRR13927737#TTGCTATAGATGTTCC-1
  SRR13927737#GCTCACTAGTTAGCGG-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-13 04:54:34.812885 : Completed GeneIntegrationMatrix : SRR13927737(1 of 4), 1.029 mins elapsed.
2023-11-13 04:54:34.81344 : Reading GeneIntegrationMatrix : SRR13927738(2 of 4), 1.029 mins elapsed.

2023-11-13 04:54:34.828893 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927738 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927738.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3736
[1] "SRR13927738#TAGCATGAGGCTAAAT-1" "SRR13927738#TGTTAGGTCGGTCAGC-1" "SRR13927738#CTCTCAGCATGTCCCT-1" "SRR13927738#ACTATTCGTCCAACCG-1"
[5] "SRR13927738#TGGTCCTTCTACTGCC-1" "SRR13927738#ATTGTCTAGACTAGGC-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-13 04:54:35.259071 : featureDF SRR13927738, Class = DFrame
DataFrame with 23525 rows and 6 columns
      seqnames     start       end    strand        name       idx
         <Rle> <integer> <integer> <integer> <character> <integer>
1         chr1     69091     70008         1       OR4F5         1
2         chr1    817371    819837         1      FAM87B         2
3         chr1    825138    859446         1   LINC01128         3
4         chr1    827522    826206         2   LINC00115         4
5         chr1    876903    868071         2      FAM41C         5
...        ...       ...       ...       ...         ...       ...
23521     chrX 155264589 155258241         2      RAB39B       923
23522     chrX 155334657 155276211         2       CLIC2       924
23523     chrX 155380787 155381134         1      H2AFB1       925
23524     chrX 155382115 155383230         1        F8A2       926
23525     chrX 155669944 155490115         2       TMLHE       927

2023-11-13 04:54:35.278524 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927738.arrow, 0.007 mins elapsed.

2023-11-13 04:54:55.129955 : mat SRR13927738, Class = dgCMatrix
mat SRR13927738: nRows = 23525, nCols = 3736
mat SRR13927738: NonZeroEntries = 41898897, EntryRange = [ 9.73975135001421e-05 , 758.763031303812 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927738#TAGCATGAGGCTAAAT-1 SRR13927738#TGTTAGGTCGGTCAGC-1 SRR13927738#CTCTCAGCATGTCCCT-1 SRR13927738#ACTATTCGTCCAACCG-1
[1,]                    .                                .                             .                             .          
[2,]                    .                                .                             .                             .          
[3,]                    0.002346352                      0.0084575                     .                             0.004081077
[4,]                    0.005836413                      .                             0.02001583                    .          
[5,]                    .                                .                             .                             .          
     SRR13927738#TGGTCCTTCTACTGCC-1
[1,]                    .          
[2,]                    .          
[3,]                    0.007640877
[4,]                    0.015048796
[5,]                    .          


2023-11-13 04:54:56.93845 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927738.arrow, 0.368 mins elapsed.

2023-11-13 04:54:57.473527 : se SRR13927738, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 3736 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(3736): SRR13927738#TAGCATGAGGCTAAAT-1 SRR13927738#TGTTAGGTCGGTCAGC-1 ... SRR13927738#GCTCCTATCCGTGCGA-1
  SRR13927738#GGCACGTAGTCGATAA-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-13 04:54:57.479964 : Completed GeneIntegrationMatrix : SRR13927738(2 of 4), 1.407 mins elapsed.
2023-11-13 04:54:57.480519 : Reading GeneIntegrationMatrix : SRR13927735(3 of 4), 1.407 mins elapsed.

2023-11-13 04:54:57.502318 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927735 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927735.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3813
[1] "SRR13927735#TTATGTCTCCAGGTAT-1" "SRR13927735#TATTGCTCATCAGAAA-1" "SRR13927735#TTCGATTGTAGGGTTG-1" "SRR13927735#CATTCATTCGGATGTT-1"
[5] "SRR13927735#ACGTTAGGTCAACTGT-1" "SRR13927735#AAATGCCCAGCAATGG-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-13 04:54:58.040524 : featureDF SRR13927735, Class = DFrame
DataFrame with 23525 rows and 6 columns
      seqnames     start       end    strand        name       idx
         <Rle> <integer> <integer> <integer> <character> <integer>
1         chr1     69091     70008         1       OR4F5         1
2         chr1    817371    819837         1      FAM87B         2
3         chr1    825138    859446         1   LINC01128         3
4         chr1    827522    826206         2   LINC00115         4
5         chr1    876903    868071         2      FAM41C         5
...        ...       ...       ...       ...         ...       ...
23521     chrX 155264589 155258241         2      RAB39B       923
23522     chrX 155334657 155276211         2       CLIC2       924
23523     chrX 155380787 155381134         1      H2AFB1       925
23524     chrX 155382115 155383230         1        F8A2       926
23525     chrX 155669944 155490115         2       TMLHE       927

2023-11-13 04:54:58.061785 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927735.arrow, 0.009 mins elapsed.

2023-11-13 04:55:19.162702 : mat SRR13927735, Class = dgCMatrix
mat SRR13927735: nRows = 23525, nCols = 3813
mat SRR13927735: NonZeroEntries = 45535634, EntryRange = [ 7.20715653308254e-05 , 673.990433360078 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927735#TTATGTCTCCAGGTAT-1 SRR13927735#TATTGCTCATCAGAAA-1 SRR13927735#TTCGATTGTAGGGTTG-1 SRR13927735#CATTCATTCGGATGTT-1
[1,]                    .                              .                               .                               .        
[2,]                    .                              .                               .                               .        
[3,]                    0.070892732                    0.057009268                     0.05523011                      0.0248301
[4,]                    0.008800039                    0.002978408                     .                               .        
[5,]                    .                              .                               .                               .        
     SRR13927735#ACGTTAGGTCAACTGT-1
[1,]                     .         
[2,]                     .         
[3,]                     0.05384571
[4,]                     0.01228688
[5,]                     .         


2023-11-13 04:55:19.802354 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927735.arrow, 0.372 mins elapsed.

2023-11-13 04:55:20.562409 : se SRR13927735, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 3813 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(3813): SRR13927735#TTATGTCTCCAGGTAT-1 SRR13927735#TATTGCTCATCAGAAA-1 ... SRR13927735#GATCATGGTGTTGTTG-1
  SRR13927735#GTCCATCCATTTCACT-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-13 04:55:20.570971 : Completed GeneIntegrationMatrix : SRR13927735(3 of 4), 1.792 mins elapsed.
2023-11-13 04:55:20.571712 : Reading GeneIntegrationMatrix : SRR13927736(4 of 4), 1.792 mins elapsed.

2023-11-13 04:55:20.598906 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927736 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927736.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3534
[1] "SRR13927736#GCTGTTCTCCTTGACC-1" "SRR13927736#AGTTACGAGTAGGTCG-1" "SRR13927736#AGTTTGGAGCATTCCA-1" "SRR13927736#GCACGGTTCGGCAATT-1"
[5] "SRR13927736#TCAGCTCTCTCTATTG-1" "SRR13927736#CCTTAATTCTACTTTG-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-13 04:55:21.344563 : featureDF SRR13927736, Class = DFrame
DataFrame with 23525 rows and 6 columns
      seqnames     start       end    strand        name       idx
         <Rle> <integer> <integer> <integer> <character> <integer>
1         chr1     69091     70008         1       OR4F5         1
2         chr1    817371    819837         1      FAM87B         2
3         chr1    825138    859446         1   LINC01128         3
4         chr1    827522    826206         2   LINC00115         4
5         chr1    876903    868071         2      FAM41C         5
...        ...       ...       ...       ...         ...       ...
23521     chrX 155264589 155258241         2      RAB39B       923
23522     chrX 155334657 155276211         2       CLIC2       924
23523     chrX 155380787 155381134         1      H2AFB1       925
23524     chrX 155382115 155383230         1        F8A2       926
23525     chrX 155669944 155490115         2       TMLHE       927

2023-11-13 04:55:21.369198 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927736.arrow, 0.013 mins elapsed.

2023-11-13 04:55:38.648602 : mat SRR13927736, Class = dgCMatrix
mat SRR13927736: nRows = 23525, nCols = 3534
mat SRR13927736: NonZeroEntries = 39003616, EntryRange = [ 8.77627031078293e-05 , 770.387276519725 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927736#GCTGTTCTCCTTGACC-1 SRR13927736#AGTTACGAGTAGGTCG-1 SRR13927736#AGTTTGGAGCATTCCA-1 SRR13927736#GCACGGTTCGGCAATT-1
[1,]                     .                             .                               .                              .         
[2,]                     .                             .                               .                              .         
[3,]                     0.07845381                    0.028344539                     0.06535130                     0.06576222
[4,]                     0.01395289                    0.008499189                     0.01810177                     0.01315728
[5,]                     .                             .                               .                              .         
     SRR13927736#TCAGCTCTCTCTATTG-1
[1,]                    .          
[2,]                    .          
[3,]                    0.016229526
[4,]                    0.003862592
[5,]                    .          


2023-11-13 04:55:39.107729 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927736.arrow, 0.308 mins elapsed.

2023-11-13 04:55:39.558079 : se SRR13927736, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 3534 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(3534): SRR13927736#GCTGTTCTCCTTGACC-1 SRR13927736#AGTTACGAGTAGGTCG-1 ... SRR13927736#TAAACCGCAAGTCTGT-1
  SRR13927736#TGCATTTGTTCTACCC-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-13 04:55:39.564537 : Completed GeneIntegrationMatrix : SRR13927736(4 of 4), 2.109 mins elapsed.
2023-11-13 04:55:39.569239 : Organizing colData, 2.109 mins elapsed.
2023-11-13 04:55:39.90348 : Organizing rowData, 2.114 mins elapsed.
2023-11-13 04:55:39.907665 : Organizing rowRanges, 2.114 mins elapsed.
2023-11-13 04:55:39.916525 : Organizing Assays (1 of 1), 2.114 mins elapsed.
2023-11-13 04:55:48.165076 : Constructing SummarizedExperiment, 2.252 mins elapsed.
2023-11-13 04:55:49.944937 : Finished Matrix Creation, 2.276 mins elapsed.

2023-11-13 04:59:59.394513 : getMatrixFromProject Input-Parameters, Class = list

getMatrixFromProject Input-Parameters$ArchRProj: length = 1

getMatrixFromProject Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromProject Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromProject Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromProject Input-Parameters$verbose: length = 1
[1] TRUE


getMatrixFromProject Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromProject Input-Parameters$threads: length = 1
[1] 1


getMatrixFromProject Input-Parameters$logFile: length = 1
[1] "x"


2023-11-13 05:00:00.288777 : Reading PeakMatrix : SRR13927737(1 of 4), 0.015 mins elapsed.

2023-11-13 05:00:00.338587 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927737 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927737.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 4439
[1] "SRR13927737#GTCACGGGTCAACAGG-1" "SRR13927737#TGCCTGTTCTTACTCA-1" "SRR13927737#TTGCGGGTCAGGTCTA-1" "SRR13927737#GGTTGCGCAAAGCTTC-1"
[5] "SRR13927737#CTTAATCTCGCGCTGA-1" "SRR13927737#TGTACGATCGGAAAGC-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-13 05:00:01.134658 : featureDF SRR13927737, Class = DFrame
DataFrame with 210721 rows and 4 columns
       seqnames       idx     start       end
          <Rle> <integer> <integer> <integer>
1          chr1         1    804677    805177
2          chr1         2    811866    812366
3          chr1         3    817090    817590
4          chr1         4    818831    819331
5          chr1         5    826543    827043
...         ...       ...       ...       ...
210717     chrX      5777 155898525 155899025
210718     chrX      5778 155899804 155900304
210719     chrX      5779 155956608 155957108
210720     chrX      5780 155966827 155967327
210721     chrX      5781 155971625 155972125

2023-11-13 05:00:01.349745 : Getting PeakMatrix from ArrowFile : SRR13927737.arrow, 0.017 mins elapsed.

2023-11-13 05:00:13.604556 : mat SRR13927737, Class = dgCMatrix
mat SRR13927737: nRows = 210721, nCols = 4439
mat SRR13927737: NonZeroEntries = 29963541, EntryRange = [ 1 , 4 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927737#GTCACGGGTCAACAGG-1 SRR13927737#TGCCTGTTCTTACTCA-1 SRR13927737#TTGCGGGTCAGGTCTA-1 SRR13927737#GGTTGCGCAAAGCTTC-1
[1,]                              .                              .                              .                              .
[2,]                              .                              2                              .                              .
[3,]                              .                              .                              2                              .
[4,]                              .                              .                              .                              .
[5,]                              .                              .                              4                              .
     SRR13927737#CTTAATCTCGCGCTGA-1
[1,]                              .
[2,]                              .
[3,]                              .
[4,]                              .
[5,]                              .


2023-11-13 05:00:13.956895 : Organizing SE PeakMatrix from ArrowFile : SRR13927737.arrow, 0.227 mins elapsed.

2023-11-13 05:00:14.555851 : se SRR13927737, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 4439 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(4439): SRR13927737#GTCACGGGTCAACAGG-1 SRR13927737#TGCCTGTTCTTACTCA-1 ... SRR13927737#TTGCTATAGATGTTCC-1
  SRR13927737#GCTCACTAGTTAGCGG-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-13 05:00:14.569304 : Completed PeakMatrix : SRR13927737(1 of 4), 0.253 mins elapsed.
2023-11-13 05:00:14.569913 : Reading PeakMatrix : SRR13927738(2 of 4), 0.253 mins elapsed.

2023-11-13 05:00:14.592583 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927738 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927738.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3736
[1] "SRR13927738#TAGCATGAGGCTAAAT-1" "SRR13927738#TGTTAGGTCGGTCAGC-1" "SRR13927738#CTCTCAGCATGTCCCT-1" "SRR13927738#ACTATTCGTCCAACCG-1"
[5] "SRR13927738#TGGTCCTTCTACTGCC-1" "SRR13927738#ATTGTCTAGACTAGGC-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-13 05:00:15.202448 : featureDF SRR13927738, Class = DFrame
DataFrame with 210721 rows and 4 columns
       seqnames       idx     start       end
          <Rle> <integer> <integer> <integer>
1          chr1         1    804677    805177
2          chr1         2    811866    812366
3          chr1         3    817090    817590
4          chr1         4    818831    819331
5          chr1         5    826543    827043
...         ...       ...       ...       ...
210717     chrX      5777 155898525 155899025
210718     chrX      5778 155899804 155900304
210719     chrX      5779 155956608 155957108
210720     chrX      5780 155966827 155967327
210721     chrX      5781 155971625 155972125

2023-11-13 05:00:15.421838 : Getting PeakMatrix from ArrowFile : SRR13927738.arrow, 0.014 mins elapsed.

2023-11-13 05:00:27.230733 : mat SRR13927738, Class = dgCMatrix
mat SRR13927738: nRows = 210721, nCols = 3736
mat SRR13927738: NonZeroEntries = 30165921, EntryRange = [ 1 , 4 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927738#TAGCATGAGGCTAAAT-1 SRR13927738#TGTTAGGTCGGTCAGC-1 SRR13927738#CTCTCAGCATGTCCCT-1 SRR13927738#ACTATTCGTCCAACCG-1
[1,]                              2                              .                              .                              .
[2,]                              .                              .                              .                              .
[3,]                              .                              .                              .                              .
[4,]                              .                              .                              .                              .
[5,]                              .                              .                              2                              .
     SRR13927738#TGGTCCTTCTACTGCC-1
[1,]                              .
[2,]                              .
[3,]                              2
[4,]                              .
[5,]                              .


2023-11-13 05:00:27.803801 : Organizing SE PeakMatrix from ArrowFile : SRR13927738.arrow, 0.22 mins elapsed.

2023-11-13 05:00:28.382213 : se SRR13927738, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 3736 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(3736): SRR13927738#TAGCATGAGGCTAAAT-1 SRR13927738#TGTTAGGTCGGTCAGC-1 ... SRR13927738#GCTCCTATCCGTGCGA-1
  SRR13927738#GGCACGTAGTCGATAA-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-13 05:00:28.391443 : Completed PeakMatrix : SRR13927738(2 of 4), 0.483 mins elapsed.
2023-11-13 05:00:28.392034 : Reading PeakMatrix : SRR13927735(3 of 4), 0.483 mins elapsed.

2023-11-13 05:00:28.417381 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927735 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927735.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3813
[1] "SRR13927735#TTATGTCTCCAGGTAT-1" "SRR13927735#TATTGCTCATCAGAAA-1" "SRR13927735#TTCGATTGTAGGGTTG-1" "SRR13927735#CATTCATTCGGATGTT-1"
[5] "SRR13927735#ACGTTAGGTCAACTGT-1" "SRR13927735#AAATGCCCAGCAATGG-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-13 05:00:29.294695 : featureDF SRR13927735, Class = DFrame
DataFrame with 210721 rows and 4 columns
       seqnames       idx     start       end
          <Rle> <integer> <integer> <integer>
1          chr1         1    804677    805177
2          chr1         2    811866    812366
3          chr1         3    817090    817590
4          chr1         4    818831    819331
5          chr1         5    826543    827043
...         ...       ...       ...       ...
210717     chrX      5777 155898525 155899025
210718     chrX      5778 155899804 155900304
210719     chrX      5779 155956608 155957108
210720     chrX      5780 155966827 155967327
210721     chrX      5781 155971625 155972125

2023-11-13 05:00:29.388529 : Getting PeakMatrix from ArrowFile : SRR13927735.arrow, 0.016 mins elapsed.

2023-11-13 05:00:40.114232 : mat SRR13927735, Class = dgCMatrix
mat SRR13927735: nRows = 210721, nCols = 3813
mat SRR13927735: NonZeroEntries = 26038084, EntryRange = [ 1 , 4 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927735#TTATGTCTCCAGGTAT-1 SRR13927735#TATTGCTCATCAGAAA-1 SRR13927735#TTCGATTGTAGGGTTG-1 SRR13927735#CATTCATTCGGATGTT-1
[1,]                              .                              .                              .                              .
[2,]                              .                              .                              .                              .
[3,]                              1                              .                              .                              .
[4,]                              .                              .                              .                              .
[5,]                              .                              .                              .                              .
     SRR13927735#ACGTTAGGTCAACTGT-1
[1,]                              2
[2,]                              .
[3,]                              .
[4,]                              .
[5,]                              .


2023-11-13 05:00:40.27785 : Organizing SE PeakMatrix from ArrowFile : SRR13927735.arrow, 0.198 mins elapsed.

2023-11-13 05:00:40.918471 : se SRR13927735, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 3813 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(3813): SRR13927735#TTATGTCTCCAGGTAT-1 SRR13927735#TATTGCTCATCAGAAA-1 ... SRR13927735#GATCATGGTGTTGTTG-1
  SRR13927735#GTCCATCCATTTCACT-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-13 05:00:40.928127 : Completed PeakMatrix : SRR13927735(3 of 4), 0.692 mins elapsed.
2023-11-13 05:00:40.928679 : Reading PeakMatrix : SRR13927736(4 of 4), 0.692 mins elapsed.

2023-11-13 05:00:40.945851 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927736 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927736.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3534
[1] "SRR13927736#GCTGTTCTCCTTGACC-1" "SRR13927736#AGTTACGAGTAGGTCG-1" "SRR13927736#AGTTTGGAGCATTCCA-1" "SRR13927736#GCACGGTTCGGCAATT-1"
[5] "SRR13927736#TCAGCTCTCTCTATTG-1" "SRR13927736#CCTTAATTCTACTTTG-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-13 05:00:41.780992 : featureDF SRR13927736, Class = DFrame
DataFrame with 210721 rows and 4 columns
       seqnames       idx     start       end
          <Rle> <integer> <integer> <integer>
1          chr1         1    804677    805177
2          chr1         2    811866    812366
3          chr1         3    817090    817590
4          chr1         4    818831    819331
5          chr1         5    826543    827043
...         ...       ...       ...       ...
210717     chrX      5777 155898525 155899025
210718     chrX      5778 155899804 155900304
210719     chrX      5779 155956608 155957108
210720     chrX      5780 155966827 155967327
210721     chrX      5781 155971625 155972125

2023-11-13 05:00:41.863605 : Getting PeakMatrix from ArrowFile : SRR13927736.arrow, 0.015 mins elapsed.

2023-11-13 05:00:52.052307 : mat SRR13927736, Class = dgCMatrix
mat SRR13927736: nRows = 210721, nCols = 3534
mat SRR13927736: NonZeroEntries = 25583439, EntryRange = [ 1 , 4 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927736#GCTGTTCTCCTTGACC-1 SRR13927736#AGTTACGAGTAGGTCG-1 SRR13927736#AGTTTGGAGCATTCCA-1 SRR13927736#GCACGGTTCGGCAATT-1
[1,]                              .                              .                              .                              .
[2,]                              .                              .                              .                              .
[3,]                              .                              .                              1                              .
[4,]                              .                              .                              .                              .
[5,]                              .                              .                              .                              .
     SRR13927736#TCAGCTCTCTCTATTG-1
[1,]                              .
[2,]                              .
[3,]                              .
[4,]                              .
[5,]                              .


2023-11-13 05:00:52.24582 : Organizing SE PeakMatrix from ArrowFile : SRR13927736.arrow, 0.188 mins elapsed.

2023-11-13 05:00:52.718909 : se SRR13927736, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 3534 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(3534): SRR13927736#GCTGTTCTCCTTGACC-1 SRR13927736#AGTTACGAGTAGGTCG-1 ... SRR13927736#TAAACCGCAAGTCTGT-1
  SRR13927736#TGCATTTGTTCTACCC-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-13 05:00:52.726352 : Completed PeakMatrix : SRR13927736(4 of 4), 0.889 mins elapsed.
2023-11-13 05:00:52.727378 : Organizing colData, 0.889 mins elapsed.
2023-11-13 05:00:53.075891 : Organizing rowData, 0.895 mins elapsed.
2023-11-13 05:00:53.087833 : Organizing rowRanges, 0.895 mins elapsed.
2023-11-13 05:00:53.102613 : Organizing Assays (1 of 1), 0.895 mins elapsed.
2023-11-13 05:00:56.579685 : Constructing SummarizedExperiment, 0.953 mins elapsed.
2023-11-13 05:01:15.527075 : Finished Matrix Creation, 1.269 mins elapsed.

2023-11-13 05:04:59.171606 : getMatrixFromProject Input-Parameters, Class = list

getMatrixFromProject Input-Parameters$ArchRProj: length = 1

getMatrixFromProject Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromProject Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromProject Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromProject Input-Parameters$verbose: length = 1
[1] TRUE


getMatrixFromProject Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromProject Input-Parameters$threads: length = 1
[1] 1


getMatrixFromProject Input-Parameters$logFile: length = 1
[1] "x"


2023-11-13 05:05:00.739054 : Reading GeneIntegrationMatrix : SRR13927737(1 of 4), 0.026 mins elapsed.

2023-11-13 05:05:00.770006 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927737 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927737.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 4439
[1] "SRR13927737#GTCACGGGTCAACAGG-1" "SRR13927737#TGCCTGTTCTTACTCA-1"
[3] "SRR13927737#TTGCGGGTCAGGTCTA-1" "SRR13927737#GGTTGCGCAAAGCTTC-1"
[5] "SRR13927737#CTTAATCTCGCGCTGA-1" "SRR13927737#TGTACGATCGGAAAGC-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-13 05:05:01.472432 : featureDF SRR13927737, Class = DFrame
DataFrame with 23525 rows and 6 columns
      seqnames     start       end    strand        name       idx
         <Rle> <integer> <integer> <integer> <character> <integer>
1         chr1     69091     70008         1       OR4F5         1
2         chr1    817371    819837         1      FAM87B         2
3         chr1    825138    859446         1   LINC01128         3
4         chr1    827522    826206         2   LINC00115         4
5         chr1    876903    868071         2      FAM41C         5
...        ...       ...       ...       ...         ...       ...
23521     chrX 155264589 155258241         2      RAB39B       923
23522     chrX 155334657 155276211         2       CLIC2       924
23523     chrX 155380787 155381134         1      H2AFB1       925
23524     chrX 155382115 155383230         1        F8A2       926
23525     chrX 155669944 155490115         2       TMLHE       927

2023-11-13 05:05:01.519482 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927737.arrow, 0.012 mins elapsed.

2023-11-13 05:05:40.089989 : mat SRR13927737, Class = dgCMatrix
mat SRR13927737: nRows = 23525, nCols = 4439
mat SRR13927737: NonZeroEntries = 50277774, EntryRange = [ 0.000148130054904683 , 376.345870564644 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927737#GTCACGGGTCAACAGG-1 SRR13927737#TGCCTGTTCTTACTCA-1
[1,]                      .                             .         
[2,]                      .                             .         
[3,]                      0.1057576                     0.09630696
[4,]                      .                             0.01806508
[5,]                      .                             .         
     SRR13927737#TTGCGGGTCAGGTCTA-1 SRR13927737#GGTTGCGCAAAGCTTC-1
[1,]                     .                              .         
[2,]                     .                              .         
[3,]                     0.09113499                     0.08396973
[4,]                     0.02096304                     0.01727603
[5,]                     .                              .         
     SRR13927737#CTTAATCTCGCGCTGA-1
[1,]                     .         
[2,]                     .         
[3,]                     0.08327649
[4,]                     0.01368917
[5,]                     .         


2023-11-13 05:05:40.557543 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927737.arrow, 0.663 mins elapsed.

2023-11-13 05:05:40.965064 : se SRR13927737, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 4439 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(4439): SRR13927737#GTCACGGGTCAACAGG-1
  SRR13927737#TGCCTGTTCTTACTCA-1 ... SRR13927737#TTGCTATAGATGTTCC-1
  SRR13927737#GCTCACTAGTTAGCGG-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-13 05:05:40.971172 : Completed GeneIntegrationMatrix : SRR13927737(1 of 4), 0.697 mins elapsed.
2023-11-13 05:05:40.971567 : Reading GeneIntegrationMatrix : SRR13927738(2 of 4), 0.697 mins elapsed.

2023-11-13 05:05:40.990647 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927738 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927738.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3736
[1] "SRR13927738#TAGCATGAGGCTAAAT-1" "SRR13927738#TGTTAGGTCGGTCAGC-1"
[3] "SRR13927738#CTCTCAGCATGTCCCT-1" "SRR13927738#ACTATTCGTCCAACCG-1"
[5] "SRR13927738#TGGTCCTTCTACTGCC-1" "SRR13927738#ATTGTCTAGACTAGGC-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-13 05:05:41.728609 : featureDF SRR13927738, Class = DFrame
DataFrame with 23525 rows and 6 columns
      seqnames     start       end    strand        name       idx
         <Rle> <integer> <integer> <integer> <character> <integer>
1         chr1     69091     70008         1       OR4F5         1
2         chr1    817371    819837         1      FAM87B         2
3         chr1    825138    859446         1   LINC01128         3
4         chr1    827522    826206         2   LINC00115         4
5         chr1    876903    868071         2      FAM41C         5
...        ...       ...       ...       ...         ...       ...
23521     chrX 155264589 155258241         2      RAB39B       923
23522     chrX 155334657 155276211         2       CLIC2       924
23523     chrX 155380787 155381134         1      H2AFB1       925
23524     chrX 155382115 155383230         1        F8A2       926
23525     chrX 155669944 155490115         2       TMLHE       927

2023-11-13 05:05:41.76208 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927738.arrow, 0.013 mins elapsed.

2023-11-13 05:06:01.091066 : mat SRR13927738, Class = dgCMatrix
mat SRR13927738: nRows = 23525, nCols = 3736
mat SRR13927738: NonZeroEntries = 41898897, EntryRange = [ 9.73975135001421e-05 , 758.763031303812 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927738#TAGCATGAGGCTAAAT-1 SRR13927738#TGTTAGGTCGGTCAGC-1
[1,]                    .                                .        
[2,]                    .                                .        
[3,]                    0.002346352                      0.0084575
[4,]                    0.005836413                      .        
[5,]                    .                                .        
     SRR13927738#CTCTCAGCATGTCCCT-1 SRR13927738#ACTATTCGTCCAACCG-1
[1,]                     .                             .          
[2,]                     .                             .          
[3,]                     .                             0.004081077
[4,]                     0.02001583                    .          
[5,]                     .                             .          
     SRR13927738#TGGTCCTTCTACTGCC-1
[1,]                    .          
[2,]                    .          
[3,]                    0.007640877
[4,]                    0.015048796
[5,]                    .          


2023-11-13 05:06:01.444391 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927738.arrow, 0.341 mins elapsed.

2023-11-13 05:06:01.781621 : se SRR13927738, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 3736 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(3736): SRR13927738#TAGCATGAGGCTAAAT-1
  SRR13927738#TGTTAGGTCGGTCAGC-1 ... SRR13927738#GCTCCTATCCGTGCGA-1
  SRR13927738#GGCACGTAGTCGATAA-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-13 05:06:01.786466 : Completed GeneIntegrationMatrix : SRR13927738(2 of 4), 1.044 mins elapsed.
2023-11-13 05:06:01.786863 : Reading GeneIntegrationMatrix : SRR13927735(3 of 4), 1.044 mins elapsed.

2023-11-13 05:06:01.805393 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927735 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927735.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3813
[1] "SRR13927735#TTATGTCTCCAGGTAT-1" "SRR13927735#TATTGCTCATCAGAAA-1"
[3] "SRR13927735#TTCGATTGTAGGGTTG-1" "SRR13927735#CATTCATTCGGATGTT-1"
[5] "SRR13927735#ACGTTAGGTCAACTGT-1" "SRR13927735#AAATGCCCAGCAATGG-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-13 05:06:02.24446 : featureDF SRR13927735, Class = DFrame
DataFrame with 23525 rows and 6 columns
      seqnames     start       end    strand        name       idx
         <Rle> <integer> <integer> <integer> <character> <integer>
1         chr1     69091     70008         1       OR4F5         1
2         chr1    817371    819837         1      FAM87B         2
3         chr1    825138    859446         1   LINC01128         3
4         chr1    827522    826206         2   LINC00115         4
5         chr1    876903    868071         2      FAM41C         5
...        ...       ...       ...       ...         ...       ...
23521     chrX 155264589 155258241         2      RAB39B       923
23522     chrX 155334657 155276211         2       CLIC2       924
23523     chrX 155380787 155381134         1      H2AFB1       925
23524     chrX 155382115 155383230         1        F8A2       926
23525     chrX 155669944 155490115         2       TMLHE       927

2023-11-13 05:06:02.260736 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927735.arrow, 0.008 mins elapsed.

2023-11-13 05:06:22.169518 : mat SRR13927735, Class = dgCMatrix
mat SRR13927735: nRows = 23525, nCols = 3813
mat SRR13927735: NonZeroEntries = 45535634, EntryRange = [ 7.20715653308254e-05 , 673.990433360078 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927735#TTATGTCTCCAGGTAT-1 SRR13927735#TATTGCTCATCAGAAA-1
[1,]                    .                              .          
[2,]                    .                              .          
[3,]                    0.070892732                    0.057009268
[4,]                    0.008800039                    0.002978408
[5,]                    .                              .          
     SRR13927735#TTCGATTGTAGGGTTG-1 SRR13927735#CATTCATTCGGATGTT-1
[1,]                     .                               .        
[2,]                     .                               .        
[3,]                     0.05523011                      0.0248301
[4,]                     .                               .        
[5,]                     .                               .        
     SRR13927735#ACGTTAGGTCAACTGT-1
[1,]                     .         
[2,]                     .         
[3,]                     0.05384571
[4,]                     0.01228688
[5,]                     .         


2023-11-13 05:06:22.518125 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927735.arrow, 0.345 mins elapsed.

2023-11-13 05:06:22.80568 : se SRR13927735, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 3813 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(3813): SRR13927735#TTATGTCTCCAGGTAT-1
  SRR13927735#TATTGCTCATCAGAAA-1 ... SRR13927735#GATCATGGTGTTGTTG-1
  SRR13927735#GTCCATCCATTTCACT-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-13 05:06:22.810813 : Completed GeneIntegrationMatrix : SRR13927735(3 of 4), 1.394 mins elapsed.
2023-11-13 05:06:22.811275 : Reading GeneIntegrationMatrix : SRR13927736(4 of 4), 1.394 mins elapsed.

2023-11-13 05:06:22.834556 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927736 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927736.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3534
[1] "SRR13927736#GCTGTTCTCCTTGACC-1" "SRR13927736#AGTTACGAGTAGGTCG-1"
[3] "SRR13927736#AGTTTGGAGCATTCCA-1" "SRR13927736#GCACGGTTCGGCAATT-1"
[5] "SRR13927736#TCAGCTCTCTCTATTG-1" "SRR13927736#CCTTAATTCTACTTTG-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-13 05:06:23.337555 : featureDF SRR13927736, Class = DFrame
DataFrame with 23525 rows and 6 columns
      seqnames     start       end    strand        name       idx
         <Rle> <integer> <integer> <integer> <character> <integer>
1         chr1     69091     70008         1       OR4F5         1
2         chr1    817371    819837         1      FAM87B         2
3         chr1    825138    859446         1   LINC01128         3
4         chr1    827522    826206         2   LINC00115         4
5         chr1    876903    868071         2      FAM41C         5
...        ...       ...       ...       ...         ...       ...
23521     chrX 155264589 155258241         2      RAB39B       923
23522     chrX 155334657 155276211         2       CLIC2       924
23523     chrX 155380787 155381134         1      H2AFB1       925
23524     chrX 155382115 155383230         1        F8A2       926
23525     chrX 155669944 155490115         2       TMLHE       927

2023-11-13 05:06:23.350998 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927736.arrow, 0.009 mins elapsed.

2023-11-13 05:06:39.167982 : mat SRR13927736, Class = dgCMatrix
mat SRR13927736: nRows = 23525, nCols = 3534
mat SRR13927736: NonZeroEntries = 39003616, EntryRange = [ 8.77627031078293e-05 , 770.387276519725 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927736#GCTGTTCTCCTTGACC-1 SRR13927736#AGTTACGAGTAGGTCG-1
[1,]                     .                             .          
[2,]                     .                             .          
[3,]                     0.07845381                    0.028344539
[4,]                     0.01395289                    0.008499189
[5,]                     .                             .          
     SRR13927736#AGTTTGGAGCATTCCA-1 SRR13927736#GCACGGTTCGGCAATT-1
[1,]                     .                              .         
[2,]                     .                              .         
[3,]                     0.06535130                     0.06576222
[4,]                     0.01810177                     0.01315728
[5,]                     .                              .         
     SRR13927736#TCAGCTCTCTCTATTG-1
[1,]                    .          
[2,]                    .          
[3,]                    0.016229526
[4,]                    0.003862592
[5,]                    .          


2023-11-13 05:06:40.684528 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927736.arrow, 0.298 mins elapsed.

2023-11-13 05:06:40.955923 : se SRR13927736, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 3534 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(3534): SRR13927736#GCTGTTCTCCTTGACC-1
  SRR13927736#AGTTACGAGTAGGTCG-1 ... SRR13927736#TAAACCGCAAGTCTGT-1
  SRR13927736#TGCATTTGTTCTACCC-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-13 05:06:40.96088 : Completed GeneIntegrationMatrix : SRR13927736(4 of 4), 1.697 mins elapsed.
2023-11-13 05:06:40.96181 : Organizing colData, 1.697 mins elapsed.
2023-11-13 05:06:41.109731 : Organizing rowData, 1.699 mins elapsed.
2023-11-13 05:06:41.112914 : Organizing rowRanges, 1.699 mins elapsed.
2023-11-13 05:06:41.118701 : Organizing Assays (1 of 1), 1.699 mins elapsed.
2023-11-13 05:06:47.935934 : Constructing SummarizedExperiment, 1.813 mins elapsed.
2023-11-13 05:06:49.216145 : Finished Matrix Creation, 1.834 mins elapsed.

2023-11-13 05:06:49.218526 : getMatrixFromProject Input-Parameters, Class = list

getMatrixFromProject Input-Parameters$ArchRProj: length = 1

getMatrixFromProject Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromProject Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromProject Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromProject Input-Parameters$verbose: length = 1
[1] TRUE


getMatrixFromProject Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromProject Input-Parameters$threads: length = 1
[1] 1


getMatrixFromProject Input-Parameters$logFile: length = 1
[1] "x"


2023-11-13 05:06:50.971184 : Reading PeakMatrix : SRR13927737(1 of 4), 0.029 mins elapsed.

2023-11-13 05:06:50.997987 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927737 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927737.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 4439
[1] "SRR13927737#GTCACGGGTCAACAGG-1" "SRR13927737#TGCCTGTTCTTACTCA-1"
[3] "SRR13927737#TTGCGGGTCAGGTCTA-1" "SRR13927737#GGTTGCGCAAAGCTTC-1"
[5] "SRR13927737#CTTAATCTCGCGCTGA-1" "SRR13927737#TGTACGATCGGAAAGC-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-13 05:06:51.825844 : featureDF SRR13927737, Class = DFrame
DataFrame with 210721 rows and 4 columns
       seqnames       idx     start       end
          <Rle> <integer> <integer> <integer>
1          chr1         1    804677    805177
2          chr1         2    811866    812366
3          chr1         3    817090    817590
4          chr1         4    818831    819331
5          chr1         5    826543    827043
...         ...       ...       ...       ...
210717     chrX      5777 155898525 155899025
210718     chrX      5778 155899804 155900304
210719     chrX      5779 155956608 155957108
210720     chrX      5780 155966827 155967327
210721     chrX      5781 155971625 155972125

2023-11-13 05:06:51.933326 : Getting PeakMatrix from ArrowFile : SRR13927737.arrow, 0.016 mins elapsed.

2023-11-13 05:07:04.196372 : mat SRR13927737, Class = dgCMatrix
mat SRR13927737: nRows = 210721, nCols = 4439
mat SRR13927737: NonZeroEntries = 29963541, EntryRange = [ 1 , 4 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927737#GTCACGGGTCAACAGG-1 SRR13927737#TGCCTGTTCTTACTCA-1
[1,]                              .                              .
[2,]                              .                              2
[3,]                              .                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927737#TTGCGGGTCAGGTCTA-1 SRR13927737#GGTTGCGCAAAGCTTC-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              2                              .
[4,]                              .                              .
[5,]                              4                              .
     SRR13927737#CTTAATCTCGCGCTGA-1
[1,]                              .
[2,]                              .
[3,]                              .
[4,]                              .
[5,]                              .


2023-11-13 05:07:04.374881 : Organizing SE PeakMatrix from ArrowFile : SRR13927737.arrow, 0.223 mins elapsed.

2023-11-13 05:07:04.669884 : se SRR13927737, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 4439 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(4439): SRR13927737#GTCACGGGTCAACAGG-1
  SRR13927737#TGCCTGTTCTTACTCA-1 ... SRR13927737#TTGCTATAGATGTTCC-1
  SRR13927737#GCTCACTAGTTAGCGG-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-13 05:07:04.679376 : Completed PeakMatrix : SRR13927737(1 of 4), 0.258 mins elapsed.
2023-11-13 05:07:04.679774 : Reading PeakMatrix : SRR13927738(2 of 4), 0.258 mins elapsed.

2023-11-13 05:07:04.692099 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927738 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927738.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3736
[1] "SRR13927738#TAGCATGAGGCTAAAT-1" "SRR13927738#TGTTAGGTCGGTCAGC-1"
[3] "SRR13927738#CTCTCAGCATGTCCCT-1" "SRR13927738#ACTATTCGTCCAACCG-1"
[5] "SRR13927738#TGGTCCTTCTACTGCC-1" "SRR13927738#ATTGTCTAGACTAGGC-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-13 05:07:05.15984 : featureDF SRR13927738, Class = DFrame
DataFrame with 210721 rows and 4 columns
       seqnames       idx     start       end
          <Rle> <integer> <integer> <integer>
1          chr1         1    804677    805177
2          chr1         2    811866    812366
3          chr1         3    817090    817590
4          chr1         4    818831    819331
5          chr1         5    826543    827043
...         ...       ...       ...       ...
210717     chrX      5777 155898525 155899025
210718     chrX      5778 155899804 155900304
210719     chrX      5779 155956608 155957108
210720     chrX      5780 155966827 155967327
210721     chrX      5781 155971625 155972125

2023-11-13 05:07:05.21058 : Getting PeakMatrix from ArrowFile : SRR13927738.arrow, 0.009 mins elapsed.

2023-11-13 05:07:16.559804 : mat SRR13927738, Class = dgCMatrix
mat SRR13927738: nRows = 210721, nCols = 3736
mat SRR13927738: NonZeroEntries = 30165921, EntryRange = [ 1 , 4 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927738#TAGCATGAGGCTAAAT-1 SRR13927738#TGTTAGGTCGGTCAGC-1
[1,]                              2                              .
[2,]                              .                              .
[3,]                              .                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927738#CTCTCAGCATGTCCCT-1 SRR13927738#ACTATTCGTCCAACCG-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              .                              .
[4,]                              .                              .
[5,]                              2                              .
     SRR13927738#TGGTCCTTCTACTGCC-1
[1,]                              .
[2,]                              .
[3,]                              2
[4,]                              .
[5,]                              .


2023-11-13 05:07:16.786464 : Organizing SE PeakMatrix from ArrowFile : SRR13927738.arrow, 0.202 mins elapsed.

2023-11-13 05:07:17.074393 : se SRR13927738, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 3736 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(3736): SRR13927738#TAGCATGAGGCTAAAT-1
  SRR13927738#TGTTAGGTCGGTCAGC-1 ... SRR13927738#GCTCCTATCCGTGCGA-1
  SRR13927738#GGCACGTAGTCGATAA-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-13 05:07:17.081408 : Completed PeakMatrix : SRR13927738(2 of 4), 0.464 mins elapsed.
2023-11-13 05:07:17.081796 : Reading PeakMatrix : SRR13927735(3 of 4), 0.464 mins elapsed.

2023-11-13 05:07:17.09441 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927735 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927735.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3813
[1] "SRR13927735#TTATGTCTCCAGGTAT-1" "SRR13927735#TATTGCTCATCAGAAA-1"
[3] "SRR13927735#TTCGATTGTAGGGTTG-1" "SRR13927735#CATTCATTCGGATGTT-1"
[5] "SRR13927735#ACGTTAGGTCAACTGT-1" "SRR13927735#AAATGCCCAGCAATGG-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-13 05:07:17.586854 : featureDF SRR13927735, Class = DFrame
DataFrame with 210721 rows and 4 columns
       seqnames       idx     start       end
          <Rle> <integer> <integer> <integer>
1          chr1         1    804677    805177
2          chr1         2    811866    812366
3          chr1         3    817090    817590
4          chr1         4    818831    819331
5          chr1         5    826543    827043
...         ...       ...       ...       ...
210717     chrX      5777 155898525 155899025
210718     chrX      5778 155899804 155900304
210719     chrX      5779 155956608 155957108
210720     chrX      5780 155966827 155967327
210721     chrX      5781 155971625 155972125

2023-11-13 05:07:17.637431 : Getting PeakMatrix from ArrowFile : SRR13927735.arrow, 0.009 mins elapsed.

2023-11-13 05:07:27.484999 : mat SRR13927735, Class = dgCMatrix
mat SRR13927735: nRows = 210721, nCols = 3813
mat SRR13927735: NonZeroEntries = 26038084, EntryRange = [ 1 , 4 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927735#TTATGTCTCCAGGTAT-1 SRR13927735#TATTGCTCATCAGAAA-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              1                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927735#TTCGATTGTAGGGTTG-1 SRR13927735#CATTCATTCGGATGTT-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              .                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927735#ACGTTAGGTCAACTGT-1
[1,]                              2
[2,]                              .
[3,]                              .
[4,]                              .
[5,]                              .


2023-11-13 05:07:27.681299 : Organizing SE PeakMatrix from ArrowFile : SRR13927735.arrow, 0.176 mins elapsed.

2023-11-13 05:07:27.954411 : se SRR13927735, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 3813 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(3813): SRR13927735#TTATGTCTCCAGGTAT-1
  SRR13927735#TATTGCTCATCAGAAA-1 ... SRR13927735#GATCATGGTGTTGTTG-1
  SRR13927735#GTCCATCCATTTCACT-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-13 05:07:27.961326 : Completed PeakMatrix : SRR13927735(3 of 4), 0.646 mins elapsed.
2023-11-13 05:07:27.961729 : Reading PeakMatrix : SRR13927736(4 of 4), 0.646 mins elapsed.

2023-11-13 05:07:27.980593 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927736 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927736.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3534
[1] "SRR13927736#GCTGTTCTCCTTGACC-1" "SRR13927736#AGTTACGAGTAGGTCG-1"
[3] "SRR13927736#AGTTTGGAGCATTCCA-1" "SRR13927736#GCACGGTTCGGCAATT-1"
[5] "SRR13927736#TCAGCTCTCTCTATTG-1" "SRR13927736#CCTTAATTCTACTTTG-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-13 05:07:28.49035 : featureDF SRR13927736, Class = DFrame
DataFrame with 210721 rows and 4 columns
       seqnames       idx     start       end
          <Rle> <integer> <integer> <integer>
1          chr1         1    804677    805177
2          chr1         2    811866    812366
3          chr1         3    817090    817590
4          chr1         4    818831    819331
5          chr1         5    826543    827043
...         ...       ...       ...       ...
210717     chrX      5777 155898525 155899025
210718     chrX      5778 155899804 155900304
210719     chrX      5779 155956608 155957108
210720     chrX      5780 155966827 155967327
210721     chrX      5781 155971625 155972125

2023-11-13 05:07:28.541008 : Getting PeakMatrix from ArrowFile : SRR13927736.arrow, 0.009 mins elapsed.

2023-11-13 05:07:38.279847 : mat SRR13927736, Class = dgCMatrix
mat SRR13927736: nRows = 210721, nCols = 3534
mat SRR13927736: NonZeroEntries = 25583439, EntryRange = [ 1 , 4 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927736#GCTGTTCTCCTTGACC-1 SRR13927736#AGTTACGAGTAGGTCG-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              .                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927736#AGTTTGGAGCATTCCA-1 SRR13927736#GCACGGTTCGGCAATT-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              1                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927736#TCAGCTCTCTCTATTG-1
[1,]                              .
[2,]                              .
[3,]                              .
[4,]                              .
[5,]                              .


2023-11-13 05:07:38.474086 : Organizing SE PeakMatrix from ArrowFile : SRR13927736.arrow, 0.175 mins elapsed.

2023-11-13 05:07:38.753417 : se SRR13927736, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 3534 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(3534): SRR13927736#GCTGTTCTCCTTGACC-1
  SRR13927736#AGTTACGAGTAGGTCG-1 ... SRR13927736#TAAACCGCAAGTCTGT-1
  SRR13927736#TGCATTTGTTCTACCC-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-13 05:07:38.760337 : Completed PeakMatrix : SRR13927736(4 of 4), 0.826 mins elapsed.
2023-11-13 05:07:38.761276 : Organizing colData, 0.826 mins elapsed.
2023-11-13 05:07:38.897606 : Organizing rowData, 0.828 mins elapsed.
2023-11-13 05:07:38.905722 : Organizing rowRanges, 0.828 mins elapsed.
2023-11-13 05:07:38.917771 : Organizing Assays (1 of 1), 0.828 mins elapsed.
2023-11-13 05:07:41.967784 : Constructing SummarizedExperiment, 0.879 mins elapsed.
2023-11-13 05:07:55.572256 : Finished Matrix Creation, 1.106 mins elapsed.

2023-11-13 05:16:57.234546 : getMatrixFromProject Input-Parameters, Class = list

getMatrixFromProject Input-Parameters$ArchRProj: length = 1

getMatrixFromProject Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromProject Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromProject Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromProject Input-Parameters$verbose: length = 1
[1] TRUE


getMatrixFromProject Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromProject Input-Parameters$threads: length = 1
[1] 1


getMatrixFromProject Input-Parameters$logFile: length = 1
[1] "x"


2023-11-13 05:16:58.242043 : Reading GeneIntegrationMatrix : SRR13927737(1 of 4), 0.017 mins elapsed.

2023-11-13 05:16:58.281132 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927737 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927737.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 4439
[1] "SRR13927737#GTCACGGGTCAACAGG-1" "SRR13927737#TGCCTGTTCTTACTCA-1" "SRR13927737#TTGCGGGTCAGGTCTA-1" "SRR13927737#GGTTGCGCAAAGCTTC-1"
[5] "SRR13927737#CTTAATCTCGCGCTGA-1" "SRR13927737#TGTACGATCGGAAAGC-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-13 05:16:59.414992 : featureDF SRR13927737, Class = DFrame
DataFrame with 23525 rows and 6 columns
      seqnames     start       end    strand        name       idx
         <Rle> <integer> <integer> <integer> <character> <integer>
1         chr1     69091     70008         1       OR4F5         1
2         chr1    817371    819837         1      FAM87B         2
3         chr1    825138    859446         1   LINC01128         3
4         chr1    827522    826206         2   LINC00115         4
5         chr1    876903    868071         2      FAM41C         5
...        ...       ...       ...       ...         ...       ...
23521     chrX 155264589 155258241         2      RAB39B       923
23522     chrX 155334657 155276211         2       CLIC2       924
23523     chrX 155380787 155381134         1      H2AFB1       925
23524     chrX 155382115 155383230         1        F8A2       926
23525     chrX 155669944 155490115         2       TMLHE       927

2023-11-13 05:16:59.441739 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927737.arrow, 0.019 mins elapsed.

2023-11-13 05:17:20.699726 : mat SRR13927737, Class = dgCMatrix
mat SRR13927737: nRows = 23525, nCols = 4439
mat SRR13927737: NonZeroEntries = 50277774, EntryRange = [ 0.000148130054904683 , 376.345870564644 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927737#GTCACGGGTCAACAGG-1 SRR13927737#TGCCTGTTCTTACTCA-1 SRR13927737#TTGCGGGTCAGGTCTA-1 SRR13927737#GGTTGCGCAAAGCTTC-1
[1,]                      .                             .                              .                              .         
[2,]                      .                             .                              .                              .         
[3,]                      0.1057576                     0.09630696                     0.09113499                     0.08396973
[4,]                      .                             0.01806508                     0.02096304                     0.01727603
[5,]                      .                             .                              .                              .         
     SRR13927737#CTTAATCTCGCGCTGA-1
[1,]                     .         
[2,]                     .         
[3,]                     0.08327649
[4,]                     0.01368917
[5,]                     .         


2023-11-13 05:17:21.450299 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927737.arrow, 0.386 mins elapsed.

2023-11-13 05:17:22.081956 : se SRR13927737, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 4439 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(4439): SRR13927737#GTCACGGGTCAACAGG-1 SRR13927737#TGCCTGTTCTTACTCA-1 ... SRR13927737#TTGCTATAGATGTTCC-1
  SRR13927737#GCTCACTAGTTAGCGG-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-13 05:17:22.088972 : Completed GeneIntegrationMatrix : SRR13927737(1 of 4), 0.414 mins elapsed.
2023-11-13 05:17:22.089537 : Reading GeneIntegrationMatrix : SRR13927738(2 of 4), 0.414 mins elapsed.

2023-11-13 05:17:22.105444 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927738 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927738.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3736
[1] "SRR13927738#TAGCATGAGGCTAAAT-1" "SRR13927738#TGTTAGGTCGGTCAGC-1" "SRR13927738#CTCTCAGCATGTCCCT-1" "SRR13927738#ACTATTCGTCCAACCG-1"
[5] "SRR13927738#TGGTCCTTCTACTGCC-1" "SRR13927738#ATTGTCTAGACTAGGC-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-13 05:17:22.835716 : featureDF SRR13927738, Class = DFrame
DataFrame with 23525 rows and 6 columns
      seqnames     start       end    strand        name       idx
         <Rle> <integer> <integer> <integer> <character> <integer>
1         chr1     69091     70008         1       OR4F5         1
2         chr1    817371    819837         1      FAM87B         2
3         chr1    825138    859446         1   LINC01128         3
4         chr1    827522    826206         2   LINC00115         4
5         chr1    876903    868071         2      FAM41C         5
...        ...       ...       ...       ...         ...       ...
23521     chrX 155264589 155258241         2      RAB39B       923
23522     chrX 155334657 155276211         2       CLIC2       924
23523     chrX 155380787 155381134         1      H2AFB1       925
23524     chrX 155382115 155383230         1        F8A2       926
23525     chrX 155669944 155490115         2       TMLHE       927

2023-11-13 05:17:22.854466 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927738.arrow, 0.012 mins elapsed.

2023-11-13 05:17:39.09653 : mat SRR13927738, Class = dgCMatrix
mat SRR13927738: nRows = 23525, nCols = 3736
mat SRR13927738: NonZeroEntries = 41898897, EntryRange = [ 9.73975135001421e-05 , 758.763031303812 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927738#TAGCATGAGGCTAAAT-1 SRR13927738#TGTTAGGTCGGTCAGC-1 SRR13927738#CTCTCAGCATGTCCCT-1 SRR13927738#ACTATTCGTCCAACCG-1
[1,]                    .                                .                             .                             .          
[2,]                    .                                .                             .                             .          
[3,]                    0.002346352                      0.0084575                     .                             0.004081077
[4,]                    0.005836413                      .                             0.02001583                    .          
[5,]                    .                                .                             .                             .          
     SRR13927738#TGGTCCTTCTACTGCC-1
[1,]                    .          
[2,]                    .          
[3,]                    0.007640877
[4,]                    0.015048796
[5,]                    .          


2023-11-13 05:17:39.778509 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927738.arrow, 0.295 mins elapsed.

2023-11-13 05:17:40.32738 : se SRR13927738, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 3736 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(3736): SRR13927738#TAGCATGAGGCTAAAT-1 SRR13927738#TGTTAGGTCGGTCAGC-1 ... SRR13927738#GCTCCTATCCGTGCGA-1
  SRR13927738#GGCACGTAGTCGATAA-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-13 05:17:40.334412 : Completed GeneIntegrationMatrix : SRR13927738(2 of 4), 0.718 mins elapsed.
2023-11-13 05:17:40.335023 : Reading GeneIntegrationMatrix : SRR13927735(3 of 4), 0.718 mins elapsed.

2023-11-13 05:17:40.358687 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927735 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927735.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3813
[1] "SRR13927735#TTATGTCTCCAGGTAT-1" "SRR13927735#TATTGCTCATCAGAAA-1" "SRR13927735#TTCGATTGTAGGGTTG-1" "SRR13927735#CATTCATTCGGATGTT-1"
[5] "SRR13927735#ACGTTAGGTCAACTGT-1" "SRR13927735#AAATGCCCAGCAATGG-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-13 05:17:40.961762 : featureDF SRR13927735, Class = DFrame
DataFrame with 23525 rows and 6 columns
      seqnames     start       end    strand        name       idx
         <Rle> <integer> <integer> <integer> <character> <integer>
1         chr1     69091     70008         1       OR4F5         1
2         chr1    817371    819837         1      FAM87B         2
3         chr1    825138    859446         1   LINC01128         3
4         chr1    827522    826206         2   LINC00115         4
5         chr1    876903    868071         2      FAM41C         5
...        ...       ...       ...       ...         ...       ...
23521     chrX 155264589 155258241         2      RAB39B       923
23522     chrX 155334657 155276211         2       CLIC2       924
23523     chrX 155380787 155381134         1      H2AFB1       925
23524     chrX 155382115 155383230         1        F8A2       926
23525     chrX 155669944 155490115         2       TMLHE       927

2023-11-13 05:17:40.984832 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927735.arrow, 0.01 mins elapsed.

2023-11-13 05:18:00.610753 : mat SRR13927735, Class = dgCMatrix
mat SRR13927735: nRows = 23525, nCols = 3813
mat SRR13927735: NonZeroEntries = 45535634, EntryRange = [ 7.20715653308254e-05 , 673.990433360078 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927735#TTATGTCTCCAGGTAT-1 SRR13927735#TATTGCTCATCAGAAA-1 SRR13927735#TTCGATTGTAGGGTTG-1 SRR13927735#CATTCATTCGGATGTT-1
[1,]                    .                              .                               .                               .        
[2,]                    .                              .                               .                               .        
[3,]                    0.070892732                    0.057009268                     0.05523011                      0.0248301
[4,]                    0.008800039                    0.002978408                     .                               .        
[5,]                    .                              .                               .                               .        
     SRR13927735#ACGTTAGGTCAACTGT-1
[1,]                     .         
[2,]                     .         
[3,]                     0.05384571
[4,]                     0.01228688
[5,]                     .         


2023-11-13 05:18:01.352928 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927735.arrow, 0.35 mins elapsed.

2023-11-13 05:18:01.857109 : se SRR13927735, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 3813 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(3813): SRR13927735#TTATGTCTCCAGGTAT-1 SRR13927735#TATTGCTCATCAGAAA-1 ... SRR13927735#GATCATGGTGTTGTTG-1
  SRR13927735#GTCCATCCATTTCACT-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-13 05:18:01.863954 : Completed GeneIntegrationMatrix : SRR13927735(3 of 4), 1.077 mins elapsed.
2023-11-13 05:18:01.86454 : Reading GeneIntegrationMatrix : SRR13927736(4 of 4), 1.077 mins elapsed.

2023-11-13 05:18:01.881019 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927736 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927736.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3534
[1] "SRR13927736#GCTGTTCTCCTTGACC-1" "SRR13927736#AGTTACGAGTAGGTCG-1" "SRR13927736#AGTTTGGAGCATTCCA-1" "SRR13927736#GCACGGTTCGGCAATT-1"
[5] "SRR13927736#TCAGCTCTCTCTATTG-1" "SRR13927736#CCTTAATTCTACTTTG-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-13 05:18:02.59556 : featureDF SRR13927736, Class = DFrame
DataFrame with 23525 rows and 6 columns
      seqnames     start       end    strand        name       idx
         <Rle> <integer> <integer> <integer> <character> <integer>
1         chr1     69091     70008         1       OR4F5         1
2         chr1    817371    819837         1      FAM87B         2
3         chr1    825138    859446         1   LINC01128         3
4         chr1    827522    826206         2   LINC00115         4
5         chr1    876903    868071         2      FAM41C         5
...        ...       ...       ...       ...         ...       ...
23521     chrX 155264589 155258241         2      RAB39B       923
23522     chrX 155334657 155276211         2       CLIC2       924
23523     chrX 155380787 155381134         1      H2AFB1       925
23524     chrX 155382115 155383230         1        F8A2       926
23525     chrX 155669944 155490115         2       TMLHE       927

2023-11-13 05:18:02.625646 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927736.arrow, 0.012 mins elapsed.

2023-11-13 05:18:17.426019 : mat SRR13927736, Class = dgCMatrix
mat SRR13927736: nRows = 23525, nCols = 3534
mat SRR13927736: NonZeroEntries = 39003616, EntryRange = [ 8.77627031078293e-05 , 770.387276519725 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927736#GCTGTTCTCCTTGACC-1 SRR13927736#AGTTACGAGTAGGTCG-1 SRR13927736#AGTTTGGAGCATTCCA-1 SRR13927736#GCACGGTTCGGCAATT-1
[1,]                     .                             .                               .                              .         
[2,]                     .                             .                               .                              .         
[3,]                     0.07845381                    0.028344539                     0.06535130                     0.06576222
[4,]                     0.01395289                    0.008499189                     0.01810177                     0.01315728
[5,]                     .                             .                               .                              .         
     SRR13927736#TCAGCTCTCTCTATTG-1
[1,]                    .          
[2,]                    .          
[3,]                    0.016229526
[4,]                    0.003862592
[5,]                    .          


2023-11-13 05:18:18.124168 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927736.arrow, 0.271 mins elapsed.

2023-11-13 05:18:18.638398 : se SRR13927736, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 3534 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(3534): SRR13927736#GCTGTTCTCCTTGACC-1 SRR13927736#AGTTACGAGTAGGTCG-1 ... SRR13927736#TAAACCGCAAGTCTGT-1
  SRR13927736#TGCATTTGTTCTACCC-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-13 05:18:18.645299 : Completed GeneIntegrationMatrix : SRR13927736(4 of 4), 1.357 mins elapsed.
2023-11-13 05:18:18.646626 : Organizing colData, 1.357 mins elapsed.
2023-11-13 05:18:19.010329 : Organizing rowData, 1.363 mins elapsed.
2023-11-13 05:18:19.014758 : Organizing rowRanges, 1.363 mins elapsed.
2023-11-13 05:18:19.022414 : Organizing Assays (1 of 1), 1.363 mins elapsed.
2023-11-13 05:18:26.200759 : Constructing SummarizedExperiment, 1.483 mins elapsed.
2023-11-13 05:18:27.798652 : Finished Matrix Creation, 1.507 mins elapsed.

2023-11-13 05:23:23.576157 : getMatrixFromProject Input-Parameters, Class = list

getMatrixFromProject Input-Parameters$ArchRProj: length = 1

getMatrixFromProject Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromProject Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromProject Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromProject Input-Parameters$verbose: length = 1
[1] TRUE


getMatrixFromProject Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromProject Input-Parameters$threads: length = 1
[1] 1


getMatrixFromProject Input-Parameters$logFile: length = 1
[1] "x"


2023-11-13 05:23:25.277605 : Reading GeneIntegrationMatrix : SRR13927737(1 of 4), 0.028 mins elapsed.

2023-11-13 05:23:25.326798 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927737 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927737.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 4439
[1] "SRR13927737#GTCACGGGTCAACAGG-1" "SRR13927737#TGCCTGTTCTTACTCA-1"
[3] "SRR13927737#TTGCGGGTCAGGTCTA-1" "SRR13927737#GGTTGCGCAAAGCTTC-1"
[5] "SRR13927737#CTTAATCTCGCGCTGA-1" "SRR13927737#TGTACGATCGGAAAGC-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-13 05:23:26.089832 : featureDF SRR13927737, Class = DFrame
DataFrame with 23525 rows and 6 columns
      seqnames     start       end    strand        name       idx
         <Rle> <integer> <integer> <integer> <character> <integer>
1         chr1     69091     70008         1       OR4F5         1
2         chr1    817371    819837         1      FAM87B         2
3         chr1    825138    859446         1   LINC01128         3
4         chr1    827522    826206         2   LINC00115         4
5         chr1    876903    868071         2      FAM41C         5
...        ...       ...       ...       ...         ...       ...
23521     chrX 155264589 155258241         2      RAB39B       923
23522     chrX 155334657 155276211         2       CLIC2       924
23523     chrX 155380787 155381134         1      H2AFB1       925
23524     chrX 155382115 155383230         1        F8A2       926
23525     chrX 155669944 155490115         2       TMLHE       927

2023-11-13 05:23:26.154579 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927737.arrow, 0.014 mins elapsed.

2023-11-13 05:24:04.750722 : mat SRR13927737, Class = dgCMatrix
mat SRR13927737: nRows = 23525, nCols = 4439
mat SRR13927737: NonZeroEntries = 50277774, EntryRange = [ 0.000148130054904683 , 376.345870564644 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927737#GTCACGGGTCAACAGG-1 SRR13927737#TGCCTGTTCTTACTCA-1
[1,]                      .                             .         
[2,]                      .                             .         
[3,]                      0.1057576                     0.09630696
[4,]                      .                             0.01806508
[5,]                      .                             .         
     SRR13927737#TTGCGGGTCAGGTCTA-1 SRR13927737#GGTTGCGCAAAGCTTC-1
[1,]                     .                              .         
[2,]                     .                              .         
[3,]                     0.09113499                     0.08396973
[4,]                     0.02096304                     0.01727603
[5,]                     .                              .         
     SRR13927737#CTTAATCTCGCGCTGA-1
[1,]                     .         
[2,]                     .         
[3,]                     0.08327649
[4,]                     0.01368917
[5,]                     .         


2023-11-13 05:24:05.220294 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927737.arrow, 0.665 mins elapsed.

2023-11-13 05:24:05.656834 : se SRR13927737, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 4439 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(4439): SRR13927737#GTCACGGGTCAACAGG-1
  SRR13927737#TGCCTGTTCTTACTCA-1 ... SRR13927737#TTGCTATAGATGTTCC-1
  SRR13927737#GCTCACTAGTTAGCGG-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-13 05:24:05.664016 : Completed GeneIntegrationMatrix : SRR13927737(1 of 4), 0.702 mins elapsed.
2023-11-13 05:24:05.664496 : Reading GeneIntegrationMatrix : SRR13927738(2 of 4), 0.702 mins elapsed.

2023-11-13 05:24:05.67886 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927738 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927738.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3736
[1] "SRR13927738#TAGCATGAGGCTAAAT-1" "SRR13927738#TGTTAGGTCGGTCAGC-1"
[3] "SRR13927738#CTCTCAGCATGTCCCT-1" "SRR13927738#ACTATTCGTCCAACCG-1"
[5] "SRR13927738#TGGTCCTTCTACTGCC-1" "SRR13927738#ATTGTCTAGACTAGGC-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-13 05:24:06.046503 : featureDF SRR13927738, Class = DFrame
DataFrame with 23525 rows and 6 columns
      seqnames     start       end    strand        name       idx
         <Rle> <integer> <integer> <integer> <character> <integer>
1         chr1     69091     70008         1       OR4F5         1
2         chr1    817371    819837         1      FAM87B         2
3         chr1    825138    859446         1   LINC01128         3
4         chr1    827522    826206         2   LINC00115         4
5         chr1    876903    868071         2      FAM41C         5
...        ...       ...       ...       ...         ...       ...
23521     chrX 155264589 155258241         2      RAB39B       923
23522     chrX 155334657 155276211         2       CLIC2       924
23523     chrX 155380787 155381134         1      H2AFB1       925
23524     chrX 155382115 155383230         1        F8A2       926
23525     chrX 155669944 155490115         2       TMLHE       927

2023-11-13 05:24:06.065006 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927738.arrow, 0.006 mins elapsed.

2023-11-13 05:24:25.049987 : mat SRR13927738, Class = dgCMatrix
mat SRR13927738: nRows = 23525, nCols = 3736
mat SRR13927738: NonZeroEntries = 41898897, EntryRange = [ 9.73975135001421e-05 , 758.763031303812 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927738#TAGCATGAGGCTAAAT-1 SRR13927738#TGTTAGGTCGGTCAGC-1
[1,]                    .                                .        
[2,]                    .                                .        
[3,]                    0.002346352                      0.0084575
[4,]                    0.005836413                      .        
[5,]                    .                                .        
     SRR13927738#CTCTCAGCATGTCCCT-1 SRR13927738#ACTATTCGTCCAACCG-1
[1,]                     .                             .          
[2,]                     .                             .          
[3,]                     .                             0.004081077
[4,]                     0.02001583                    .          
[5,]                     .                             .          
     SRR13927738#TGGTCCTTCTACTGCC-1
[1,]                    .          
[2,]                    .          
[3,]                    0.007640877
[4,]                    0.015048796
[5,]                    .          


2023-11-13 05:24:25.412877 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927738.arrow, 0.329 mins elapsed.

2023-11-13 05:24:25.763298 : se SRR13927738, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 3736 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(3736): SRR13927738#TAGCATGAGGCTAAAT-1
  SRR13927738#TGTTAGGTCGGTCAGC-1 ... SRR13927738#GCTCCTATCCGTGCGA-1
  SRR13927738#GGCACGTAGTCGATAA-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-13 05:24:25.768244 : Completed GeneIntegrationMatrix : SRR13927738(2 of 4), 1.037 mins elapsed.
2023-11-13 05:24:25.76864 : Reading GeneIntegrationMatrix : SRR13927735(3 of 4), 1.037 mins elapsed.

2023-11-13 05:24:25.787712 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927735 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927735.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3813
[1] "SRR13927735#TTATGTCTCCAGGTAT-1" "SRR13927735#TATTGCTCATCAGAAA-1"
[3] "SRR13927735#TTCGATTGTAGGGTTG-1" "SRR13927735#CATTCATTCGGATGTT-1"
[5] "SRR13927735#ACGTTAGGTCAACTGT-1" "SRR13927735#AAATGCCCAGCAATGG-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-13 05:24:26.282065 : featureDF SRR13927735, Class = DFrame
DataFrame with 23525 rows and 6 columns
      seqnames     start       end    strand        name       idx
         <Rle> <integer> <integer> <integer> <character> <integer>
1         chr1     69091     70008         1       OR4F5         1
2         chr1    817371    819837         1      FAM87B         2
3         chr1    825138    859446         1   LINC01128         3
4         chr1    827522    826206         2   LINC00115         4
5         chr1    876903    868071         2      FAM41C         5
...        ...       ...       ...       ...         ...       ...
23521     chrX 155264589 155258241         2      RAB39B       923
23522     chrX 155334657 155276211         2       CLIC2       924
23523     chrX 155380787 155381134         1      H2AFB1       925
23524     chrX 155382115 155383230         1        F8A2       926
23525     chrX 155669944 155490115         2       TMLHE       927

2023-11-13 05:24:26.312454 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927735.arrow, 0.009 mins elapsed.

2023-11-13 05:24:46.471365 : mat SRR13927735, Class = dgCMatrix
mat SRR13927735: nRows = 23525, nCols = 3813
mat SRR13927735: NonZeroEntries = 45535634, EntryRange = [ 7.20715653308254e-05 , 673.990433360078 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927735#TTATGTCTCCAGGTAT-1 SRR13927735#TATTGCTCATCAGAAA-1
[1,]                    .                              .          
[2,]                    .                              .          
[3,]                    0.070892732                    0.057009268
[4,]                    0.008800039                    0.002978408
[5,]                    .                              .          
     SRR13927735#TTCGATTGTAGGGTTG-1 SRR13927735#CATTCATTCGGATGTT-1
[1,]                     .                               .        
[2,]                     .                               .        
[3,]                     0.05523011                      0.0248301
[4,]                     .                               .        
[5,]                     .                               .        
     SRR13927735#ACGTTAGGTCAACTGT-1
[1,]                     .         
[2,]                     .         
[3,]                     0.05384571
[4,]                     0.01228688
[5,]                     .         


2023-11-13 05:24:46.821856 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927735.arrow, 0.351 mins elapsed.

2023-11-13 05:24:47.098036 : se SRR13927735, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 3813 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(3813): SRR13927735#TTATGTCTCCAGGTAT-1
  SRR13927735#TATTGCTCATCAGAAA-1 ... SRR13927735#GATCATGGTGTTGTTG-1
  SRR13927735#GTCCATCCATTTCACT-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-13 05:24:47.102941 : Completed GeneIntegrationMatrix : SRR13927735(3 of 4), 1.392 mins elapsed.
2023-11-13 05:24:47.103347 : Reading GeneIntegrationMatrix : SRR13927736(4 of 4), 1.392 mins elapsed.

2023-11-13 05:24:47.116586 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927736 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927736.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3534
[1] "SRR13927736#GCTGTTCTCCTTGACC-1" "SRR13927736#AGTTACGAGTAGGTCG-1"
[3] "SRR13927736#AGTTTGGAGCATTCCA-1" "SRR13927736#GCACGGTTCGGCAATT-1"
[5] "SRR13927736#TCAGCTCTCTCTATTG-1" "SRR13927736#CCTTAATTCTACTTTG-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-13 05:24:47.732725 : featureDF SRR13927736, Class = DFrame
DataFrame with 23525 rows and 6 columns
      seqnames     start       end    strand        name       idx
         <Rle> <integer> <integer> <integer> <character> <integer>
1         chr1     69091     70008         1       OR4F5         1
2         chr1    817371    819837         1      FAM87B         2
3         chr1    825138    859446         1   LINC01128         3
4         chr1    827522    826206         2   LINC00115         4
5         chr1    876903    868071         2      FAM41C         5
...        ...       ...       ...       ...         ...       ...
23521     chrX 155264589 155258241         2      RAB39B       923
23522     chrX 155334657 155276211         2       CLIC2       924
23523     chrX 155380787 155381134         1      H2AFB1       925
23524     chrX 155382115 155383230         1        F8A2       926
23525     chrX 155669944 155490115         2       TMLHE       927

2023-11-13 05:24:47.764773 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927736.arrow, 0.011 mins elapsed.

2023-11-13 05:25:04.044185 : mat SRR13927736, Class = dgCMatrix
mat SRR13927736: nRows = 23525, nCols = 3534
mat SRR13927736: NonZeroEntries = 39003616, EntryRange = [ 8.77627031078293e-05 , 770.387276519725 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927736#GCTGTTCTCCTTGACC-1 SRR13927736#AGTTACGAGTAGGTCG-1
[1,]                     .                             .          
[2,]                     .                             .          
[3,]                     0.07845381                    0.028344539
[4,]                     0.01395289                    0.008499189
[5,]                     .                             .          
     SRR13927736#AGTTTGGAGCATTCCA-1 SRR13927736#GCACGGTTCGGCAATT-1
[1,]                     .                              .         
[2,]                     .                              .         
[3,]                     0.06535130                     0.06576222
[4,]                     0.01810177                     0.01315728
[5,]                     .                              .         
     SRR13927736#TCAGCTCTCTCTATTG-1
[1,]                    .          
[2,]                    .          
[3,]                    0.016229526
[4,]                    0.003862592
[5,]                    .          


2023-11-13 05:25:05.555387 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927736.arrow, 0.307 mins elapsed.

2023-11-13 05:25:05.827184 : se SRR13927736, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 3534 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(3534): SRR13927736#GCTGTTCTCCTTGACC-1
  SRR13927736#AGTTACGAGTAGGTCG-1 ... SRR13927736#TAAACCGCAAGTCTGT-1
  SRR13927736#TGCATTTGTTCTACCC-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-13 05:25:05.832487 : Completed GeneIntegrationMatrix : SRR13927736(4 of 4), 1.704 mins elapsed.
2023-11-13 05:25:05.833478 : Organizing colData, 1.704 mins elapsed.
2023-11-13 05:25:05.990471 : Organizing rowData, 1.707 mins elapsed.
2023-11-13 05:25:05.993644 : Organizing rowRanges, 1.707 mins elapsed.
2023-11-13 05:25:05.999488 : Organizing Assays (1 of 1), 1.707 mins elapsed.
2023-11-13 05:25:13.178097 : Constructing SummarizedExperiment, 1.827 mins elapsed.
2023-11-13 05:25:14.452706 : Finished Matrix Creation, 1.848 mins elapsed.

2023-11-13 05:25:14.455095 : getMatrixFromProject Input-Parameters, Class = list

getMatrixFromProject Input-Parameters$ArchRProj: length = 1

getMatrixFromProject Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromProject Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromProject Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromProject Input-Parameters$verbose: length = 1
[1] TRUE


getMatrixFromProject Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromProject Input-Parameters$threads: length = 1
[1] 1


getMatrixFromProject Input-Parameters$logFile: length = 1
[1] "x"


2023-11-13 05:25:15.852462 : Reading PeakMatrix : SRR13927737(1 of 4), 0.023 mins elapsed.

2023-11-13 05:25:15.887017 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927737 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927737.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 4439
[1] "SRR13927737#GTCACGGGTCAACAGG-1" "SRR13927737#TGCCTGTTCTTACTCA-1"
[3] "SRR13927737#TTGCGGGTCAGGTCTA-1" "SRR13927737#GGTTGCGCAAAGCTTC-1"
[5] "SRR13927737#CTTAATCTCGCGCTGA-1" "SRR13927737#TGTACGATCGGAAAGC-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-13 05:25:16.734446 : featureDF SRR13927737, Class = DFrame
DataFrame with 210721 rows and 4 columns
       seqnames       idx     start       end
          <Rle> <integer> <integer> <integer>
1          chr1         1    804677    805177
2          chr1         2    811866    812366
3          chr1         3    817090    817590
4          chr1         4    818831    819331
5          chr1         5    826543    827043
...         ...       ...       ...       ...
210717     chrX      5777 155898525 155899025
210718     chrX      5778 155899804 155900304
210719     chrX      5779 155956608 155957108
210720     chrX      5780 155966827 155967327
210721     chrX      5781 155971625 155972125

2023-11-13 05:25:16.837459 : Getting PeakMatrix from ArrowFile : SRR13927737.arrow, 0.016 mins elapsed.

2023-11-13 05:25:29.277086 : mat SRR13927737, Class = dgCMatrix
mat SRR13927737: nRows = 210721, nCols = 4439
mat SRR13927737: NonZeroEntries = 29963541, EntryRange = [ 1 , 4 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927737#GTCACGGGTCAACAGG-1 SRR13927737#TGCCTGTTCTTACTCA-1
[1,]                              .                              .
[2,]                              .                              2
[3,]                              .                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927737#TTGCGGGTCAGGTCTA-1 SRR13927737#GGTTGCGCAAAGCTTC-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              2                              .
[4,]                              .                              .
[5,]                              4                              .
     SRR13927737#CTTAATCTCGCGCTGA-1
[1,]                              .
[2,]                              .
[3,]                              .
[4,]                              .
[5,]                              .


2023-11-13 05:25:29.461382 : Organizing SE PeakMatrix from ArrowFile : SRR13927737.arrow, 0.226 mins elapsed.

2023-11-13 05:25:29.767586 : se SRR13927737, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 4439 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(4439): SRR13927737#GTCACGGGTCAACAGG-1
  SRR13927737#TGCCTGTTCTTACTCA-1 ... SRR13927737#TTGCTATAGATGTTCC-1
  SRR13927737#GCTCACTAGTTAGCGG-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-13 05:25:29.777425 : Completed PeakMatrix : SRR13927737(1 of 4), 0.255 mins elapsed.
2023-11-13 05:25:29.777854 : Reading PeakMatrix : SRR13927738(2 of 4), 0.255 mins elapsed.

2023-11-13 05:25:29.790448 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927738 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927738.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3736
[1] "SRR13927738#TAGCATGAGGCTAAAT-1" "SRR13927738#TGTTAGGTCGGTCAGC-1"
[3] "SRR13927738#CTCTCAGCATGTCCCT-1" "SRR13927738#ACTATTCGTCCAACCG-1"
[5] "SRR13927738#TGGTCCTTCTACTGCC-1" "SRR13927738#ATTGTCTAGACTAGGC-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-13 05:25:30.588141 : featureDF SRR13927738, Class = DFrame
DataFrame with 210721 rows and 4 columns
       seqnames       idx     start       end
          <Rle> <integer> <integer> <integer>
1          chr1         1    804677    805177
2          chr1         2    811866    812366
3          chr1         3    817090    817590
4          chr1         4    818831    819331
5          chr1         5    826543    827043
...         ...       ...       ...       ...
210717     chrX      5777 155898525 155899025
210718     chrX      5778 155899804 155900304
210719     chrX      5779 155956608 155957108
210720     chrX      5780 155966827 155967327
210721     chrX      5781 155971625 155972125

2023-11-13 05:25:30.690538 : Getting PeakMatrix from ArrowFile : SRR13927738.arrow, 0.015 mins elapsed.

2023-11-13 05:25:42.520498 : mat SRR13927738, Class = dgCMatrix
mat SRR13927738: nRows = 210721, nCols = 3736
mat SRR13927738: NonZeroEntries = 30165921, EntryRange = [ 1 , 4 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927738#TAGCATGAGGCTAAAT-1 SRR13927738#TGTTAGGTCGGTCAGC-1
[1,]                              2                              .
[2,]                              .                              .
[3,]                              .                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927738#CTCTCAGCATGTCCCT-1 SRR13927738#ACTATTCGTCCAACCG-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              .                              .
[4,]                              .                              .
[5,]                              2                              .
     SRR13927738#TGGTCCTTCTACTGCC-1
[1,]                              .
[2,]                              .
[3,]                              2
[4,]                              .
[5,]                              .


2023-11-13 05:25:42.749537 : Organizing SE PeakMatrix from ArrowFile : SRR13927738.arrow, 0.216 mins elapsed.

2023-11-13 05:25:43.033108 : se SRR13927738, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 3736 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(3736): SRR13927738#TAGCATGAGGCTAAAT-1
  SRR13927738#TGTTAGGTCGGTCAGC-1 ... SRR13927738#GCTCCTATCCGTGCGA-1
  SRR13927738#GGCACGTAGTCGATAA-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-13 05:25:43.040135 : Completed PeakMatrix : SRR13927738(2 of 4), 0.476 mins elapsed.
2023-11-13 05:25:43.040558 : Reading PeakMatrix : SRR13927735(3 of 4), 0.476 mins elapsed.

2023-11-13 05:25:43.053016 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927735 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927735.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3813
[1] "SRR13927735#TTATGTCTCCAGGTAT-1" "SRR13927735#TATTGCTCATCAGAAA-1"
[3] "SRR13927735#TTCGATTGTAGGGTTG-1" "SRR13927735#CATTCATTCGGATGTT-1"
[5] "SRR13927735#ACGTTAGGTCAACTGT-1" "SRR13927735#AAATGCCCAGCAATGG-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-13 05:25:43.663469 : featureDF SRR13927735, Class = DFrame
DataFrame with 210721 rows and 4 columns
       seqnames       idx     start       end
          <Rle> <integer> <integer> <integer>
1          chr1         1    804677    805177
2          chr1         2    811866    812366
3          chr1         3    817090    817590
4          chr1         4    818831    819331
5          chr1         5    826543    827043
...         ...       ...       ...       ...
210717     chrX      5777 155898525 155899025
210718     chrX      5778 155899804 155900304
210719     chrX      5779 155956608 155957108
210720     chrX      5780 155966827 155967327
210721     chrX      5781 155971625 155972125

2023-11-13 05:25:43.715048 : Getting PeakMatrix from ArrowFile : SRR13927735.arrow, 0.011 mins elapsed.

2023-11-13 05:25:53.716733 : mat SRR13927735, Class = dgCMatrix
mat SRR13927735: nRows = 210721, nCols = 3813
mat SRR13927735: NonZeroEntries = 26038084, EntryRange = [ 1 , 4 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927735#TTATGTCTCCAGGTAT-1 SRR13927735#TATTGCTCATCAGAAA-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              1                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927735#TTCGATTGTAGGGTTG-1 SRR13927735#CATTCATTCGGATGTT-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              .                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927735#ACGTTAGGTCAACTGT-1
[1,]                              2
[2,]                              .
[3,]                              .
[4,]                              .
[5,]                              .


2023-11-13 05:25:53.916479 : Organizing SE PeakMatrix from ArrowFile : SRR13927735.arrow, 0.181 mins elapsed.

2023-11-13 05:25:54.193991 : se SRR13927735, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 3813 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(3813): SRR13927735#TTATGTCTCCAGGTAT-1
  SRR13927735#TATTGCTCATCAGAAA-1 ... SRR13927735#GATCATGGTGTTGTTG-1
  SRR13927735#GTCCATCCATTTCACT-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-13 05:25:54.201056 : Completed PeakMatrix : SRR13927735(3 of 4), 0.662 mins elapsed.
2023-11-13 05:25:54.201489 : Reading PeakMatrix : SRR13927736(4 of 4), 0.662 mins elapsed.

2023-11-13 05:25:54.221851 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927736 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927736.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3534
[1] "SRR13927736#GCTGTTCTCCTTGACC-1" "SRR13927736#AGTTACGAGTAGGTCG-1"
[3] "SRR13927736#AGTTTGGAGCATTCCA-1" "SRR13927736#GCACGGTTCGGCAATT-1"
[5] "SRR13927736#TCAGCTCTCTCTATTG-1" "SRR13927736#CCTTAATTCTACTTTG-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-13 05:25:54.785621 : featureDF SRR13927736, Class = DFrame
DataFrame with 210721 rows and 4 columns
       seqnames       idx     start       end
          <Rle> <integer> <integer> <integer>
1          chr1         1    804677    805177
2          chr1         2    811866    812366
3          chr1         3    817090    817590
4          chr1         4    818831    819331
5          chr1         5    826543    827043
...         ...       ...       ...       ...
210717     chrX      5777 155898525 155899025
210718     chrX      5778 155899804 155900304
210719     chrX      5779 155956608 155957108
210720     chrX      5780 155966827 155967327
210721     chrX      5781 155971625 155972125

2023-11-13 05:25:54.837571 : Getting PeakMatrix from ArrowFile : SRR13927736.arrow, 0.01 mins elapsed.

2023-11-13 05:26:04.913018 : mat SRR13927736, Class = dgCMatrix
mat SRR13927736: nRows = 210721, nCols = 3534
mat SRR13927736: NonZeroEntries = 25583439, EntryRange = [ 1 , 4 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927736#GCTGTTCTCCTTGACC-1 SRR13927736#AGTTACGAGTAGGTCG-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              .                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927736#AGTTTGGAGCATTCCA-1 SRR13927736#GCACGGTTCGGCAATT-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              1                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927736#TCAGCTCTCTCTATTG-1
[1,]                              .
[2,]                              .
[3,]                              .
[4,]                              .
[5,]                              .


2023-11-13 05:26:05.109058 : Organizing SE PeakMatrix from ArrowFile : SRR13927736.arrow, 0.181 mins elapsed.

2023-11-13 05:26:05.382068 : se SRR13927736, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 3534 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(3534): SRR13927736#GCTGTTCTCCTTGACC-1
  SRR13927736#AGTTACGAGTAGGTCG-1 ... SRR13927736#TAAACCGCAAGTCTGT-1
  SRR13927736#TGCATTTGTTCTACCC-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-13 05:26:05.389011 : Completed PeakMatrix : SRR13927736(4 of 4), 0.849 mins elapsed.
2023-11-13 05:26:05.389911 : Organizing colData, 0.849 mins elapsed.
2023-11-13 05:26:05.523402 : Organizing rowData, 0.851 mins elapsed.
2023-11-13 05:26:05.531371 : Organizing rowRanges, 0.851 mins elapsed.
2023-11-13 05:26:05.542171 : Organizing Assays (1 of 1), 0.851 mins elapsed.
2023-11-13 05:26:08.613307 : Constructing SummarizedExperiment, 0.903 mins elapsed.
2023-11-13 05:26:21.986284 : Finished Matrix Creation, 1.126 mins elapsed.

2023-11-14 04:34:54.042194 : getMatrixFromProject Input-Parameters, Class = list

getMatrixFromProject Input-Parameters$ArchRProj: length = 1

getMatrixFromProject Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromProject Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromProject Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromProject Input-Parameters$verbose: length = 1
[1] TRUE


getMatrixFromProject Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromProject Input-Parameters$threads: length = 1
[1] 1


getMatrixFromProject Input-Parameters$logFile: length = 1
[1] "x"


2023-11-14 04:34:58.670851 : Reading GeneIntegrationMatrix : SRR13927737(1 of 4), 0.077 mins elapsed.

2023-11-14 04:34:58.712178 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927737 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927737.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 4439
[1] "SRR13927737#GTCACGGGTCAACAGG-1" "SRR13927737#TGCCTGTTCTTACTCA-1"
[3] "SRR13927737#TTGCGGGTCAGGTCTA-1" "SRR13927737#GGTTGCGCAAAGCTTC-1"
[5] "SRR13927737#CTTAATCTCGCGCTGA-1" "SRR13927737#TGTACGATCGGAAAGC-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-14 04:34:59.677593 : featureDF SRR13927737, Class = DFrame
DataFrame with 23525 rows and 6 columns
      seqnames     start       end    strand        name       idx
         <Rle> <integer> <integer> <integer> <character> <integer>
1         chr1     69091     70008         1       OR4F5         1
2         chr1    817371    819837         1      FAM87B         2
3         chr1    825138    859446         1   LINC01128         3
4         chr1    827522    826206         2   LINC00115         4
5         chr1    876903    868071         2      FAM41C         5
...        ...       ...       ...       ...         ...       ...
23521     chrX 155264589 155258241         2      RAB39B       923
23522     chrX 155334657 155276211         2       CLIC2       924
23523     chrX 155380787 155381134         1      H2AFB1       925
23524     chrX 155382115 155383230         1        F8A2       926
23525     chrX 155669944 155490115         2       TMLHE       927

2023-11-14 04:34:59.70981 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927737.arrow, 0.017 mins elapsed.

2023-11-14 04:35:20.648855 : mat SRR13927737, Class = dgCMatrix
mat SRR13927737: nRows = 23525, nCols = 4439
mat SRR13927737: NonZeroEntries = 50277774, EntryRange = [ 0.000148130054904683 , 376.345870564644 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927737#GTCACGGGTCAACAGG-1 SRR13927737#TGCCTGTTCTTACTCA-1
[1,]                      .                             .         
[2,]                      .                             .         
[3,]                      0.1057576                     0.09630696
[4,]                      .                             0.01806508
[5,]                      .                             .         
     SRR13927737#TTGCGGGTCAGGTCTA-1 SRR13927737#GGTTGCGCAAAGCTTC-1
[1,]                     .                              .         
[2,]                     .                              .         
[3,]                     0.09113499                     0.08396973
[4,]                     0.02096304                     0.01727603
[5,]                     .                              .         
     SRR13927737#CTTAATCTCGCGCTGA-1
[1,]                     .         
[2,]                     .         
[3,]                     0.08327649
[4,]                     0.01368917
[5,]                     .         


2023-11-14 04:35:21.129591 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927737.arrow, 0.374 mins elapsed.

2023-11-14 04:35:21.45731 : se SRR13927737, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 4439 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(4439): SRR13927737#GTCACGGGTCAACAGG-1
  SRR13927737#TGCCTGTTCTTACTCA-1 ... SRR13927737#TTGCTATAGATGTTCC-1
  SRR13927737#GCTCACTAGTTAGCGG-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-14 04:35:21.464407 : Completed GeneIntegrationMatrix : SRR13927737(1 of 4), 0.457 mins elapsed.
2023-11-14 04:35:21.464982 : Reading GeneIntegrationMatrix : SRR13927738(2 of 4), 0.457 mins elapsed.

2023-11-14 04:35:21.489317 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927738 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927738.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3736
[1] "SRR13927738#TAGCATGAGGCTAAAT-1" "SRR13927738#TGTTAGGTCGGTCAGC-1"
[3] "SRR13927738#CTCTCAGCATGTCCCT-1" "SRR13927738#ACTATTCGTCCAACCG-1"
[5] "SRR13927738#TGGTCCTTCTACTGCC-1" "SRR13927738#ATTGTCTAGACTAGGC-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-14 04:35:22.657015 : featureDF SRR13927738, Class = DFrame
DataFrame with 23525 rows and 6 columns
      seqnames     start       end    strand        name       idx
         <Rle> <integer> <integer> <integer> <character> <integer>
1         chr1     69091     70008         1       OR4F5         1
2         chr1    817371    819837         1      FAM87B         2
3         chr1    825138    859446         1   LINC01128         3
4         chr1    827522    826206         2   LINC00115         4
5         chr1    876903    868071         2      FAM41C         5
...        ...       ...       ...       ...         ...       ...
23521     chrX 155264589 155258241         2      RAB39B       923
23522     chrX 155334657 155276211         2       CLIC2       924
23523     chrX 155380787 155381134         1      H2AFB1       925
23524     chrX 155382115 155383230         1        F8A2       926
23525     chrX 155669944 155490115         2       TMLHE       927

2023-11-14 04:35:22.690381 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927738.arrow, 0.02 mins elapsed.

2023-11-14 04:35:42.787937 : mat SRR13927738, Class = dgCMatrix
mat SRR13927738: nRows = 23525, nCols = 3736
mat SRR13927738: NonZeroEntries = 41898897, EntryRange = [ 9.73975135001421e-05 , 758.763031303812 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927738#TAGCATGAGGCTAAAT-1 SRR13927738#TGTTAGGTCGGTCAGC-1
[1,]                    .                                .        
[2,]                    .                                .        
[3,]                    0.002346352                      0.0084575
[4,]                    0.005836413                      .        
[5,]                    .                                .        
     SRR13927738#CTCTCAGCATGTCCCT-1 SRR13927738#ACTATTCGTCCAACCG-1
[1,]                     .                             .          
[2,]                     .                             .          
[3,]                     .                             0.004081077
[4,]                     0.02001583                    .          
[5,]                     .                             .          
     SRR13927738#TGGTCCTTCTACTGCC-1
[1,]                    .          
[2,]                    .          
[3,]                    0.007640877
[4,]                    0.015048796
[5,]                    .          


2023-11-14 04:35:43.1193 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927738.arrow, 0.361 mins elapsed.

2023-11-14 04:35:43.430782 : se SRR13927738, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 3736 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(3736): SRR13927738#TAGCATGAGGCTAAAT-1
  SRR13927738#TGTTAGGTCGGTCAGC-1 ... SRR13927738#GCTCCTATCCGTGCGA-1
  SRR13927738#GGCACGTAGTCGATAA-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-14 04:35:43.436557 : Completed GeneIntegrationMatrix : SRR13927738(2 of 4), 0.823 mins elapsed.
2023-11-14 04:35:43.437079 : Reading GeneIntegrationMatrix : SRR13927735(3 of 4), 0.823 mins elapsed.

2023-11-14 04:35:43.470225 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927735 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927735.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3813
[1] "SRR13927735#TTATGTCTCCAGGTAT-1" "SRR13927735#TATTGCTCATCAGAAA-1"
[3] "SRR13927735#TTCGATTGTAGGGTTG-1" "SRR13927735#CATTCATTCGGATGTT-1"
[5] "SRR13927735#ACGTTAGGTCAACTGT-1" "SRR13927735#AAATGCCCAGCAATGG-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-14 04:35:44.64561 : featureDF SRR13927735, Class = DFrame
DataFrame with 23525 rows and 6 columns
      seqnames     start       end    strand        name       idx
         <Rle> <integer> <integer> <integer> <character> <integer>
1         chr1     69091     70008         1       OR4F5         1
2         chr1    817371    819837         1      FAM87B         2
3         chr1    825138    859446         1   LINC01128         3
4         chr1    827522    826206         2   LINC00115         4
5         chr1    876903    868071         2      FAM41C         5
...        ...       ...       ...       ...         ...       ...
23521     chrX 155264589 155258241         2      RAB39B       923
23522     chrX 155334657 155276211         2       CLIC2       924
23523     chrX 155380787 155381134         1      H2AFB1       925
23524     chrX 155382115 155383230         1        F8A2       926
23525     chrX 155669944 155490115         2       TMLHE       927

2023-11-14 04:35:44.679237 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927735.arrow, 0.02 mins elapsed.

2023-11-14 04:36:04.718519 : mat SRR13927735, Class = dgCMatrix
mat SRR13927735: nRows = 23525, nCols = 3813
mat SRR13927735: NonZeroEntries = 45535634, EntryRange = [ 7.20715653308254e-05 , 673.990433360078 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927735#TTATGTCTCCAGGTAT-1 SRR13927735#TATTGCTCATCAGAAA-1
[1,]                    .                              .          
[2,]                    .                              .          
[3,]                    0.070892732                    0.057009268
[4,]                    0.008800039                    0.002978408
[5,]                    .                              .          
     SRR13927735#TTCGATTGTAGGGTTG-1 SRR13927735#CATTCATTCGGATGTT-1
[1,]                     .                               .        
[2,]                     .                               .        
[3,]                     0.05523011                      0.0248301
[4,]                     .                               .        
[5,]                     .                               .        
     SRR13927735#ACGTTAGGTCAACTGT-1
[1,]                     .         
[2,]                     .         
[3,]                     0.05384571
[4,]                     0.01228688
[5,]                     .         


2023-11-14 04:36:05.108048 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927735.arrow, 0.361 mins elapsed.

2023-11-14 04:36:05.400267 : se SRR13927735, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 3813 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(3813): SRR13927735#TTATGTCTCCAGGTAT-1
  SRR13927735#TATTGCTCATCAGAAA-1 ... SRR13927735#GATCATGGTGTTGTTG-1
  SRR13927735#GTCCATCCATTTCACT-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-14 04:36:05.406029 : Completed GeneIntegrationMatrix : SRR13927735(3 of 4), 1.189 mins elapsed.
2023-11-14 04:36:05.406534 : Reading GeneIntegrationMatrix : SRR13927736(4 of 4), 1.189 mins elapsed.

2023-11-14 04:36:05.429178 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927736 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927736.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3534
[1] "SRR13927736#GCTGTTCTCCTTGACC-1" "SRR13927736#AGTTACGAGTAGGTCG-1"
[3] "SRR13927736#AGTTTGGAGCATTCCA-1" "SRR13927736#GCACGGTTCGGCAATT-1"
[5] "SRR13927736#TCAGCTCTCTCTATTG-1" "SRR13927736#CCTTAATTCTACTTTG-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-14 04:36:06.44368 : featureDF SRR13927736, Class = DFrame
DataFrame with 23525 rows and 6 columns
      seqnames     start       end    strand        name       idx
         <Rle> <integer> <integer> <integer> <character> <integer>
1         chr1     69091     70008         1       OR4F5         1
2         chr1    817371    819837         1      FAM87B         2
3         chr1    825138    859446         1   LINC01128         3
4         chr1    827522    826206         2   LINC00115         4
5         chr1    876903    868071         2      FAM41C         5
...        ...       ...       ...       ...         ...       ...
23521     chrX 155264589 155258241         2      RAB39B       923
23522     chrX 155334657 155276211         2       CLIC2       924
23523     chrX 155380787 155381134         1      H2AFB1       925
23524     chrX 155382115 155383230         1        F8A2       926
23525     chrX 155669944 155490115         2       TMLHE       927

2023-11-14 04:36:06.47891 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927736.arrow, 0.017 mins elapsed.

2023-11-14 04:36:23.242969 : mat SRR13927736, Class = dgCMatrix
mat SRR13927736: nRows = 23525, nCols = 3534
mat SRR13927736: NonZeroEntries = 39003616, EntryRange = [ 8.77627031078293e-05 , 770.387276519725 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927736#GCTGTTCTCCTTGACC-1 SRR13927736#AGTTACGAGTAGGTCG-1
[1,]                     .                             .          
[2,]                     .                             .          
[3,]                     0.07845381                    0.028344539
[4,]                     0.01395289                    0.008499189
[5,]                     .                             .          
     SRR13927736#AGTTTGGAGCATTCCA-1 SRR13927736#GCACGGTTCGGCAATT-1
[1,]                     .                              .         
[2,]                     .                              .         
[3,]                     0.06535130                     0.06576222
[4,]                     0.01810177                     0.01315728
[5,]                     .                              .         
     SRR13927736#TCAGCTCTCTCTATTG-1
[1,]                    .          
[2,]                    .          
[3,]                    0.016229526
[4,]                    0.003862592
[5,]                    .          


2023-11-14 04:36:23.567478 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927736.arrow, 0.302 mins elapsed.

2023-11-14 04:36:23.832512 : se SRR13927736, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 3534 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(3534): SRR13927736#GCTGTTCTCCTTGACC-1
  SRR13927736#AGTTACGAGTAGGTCG-1 ... SRR13927736#TAAACCGCAAGTCTGT-1
  SRR13927736#TGCATTTGTTCTACCC-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-14 04:36:23.837846 : Completed GeneIntegrationMatrix : SRR13927736(4 of 4), 1.497 mins elapsed.
2023-11-14 04:36:23.839036 : Organizing colData, 1.497 mins elapsed.
2023-11-14 04:36:23.986064 : Organizing rowData, 1.499 mins elapsed.
2023-11-14 04:36:23.989325 : Organizing rowRanges, 1.499 mins elapsed.
2023-11-14 04:36:23.994774 : Organizing Assays (1 of 1), 1.499 mins elapsed.
2023-11-14 04:36:33.433822 : Constructing SummarizedExperiment, 1.657 mins elapsed.
2023-11-14 04:36:36.218726 : Finished Matrix Creation, 1.703 mins elapsed.

2023-11-14 04:36:36.221188 : getMatrixFromProject Input-Parameters, Class = list

getMatrixFromProject Input-Parameters$ArchRProj: length = 1

getMatrixFromProject Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromProject Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromProject Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromProject Input-Parameters$verbose: length = 1
[1] TRUE


getMatrixFromProject Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromProject Input-Parameters$threads: length = 1
[1] 1


getMatrixFromProject Input-Parameters$logFile: length = 1
[1] "x"


2023-11-14 04:36:38.846849 : Reading PeakMatrix : SRR13927737(1 of 4), 0.044 mins elapsed.

2023-11-14 04:36:38.881399 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927737 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927737.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 4439
[1] "SRR13927737#GTCACGGGTCAACAGG-1" "SRR13927737#TGCCTGTTCTTACTCA-1"
[3] "SRR13927737#TTGCGGGTCAGGTCTA-1" "SRR13927737#GGTTGCGCAAAGCTTC-1"
[5] "SRR13927737#CTTAATCTCGCGCTGA-1" "SRR13927737#TGTACGATCGGAAAGC-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-14 04:36:39.877393 : featureDF SRR13927737, Class = DFrame
DataFrame with 210721 rows and 4 columns
       seqnames       idx     start       end
          <Rle> <integer> <integer> <integer>
1          chr1         1    804677    805177
2          chr1         2    811866    812366
3          chr1         3    817090    817590
4          chr1         4    818831    819331
5          chr1         5    826543    827043
...         ...       ...       ...       ...
210717     chrX      5777 155898525 155899025
210718     chrX      5778 155899804 155900304
210719     chrX      5779 155956608 155957108
210720     chrX      5780 155966827 155967327
210721     chrX      5781 155971625 155972125

2023-11-14 04:36:39.980603 : Getting PeakMatrix from ArrowFile : SRR13927737.arrow, 0.018 mins elapsed.

2023-11-14 04:36:53.249179 : mat SRR13927737, Class = dgCMatrix
mat SRR13927737: nRows = 210721, nCols = 4439
mat SRR13927737: NonZeroEntries = 29963541, EntryRange = [ 1 , 4 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927737#GTCACGGGTCAACAGG-1 SRR13927737#TGCCTGTTCTTACTCA-1
[1,]                              .                              .
[2,]                              .                              2
[3,]                              .                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927737#TTGCGGGTCAGGTCTA-1 SRR13927737#GGTTGCGCAAAGCTTC-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              2                              .
[4,]                              .                              .
[5,]                              4                              .
     SRR13927737#CTTAATCTCGCGCTGA-1
[1,]                              .
[2,]                              .
[3,]                              .
[4,]                              .
[5,]                              .


2023-11-14 04:36:53.481856 : Organizing SE PeakMatrix from ArrowFile : SRR13927737.arrow, 0.243 mins elapsed.

2023-11-14 04:36:53.79917 : se SRR13927737, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 4439 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(4439): SRR13927737#GTCACGGGTCAACAGG-1
  SRR13927737#TGCCTGTTCTTACTCA-1 ... SRR13927737#TTGCTATAGATGTTCC-1
  SRR13927737#GCTCACTAGTTAGCGG-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-14 04:36:53.80931 : Completed PeakMatrix : SRR13927737(1 of 4), 0.293 mins elapsed.
2023-11-14 04:36:53.809827 : Reading PeakMatrix : SRR13927738(2 of 4), 0.293 mins elapsed.

2023-11-14 04:36:53.826733 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927738 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927738.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3736
[1] "SRR13927738#TAGCATGAGGCTAAAT-1" "SRR13927738#TGTTAGGTCGGTCAGC-1"
[3] "SRR13927738#CTCTCAGCATGTCCCT-1" "SRR13927738#ACTATTCGTCCAACCG-1"
[5] "SRR13927738#TGGTCCTTCTACTGCC-1" "SRR13927738#ATTGTCTAGACTAGGC-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-14 04:36:54.592878 : featureDF SRR13927738, Class = DFrame
DataFrame with 210721 rows and 4 columns
       seqnames       idx     start       end
          <Rle> <integer> <integer> <integer>
1          chr1         1    804677    805177
2          chr1         2    811866    812366
3          chr1         3    817090    817590
4          chr1         4    818831    819331
5          chr1         5    826543    827043
...         ...       ...       ...       ...
210717     chrX      5777 155898525 155899025
210718     chrX      5778 155899804 155900304
210719     chrX      5779 155956608 155957108
210720     chrX      5780 155966827 155967327
210721     chrX      5781 155971625 155972125

2023-11-14 04:36:54.646801 : Getting PeakMatrix from ArrowFile : SRR13927738.arrow, 0.014 mins elapsed.

2023-11-14 04:37:07.709293 : mat SRR13927738, Class = dgCMatrix
mat SRR13927738: nRows = 210721, nCols = 3736
mat SRR13927738: NonZeroEntries = 30165921, EntryRange = [ 1 , 4 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927738#TAGCATGAGGCTAAAT-1 SRR13927738#TGTTAGGTCGGTCAGC-1
[1,]                              2                              .
[2,]                              .                              .
[3,]                              .                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927738#CTCTCAGCATGTCCCT-1 SRR13927738#ACTATTCGTCCAACCG-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              .                              .
[4,]                              .                              .
[5,]                              2                              .
     SRR13927738#TGGTCCTTCTACTGCC-1
[1,]                              .
[2,]                              .
[3,]                              2
[4,]                              .
[5,]                              .


2023-11-14 04:37:07.903662 : Organizing SE PeakMatrix from ArrowFile : SRR13927738.arrow, 0.235 mins elapsed.

2023-11-14 04:37:08.19677 : se SRR13927738, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 3736 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(3736): SRR13927738#TAGCATGAGGCTAAAT-1
  SRR13927738#TGTTAGGTCGGTCAGC-1 ... SRR13927738#GCTCCTATCCGTGCGA-1
  SRR13927738#GGCACGTAGTCGATAA-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-14 04:37:08.204336 : Completed PeakMatrix : SRR13927738(2 of 4), 0.533 mins elapsed.
2023-11-14 04:37:08.204883 : Reading PeakMatrix : SRR13927735(3 of 4), 0.533 mins elapsed.

2023-11-14 04:37:08.259334 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927735 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927735.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3813
[1] "SRR13927735#TTATGTCTCCAGGTAT-1" "SRR13927735#TATTGCTCATCAGAAA-1"
[3] "SRR13927735#TTCGATTGTAGGGTTG-1" "SRR13927735#CATTCATTCGGATGTT-1"
[5] "SRR13927735#ACGTTAGGTCAACTGT-1" "SRR13927735#AAATGCCCAGCAATGG-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-14 04:37:09.328358 : featureDF SRR13927735, Class = DFrame
DataFrame with 210721 rows and 4 columns
       seqnames       idx     start       end
          <Rle> <integer> <integer> <integer>
1          chr1         1    804677    805177
2          chr1         2    811866    812366
3          chr1         3    817090    817590
4          chr1         4    818831    819331
5          chr1         5    826543    827043
...         ...       ...       ...       ...
210717     chrX      5777 155898525 155899025
210718     chrX      5778 155899804 155900304
210719     chrX      5779 155956608 155957108
210720     chrX      5780 155966827 155967327
210721     chrX      5781 155971625 155972125

2023-11-14 04:37:09.432828 : Getting PeakMatrix from ArrowFile : SRR13927735.arrow, 0.02 mins elapsed.

2023-11-14 04:37:21.096422 : mat SRR13927735, Class = dgCMatrix
mat SRR13927735: nRows = 210721, nCols = 3813
mat SRR13927735: NonZeroEntries = 26038084, EntryRange = [ 1 , 4 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927735#TTATGTCTCCAGGTAT-1 SRR13927735#TATTGCTCATCAGAAA-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              1                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927735#TTCGATTGTAGGGTTG-1 SRR13927735#CATTCATTCGGATGTT-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              .                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927735#ACGTTAGGTCAACTGT-1
[1,]                              2
[2,]                              .
[3,]                              .
[4,]                              .
[5,]                              .


2023-11-14 04:37:21.276463 : Organizing SE PeakMatrix from ArrowFile : SRR13927735.arrow, 0.217 mins elapsed.

2023-11-14 04:37:21.582537 : se SRR13927735, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 3813 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(3813): SRR13927735#TTATGTCTCCAGGTAT-1
  SRR13927735#TATTGCTCATCAGAAA-1 ... SRR13927735#GATCATGGTGTTGTTG-1
  SRR13927735#GTCCATCCATTTCACT-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-14 04:37:21.59084 : Completed PeakMatrix : SRR13927735(3 of 4), 0.756 mins elapsed.
2023-11-14 04:37:21.591398 : Reading PeakMatrix : SRR13927736(4 of 4), 0.756 mins elapsed.

2023-11-14 04:37:21.610724 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927736 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927736.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3534
[1] "SRR13927736#GCTGTTCTCCTTGACC-1" "SRR13927736#AGTTACGAGTAGGTCG-1"
[3] "SRR13927736#AGTTTGGAGCATTCCA-1" "SRR13927736#GCACGGTTCGGCAATT-1"
[5] "SRR13927736#TCAGCTCTCTCTATTG-1" "SRR13927736#CCTTAATTCTACTTTG-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-14 04:37:22.380122 : featureDF SRR13927736, Class = DFrame
DataFrame with 210721 rows and 4 columns
       seqnames       idx     start       end
          <Rle> <integer> <integer> <integer>
1          chr1         1    804677    805177
2          chr1         2    811866    812366
3          chr1         3    817090    817590
4          chr1         4    818831    819331
5          chr1         5    826543    827043
...         ...       ...       ...       ...
210717     chrX      5777 155898525 155899025
210718     chrX      5778 155899804 155900304
210719     chrX      5779 155956608 155957108
210720     chrX      5780 155966827 155967327
210721     chrX      5781 155971625 155972125

2023-11-14 04:37:22.432963 : Getting PeakMatrix from ArrowFile : SRR13927736.arrow, 0.014 mins elapsed.

2023-11-14 04:37:32.974291 : mat SRR13927736, Class = dgCMatrix
mat SRR13927736: nRows = 210721, nCols = 3534
mat SRR13927736: NonZeroEntries = 25583439, EntryRange = [ 1 , 4 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927736#GCTGTTCTCCTTGACC-1 SRR13927736#AGTTACGAGTAGGTCG-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              .                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927736#AGTTTGGAGCATTCCA-1 SRR13927736#GCACGGTTCGGCAATT-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              1                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927736#TCAGCTCTCTCTATTG-1
[1,]                              .
[2,]                              .
[3,]                              .
[4,]                              .
[5,]                              .


2023-11-14 04:37:33.18388 : Organizing SE PeakMatrix from ArrowFile : SRR13927736.arrow, 0.193 mins elapsed.

2023-11-14 04:37:33.487366 : se SRR13927736, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 3534 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(3534): SRR13927736#GCTGTTCTCCTTGACC-1
  SRR13927736#AGTTACGAGTAGGTCG-1 ... SRR13927736#TAAACCGCAAGTCTGT-1
  SRR13927736#TGCATTTGTTCTACCC-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-14 04:37:33.494921 : Completed PeakMatrix : SRR13927736(4 of 4), 0.955 mins elapsed.
2023-11-14 04:37:33.496065 : Organizing colData, 0.955 mins elapsed.
2023-11-14 04:37:33.647049 : Organizing rowData, 0.957 mins elapsed.
2023-11-14 04:37:33.655823 : Organizing rowRanges, 0.957 mins elapsed.
2023-11-14 04:37:33.66738 : Organizing Assays (1 of 1), 0.957 mins elapsed.
2023-11-14 04:37:37.319184 : Constructing SummarizedExperiment, 1.018 mins elapsed.
2023-11-14 04:37:53.286671 : Finished Matrix Creation, 1.284 mins elapsed.

2023-11-14 06:13:52.95926 : getMatrixFromProject Input-Parameters, Class = list

getMatrixFromProject Input-Parameters$ArchRProj: length = 1

getMatrixFromProject Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromProject Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromProject Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromProject Input-Parameters$verbose: length = 1
[1] TRUE


getMatrixFromProject Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromProject Input-Parameters$threads: length = 1
[1] 1


getMatrixFromProject Input-Parameters$logFile: length = 1
[1] "x"


2023-11-14 06:13:53.965478 : Reading GeneIntegrationMatrix : SRR13927737(1 of 4), 0.017 mins elapsed.

2023-11-14 06:13:53.986277 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927737 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927737.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 4439
[1] "SRR13927737#GTCACGGGTCAACAGG-1" "SRR13927737#TGCCTGTTCTTACTCA-1"
[3] "SRR13927737#TTGCGGGTCAGGTCTA-1" "SRR13927737#GGTTGCGCAAAGCTTC-1"
[5] "SRR13927737#CTTAATCTCGCGCTGA-1" "SRR13927737#TGTACGATCGGAAAGC-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-14 06:13:54.309119 : featureDF SRR13927737, Class = DFrame
DataFrame with 23525 rows and 6 columns
      seqnames     start       end    strand        name       idx
         <Rle> <integer> <integer> <integer> <character> <integer>
1         chr1     69091     70008         1       OR4F5         1
2         chr1    817371    819837         1      FAM87B         2
3         chr1    825138    859446         1   LINC01128         3
4         chr1    827522    826206         2   LINC00115         4
5         chr1    876903    868071         2      FAM41C         5
...        ...       ...       ...       ...         ...       ...
23521     chrX 155264589 155258241         2      RAB39B       923
23522     chrX 155334657 155276211         2       CLIC2       924
23523     chrX 155380787 155381134         1      H2AFB1       925
23524     chrX 155382115 155383230         1        F8A2       926
23525     chrX 155669944 155490115         2       TMLHE       927

2023-11-14 06:13:54.356765 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927737.arrow, 0.006 mins elapsed.

2023-11-14 06:14:49.362983 : getMatrixFromProject Input-Parameters, Class = list

getMatrixFromProject Input-Parameters$ArchRProj: length = 1

getMatrixFromProject Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromProject Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromProject Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromProject Input-Parameters$verbose: length = 1
[1] TRUE


getMatrixFromProject Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromProject Input-Parameters$threads: length = 1
[1] 1


getMatrixFromProject Input-Parameters$logFile: length = 1
[1] "x"


2023-11-14 06:14:49.713162 : Reading GeneIntegrationMatrix : SRR13927737(1 of 4), 0.006 mins elapsed.

2023-11-14 06:14:49.728871 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927737 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927737.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 4439
[1] "SRR13927737#GTCACGGGTCAACAGG-1" "SRR13927737#TGCCTGTTCTTACTCA-1"
[3] "SRR13927737#TTGCGGGTCAGGTCTA-1" "SRR13927737#GGTTGCGCAAAGCTTC-1"
[5] "SRR13927737#CTTAATCTCGCGCTGA-1" "SRR13927737#TGTACGATCGGAAAGC-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-14 06:14:49.958341 : featureDF SRR13927737, Class = DFrame
DataFrame with 23525 rows and 6 columns
      seqnames     start       end    strand        name       idx
         <Rle> <integer> <integer> <integer> <character> <integer>
1         chr1     69091     70008         1       OR4F5         1
2         chr1    817371    819837         1      FAM87B         2
3         chr1    825138    859446         1   LINC01128         3
4         chr1    827522    826206         2   LINC00115         4
5         chr1    876903    868071         2      FAM41C         5
...        ...       ...       ...       ...         ...       ...
23521     chrX 155264589 155258241         2      RAB39B       923
23522     chrX 155334657 155276211         2       CLIC2       924
23523     chrX 155380787 155381134         1      H2AFB1       925
23524     chrX 155382115 155383230         1        F8A2       926
23525     chrX 155669944 155490115         2       TMLHE       927

2023-11-14 06:14:49.991564 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927737.arrow, 0.004 mins elapsed.

2023-11-14 06:15:23.264313 : mat SRR13927737, Class = dgCMatrix
mat SRR13927737: nRows = 23525, nCols = 4439
mat SRR13927737: NonZeroEntries = 50277774, EntryRange = [ 0.000148130054904683 , 376.345870564644 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927737#GTCACGGGTCAACAGG-1 SRR13927737#TGCCTGTTCTTACTCA-1
[1,]                      .                             .         
[2,]                      .                             .         
[3,]                      0.1057576                     0.09630696
[4,]                      .                             0.01806508
[5,]                      .                             .         
     SRR13927737#TTGCGGGTCAGGTCTA-1 SRR13927737#GGTTGCGCAAAGCTTC-1
[1,]                     .                              .         
[2,]                     .                              .         
[3,]                     0.09113499                     0.08396973
[4,]                     0.02096304                     0.01727603
[5,]                     .                              .         
     SRR13927737#CTTAATCTCGCGCTGA-1
[1,]                     .         
[2,]                     .         
[3,]                     0.08327649
[4,]                     0.01368917
[5,]                     .         


2023-11-14 06:15:23.738251 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927737.arrow, 0.567 mins elapsed.

2023-11-14 06:15:24.096991 : se SRR13927737, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 4439 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(4439): SRR13927737#GTCACGGGTCAACAGG-1
  SRR13927737#TGCCTGTTCTTACTCA-1 ... SRR13927737#TTGCTATAGATGTTCC-1
  SRR13927737#GCTCACTAGTTAGCGG-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-14 06:15:24.103542 : Completed GeneIntegrationMatrix : SRR13927737(1 of 4), 0.579 mins elapsed.
2023-11-14 06:15:24.104048 : Reading GeneIntegrationMatrix : SRR13927738(2 of 4), 0.579 mins elapsed.

2023-11-14 06:15:24.11641 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927738 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927738.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3736
[1] "SRR13927738#TAGCATGAGGCTAAAT-1" "SRR13927738#TGTTAGGTCGGTCAGC-1"
[3] "SRR13927738#CTCTCAGCATGTCCCT-1" "SRR13927738#ACTATTCGTCCAACCG-1"
[5] "SRR13927738#TGGTCCTTCTACTGCC-1" "SRR13927738#ATTGTCTAGACTAGGC-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-14 06:15:24.297146 : featureDF SRR13927738, Class = DFrame
DataFrame with 23525 rows and 6 columns
      seqnames     start       end    strand        name       idx
         <Rle> <integer> <integer> <integer> <character> <integer>
1         chr1     69091     70008         1       OR4F5         1
2         chr1    817371    819837         1      FAM87B         2
3         chr1    825138    859446         1   LINC01128         3
4         chr1    827522    826206         2   LINC00115         4
5         chr1    876903    868071         2      FAM41C         5
...        ...       ...       ...       ...         ...       ...
23521     chrX 155264589 155258241         2      RAB39B       923
23522     chrX 155334657 155276211         2       CLIC2       924
23523     chrX 155380787 155381134         1      H2AFB1       925
23524     chrX 155382115 155383230         1        F8A2       926
23525     chrX 155669944 155490115         2       TMLHE       927

2023-11-14 06:15:24.313689 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927738.arrow, 0.003 mins elapsed.

2023-11-14 06:15:42.679327 : mat SRR13927738, Class = dgCMatrix
mat SRR13927738: nRows = 23525, nCols = 3736
mat SRR13927738: NonZeroEntries = 41898897, EntryRange = [ 9.73975135001421e-05 , 758.763031303812 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927738#TAGCATGAGGCTAAAT-1 SRR13927738#TGTTAGGTCGGTCAGC-1
[1,]                    .                                .        
[2,]                    .                                .        
[3,]                    0.002346352                      0.0084575
[4,]                    0.005836413                      .        
[5,]                    .                                .        
     SRR13927738#CTCTCAGCATGTCCCT-1 SRR13927738#ACTATTCGTCCAACCG-1
[1,]                     .                             .          
[2,]                     .                             .          
[3,]                     .                             0.004081077
[4,]                     0.02001583                    .          
[5,]                     .                             .          
     SRR13927738#TGGTCCTTCTACTGCC-1
[1,]                    .          
[2,]                    .          
[3,]                    0.007640877
[4,]                    0.015048796
[5,]                    .          


2023-11-14 06:15:43.025897 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927738.arrow, 0.315 mins elapsed.

2023-11-14 06:15:44.02805 : se SRR13927738, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 3736 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(3736): SRR13927738#TAGCATGAGGCTAAAT-1
  SRR13927738#TGTTAGGTCGGTCAGC-1 ... SRR13927738#GCTCCTATCCGTGCGA-1
  SRR13927738#GGCACGTAGTCGATAA-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-14 06:15:44.03347 : Completed GeneIntegrationMatrix : SRR13927738(2 of 4), 0.911 mins elapsed.
2023-11-14 06:15:44.033971 : Reading GeneIntegrationMatrix : SRR13927735(3 of 4), 0.911 mins elapsed.

2023-11-14 06:15:44.046197 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927735 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927735.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3813
[1] "SRR13927735#TTATGTCTCCAGGTAT-1" "SRR13927735#TATTGCTCATCAGAAA-1"
[3] "SRR13927735#TTCGATTGTAGGGTTG-1" "SRR13927735#CATTCATTCGGATGTT-1"
[5] "SRR13927735#ACGTTAGGTCAACTGT-1" "SRR13927735#AAATGCCCAGCAATGG-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-14 06:15:44.232858 : featureDF SRR13927735, Class = DFrame
DataFrame with 23525 rows and 6 columns
      seqnames     start       end    strand        name       idx
         <Rle> <integer> <integer> <integer> <character> <integer>
1         chr1     69091     70008         1       OR4F5         1
2         chr1    817371    819837         1      FAM87B         2
3         chr1    825138    859446         1   LINC01128         3
4         chr1    827522    826206         2   LINC00115         4
5         chr1    876903    868071         2      FAM41C         5
...        ...       ...       ...       ...         ...       ...
23521     chrX 155264589 155258241         2      RAB39B       923
23522     chrX 155334657 155276211         2       CLIC2       924
23523     chrX 155380787 155381134         1      H2AFB1       925
23524     chrX 155382115 155383230         1        F8A2       926
23525     chrX 155669944 155490115         2       TMLHE       927

2023-11-14 06:15:44.247892 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927735.arrow, 0.003 mins elapsed.

2023-11-14 06:16:03.594877 : mat SRR13927735, Class = dgCMatrix
mat SRR13927735: nRows = 23525, nCols = 3813
mat SRR13927735: NonZeroEntries = 45535634, EntryRange = [ 7.20715653308254e-05 , 673.990433360078 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927735#TTATGTCTCCAGGTAT-1 SRR13927735#TATTGCTCATCAGAAA-1
[1,]                    .                              .          
[2,]                    .                              .          
[3,]                    0.070892732                    0.057009268
[4,]                    0.008800039                    0.002978408
[5,]                    .                              .          
     SRR13927735#TTCGATTGTAGGGTTG-1 SRR13927735#CATTCATTCGGATGTT-1
[1,]                     .                               .        
[2,]                     .                               .        
[3,]                     0.05523011                      0.0248301
[4,]                     .                               .        
[5,]                     .                               .        
     SRR13927735#ACGTTAGGTCAACTGT-1
[1,]                     .         
[2,]                     .         
[3,]                     0.05384571
[4,]                     0.01228688
[5,]                     .         


2023-11-14 06:16:03.957084 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927735.arrow, 0.332 mins elapsed.

2023-11-14 06:16:04.206586 : se SRR13927735, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 3813 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(3813): SRR13927735#TTATGTCTCCAGGTAT-1
  SRR13927735#TATTGCTCATCAGAAA-1 ... SRR13927735#GATCATGGTGTTGTTG-1
  SRR13927735#GTCCATCCATTTCACT-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-14 06:16:04.211783 : Completed GeneIntegrationMatrix : SRR13927735(3 of 4), 1.247 mins elapsed.
2023-11-14 06:16:04.212266 : Reading GeneIntegrationMatrix : SRR13927736(4 of 4), 1.248 mins elapsed.

2023-11-14 06:16:04.223965 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927736 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927736.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3534
[1] "SRR13927736#GCTGTTCTCCTTGACC-1" "SRR13927736#AGTTACGAGTAGGTCG-1"
[3] "SRR13927736#AGTTTGGAGCATTCCA-1" "SRR13927736#GCACGGTTCGGCAATT-1"
[5] "SRR13927736#TCAGCTCTCTCTATTG-1" "SRR13927736#CCTTAATTCTACTTTG-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-14 06:16:04.398181 : featureDF SRR13927736, Class = DFrame
DataFrame with 23525 rows and 6 columns
      seqnames     start       end    strand        name       idx
         <Rle> <integer> <integer> <integer> <character> <integer>
1         chr1     69091     70008         1       OR4F5         1
2         chr1    817371    819837         1      FAM87B         2
3         chr1    825138    859446         1   LINC01128         3
4         chr1    827522    826206         2   LINC00115         4
5         chr1    876903    868071         2      FAM41C         5
...        ...       ...       ...       ...         ...       ...
23521     chrX 155264589 155258241         2      RAB39B       923
23522     chrX 155334657 155276211         2       CLIC2       924
23523     chrX 155380787 155381134         1      H2AFB1       925
23524     chrX 155382115 155383230         1        F8A2       926
23525     chrX 155669944 155490115         2       TMLHE       927

2023-11-14 06:16:04.412576 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927736.arrow, 0.003 mins elapsed.

2023-11-14 06:16:20.64852 : mat SRR13927736, Class = dgCMatrix
mat SRR13927736: nRows = 23525, nCols = 3534
mat SRR13927736: NonZeroEntries = 39003616, EntryRange = [ 8.77627031078293e-05 , 770.387276519725 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927736#GCTGTTCTCCTTGACC-1 SRR13927736#AGTTACGAGTAGGTCG-1
[1,]                     .                             .          
[2,]                     .                             .          
[3,]                     0.07845381                    0.028344539
[4,]                     0.01395289                    0.008499189
[5,]                     .                             .          
     SRR13927736#AGTTTGGAGCATTCCA-1 SRR13927736#GCACGGTTCGGCAATT-1
[1,]                     .                              .         
[2,]                     .                              .         
[3,]                     0.06535130                     0.06576222
[4,]                     0.01810177                     0.01315728
[5,]                     .                              .         
     SRR13927736#TCAGCTCTCTCTATTG-1
[1,]                    .          
[2,]                    .          
[3,]                    0.016229526
[4,]                    0.003862592
[5,]                    .          


2023-11-14 06:16:20.978674 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927736.arrow, 0.279 mins elapsed.

2023-11-14 06:16:21.23631 : se SRR13927736, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 3534 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(3534): SRR13927736#GCTGTTCTCCTTGACC-1
  SRR13927736#AGTTACGAGTAGGTCG-1 ... SRR13927736#TAAACCGCAAGTCTGT-1
  SRR13927736#TGCATTTGTTCTACCC-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-14 06:16:21.241697 : Completed GeneIntegrationMatrix : SRR13927736(4 of 4), 1.531 mins elapsed.
2023-11-14 06:16:21.242779 : Organizing colData, 1.531 mins elapsed.
2023-11-14 06:16:21.392269 : Organizing rowData, 1.534 mins elapsed.
2023-11-14 06:16:21.395642 : Organizing rowRanges, 1.534 mins elapsed.
2023-11-14 06:16:21.40152 : Organizing Assays (1 of 1), 1.534 mins elapsed.
2023-11-14 06:16:27.818254 : Constructing SummarizedExperiment, 1.641 mins elapsed.
2023-11-14 06:16:28.549128 : Finished Matrix Creation, 1.653 mins elapsed.

2023-11-14 06:16:28.551666 : getMatrixFromProject Input-Parameters, Class = list

getMatrixFromProject Input-Parameters$ArchRProj: length = 1

getMatrixFromProject Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromProject Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromProject Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromProject Input-Parameters$verbose: length = 1
[1] TRUE


getMatrixFromProject Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromProject Input-Parameters$threads: length = 1
[1] 1


getMatrixFromProject Input-Parameters$logFile: length = 1
[1] "x"


2023-11-14 06:16:29.214335 : Reading PeakMatrix : SRR13927737(1 of 4), 0.011 mins elapsed.

2023-11-14 06:16:29.228618 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927737 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927737.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 4439
[1] "SRR13927737#GTCACGGGTCAACAGG-1" "SRR13927737#TGCCTGTTCTTACTCA-1"
[3] "SRR13927737#TTGCGGGTCAGGTCTA-1" "SRR13927737#GGTTGCGCAAAGCTTC-1"
[5] "SRR13927737#CTTAATCTCGCGCTGA-1" "SRR13927737#TGTACGATCGGAAAGC-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-14 06:16:29.509403 : featureDF SRR13927737, Class = DFrame
DataFrame with 210721 rows and 4 columns
       seqnames       idx     start       end
          <Rle> <integer> <integer> <integer>
1          chr1         1    804677    805177
2          chr1         2    811866    812366
3          chr1         3    817090    817590
4          chr1         4    818831    819331
5          chr1         5    826543    827043
...         ...       ...       ...       ...
210717     chrX      5777 155898525 155899025
210718     chrX      5778 155899804 155900304
210719     chrX      5779 155956608 155957108
210720     chrX      5780 155966827 155967327
210721     chrX      5781 155971625 155972125

2023-11-14 06:16:29.560348 : Getting PeakMatrix from ArrowFile : SRR13927737.arrow, 0.006 mins elapsed.

2023-11-14 06:16:42.052511 : mat SRR13927737, Class = dgCMatrix
mat SRR13927737: nRows = 210721, nCols = 4439
mat SRR13927737: NonZeroEntries = 29963541, EntryRange = [ 1 , 4 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927737#GTCACGGGTCAACAGG-1 SRR13927737#TGCCTGTTCTTACTCA-1
[1,]                              .                              .
[2,]                              .                              2
[3,]                              .                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927737#TTGCGGGTCAGGTCTA-1 SRR13927737#GGTTGCGCAAAGCTTC-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              2                              .
[4,]                              .                              .
[5,]                              4                              .
     SRR13927737#CTTAATCTCGCGCTGA-1
[1,]                              .
[2,]                              .
[3,]                              .
[4,]                              .
[5,]                              .


2023-11-14 06:16:42.239958 : Organizing SE PeakMatrix from ArrowFile : SRR13927737.arrow, 0.217 mins elapsed.

2023-11-14 06:16:42.558903 : se SRR13927737, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 4439 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(4439): SRR13927737#GTCACGGGTCAACAGG-1
  SRR13927737#TGCCTGTTCTTACTCA-1 ... SRR13927737#TTGCTATAGATGTTCC-1
  SRR13927737#GCTCACTAGTTAGCGG-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-14 06:16:42.569254 : Completed PeakMatrix : SRR13927737(1 of 4), 0.234 mins elapsed.
2023-11-14 06:16:42.569818 : Reading PeakMatrix : SRR13927738(2 of 4), 0.234 mins elapsed.

2023-11-14 06:16:42.582556 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927738 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927738.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3736
[1] "SRR13927738#TAGCATGAGGCTAAAT-1" "SRR13927738#TGTTAGGTCGGTCAGC-1"
[3] "SRR13927738#CTCTCAGCATGTCCCT-1" "SRR13927738#ACTATTCGTCCAACCG-1"
[5] "SRR13927738#TGGTCCTTCTACTGCC-1" "SRR13927738#ATTGTCTAGACTAGGC-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-14 06:16:42.842114 : featureDF SRR13927738, Class = DFrame
DataFrame with 210721 rows and 4 columns
       seqnames       idx     start       end
          <Rle> <integer> <integer> <integer>
1          chr1         1    804677    805177
2          chr1         2    811866    812366
3          chr1         3    817090    817590
4          chr1         4    818831    819331
5          chr1         5    826543    827043
...         ...       ...       ...       ...
210717     chrX      5777 155898525 155899025
210718     chrX      5778 155899804 155900304
210719     chrX      5779 155956608 155957108
210720     chrX      5780 155966827 155967327
210721     chrX      5781 155971625 155972125

2023-11-14 06:16:42.883088 : Getting PeakMatrix from ArrowFile : SRR13927738.arrow, 0.005 mins elapsed.

2023-11-14 06:16:54.841495 : mat SRR13927738, Class = dgCMatrix
mat SRR13927738: nRows = 210721, nCols = 3736
mat SRR13927738: NonZeroEntries = 30165921, EntryRange = [ 1 , 4 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927738#TAGCATGAGGCTAAAT-1 SRR13927738#TGTTAGGTCGGTCAGC-1
[1,]                              2                              .
[2,]                              .                              .
[3,]                              .                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927738#CTCTCAGCATGTCCCT-1 SRR13927738#ACTATTCGTCCAACCG-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              .                              .
[4,]                              .                              .
[5,]                              2                              .
     SRR13927738#TGGTCCTTCTACTGCC-1
[1,]                              .
[2,]                              .
[3,]                              2
[4,]                              .
[5,]                              .


2023-11-14 06:16:55.096251 : Organizing SE PeakMatrix from ArrowFile : SRR13927738.arrow, 0.209 mins elapsed.

2023-11-14 06:16:55.381942 : se SRR13927738, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 3736 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(3736): SRR13927738#TAGCATGAGGCTAAAT-1
  SRR13927738#TGTTAGGTCGGTCAGC-1 ... SRR13927738#GCTCCTATCCGTGCGA-1
  SRR13927738#GGCACGTAGTCGATAA-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-14 06:16:55.388936 : Completed PeakMatrix : SRR13927738(2 of 4), 0.447 mins elapsed.
2023-11-14 06:16:55.389415 : Reading PeakMatrix : SRR13927735(3 of 4), 0.447 mins elapsed.

2023-11-14 06:16:55.402107 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927735 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927735.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3813
[1] "SRR13927735#TTATGTCTCCAGGTAT-1" "SRR13927735#TATTGCTCATCAGAAA-1"
[3] "SRR13927735#TTCGATTGTAGGGTTG-1" "SRR13927735#CATTCATTCGGATGTT-1"
[5] "SRR13927735#ACGTTAGGTCAACTGT-1" "SRR13927735#AAATGCCCAGCAATGG-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-14 06:16:55.720554 : featureDF SRR13927735, Class = DFrame
DataFrame with 210721 rows and 4 columns
       seqnames       idx     start       end
          <Rle> <integer> <integer> <integer>
1          chr1         1    804677    805177
2          chr1         2    811866    812366
3          chr1         3    817090    817590
4          chr1         4    818831    819331
5          chr1         5    826543    827043
...         ...       ...       ...       ...
210717     chrX      5777 155898525 155899025
210718     chrX      5778 155899804 155900304
210719     chrX      5779 155956608 155957108
210720     chrX      5780 155966827 155967327
210721     chrX      5781 155971625 155972125

2023-11-14 06:16:55.76092 : Getting PeakMatrix from ArrowFile : SRR13927735.arrow, 0.006 mins elapsed.

2023-11-14 06:17:05.899371 : mat SRR13927735, Class = dgCMatrix
mat SRR13927735: nRows = 210721, nCols = 3813
mat SRR13927735: NonZeroEntries = 26038084, EntryRange = [ 1 , 4 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927735#TTATGTCTCCAGGTAT-1 SRR13927735#TATTGCTCATCAGAAA-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              1                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927735#TTCGATTGTAGGGTTG-1 SRR13927735#CATTCATTCGGATGTT-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              .                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927735#ACGTTAGGTCAACTGT-1
[1,]                              2
[2,]                              .
[3,]                              .
[4,]                              .
[5,]                              .


2023-11-14 06:17:06.107453 : Organizing SE PeakMatrix from ArrowFile : SRR13927735.arrow, 0.178 mins elapsed.

2023-11-14 06:17:06.384486 : se SRR13927735, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 3813 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(3813): SRR13927735#TTATGTCTCCAGGTAT-1
  SRR13927735#TATTGCTCATCAGAAA-1 ... SRR13927735#GATCATGGTGTTGTTG-1
  SRR13927735#GTCCATCCATTTCACT-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-14 06:17:06.391737 : Completed PeakMatrix : SRR13927735(3 of 4), 0.631 mins elapsed.
2023-11-14 06:17:06.392242 : Reading PeakMatrix : SRR13927736(4 of 4), 0.631 mins elapsed.

2023-11-14 06:17:06.404649 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927736 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927736.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3534
[1] "SRR13927736#GCTGTTCTCCTTGACC-1" "SRR13927736#AGTTACGAGTAGGTCG-1"
[3] "SRR13927736#AGTTTGGAGCATTCCA-1" "SRR13927736#GCACGGTTCGGCAATT-1"
[5] "SRR13927736#TCAGCTCTCTCTATTG-1" "SRR13927736#CCTTAATTCTACTTTG-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-14 06:17:06.69032 : featureDF SRR13927736, Class = DFrame
DataFrame with 210721 rows and 4 columns
       seqnames       idx     start       end
          <Rle> <integer> <integer> <integer>
1          chr1         1    804677    805177
2          chr1         2    811866    812366
3          chr1         3    817090    817590
4          chr1         4    818831    819331
5          chr1         5    826543    827043
...         ...       ...       ...       ...
210717     chrX      5777 155898525 155899025
210718     chrX      5778 155899804 155900304
210719     chrX      5779 155956608 155957108
210720     chrX      5780 155966827 155967327
210721     chrX      5781 155971625 155972125

2023-11-14 06:17:06.729396 : Getting PeakMatrix from ArrowFile : SRR13927736.arrow, 0.005 mins elapsed.

2023-11-14 06:17:17.067418 : mat SRR13927736, Class = dgCMatrix
mat SRR13927736: nRows = 210721, nCols = 3534
mat SRR13927736: NonZeroEntries = 25583439, EntryRange = [ 1 , 4 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927736#GCTGTTCTCCTTGACC-1 SRR13927736#AGTTACGAGTAGGTCG-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              .                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927736#AGTTTGGAGCATTCCA-1 SRR13927736#GCACGGTTCGGCAATT-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              1                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927736#TCAGCTCTCTCTATTG-1
[1,]                              .
[2,]                              .
[3,]                              .
[4,]                              .
[5,]                              .


2023-11-14 06:17:17.274543 : Organizing SE PeakMatrix from ArrowFile : SRR13927736.arrow, 0.181 mins elapsed.

2023-11-14 06:17:17.565583 : se SRR13927736, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 3534 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(3534): SRR13927736#GCTGTTCTCCTTGACC-1
  SRR13927736#AGTTACGAGTAGGTCG-1 ... SRR13927736#TAAACCGCAAGTCTGT-1
  SRR13927736#TGCATTTGTTCTACCC-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-14 06:17:17.57298 : Completed PeakMatrix : SRR13927736(4 of 4), 0.817 mins elapsed.
2023-11-14 06:17:17.574028 : Organizing colData, 0.817 mins elapsed.
2023-11-14 06:17:17.708702 : Organizing rowData, 0.819 mins elapsed.
2023-11-14 06:17:17.717264 : Organizing rowRanges, 0.819 mins elapsed.
2023-11-14 06:17:17.728079 : Organizing Assays (1 of 1), 0.82 mins elapsed.
2023-11-14 06:17:20.955517 : Constructing SummarizedExperiment, 0.873 mins elapsed.
2023-11-14 06:17:34.28632 : Finished Matrix Creation, 1.096 mins elapsed.

2023-11-14 06:22:34.817651 : getMatrixFromProject Input-Parameters, Class = list

getMatrixFromProject Input-Parameters$ArchRProj: length = 1

getMatrixFromProject Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromProject Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromProject Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromProject Input-Parameters$verbose: length = 1
[1] TRUE


getMatrixFromProject Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromProject Input-Parameters$threads: length = 1
[1] 1


getMatrixFromProject Input-Parameters$logFile: length = 1
[1] "x"


2023-11-14 06:22:35.261912 : Reading GeneIntegrationMatrix : SRR13927737(1 of 4), 0.007 mins elapsed.

2023-11-14 06:22:35.276826 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927737 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927737.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 4439
[1] "SRR13927737#GTCACGGGTCAACAGG-1" "SRR13927737#TGCCTGTTCTTACTCA-1"
[3] "SRR13927737#TTGCGGGTCAGGTCTA-1" "SRR13927737#GGTTGCGCAAAGCTTC-1"
[5] "SRR13927737#CTTAATCTCGCGCTGA-1" "SRR13927737#TGTACGATCGGAAAGC-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-14 06:22:35.525662 : featureDF SRR13927737, Class = DFrame
DataFrame with 23525 rows and 6 columns
      seqnames     start       end    strand        name       idx
         <Rle> <integer> <integer> <integer> <character> <integer>
1         chr1     69091     70008         1       OR4F5         1
2         chr1    817371    819837         1      FAM87B         2
3         chr1    825138    859446         1   LINC01128         3
4         chr1    827522    826206         2   LINC00115         4
5         chr1    876903    868071         2      FAM41C         5
...        ...       ...       ...       ...         ...       ...
23521     chrX 155264589 155258241         2      RAB39B       923
23522     chrX 155334657 155276211         2       CLIC2       924
23523     chrX 155380787 155381134         1      H2AFB1       925
23524     chrX 155382115 155383230         1        F8A2       926
23525     chrX 155669944 155490115         2       TMLHE       927

2023-11-14 06:22:35.579642 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927737.arrow, 0.005 mins elapsed.

2023-11-14 06:23:09.145415 : mat SRR13927737, Class = dgCMatrix
mat SRR13927737: nRows = 23525, nCols = 4439
mat SRR13927737: NonZeroEntries = 50277774, EntryRange = [ 0.000148130054904683 , 376.345870564644 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927737#GTCACGGGTCAACAGG-1 SRR13927737#TGCCTGTTCTTACTCA-1
[1,]                      .                             .         
[2,]                      .                             .         
[3,]                      0.1057576                     0.09630696
[4,]                      .                             0.01806508
[5,]                      .                             .         
     SRR13927737#TTGCGGGTCAGGTCTA-1 SRR13927737#GGTTGCGCAAAGCTTC-1
[1,]                     .                              .         
[2,]                     .                              .         
[3,]                     0.09113499                     0.08396973
[4,]                     0.02096304                     0.01727603
[5,]                     .                              .         
     SRR13927737#CTTAATCTCGCGCTGA-1
[1,]                     .         
[2,]                     .         
[3,]                     0.08327649
[4,]                     0.01368917
[5,]                     .         


2023-11-14 06:23:09.618093 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927737.arrow, 0.572 mins elapsed.

2023-11-14 06:23:09.983389 : se SRR13927737, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 4439 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(4439): SRR13927737#GTCACGGGTCAACAGG-1
  SRR13927737#TGCCTGTTCTTACTCA-1 ... SRR13927737#TTGCTATAGATGTTCC-1
  SRR13927737#GCTCACTAGTTAGCGG-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-14 06:23:09.989983 : Completed GeneIntegrationMatrix : SRR13927737(1 of 4), 0.586 mins elapsed.
2023-11-14 06:23:09.99056 : Reading GeneIntegrationMatrix : SRR13927738(2 of 4), 0.586 mins elapsed.

2023-11-14 06:23:10.003118 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927738 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927738.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3736
[1] "SRR13927738#TAGCATGAGGCTAAAT-1" "SRR13927738#TGTTAGGTCGGTCAGC-1"
[3] "SRR13927738#CTCTCAGCATGTCCCT-1" "SRR13927738#ACTATTCGTCCAACCG-1"
[5] "SRR13927738#TGGTCCTTCTACTGCC-1" "SRR13927738#ATTGTCTAGACTAGGC-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-14 06:23:10.19484 : featureDF SRR13927738, Class = DFrame
DataFrame with 23525 rows and 6 columns
      seqnames     start       end    strand        name       idx
         <Rle> <integer> <integer> <integer> <character> <integer>
1         chr1     69091     70008         1       OR4F5         1
2         chr1    817371    819837         1      FAM87B         2
3         chr1    825138    859446         1   LINC01128         3
4         chr1    827522    826206         2   LINC00115         4
5         chr1    876903    868071         2      FAM41C         5
...        ...       ...       ...       ...         ...       ...
23521     chrX 155264589 155258241         2      RAB39B       923
23522     chrX 155334657 155276211         2       CLIC2       924
23523     chrX 155380787 155381134         1      H2AFB1       925
23524     chrX 155382115 155383230         1        F8A2       926
23525     chrX 155669944 155490115         2       TMLHE       927

2023-11-14 06:23:10.211825 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927738.arrow, 0.003 mins elapsed.

2023-11-14 06:23:27.558019 : mat SRR13927738, Class = dgCMatrix
mat SRR13927738: nRows = 23525, nCols = 3736
mat SRR13927738: NonZeroEntries = 41898897, EntryRange = [ 9.73975135001421e-05 , 758.763031303812 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927738#TAGCATGAGGCTAAAT-1 SRR13927738#TGTTAGGTCGGTCAGC-1
[1,]                    .                                .        
[2,]                    .                                .        
[3,]                    0.002346352                      0.0084575
[4,]                    0.005836413                      .        
[5,]                    .                                .        
     SRR13927738#CTCTCAGCATGTCCCT-1 SRR13927738#ACTATTCGTCCAACCG-1
[1,]                     .                             .          
[2,]                     .                             .          
[3,]                     .                             0.004081077
[4,]                     0.02001583                    .          
[5,]                     .                             .          
     SRR13927738#TGGTCCTTCTACTGCC-1
[1,]                    .          
[2,]                    .          
[3,]                    0.007640877
[4,]                    0.015048796
[5,]                    .          


2023-11-14 06:23:27.895305 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927738.arrow, 0.298 mins elapsed.

2023-11-14 06:23:28.907256 : se SRR13927738, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 3736 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(3736): SRR13927738#TAGCATGAGGCTAAAT-1
  SRR13927738#TGTTAGGTCGGTCAGC-1 ... SRR13927738#GCTCCTATCCGTGCGA-1
  SRR13927738#GGCACGTAGTCGATAA-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-14 06:23:28.912612 : Completed GeneIntegrationMatrix : SRR13927738(2 of 4), 0.902 mins elapsed.
2023-11-14 06:23:28.913103 : Reading GeneIntegrationMatrix : SRR13927735(3 of 4), 0.902 mins elapsed.

2023-11-14 06:23:28.925262 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927735 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927735.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3813
[1] "SRR13927735#TTATGTCTCCAGGTAT-1" "SRR13927735#TATTGCTCATCAGAAA-1"
[3] "SRR13927735#TTCGATTGTAGGGTTG-1" "SRR13927735#CATTCATTCGGATGTT-1"
[5] "SRR13927735#ACGTTAGGTCAACTGT-1" "SRR13927735#AAATGCCCAGCAATGG-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-14 06:23:29.145978 : featureDF SRR13927735, Class = DFrame
DataFrame with 23525 rows and 6 columns
      seqnames     start       end    strand        name       idx
         <Rle> <integer> <integer> <integer> <character> <integer>
1         chr1     69091     70008         1       OR4F5         1
2         chr1    817371    819837         1      FAM87B         2
3         chr1    825138    859446         1   LINC01128         3
4         chr1    827522    826206         2   LINC00115         4
5         chr1    876903    868071         2      FAM41C         5
...        ...       ...       ...       ...         ...       ...
23521     chrX 155264589 155258241         2      RAB39B       923
23522     chrX 155334657 155276211         2       CLIC2       924
23523     chrX 155380787 155381134         1      H2AFB1       925
23524     chrX 155382115 155383230         1        F8A2       926
23525     chrX 155669944 155490115         2       TMLHE       927

2023-11-14 06:23:29.161376 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927735.arrow, 0.004 mins elapsed.

2023-11-14 06:23:46.738884 : mat SRR13927735, Class = dgCMatrix
mat SRR13927735: nRows = 23525, nCols = 3813
mat SRR13927735: NonZeroEntries = 45535634, EntryRange = [ 7.20715653308254e-05 , 673.990433360078 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927735#TTATGTCTCCAGGTAT-1 SRR13927735#TATTGCTCATCAGAAA-1
[1,]                    .                              .          
[2,]                    .                              .          
[3,]                    0.070892732                    0.057009268
[4,]                    0.008800039                    0.002978408
[5,]                    .                              .          
     SRR13927735#TTCGATTGTAGGGTTG-1 SRR13927735#CATTCATTCGGATGTT-1
[1,]                     .                               .        
[2,]                     .                               .        
[3,]                     0.05523011                      0.0248301
[4,]                     .                               .        
[5,]                     .                               .        
     SRR13927735#ACGTTAGGTCAACTGT-1
[1,]                     .         
[2,]                     .         
[3,]                     0.05384571
[4,]                     0.01228688
[5,]                     .         


2023-11-14 06:23:47.107289 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927735.arrow, 0.303 mins elapsed.

2023-11-14 06:23:47.359538 : se SRR13927735, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 3813 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(3813): SRR13927735#TTATGTCTCCAGGTAT-1
  SRR13927735#TATTGCTCATCAGAAA-1 ... SRR13927735#GATCATGGTGTTGTTG-1
  SRR13927735#GTCCATCCATTTCACT-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-14 06:23:47.364258 : Completed GeneIntegrationMatrix : SRR13927735(3 of 4), 1.209 mins elapsed.
2023-11-14 06:23:47.364642 : Reading GeneIntegrationMatrix : SRR13927736(4 of 4), 1.209 mins elapsed.

2023-11-14 06:23:47.376507 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927736 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927736.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3534
[1] "SRR13927736#GCTGTTCTCCTTGACC-1" "SRR13927736#AGTTACGAGTAGGTCG-1"
[3] "SRR13927736#AGTTTGGAGCATTCCA-1" "SRR13927736#GCACGGTTCGGCAATT-1"
[5] "SRR13927736#TCAGCTCTCTCTATTG-1" "SRR13927736#CCTTAATTCTACTTTG-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-14 06:23:47.543532 : featureDF SRR13927736, Class = DFrame
DataFrame with 23525 rows and 6 columns
      seqnames     start       end    strand        name       idx
         <Rle> <integer> <integer> <integer> <character> <integer>
1         chr1     69091     70008         1       OR4F5         1
2         chr1    817371    819837         1      FAM87B         2
3         chr1    825138    859446         1   LINC01128         3
4         chr1    827522    826206         2   LINC00115         4
5         chr1    876903    868071         2      FAM41C         5
...        ...       ...       ...       ...         ...       ...
23521     chrX 155264589 155258241         2      RAB39B       923
23522     chrX 155334657 155276211         2       CLIC2       924
23523     chrX 155380787 155381134         1      H2AFB1       925
23524     chrX 155382115 155383230         1        F8A2       926
23525     chrX 155669944 155490115         2       TMLHE       927

2023-11-14 06:23:47.557872 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927736.arrow, 0.003 mins elapsed.

2023-11-14 06:24:02.747429 : mat SRR13927736, Class = dgCMatrix
mat SRR13927736: nRows = 23525, nCols = 3534
mat SRR13927736: NonZeroEntries = 39003616, EntryRange = [ 8.77627031078293e-05 , 770.387276519725 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927736#GCTGTTCTCCTTGACC-1 SRR13927736#AGTTACGAGTAGGTCG-1
[1,]                     .                             .          
[2,]                     .                             .          
[3,]                     0.07845381                    0.028344539
[4,]                     0.01395289                    0.008499189
[5,]                     .                             .          
     SRR13927736#AGTTTGGAGCATTCCA-1 SRR13927736#GCACGGTTCGGCAATT-1
[1,]                     .                              .         
[2,]                     .                              .         
[3,]                     0.06535130                     0.06576222
[4,]                     0.01810177                     0.01315728
[5,]                     .                              .         
     SRR13927736#TCAGCTCTCTCTATTG-1
[1,]                    .          
[2,]                    .          
[3,]                    0.016229526
[4,]                    0.003862592
[5,]                    .          


2023-11-14 06:24:03.071558 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927736.arrow, 0.262 mins elapsed.

2023-11-14 06:24:03.328911 : se SRR13927736, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 3534 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(3534): SRR13927736#GCTGTTCTCCTTGACC-1
  SRR13927736#AGTTACGAGTAGGTCG-1 ... SRR13927736#TAAACCGCAAGTCTGT-1
  SRR13927736#TGCATTTGTTCTACCC-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-14 06:24:03.333871 : Completed GeneIntegrationMatrix : SRR13927736(4 of 4), 1.475 mins elapsed.
2023-11-14 06:24:03.334839 : Organizing colData, 1.475 mins elapsed.
2023-11-14 06:24:03.476155 : Organizing rowData, 1.478 mins elapsed.
2023-11-14 06:24:03.479306 : Organizing rowRanges, 1.478 mins elapsed.
2023-11-14 06:24:03.484935 : Organizing Assays (1 of 1), 1.478 mins elapsed.
2023-11-14 06:24:09.795546 : Constructing SummarizedExperiment, 1.583 mins elapsed.
2023-11-14 06:24:10.517282 : Finished Matrix Creation, 1.595 mins elapsed.

2023-11-14 06:24:10.519619 : getMatrixFromProject Input-Parameters, Class = list

getMatrixFromProject Input-Parameters$ArchRProj: length = 1

getMatrixFromProject Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromProject Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromProject Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromProject Input-Parameters$verbose: length = 1
[1] TRUE


getMatrixFromProject Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromProject Input-Parameters$threads: length = 1
[1] 1


getMatrixFromProject Input-Parameters$logFile: length = 1
[1] "x"


2023-11-14 06:24:11.047502 : Reading PeakMatrix : SRR13927737(1 of 4), 0.009 mins elapsed.

2023-11-14 06:24:11.061463 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927737 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927737.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 4439
[1] "SRR13927737#GTCACGGGTCAACAGG-1" "SRR13927737#TGCCTGTTCTTACTCA-1"
[3] "SRR13927737#TTGCGGGTCAGGTCTA-1" "SRR13927737#GGTTGCGCAAAGCTTC-1"
[5] "SRR13927737#CTTAATCTCGCGCTGA-1" "SRR13927737#TGTACGATCGGAAAGC-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-14 06:24:11.340264 : featureDF SRR13927737, Class = DFrame
DataFrame with 210721 rows and 4 columns
       seqnames       idx     start       end
          <Rle> <integer> <integer> <integer>
1          chr1         1    804677    805177
2          chr1         2    811866    812366
3          chr1         3    817090    817590
4          chr1         4    818831    819331
5          chr1         5    826543    827043
...         ...       ...       ...       ...
210717     chrX      5777 155898525 155899025
210718     chrX      5778 155899804 155900304
210719     chrX      5779 155956608 155957108
210720     chrX      5780 155966827 155967327
210721     chrX      5781 155971625 155972125

2023-11-14 06:24:11.3929 : Getting PeakMatrix from ArrowFile : SRR13927737.arrow, 0.006 mins elapsed.

2023-11-14 06:24:22.997259 : mat SRR13927737, Class = dgCMatrix
mat SRR13927737: nRows = 210721, nCols = 4439
mat SRR13927737: NonZeroEntries = 29963541, EntryRange = [ 1 , 4 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927737#GTCACGGGTCAACAGG-1 SRR13927737#TGCCTGTTCTTACTCA-1
[1,]                              .                              .
[2,]                              .                              2
[3,]                              .                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927737#TTGCGGGTCAGGTCTA-1 SRR13927737#GGTTGCGCAAAGCTTC-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              2                              .
[4,]                              .                              .
[5,]                              4                              .
     SRR13927737#CTTAATCTCGCGCTGA-1
[1,]                              .
[2,]                              .
[3,]                              .
[4,]                              .
[5,]                              .


2023-11-14 06:24:23.185595 : Organizing SE PeakMatrix from ArrowFile : SRR13927737.arrow, 0.202 mins elapsed.

2023-11-14 06:24:23.483414 : se SRR13927737, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 4439 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(4439): SRR13927737#GTCACGGGTCAACAGG-1
  SRR13927737#TGCCTGTTCTTACTCA-1 ... SRR13927737#TTGCTATAGATGTTCC-1
  SRR13927737#GCTCACTAGTTAGCGG-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-14 06:24:23.49325 : Completed PeakMatrix : SRR13927737(1 of 4), 0.216 mins elapsed.
2023-11-14 06:24:23.493677 : Reading PeakMatrix : SRR13927738(2 of 4), 0.216 mins elapsed.

2023-11-14 06:24:23.506397 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927738 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927738.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3736
[1] "SRR13927738#TAGCATGAGGCTAAAT-1" "SRR13927738#TGTTAGGTCGGTCAGC-1"
[3] "SRR13927738#CTCTCAGCATGTCCCT-1" "SRR13927738#ACTATTCGTCCAACCG-1"
[5] "SRR13927738#TGGTCCTTCTACTGCC-1" "SRR13927738#ATTGTCTAGACTAGGC-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-14 06:24:23.765626 : featureDF SRR13927738, Class = DFrame
DataFrame with 210721 rows and 4 columns
       seqnames       idx     start       end
          <Rle> <integer> <integer> <integer>
1          chr1         1    804677    805177
2          chr1         2    811866    812366
3          chr1         3    817090    817590
4          chr1         4    818831    819331
5          chr1         5    826543    827043
...         ...       ...       ...       ...
210717     chrX      5777 155898525 155899025
210718     chrX      5778 155899804 155900304
210719     chrX      5779 155956608 155957108
210720     chrX      5780 155966827 155967327
210721     chrX      5781 155971625 155972125

2023-11-14 06:24:23.810143 : Getting PeakMatrix from ArrowFile : SRR13927738.arrow, 0.005 mins elapsed.

2023-11-14 06:24:35.984252 : mat SRR13927738, Class = dgCMatrix
mat SRR13927738: nRows = 210721, nCols = 3736
mat SRR13927738: NonZeroEntries = 30165921, EntryRange = [ 1 , 4 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927738#TAGCATGAGGCTAAAT-1 SRR13927738#TGTTAGGTCGGTCAGC-1
[1,]                              2                              .
[2,]                              .                              .
[3,]                              .                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927738#CTCTCAGCATGTCCCT-1 SRR13927738#ACTATTCGTCCAACCG-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              .                              .
[4,]                              .                              .
[5,]                              2                              .
     SRR13927738#TGGTCCTTCTACTGCC-1
[1,]                              .
[2,]                              .
[3,]                              2
[4,]                              .
[5,]                              .


2023-11-14 06:24:36.219008 : Organizing SE PeakMatrix from ArrowFile : SRR13927738.arrow, 0.212 mins elapsed.

2023-11-14 06:24:36.498009 : se SRR13927738, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 3736 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(3736): SRR13927738#TAGCATGAGGCTAAAT-1
  SRR13927738#TGTTAGGTCGGTCAGC-1 ... SRR13927738#GCTCCTATCCGTGCGA-1
  SRR13927738#GGCACGTAGTCGATAA-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-14 06:24:36.504769 : Completed PeakMatrix : SRR13927738(2 of 4), 0.433 mins elapsed.
2023-11-14 06:24:36.505191 : Reading PeakMatrix : SRR13927735(3 of 4), 0.433 mins elapsed.

2023-11-14 06:24:36.517648 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927735 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927735.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3813
[1] "SRR13927735#TTATGTCTCCAGGTAT-1" "SRR13927735#TATTGCTCATCAGAAA-1"
[3] "SRR13927735#TTCGATTGTAGGGTTG-1" "SRR13927735#CATTCATTCGGATGTT-1"
[5] "SRR13927735#ACGTTAGGTCAACTGT-1" "SRR13927735#AAATGCCCAGCAATGG-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-14 06:24:36.775162 : featureDF SRR13927735, Class = DFrame
DataFrame with 210721 rows and 4 columns
       seqnames       idx     start       end
          <Rle> <integer> <integer> <integer>
1          chr1         1    804677    805177
2          chr1         2    811866    812366
3          chr1         3    817090    817590
4          chr1         4    818831    819331
5          chr1         5    826543    827043
...         ...       ...       ...       ...
210717     chrX      5777 155898525 155899025
210718     chrX      5778 155899804 155900304
210719     chrX      5779 155956608 155957108
210720     chrX      5780 155966827 155967327
210721     chrX      5781 155971625 155972125

2023-11-14 06:24:36.816079 : Getting PeakMatrix from ArrowFile : SRR13927735.arrow, 0.005 mins elapsed.

2023-11-14 06:24:46.446446 : mat SRR13927735, Class = dgCMatrix
mat SRR13927735: nRows = 210721, nCols = 3813
mat SRR13927735: NonZeroEntries = 26038084, EntryRange = [ 1 , 4 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927735#TTATGTCTCCAGGTAT-1 SRR13927735#TATTGCTCATCAGAAA-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              1                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927735#TTCGATTGTAGGGTTG-1 SRR13927735#CATTCATTCGGATGTT-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              .                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927735#ACGTTAGGTCAACTGT-1
[1,]                              2
[2,]                              .
[3,]                              .
[4,]                              .
[5,]                              .


2023-11-14 06:24:46.654355 : Organizing SE PeakMatrix from ArrowFile : SRR13927735.arrow, 0.169 mins elapsed.

2023-11-14 06:24:46.98007 : se SRR13927735, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 3813 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(3813): SRR13927735#TTATGTCTCCAGGTAT-1
  SRR13927735#TATTGCTCATCAGAAA-1 ... SRR13927735#GATCATGGTGTTGTTG-1
  SRR13927735#GTCCATCCATTTCACT-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-14 06:24:46.986727 : Completed PeakMatrix : SRR13927735(3 of 4), 0.608 mins elapsed.
2023-11-14 06:24:46.987133 : Reading PeakMatrix : SRR13927736(4 of 4), 0.608 mins elapsed.

2023-11-14 06:24:46.998779 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927736 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927736.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3534
[1] "SRR13927736#GCTGTTCTCCTTGACC-1" "SRR13927736#AGTTACGAGTAGGTCG-1"
[3] "SRR13927736#AGTTTGGAGCATTCCA-1" "SRR13927736#GCACGGTTCGGCAATT-1"
[5] "SRR13927736#TCAGCTCTCTCTATTG-1" "SRR13927736#CCTTAATTCTACTTTG-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-14 06:24:47.254298 : featureDF SRR13927736, Class = DFrame
DataFrame with 210721 rows and 4 columns
       seqnames       idx     start       end
          <Rle> <integer> <integer> <integer>
1          chr1         1    804677    805177
2          chr1         2    811866    812366
3          chr1         3    817090    817590
4          chr1         4    818831    819331
5          chr1         5    826543    827043
...         ...       ...       ...       ...
210717     chrX      5777 155898525 155899025
210718     chrX      5778 155899804 155900304
210719     chrX      5779 155956608 155957108
210720     chrX      5780 155966827 155967327
210721     chrX      5781 155971625 155972125

2023-11-14 06:24:47.293727 : Getting PeakMatrix from ArrowFile : SRR13927736.arrow, 0.005 mins elapsed.

2023-11-14 06:24:57.019682 : mat SRR13927736, Class = dgCMatrix
mat SRR13927736: nRows = 210721, nCols = 3534
mat SRR13927736: NonZeroEntries = 25583439, EntryRange = [ 1 , 4 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927736#GCTGTTCTCCTTGACC-1 SRR13927736#AGTTACGAGTAGGTCG-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              .                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927736#AGTTTGGAGCATTCCA-1 SRR13927736#GCACGGTTCGGCAATT-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              1                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927736#TCAGCTCTCTCTATTG-1
[1,]                              .
[2,]                              .
[3,]                              .
[4,]                              .
[5,]                              .


2023-11-14 06:24:57.222024 : Organizing SE PeakMatrix from ArrowFile : SRR13927736.arrow, 0.17 mins elapsed.

2023-11-14 06:24:57.492102 : se SRR13927736, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 3534 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(3534): SRR13927736#GCTGTTCTCCTTGACC-1
  SRR13927736#AGTTACGAGTAGGTCG-1 ... SRR13927736#TAAACCGCAAGTCTGT-1
  SRR13927736#TGCATTTGTTCTACCC-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-14 06:24:57.49905 : Completed PeakMatrix : SRR13927736(4 of 4), 0.783 mins elapsed.
2023-11-14 06:24:57.49993 : Organizing colData, 0.783 mins elapsed.
2023-11-14 06:24:57.639328 : Organizing rowData, 0.785 mins elapsed.
2023-11-14 06:24:57.647971 : Organizing rowRanges, 0.785 mins elapsed.
2023-11-14 06:24:57.658803 : Organizing Assays (1 of 1), 0.786 mins elapsed.
2023-11-14 06:25:00.849687 : Constructing SummarizedExperiment, 0.839 mins elapsed.
2023-11-14 06:25:13.839966 : Finished Matrix Creation, 1.055 mins elapsed.

2023-11-14 06:35:15.190466 : getMatrixFromProject Input-Parameters, Class = list

getMatrixFromProject Input-Parameters$ArchRProj: length = 1

getMatrixFromProject Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromProject Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromProject Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromProject Input-Parameters$verbose: length = 1
[1] TRUE


getMatrixFromProject Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromProject Input-Parameters$threads: length = 1
[1] 1


getMatrixFromProject Input-Parameters$logFile: length = 1
[1] "x"


2023-11-14 06:35:15.527953 : Reading GeneIntegrationMatrix : SRR13927737(1 of 4), 0.006 mins elapsed.

2023-11-14 06:35:15.543221 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927737 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927737.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 4439
[1] "SRR13927737#GTCACGGGTCAACAGG-1" "SRR13927737#TGCCTGTTCTTACTCA-1"
[3] "SRR13927737#TTGCGGGTCAGGTCTA-1" "SRR13927737#GGTTGCGCAAAGCTTC-1"
[5] "SRR13927737#CTTAATCTCGCGCTGA-1" "SRR13927737#TGTACGATCGGAAAGC-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-14 06:35:15.770224 : featureDF SRR13927737, Class = DFrame
DataFrame with 23525 rows and 6 columns
      seqnames     start       end    strand        name       idx
         <Rle> <integer> <integer> <integer> <character> <integer>
1         chr1     69091     70008         1       OR4F5         1
2         chr1    817371    819837         1      FAM87B         2
3         chr1    825138    859446         1   LINC01128         3
4         chr1    827522    826206         2   LINC00115         4
5         chr1    876903    868071         2      FAM41C         5
...        ...       ...       ...       ...         ...       ...
23521     chrX 155264589 155258241         2      RAB39B       923
23522     chrX 155334657 155276211         2       CLIC2       924
23523     chrX 155380787 155381134         1      H2AFB1       925
23524     chrX 155382115 155383230         1        F8A2       926
23525     chrX 155669944 155490115         2       TMLHE       927

2023-11-14 06:35:15.803408 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927737.arrow, 0.004 mins elapsed.

2023-11-14 06:35:48.267537 : mat SRR13927737, Class = dgCMatrix
mat SRR13927737: nRows = 23525, nCols = 4439
mat SRR13927737: NonZeroEntries = 50277774, EntryRange = [ 0.000148130054904683 , 376.345870564644 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927737#GTCACGGGTCAACAGG-1 SRR13927737#TGCCTGTTCTTACTCA-1
[1,]                      .                             .         
[2,]                      .                             .         
[3,]                      0.1057576                     0.09630696
[4,]                      .                             0.01806508
[5,]                      .                             .         
     SRR13927737#TTGCGGGTCAGGTCTA-1 SRR13927737#GGTTGCGCAAAGCTTC-1
[1,]                     .                              .         
[2,]                     .                              .         
[3,]                     0.09113499                     0.08396973
[4,]                     0.02096304                     0.01727603
[5,]                     .                              .         
     SRR13927737#CTTAATCTCGCGCTGA-1
[1,]                     .         
[2,]                     .         
[3,]                     0.08327649
[4,]                     0.01368917
[5,]                     .         


2023-11-14 06:35:48.737962 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927737.arrow, 0.553 mins elapsed.

2023-11-14 06:35:49.074943 : se SRR13927737, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 4439 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(4439): SRR13927737#GTCACGGGTCAACAGG-1
  SRR13927737#TGCCTGTTCTTACTCA-1 ... SRR13927737#TTGCTATAGATGTTCC-1
  SRR13927737#GCTCACTAGTTAGCGG-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-14 06:35:49.081278 : Completed GeneIntegrationMatrix : SRR13927737(1 of 4), 0.565 mins elapsed.
2023-11-14 06:35:49.081712 : Reading GeneIntegrationMatrix : SRR13927738(2 of 4), 0.565 mins elapsed.

2023-11-14 06:35:49.094614 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927738 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927738.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3736
[1] "SRR13927738#TAGCATGAGGCTAAAT-1" "SRR13927738#TGTTAGGTCGGTCAGC-1"
[3] "SRR13927738#CTCTCAGCATGTCCCT-1" "SRR13927738#ACTATTCGTCCAACCG-1"
[5] "SRR13927738#TGGTCCTTCTACTGCC-1" "SRR13927738#ATTGTCTAGACTAGGC-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-14 06:35:49.281341 : featureDF SRR13927738, Class = DFrame
DataFrame with 23525 rows and 6 columns
      seqnames     start       end    strand        name       idx
         <Rle> <integer> <integer> <integer> <character> <integer>
1         chr1     69091     70008         1       OR4F5         1
2         chr1    817371    819837         1      FAM87B         2
3         chr1    825138    859446         1   LINC01128         3
4         chr1    827522    826206         2   LINC00115         4
5         chr1    876903    868071         2      FAM41C         5
...        ...       ...       ...       ...         ...       ...
23521     chrX 155264589 155258241         2      RAB39B       923
23522     chrX 155334657 155276211         2       CLIC2       924
23523     chrX 155380787 155381134         1      H2AFB1       925
23524     chrX 155382115 155383230         1        F8A2       926
23525     chrX 155669944 155490115         2       TMLHE       927

2023-11-14 06:35:49.29518 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927738.arrow, 0.003 mins elapsed.

2023-11-14 06:36:06.600244 : mat SRR13927738, Class = dgCMatrix
mat SRR13927738: nRows = 23525, nCols = 3736
mat SRR13927738: NonZeroEntries = 41898897, EntryRange = [ 9.73975135001421e-05 , 758.763031303812 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927738#TAGCATGAGGCTAAAT-1 SRR13927738#TGTTAGGTCGGTCAGC-1
[1,]                    .                                .        
[2,]                    .                                .        
[3,]                    0.002346352                      0.0084575
[4,]                    0.005836413                      .        
[5,]                    .                                .        
     SRR13927738#CTCTCAGCATGTCCCT-1 SRR13927738#ACTATTCGTCCAACCG-1
[1,]                     .                             .          
[2,]                     .                             .          
[3,]                     .                             0.004081077
[4,]                     0.02001583                    .          
[5,]                     .                             .          
     SRR13927738#TGGTCCTTCTACTGCC-1
[1,]                    .          
[2,]                    .          
[3,]                    0.007640877
[4,]                    0.015048796
[5,]                    .          


2023-11-14 06:36:06.938364 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927738.arrow, 0.297 mins elapsed.

2023-11-14 06:36:07.941979 : se SRR13927738, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 3736 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(3736): SRR13927738#TAGCATGAGGCTAAAT-1
  SRR13927738#TGTTAGGTCGGTCAGC-1 ... SRR13927738#GCTCCTATCCGTGCGA-1
  SRR13927738#GGCACGTAGTCGATAA-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-14 06:36:07.947145 : Completed GeneIntegrationMatrix : SRR13927738(2 of 4), 0.879 mins elapsed.
2023-11-14 06:36:07.947985 : Reading GeneIntegrationMatrix : SRR13927735(3 of 4), 0.879 mins elapsed.

2023-11-14 06:36:07.962063 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927735 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927735.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3813
[1] "SRR13927735#TTATGTCTCCAGGTAT-1" "SRR13927735#TATTGCTCATCAGAAA-1"
[3] "SRR13927735#TTCGATTGTAGGGTTG-1" "SRR13927735#CATTCATTCGGATGTT-1"
[5] "SRR13927735#ACGTTAGGTCAACTGT-1" "SRR13927735#AAATGCCCAGCAATGG-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-14 06:36:08.147994 : featureDF SRR13927735, Class = DFrame
DataFrame with 23525 rows and 6 columns
      seqnames     start       end    strand        name       idx
         <Rle> <integer> <integer> <integer> <character> <integer>
1         chr1     69091     70008         1       OR4F5         1
2         chr1    817371    819837         1      FAM87B         2
3         chr1    825138    859446         1   LINC01128         3
4         chr1    827522    826206         2   LINC00115         4
5         chr1    876903    868071         2      FAM41C         5
...        ...       ...       ...       ...         ...       ...
23521     chrX 155264589 155258241         2      RAB39B       923
23522     chrX 155334657 155276211         2       CLIC2       924
23523     chrX 155380787 155381134         1      H2AFB1       925
23524     chrX 155382115 155383230         1        F8A2       926
23525     chrX 155669944 155490115         2       TMLHE       927

2023-11-14 06:36:08.162743 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927735.arrow, 0.003 mins elapsed.

2023-11-14 06:36:25.662636 : mat SRR13927735, Class = dgCMatrix
mat SRR13927735: nRows = 23525, nCols = 3813
mat SRR13927735: NonZeroEntries = 45535634, EntryRange = [ 7.20715653308254e-05 , 673.990433360078 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927735#TTATGTCTCCAGGTAT-1 SRR13927735#TATTGCTCATCAGAAA-1
[1,]                    .                              .          
[2,]                    .                              .          
[3,]                    0.070892732                    0.057009268
[4,]                    0.008800039                    0.002978408
[5,]                    .                              .          
     SRR13927735#TTCGATTGTAGGGTTG-1 SRR13927735#CATTCATTCGGATGTT-1
[1,]                     .                               .        
[2,]                     .                               .        
[3,]                     0.05523011                      0.0248301
[4,]                     .                               .        
[5,]                     .                               .        
     SRR13927735#ACGTTAGGTCAACTGT-1
[1,]                     .         
[2,]                     .         
[3,]                     0.05384571
[4,]                     0.01228688
[5,]                     .         


2023-11-14 06:36:26.017622 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927735.arrow, 0.301 mins elapsed.

2023-11-14 06:36:26.259511 : se SRR13927735, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 3813 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(3813): SRR13927735#TTATGTCTCCAGGTAT-1
  SRR13927735#TATTGCTCATCAGAAA-1 ... SRR13927735#GATCATGGTGTTGTTG-1
  SRR13927735#GTCCATCCATTTCACT-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-14 06:36:26.264706 : Completed GeneIntegrationMatrix : SRR13927735(3 of 4), 1.185 mins elapsed.
2023-11-14 06:36:26.265116 : Reading GeneIntegrationMatrix : SRR13927736(4 of 4), 1.185 mins elapsed.

2023-11-14 06:36:26.27723 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927736 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927736.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3534
[1] "SRR13927736#GCTGTTCTCCTTGACC-1" "SRR13927736#AGTTACGAGTAGGTCG-1"
[3] "SRR13927736#AGTTTGGAGCATTCCA-1" "SRR13927736#GCACGGTTCGGCAATT-1"
[5] "SRR13927736#TCAGCTCTCTCTATTG-1" "SRR13927736#CCTTAATTCTACTTTG-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-14 06:36:26.470202 : featureDF SRR13927736, Class = DFrame
DataFrame with 23525 rows and 6 columns
      seqnames     start       end    strand        name       idx
         <Rle> <integer> <integer> <integer> <character> <integer>
1         chr1     69091     70008         1       OR4F5         1
2         chr1    817371    819837         1      FAM87B         2
3         chr1    825138    859446         1   LINC01128         3
4         chr1    827522    826206         2   LINC00115         4
5         chr1    876903    868071         2      FAM41C         5
...        ...       ...       ...       ...         ...       ...
23521     chrX 155264589 155258241         2      RAB39B       923
23522     chrX 155334657 155276211         2       CLIC2       924
23523     chrX 155380787 155381134         1      H2AFB1       925
23524     chrX 155382115 155383230         1        F8A2       926
23525     chrX 155669944 155490115         2       TMLHE       927

2023-11-14 06:36:26.483498 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927736.arrow, 0.003 mins elapsed.

2023-11-14 06:36:41.58275 : mat SRR13927736, Class = dgCMatrix
mat SRR13927736: nRows = 23525, nCols = 3534
mat SRR13927736: NonZeroEntries = 39003616, EntryRange = [ 8.77627031078293e-05 , 770.387276519725 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927736#GCTGTTCTCCTTGACC-1 SRR13927736#AGTTACGAGTAGGTCG-1
[1,]                     .                             .          
[2,]                     .                             .          
[3,]                     0.07845381                    0.028344539
[4,]                     0.01395289                    0.008499189
[5,]                     .                             .          
     SRR13927736#AGTTTGGAGCATTCCA-1 SRR13927736#GCACGGTTCGGCAATT-1
[1,]                     .                              .         
[2,]                     .                              .         
[3,]                     0.06535130                     0.06576222
[4,]                     0.01810177                     0.01315728
[5,]                     .                              .         
     SRR13927736#TCAGCTCTCTCTATTG-1
[1,]                    .          
[2,]                    .          
[3,]                    0.016229526
[4,]                    0.003862592
[5,]                    .          


2023-11-14 06:36:41.885075 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927736.arrow, 0.26 mins elapsed.

2023-11-14 06:36:42.13638 : se SRR13927736, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 3534 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(3534): SRR13927736#GCTGTTCTCCTTGACC-1
  SRR13927736#AGTTACGAGTAGGTCG-1 ... SRR13927736#TAAACCGCAAGTCTGT-1
  SRR13927736#TGCATTTGTTCTACCC-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-14 06:36:42.141623 : Completed GeneIntegrationMatrix : SRR13927736(4 of 4), 1.449 mins elapsed.
2023-11-14 06:36:42.142594 : Organizing colData, 1.449 mins elapsed.
2023-11-14 06:36:42.286047 : Organizing rowData, 1.452 mins elapsed.
2023-11-14 06:36:42.289285 : Organizing rowRanges, 1.452 mins elapsed.
2023-11-14 06:36:42.295074 : Organizing Assays (1 of 1), 1.452 mins elapsed.
2023-11-14 06:36:48.605098 : Constructing SummarizedExperiment, 1.557 mins elapsed.
2023-11-14 06:36:49.324662 : Finished Matrix Creation, 1.569 mins elapsed.

2023-11-14 06:36:49.327018 : getMatrixFromProject Input-Parameters, Class = list

getMatrixFromProject Input-Parameters$ArchRProj: length = 1

getMatrixFromProject Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromProject Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromProject Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromProject Input-Parameters$verbose: length = 1
[1] TRUE


getMatrixFromProject Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromProject Input-Parameters$threads: length = 1
[1] 1


getMatrixFromProject Input-Parameters$logFile: length = 1
[1] "x"


2023-11-14 06:36:49.777542 : Reading PeakMatrix : SRR13927737(1 of 4), 0.008 mins elapsed.

2023-11-14 06:36:49.790822 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927737 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927737.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 4439
[1] "SRR13927737#GTCACGGGTCAACAGG-1" "SRR13927737#TGCCTGTTCTTACTCA-1"
[3] "SRR13927737#TTGCGGGTCAGGTCTA-1" "SRR13927737#GGTTGCGCAAAGCTTC-1"
[5] "SRR13927737#CTTAATCTCGCGCTGA-1" "SRR13927737#TGTACGATCGGAAAGC-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-14 06:36:50.080656 : featureDF SRR13927737, Class = DFrame
DataFrame with 210721 rows and 4 columns
       seqnames       idx     start       end
          <Rle> <integer> <integer> <integer>
1          chr1         1    804677    805177
2          chr1         2    811866    812366
3          chr1         3    817090    817590
4          chr1         4    818831    819331
5          chr1         5    826543    827043
...         ...       ...       ...       ...
210717     chrX      5777 155898525 155899025
210718     chrX      5778 155899804 155900304
210719     chrX      5779 155956608 155957108
210720     chrX      5780 155966827 155967327
210721     chrX      5781 155971625 155972125

2023-11-14 06:36:50.132908 : Getting PeakMatrix from ArrowFile : SRR13927737.arrow, 0.006 mins elapsed.

2023-11-14 06:37:01.908879 : mat SRR13927737, Class = dgCMatrix
mat SRR13927737: nRows = 210721, nCols = 4439
mat SRR13927737: NonZeroEntries = 29963541, EntryRange = [ 1 , 4 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927737#GTCACGGGTCAACAGG-1 SRR13927737#TGCCTGTTCTTACTCA-1
[1,]                              .                              .
[2,]                              .                              2
[3,]                              .                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927737#TTGCGGGTCAGGTCTA-1 SRR13927737#GGTTGCGCAAAGCTTC-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              2                              .
[4,]                              .                              .
[5,]                              4                              .
     SRR13927737#CTTAATCTCGCGCTGA-1
[1,]                              .
[2,]                              .
[3,]                              .
[4,]                              .
[5,]                              .


2023-11-14 06:37:02.094123 : Organizing SE PeakMatrix from ArrowFile : SRR13927737.arrow, 0.205 mins elapsed.

2023-11-14 06:37:02.397633 : se SRR13927737, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 4439 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(4439): SRR13927737#GTCACGGGTCAACAGG-1
  SRR13927737#TGCCTGTTCTTACTCA-1 ... SRR13927737#TTGCTATAGATGTTCC-1
  SRR13927737#GCTCACTAGTTAGCGG-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-14 06:37:02.406753 : Completed PeakMatrix : SRR13927737(1 of 4), 0.218 mins elapsed.
2023-11-14 06:37:02.407127 : Reading PeakMatrix : SRR13927738(2 of 4), 0.218 mins elapsed.

2023-11-14 06:37:02.418982 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927738 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927738.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3736
[1] "SRR13927738#TAGCATGAGGCTAAAT-1" "SRR13927738#TGTTAGGTCGGTCAGC-1"
[3] "SRR13927738#CTCTCAGCATGTCCCT-1" "SRR13927738#ACTATTCGTCCAACCG-1"
[5] "SRR13927738#TGGTCCTTCTACTGCC-1" "SRR13927738#ATTGTCTAGACTAGGC-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-14 06:37:02.663679 : featureDF SRR13927738, Class = DFrame
DataFrame with 210721 rows and 4 columns
       seqnames       idx     start       end
          <Rle> <integer> <integer> <integer>
1          chr1         1    804677    805177
2          chr1         2    811866    812366
3          chr1         3    817090    817590
4          chr1         4    818831    819331
5          chr1         5    826543    827043
...         ...       ...       ...       ...
210717     chrX      5777 155898525 155899025
210718     chrX      5778 155899804 155900304
210719     chrX      5779 155956608 155957108
210720     chrX      5780 155966827 155967327
210721     chrX      5781 155971625 155972125

2023-11-14 06:37:02.703845 : Getting PeakMatrix from ArrowFile : SRR13927738.arrow, 0.005 mins elapsed.

2023-11-14 06:37:14.038735 : mat SRR13927738, Class = dgCMatrix
mat SRR13927738: nRows = 210721, nCols = 3736
mat SRR13927738: NonZeroEntries = 30165921, EntryRange = [ 1 , 4 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927738#TAGCATGAGGCTAAAT-1 SRR13927738#TGTTAGGTCGGTCAGC-1
[1,]                              2                              .
[2,]                              .                              .
[3,]                              .                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927738#CTCTCAGCATGTCCCT-1 SRR13927738#ACTATTCGTCCAACCG-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              .                              .
[4,]                              .                              .
[5,]                              2                              .
     SRR13927738#TGGTCCTTCTACTGCC-1
[1,]                              .
[2,]                              .
[3,]                              2
[4,]                              .
[5,]                              .


2023-11-14 06:37:14.264364 : Organizing SE PeakMatrix from ArrowFile : SRR13927738.arrow, 0.197 mins elapsed.

2023-11-14 06:37:14.554611 : se SRR13927738, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 3736 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(3736): SRR13927738#TAGCATGAGGCTAAAT-1
  SRR13927738#TGTTAGGTCGGTCAGC-1 ... SRR13927738#GCTCCTATCCGTGCGA-1
  SRR13927738#GGCACGTAGTCGATAA-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-14 06:37:14.561294 : Completed PeakMatrix : SRR13927738(2 of 4), 0.421 mins elapsed.
2023-11-14 06:37:14.561702 : Reading PeakMatrix : SRR13927735(3 of 4), 0.421 mins elapsed.

2023-11-14 06:37:14.573953 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927735 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927735.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3813
[1] "SRR13927735#TTATGTCTCCAGGTAT-1" "SRR13927735#TATTGCTCATCAGAAA-1"
[3] "SRR13927735#TTCGATTGTAGGGTTG-1" "SRR13927735#CATTCATTCGGATGTT-1"
[5] "SRR13927735#ACGTTAGGTCAACTGT-1" "SRR13927735#AAATGCCCAGCAATGG-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-14 06:37:14.837072 : featureDF SRR13927735, Class = DFrame
DataFrame with 210721 rows and 4 columns
       seqnames       idx     start       end
          <Rle> <integer> <integer> <integer>
1          chr1         1    804677    805177
2          chr1         2    811866    812366
3          chr1         3    817090    817590
4          chr1         4    818831    819331
5          chr1         5    826543    827043
...         ...       ...       ...       ...
210717     chrX      5777 155898525 155899025
210718     chrX      5778 155899804 155900304
210719     chrX      5779 155956608 155957108
210720     chrX      5780 155966827 155967327
210721     chrX      5781 155971625 155972125

2023-11-14 06:37:14.88193 : Getting PeakMatrix from ArrowFile : SRR13927735.arrow, 0.005 mins elapsed.

2023-11-14 06:37:25.245642 : mat SRR13927735, Class = dgCMatrix
mat SRR13927735: nRows = 210721, nCols = 3813
mat SRR13927735: NonZeroEntries = 26038084, EntryRange = [ 1 , 4 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927735#TTATGTCTCCAGGTAT-1 SRR13927735#TATTGCTCATCAGAAA-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              1                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927735#TTCGATTGTAGGGTTG-1 SRR13927735#CATTCATTCGGATGTT-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              .                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927735#ACGTTAGGTCAACTGT-1
[1,]                              2
[2,]                              .
[3,]                              .
[4,]                              .
[5,]                              .


2023-11-14 06:37:25.452072 : Organizing SE PeakMatrix from ArrowFile : SRR13927735.arrow, 0.181 mins elapsed.

2023-11-14 06:37:25.751175 : se SRR13927735, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 3813 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(3813): SRR13927735#TTATGTCTCCAGGTAT-1
  SRR13927735#TATTGCTCATCAGAAA-1 ... SRR13927735#GATCATGGTGTTGTTG-1
  SRR13927735#GTCCATCCATTTCACT-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-14 06:37:25.757921 : Completed PeakMatrix : SRR13927735(3 of 4), 0.607 mins elapsed.
2023-11-14 06:37:25.758317 : Reading PeakMatrix : SRR13927736(4 of 4), 0.607 mins elapsed.

2023-11-14 06:37:25.769953 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927736 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927736.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3534
[1] "SRR13927736#GCTGTTCTCCTTGACC-1" "SRR13927736#AGTTACGAGTAGGTCG-1"
[3] "SRR13927736#AGTTTGGAGCATTCCA-1" "SRR13927736#GCACGGTTCGGCAATT-1"
[5] "SRR13927736#TCAGCTCTCTCTATTG-1" "SRR13927736#CCTTAATTCTACTTTG-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-14 06:37:26.02523 : featureDF SRR13927736, Class = DFrame
DataFrame with 210721 rows and 4 columns
       seqnames       idx     start       end
          <Rle> <integer> <integer> <integer>
1          chr1         1    804677    805177
2          chr1         2    811866    812366
3          chr1         3    817090    817590
4          chr1         4    818831    819331
5          chr1         5    826543    827043
...         ...       ...       ...       ...
210717     chrX      5777 155898525 155899025
210718     chrX      5778 155899804 155900304
210719     chrX      5779 155956608 155957108
210720     chrX      5780 155966827 155967327
210721     chrX      5781 155971625 155972125

2023-11-14 06:37:26.065293 : Getting PeakMatrix from ArrowFile : SRR13927736.arrow, 0.005 mins elapsed.

2023-11-14 06:37:35.834568 : mat SRR13927736, Class = dgCMatrix
mat SRR13927736: nRows = 210721, nCols = 3534
mat SRR13927736: NonZeroEntries = 25583439, EntryRange = [ 1 , 4 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927736#GCTGTTCTCCTTGACC-1 SRR13927736#AGTTACGAGTAGGTCG-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              .                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927736#AGTTTGGAGCATTCCA-1 SRR13927736#GCACGGTTCGGCAATT-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              1                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927736#TCAGCTCTCTCTATTG-1
[1,]                              .
[2,]                              .
[3,]                              .
[4,]                              .
[5,]                              .


2023-11-14 06:37:36.033768 : Organizing SE PeakMatrix from ArrowFile : SRR13927736.arrow, 0.171 mins elapsed.

2023-11-14 06:37:36.303172 : se SRR13927736, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 3534 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(3534): SRR13927736#GCTGTTCTCCTTGACC-1
  SRR13927736#AGTTACGAGTAGGTCG-1 ... SRR13927736#TAAACCGCAAGTCTGT-1
  SRR13927736#TGCATTTGTTCTACCC-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-14 06:37:36.309711 : Completed PeakMatrix : SRR13927736(4 of 4), 0.783 mins elapsed.
2023-11-14 06:37:36.31054 : Organizing colData, 0.783 mins elapsed.
2023-11-14 06:37:36.446 : Organizing rowData, 0.785 mins elapsed.
2023-11-14 06:37:36.453816 : Organizing rowRanges, 0.785 mins elapsed.
2023-11-14 06:37:36.464189 : Organizing Assays (1 of 1), 0.786 mins elapsed.
2023-11-14 06:37:39.490331 : Constructing SummarizedExperiment, 0.836 mins elapsed.
2023-11-14 06:37:52.395633 : Finished Matrix Creation, 1.051 mins elapsed.

2023-11-14 06:58:06.962048 : getMatrixFromProject Input-Parameters, Class = list

getMatrixFromProject Input-Parameters$ArchRProj: length = 1

getMatrixFromProject Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromProject Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromProject Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromProject Input-Parameters$verbose: length = 1
[1] TRUE


getMatrixFromProject Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromProject Input-Parameters$threads: length = 1
[1] 1


getMatrixFromProject Input-Parameters$logFile: length = 1
[1] "x"


2023-11-14 06:58:07.414514 : Reading GeneIntegrationMatrix : SRR13927737(1 of 4), 0.008 mins elapsed.

2023-11-14 06:58:07.431345 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927737 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927737.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 4439
[1] "SRR13927737#GTCACGGGTCAACAGG-1" "SRR13927737#TGCCTGTTCTTACTCA-1"
[3] "SRR13927737#TTGCGGGTCAGGTCTA-1" "SRR13927737#GGTTGCGCAAAGCTTC-1"
[5] "SRR13927737#CTTAATCTCGCGCTGA-1" "SRR13927737#TGTACGATCGGAAAGC-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-14 06:58:07.689675 : featureDF SRR13927737, Class = DFrame
DataFrame with 23525 rows and 6 columns
      seqnames     start       end    strand        name       idx
         <Rle> <integer> <integer> <integer> <character> <integer>
1         chr1     69091     70008         1       OR4F5         1
2         chr1    817371    819837         1      FAM87B         2
3         chr1    825138    859446         1   LINC01128         3
4         chr1    827522    826206         2   LINC00115         4
5         chr1    876903    868071         2      FAM41C         5
...        ...       ...       ...       ...         ...       ...
23521     chrX 155264589 155258241         2      RAB39B       923
23522     chrX 155334657 155276211         2       CLIC2       924
23523     chrX 155380787 155381134         1      H2AFB1       925
23524     chrX 155382115 155383230         1        F8A2       926
23525     chrX 155669944 155490115         2       TMLHE       927

2023-11-14 06:58:07.728291 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927737.arrow, 0.005 mins elapsed.

2023-11-14 06:58:40.570216 : mat SRR13927737, Class = dgCMatrix
mat SRR13927737: nRows = 23525, nCols = 4439
mat SRR13927737: NonZeroEntries = 50277774, EntryRange = [ 0.000148130054904683 , 376.345870564644 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927737#GTCACGGGTCAACAGG-1 SRR13927737#TGCCTGTTCTTACTCA-1
[1,]                      .                             .         
[2,]                      .                             .         
[3,]                      0.1057576                     0.09630696
[4,]                      .                             0.01806508
[5,]                      .                             .         
     SRR13927737#TTGCGGGTCAGGTCTA-1 SRR13927737#GGTTGCGCAAAGCTTC-1
[1,]                     .                              .         
[2,]                     .                              .         
[3,]                     0.09113499                     0.08396973
[4,]                     0.02096304                     0.01727603
[5,]                     .                              .         
     SRR13927737#CTTAATCTCGCGCTGA-1
[1,]                     .         
[2,]                     .         
[3,]                     0.08327649
[4,]                     0.01368917
[5,]                     .         


2023-11-14 06:58:41.042058 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927737.arrow, 0.56 mins elapsed.

2023-11-14 06:58:41.606438 : se SRR13927737, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 4439 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(4439): SRR13927737#GTCACGGGTCAACAGG-1
  SRR13927737#TGCCTGTTCTTACTCA-1 ... SRR13927737#TTGCTATAGATGTTCC-1
  SRR13927737#GCTCACTAGTTAGCGG-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-14 06:58:41.613993 : Completed GeneIntegrationMatrix : SRR13927737(1 of 4), 0.578 mins elapsed.
2023-11-14 06:58:41.614496 : Reading GeneIntegrationMatrix : SRR13927738(2 of 4), 0.578 mins elapsed.

2023-11-14 06:58:41.628477 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927738 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927738.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3736
[1] "SRR13927738#TAGCATGAGGCTAAAT-1" "SRR13927738#TGTTAGGTCGGTCAGC-1"
[3] "SRR13927738#CTCTCAGCATGTCCCT-1" "SRR13927738#ACTATTCGTCCAACCG-1"
[5] "SRR13927738#TGGTCCTTCTACTGCC-1" "SRR13927738#ATTGTCTAGACTAGGC-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-14 06:58:41.850234 : featureDF SRR13927738, Class = DFrame
DataFrame with 23525 rows and 6 columns
      seqnames     start       end    strand        name       idx
         <Rle> <integer> <integer> <integer> <character> <integer>
1         chr1     69091     70008         1       OR4F5         1
2         chr1    817371    819837         1      FAM87B         2
3         chr1    825138    859446         1   LINC01128         3
4         chr1    827522    826206         2   LINC00115         4
5         chr1    876903    868071         2      FAM41C         5
...        ...       ...       ...       ...         ...       ...
23521     chrX 155264589 155258241         2      RAB39B       923
23522     chrX 155334657 155276211         2       CLIC2       924
23523     chrX 155380787 155381134         1      H2AFB1       925
23524     chrX 155382115 155383230         1        F8A2       926
23525     chrX 155669944 155490115         2       TMLHE       927

2023-11-14 06:58:41.869217 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927738.arrow, 0.004 mins elapsed.

2023-11-14 06:58:59.45431 : mat SRR13927738, Class = dgCMatrix
mat SRR13927738: nRows = 23525, nCols = 3736
mat SRR13927738: NonZeroEntries = 41898897, EntryRange = [ 9.73975135001421e-05 , 758.763031303812 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927738#TAGCATGAGGCTAAAT-1 SRR13927738#TGTTAGGTCGGTCAGC-1
[1,]                    .                                .        
[2,]                    .                                .        
[3,]                    0.002346352                      0.0084575
[4,]                    0.005836413                      .        
[5,]                    .                                .        
     SRR13927738#CTCTCAGCATGTCCCT-1 SRR13927738#ACTATTCGTCCAACCG-1
[1,]                     .                             .          
[2,]                     .                             .          
[3,]                     .                             0.004081077
[4,]                     0.02001583                    .          
[5,]                     .                             .          
     SRR13927738#TGGTCCTTCTACTGCC-1
[1,]                    .          
[2,]                    .          
[3,]                    0.007640877
[4,]                    0.015048796
[5,]                    .          


2023-11-14 06:58:59.785673 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927738.arrow, 0.303 mins elapsed.

2023-11-14 06:59:00.796656 : se SRR13927738, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 3736 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(3736): SRR13927738#TAGCATGAGGCTAAAT-1
  SRR13927738#TGTTAGGTCGGTCAGC-1 ... SRR13927738#GCTCCTATCCGTGCGA-1
  SRR13927738#GGCACGTAGTCGATAA-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-14 06:59:00.801816 : Completed GeneIntegrationMatrix : SRR13927738(2 of 4), 0.897 mins elapsed.
2023-11-14 06:59:00.802229 : Reading GeneIntegrationMatrix : SRR13927735(3 of 4), 0.897 mins elapsed.

2023-11-14 06:59:00.814285 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927735 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927735.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3813
[1] "SRR13927735#TTATGTCTCCAGGTAT-1" "SRR13927735#TATTGCTCATCAGAAA-1"
[3] "SRR13927735#TTCGATTGTAGGGTTG-1" "SRR13927735#CATTCATTCGGATGTT-1"
[5] "SRR13927735#ACGTTAGGTCAACTGT-1" "SRR13927735#AAATGCCCAGCAATGG-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-14 06:59:01.006088 : featureDF SRR13927735, Class = DFrame
DataFrame with 23525 rows and 6 columns
      seqnames     start       end    strand        name       idx
         <Rle> <integer> <integer> <integer> <character> <integer>
1         chr1     69091     70008         1       OR4F5         1
2         chr1    817371    819837         1      FAM87B         2
3         chr1    825138    859446         1   LINC01128         3
4         chr1    827522    826206         2   LINC00115         4
5         chr1    876903    868071         2      FAM41C         5
...        ...       ...       ...       ...         ...       ...
23521     chrX 155264589 155258241         2      RAB39B       923
23522     chrX 155334657 155276211         2       CLIC2       924
23523     chrX 155380787 155381134         1      H2AFB1       925
23524     chrX 155382115 155383230         1        F8A2       926
23525     chrX 155669944 155490115         2       TMLHE       927

2023-11-14 06:59:01.020396 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927735.arrow, 0.003 mins elapsed.

2023-11-14 06:59:18.550255 : mat SRR13927735, Class = dgCMatrix
mat SRR13927735: nRows = 23525, nCols = 3813
mat SRR13927735: NonZeroEntries = 45535634, EntryRange = [ 7.20715653308254e-05 , 673.990433360078 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927735#TTATGTCTCCAGGTAT-1 SRR13927735#TATTGCTCATCAGAAA-1
[1,]                    .                              .          
[2,]                    .                              .          
[3,]                    0.070892732                    0.057009268
[4,]                    0.008800039                    0.002978408
[5,]                    .                              .          
     SRR13927735#TTCGATTGTAGGGTTG-1 SRR13927735#CATTCATTCGGATGTT-1
[1,]                     .                               .        
[2,]                     .                               .        
[3,]                     0.05523011                      0.0248301
[4,]                     .                               .        
[5,]                     .                               .        
     SRR13927735#ACGTTAGGTCAACTGT-1
[1,]                     .         
[2,]                     .         
[3,]                     0.05384571
[4,]                     0.01228688
[5,]                     .         


2023-11-14 06:59:18.910704 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927735.arrow, 0.302 mins elapsed.

2023-11-14 06:59:19.16346 : se SRR13927735, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 3813 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(3813): SRR13927735#TTATGTCTCCAGGTAT-1
  SRR13927735#TATTGCTCATCAGAAA-1 ... SRR13927735#GATCATGGTGTTGTTG-1
  SRR13927735#GTCCATCCATTTCACT-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-14 06:59:19.168503 : Completed GeneIntegrationMatrix : SRR13927735(3 of 4), 1.203 mins elapsed.
2023-11-14 06:59:19.168909 : Reading GeneIntegrationMatrix : SRR13927736(4 of 4), 1.203 mins elapsed.

2023-11-14 06:59:19.180782 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927736 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927736.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3534
[1] "SRR13927736#GCTGTTCTCCTTGACC-1" "SRR13927736#AGTTACGAGTAGGTCG-1"
[3] "SRR13927736#AGTTTGGAGCATTCCA-1" "SRR13927736#GCACGGTTCGGCAATT-1"
[5] "SRR13927736#TCAGCTCTCTCTATTG-1" "SRR13927736#CCTTAATTCTACTTTG-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-14 06:59:19.466416 : featureDF SRR13927736, Class = DFrame
DataFrame with 23525 rows and 6 columns
      seqnames     start       end    strand        name       idx
         <Rle> <integer> <integer> <integer> <character> <integer>
1         chr1     69091     70008         1       OR4F5         1
2         chr1    817371    819837         1      FAM87B         2
3         chr1    825138    859446         1   LINC01128         3
4         chr1    827522    826206         2   LINC00115         4
5         chr1    876903    868071         2      FAM41C         5
...        ...       ...       ...       ...         ...       ...
23521     chrX 155264589 155258241         2      RAB39B       923
23522     chrX 155334657 155276211         2       CLIC2       924
23523     chrX 155380787 155381134         1      H2AFB1       925
23524     chrX 155382115 155383230         1        F8A2       926
23525     chrX 155669944 155490115         2       TMLHE       927

2023-11-14 06:59:19.48014 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927736.arrow, 0.005 mins elapsed.

2023-11-14 06:59:34.923305 : mat SRR13927736, Class = dgCMatrix
mat SRR13927736: nRows = 23525, nCols = 3534
mat SRR13927736: NonZeroEntries = 39003616, EntryRange = [ 8.77627031078293e-05 , 770.387276519725 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927736#GCTGTTCTCCTTGACC-1 SRR13927736#AGTTACGAGTAGGTCG-1
[1,]                     .                             .          
[2,]                     .                             .          
[3,]                     0.07845381                    0.028344539
[4,]                     0.01395289                    0.008499189
[5,]                     .                             .          
     SRR13927736#AGTTTGGAGCATTCCA-1 SRR13927736#GCACGGTTCGGCAATT-1
[1,]                     .                              .         
[2,]                     .                              .         
[3,]                     0.06535130                     0.06576222
[4,]                     0.01810177                     0.01315728
[5,]                     .                              .         
     SRR13927736#TCAGCTCTCTCTATTG-1
[1,]                    .          
[2,]                    .          
[3,]                    0.016229526
[4,]                    0.003862592
[5,]                    .          


2023-11-14 06:59:35.239313 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927736.arrow, 0.268 mins elapsed.

2023-11-14 06:59:35.494564 : se SRR13927736, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 3534 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(3534): SRR13927736#GCTGTTCTCCTTGACC-1
  SRR13927736#AGTTACGAGTAGGTCG-1 ... SRR13927736#TAAACCGCAAGTCTGT-1
  SRR13927736#TGCATTTGTTCTACCC-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-14 06:59:35.499233 : Completed GeneIntegrationMatrix : SRR13927736(4 of 4), 1.476 mins elapsed.
2023-11-14 06:59:35.500133 : Organizing colData, 1.476 mins elapsed.
2023-11-14 06:59:35.640486 : Organizing rowData, 1.478 mins elapsed.
2023-11-14 06:59:35.643708 : Organizing rowRanges, 1.478 mins elapsed.
2023-11-14 06:59:35.64958 : Organizing Assays (1 of 1), 1.478 mins elapsed.
2023-11-14 06:59:41.920555 : Constructing SummarizedExperiment, 1.583 mins elapsed.
2023-11-14 06:59:42.635093 : Finished Matrix Creation, 1.595 mins elapsed.

2023-11-14 06:59:42.637406 : getMatrixFromProject Input-Parameters, Class = list

getMatrixFromProject Input-Parameters$ArchRProj: length = 1

getMatrixFromProject Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromProject Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromProject Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromProject Input-Parameters$verbose: length = 1
[1] TRUE


getMatrixFromProject Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromProject Input-Parameters$threads: length = 1
[1] 1


getMatrixFromProject Input-Parameters$logFile: length = 1
[1] "x"


2023-11-14 06:59:43.608027 : Reading PeakMatrix : SRR13927737(1 of 4), 0.016 mins elapsed.

2023-11-14 06:59:43.634446 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927737 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927737.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 4439
[1] "SRR13927737#GTCACGGGTCAACAGG-1" "SRR13927737#TGCCTGTTCTTACTCA-1"
[3] "SRR13927737#TTGCGGGTCAGGTCTA-1" "SRR13927737#GGTTGCGCAAAGCTTC-1"
[5] "SRR13927737#CTTAATCTCGCGCTGA-1" "SRR13927737#TGTACGATCGGAAAGC-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-14 06:59:44.062055 : featureDF SRR13927737, Class = DFrame
DataFrame with 210721 rows and 4 columns
       seqnames       idx     start       end
          <Rle> <integer> <integer> <integer>
1          chr1         1    804677    805177
2          chr1         2    811866    812366
3          chr1         3    817090    817590
4          chr1         4    818831    819331
5          chr1         5    826543    827043
...         ...       ...       ...       ...
210717     chrX      5777 155898525 155899025
210718     chrX      5778 155899804 155900304
210719     chrX      5779 155956608 155957108
210720     chrX      5780 155966827 155967327
210721     chrX      5781 155971625 155972125

2023-11-14 06:59:44.128192 : Getting PeakMatrix from ArrowFile : SRR13927737.arrow, 0.008 mins elapsed.

2023-11-14 06:59:56.385721 : mat SRR13927737, Class = dgCMatrix
mat SRR13927737: nRows = 210721, nCols = 4439
mat SRR13927737: NonZeroEntries = 29963541, EntryRange = [ 1 , 4 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927737#GTCACGGGTCAACAGG-1 SRR13927737#TGCCTGTTCTTACTCA-1
[1,]                              .                              .
[2,]                              .                              2
[3,]                              .                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927737#TTGCGGGTCAGGTCTA-1 SRR13927737#GGTTGCGCAAAGCTTC-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              2                              .
[4,]                              .                              .
[5,]                              4                              .
     SRR13927737#CTTAATCTCGCGCTGA-1
[1,]                              .
[2,]                              .
[3,]                              .
[4,]                              .
[5,]                              .


2023-11-14 06:59:56.577961 : Organizing SE PeakMatrix from ArrowFile : SRR13927737.arrow, 0.216 mins elapsed.

2023-11-14 06:59:56.878696 : se SRR13927737, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 4439 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(4439): SRR13927737#GTCACGGGTCAACAGG-1
  SRR13927737#TGCCTGTTCTTACTCA-1 ... SRR13927737#TTGCTATAGATGTTCC-1
  SRR13927737#GCTCACTAGTTAGCGG-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-14 06:59:56.887922 : Completed PeakMatrix : SRR13927737(1 of 4), 0.238 mins elapsed.
2023-11-14 06:59:56.888318 : Reading PeakMatrix : SRR13927738(2 of 4), 0.238 mins elapsed.

2023-11-14 06:59:56.900238 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927738 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927738.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3736
[1] "SRR13927738#TAGCATGAGGCTAAAT-1" "SRR13927738#TGTTAGGTCGGTCAGC-1"
[3] "SRR13927738#CTCTCAGCATGTCCCT-1" "SRR13927738#ACTATTCGTCCAACCG-1"
[5] "SRR13927738#TGGTCCTTCTACTGCC-1" "SRR13927738#ATTGTCTAGACTAGGC-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-14 06:59:57.136898 : featureDF SRR13927738, Class = DFrame
DataFrame with 210721 rows and 4 columns
       seqnames       idx     start       end
          <Rle> <integer> <integer> <integer>
1          chr1         1    804677    805177
2          chr1         2    811866    812366
3          chr1         3    817090    817590
4          chr1         4    818831    819331
5          chr1         5    826543    827043
...         ...       ...       ...       ...
210717     chrX      5777 155898525 155899025
210718     chrX      5778 155899804 155900304
210719     chrX      5779 155956608 155957108
210720     chrX      5780 155966827 155967327
210721     chrX      5781 155971625 155972125

2023-11-14 06:59:57.179601 : Getting PeakMatrix from ArrowFile : SRR13927738.arrow, 0.005 mins elapsed.

2023-11-14 07:00:08.485265 : mat SRR13927738, Class = dgCMatrix
mat SRR13927738: nRows = 210721, nCols = 3736
mat SRR13927738: NonZeroEntries = 30165921, EntryRange = [ 1 , 4 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927738#TAGCATGAGGCTAAAT-1 SRR13927738#TGTTAGGTCGGTCAGC-1
[1,]                              2                              .
[2,]                              .                              .
[3,]                              .                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927738#CTCTCAGCATGTCCCT-1 SRR13927738#ACTATTCGTCCAACCG-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              .                              .
[4,]                              .                              .
[5,]                              2                              .
     SRR13927738#TGGTCCTTCTACTGCC-1
[1,]                              .
[2,]                              .
[3,]                              2
[4,]                              .
[5,]                              .


2023-11-14 07:00:08.711497 : Organizing SE PeakMatrix from ArrowFile : SRR13927738.arrow, 0.197 mins elapsed.

2023-11-14 07:00:08.992138 : se SRR13927738, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 3736 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(3736): SRR13927738#TAGCATGAGGCTAAAT-1
  SRR13927738#TGTTAGGTCGGTCAGC-1 ... SRR13927738#GCTCCTATCCGTGCGA-1
  SRR13927738#GGCACGTAGTCGATAA-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-14 07:00:08.999252 : Completed PeakMatrix : SRR13927738(2 of 4), 0.439 mins elapsed.
2023-11-14 07:00:08.999668 : Reading PeakMatrix : SRR13927735(3 of 4), 0.439 mins elapsed.

2023-11-14 07:00:09.012646 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927735 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927735.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3813
[1] "SRR13927735#TTATGTCTCCAGGTAT-1" "SRR13927735#TATTGCTCATCAGAAA-1"
[3] "SRR13927735#TTCGATTGTAGGGTTG-1" "SRR13927735#CATTCATTCGGATGTT-1"
[5] "SRR13927735#ACGTTAGGTCAACTGT-1" "SRR13927735#AAATGCCCAGCAATGG-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-14 07:00:09.259805 : featureDF SRR13927735, Class = DFrame
DataFrame with 210721 rows and 4 columns
       seqnames       idx     start       end
          <Rle> <integer> <integer> <integer>
1          chr1         1    804677    805177
2          chr1         2    811866    812366
3          chr1         3    817090    817590
4          chr1         4    818831    819331
5          chr1         5    826543    827043
...         ...       ...       ...       ...
210717     chrX      5777 155898525 155899025
210718     chrX      5778 155899804 155900304
210719     chrX      5779 155956608 155957108
210720     chrX      5780 155966827 155967327
210721     chrX      5781 155971625 155972125

2023-11-14 07:00:09.300091 : Getting PeakMatrix from ArrowFile : SRR13927735.arrow, 0.005 mins elapsed.

2023-11-14 07:00:19.063619 : mat SRR13927735, Class = dgCMatrix
mat SRR13927735: nRows = 210721, nCols = 3813
mat SRR13927735: NonZeroEntries = 26038084, EntryRange = [ 1 , 4 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927735#TTATGTCTCCAGGTAT-1 SRR13927735#TATTGCTCATCAGAAA-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              1                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927735#TTCGATTGTAGGGTTG-1 SRR13927735#CATTCATTCGGATGTT-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              .                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927735#ACGTTAGGTCAACTGT-1
[1,]                              2
[2,]                              .
[3,]                              .
[4,]                              .
[5,]                              .


2023-11-14 07:00:19.266016 : Organizing SE PeakMatrix from ArrowFile : SRR13927735.arrow, 0.171 mins elapsed.

2023-11-14 07:00:19.536166 : se SRR13927735, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 3813 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(3813): SRR13927735#TTATGTCTCCAGGTAT-1
  SRR13927735#TATTGCTCATCAGAAA-1 ... SRR13927735#GATCATGGTGTTGTTG-1
  SRR13927735#GTCCATCCATTTCACT-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-14 07:00:19.543207 : Completed PeakMatrix : SRR13927735(3 of 4), 0.615 mins elapsed.
2023-11-14 07:00:19.543627 : Reading PeakMatrix : SRR13927736(4 of 4), 0.615 mins elapsed.

2023-11-14 07:00:19.555483 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927736 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927736.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3534
[1] "SRR13927736#GCTGTTCTCCTTGACC-1" "SRR13927736#AGTTACGAGTAGGTCG-1"
[3] "SRR13927736#AGTTTGGAGCATTCCA-1" "SRR13927736#GCACGGTTCGGCAATT-1"
[5] "SRR13927736#TCAGCTCTCTCTATTG-1" "SRR13927736#CCTTAATTCTACTTTG-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-14 07:00:19.869977 : featureDF SRR13927736, Class = DFrame
DataFrame with 210721 rows and 4 columns
       seqnames       idx     start       end
          <Rle> <integer> <integer> <integer>
1          chr1         1    804677    805177
2          chr1         2    811866    812366
3          chr1         3    817090    817590
4          chr1         4    818831    819331
5          chr1         5    826543    827043
...         ...       ...       ...       ...
210717     chrX      5777 155898525 155899025
210718     chrX      5778 155899804 155900304
210719     chrX      5779 155956608 155957108
210720     chrX      5780 155966827 155967327
210721     chrX      5781 155971625 155972125

2023-11-14 07:00:19.908096 : Getting PeakMatrix from ArrowFile : SRR13927736.arrow, 0.006 mins elapsed.

2023-11-14 07:00:29.731047 : mat SRR13927736, Class = dgCMatrix
mat SRR13927736: nRows = 210721, nCols = 3534
mat SRR13927736: NonZeroEntries = 25583439, EntryRange = [ 1 , 4 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927736#GCTGTTCTCCTTGACC-1 SRR13927736#AGTTACGAGTAGGTCG-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              .                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927736#AGTTTGGAGCATTCCA-1 SRR13927736#GCACGGTTCGGCAATT-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              1                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927736#TCAGCTCTCTCTATTG-1
[1,]                              .
[2,]                              .
[3,]                              .
[4,]                              .
[5,]                              .


2023-11-14 07:00:29.938278 : Organizing SE PeakMatrix from ArrowFile : SRR13927736.arrow, 0.173 mins elapsed.

2023-11-14 07:00:30.214436 : se SRR13927736, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 3534 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(3534): SRR13927736#GCTGTTCTCCTTGACC-1
  SRR13927736#AGTTACGAGTAGGTCG-1 ... SRR13927736#TAAACCGCAAGTCTGT-1
  SRR13927736#TGCATTTGTTCTACCC-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-14 07:00:30.221404 : Completed PeakMatrix : SRR13927736(4 of 4), 0.793 mins elapsed.
2023-11-14 07:00:30.222328 : Organizing colData, 0.793 mins elapsed.
2023-11-14 07:00:30.359165 : Organizing rowData, 0.795 mins elapsed.
2023-11-14 07:00:30.367418 : Organizing rowRanges, 0.795 mins elapsed.
2023-11-14 07:00:30.377984 : Organizing Assays (1 of 1), 0.796 mins elapsed.
2023-11-14 07:00:33.563509 : Constructing SummarizedExperiment, 0.849 mins elapsed.
2023-11-14 07:00:46.585763 : Finished Matrix Creation, 1.066 mins elapsed.

2023-11-14 08:19:05.970704 : getMatrixFromProject Input-Parameters, Class = list

getMatrixFromProject Input-Parameters$ArchRProj: length = 1

getMatrixFromProject Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromProject Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromProject Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromProject Input-Parameters$verbose: length = 1
[1] TRUE


getMatrixFromProject Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromProject Input-Parameters$threads: length = 1
[1] 1


getMatrixFromProject Input-Parameters$logFile: length = 1
[1] "x"


2023-11-14 08:19:10.016611 : Reading GeneIntegrationMatrix : SRR13927737(1 of 4), 0.07 mins elapsed.

2023-11-14 08:19:10.047739 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927737 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927737.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 4439
[1] "SRR13927737#GTCACGGGTCAACAGG-1" "SRR13927737#TGCCTGTTCTTACTCA-1"
[3] "SRR13927737#TTGCGGGTCAGGTCTA-1" "SRR13927737#GGTTGCGCAAAGCTTC-1"
[5] "SRR13927737#CTTAATCTCGCGCTGA-1" "SRR13927737#TGTACGATCGGAAAGC-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-14 08:19:10.462971 : featureDF SRR13927737, Class = DFrame
DataFrame with 23525 rows and 6 columns
      seqnames     start       end    strand        name       idx
         <Rle> <integer> <integer> <integer> <character> <integer>
1         chr1     69091     70008         1       OR4F5         1
2         chr1    817371    819837         1      FAM87B         2
3         chr1    825138    859446         1   LINC01128         3
4         chr1    827522    826206         2   LINC00115         4
5         chr1    876903    868071         2      FAM41C         5
...        ...       ...       ...       ...         ...       ...
23521     chrX 155264589 155258241         2      RAB39B       923
23522     chrX 155334657 155276211         2       CLIC2       924
23523     chrX 155380787 155381134         1      H2AFB1       925
23524     chrX 155382115 155383230         1        F8A2       926
23525     chrX 155669944 155490115         2       TMLHE       927

2023-11-14 08:19:10.485398 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927737.arrow, 0.007 mins elapsed.

2023-11-14 08:19:31.321422 : mat SRR13927737, Class = dgCMatrix
mat SRR13927737: nRows = 23525, nCols = 4439
mat SRR13927737: NonZeroEntries = 50277774, EntryRange = [ 0.000148130054904683 , 376.345870564644 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927737#GTCACGGGTCAACAGG-1 SRR13927737#TGCCTGTTCTTACTCA-1
[1,]                      .                             .         
[2,]                      .                             .         
[3,]                      0.1057576                     0.09630696
[4,]                      .                             0.01806508
[5,]                      .                             .         
     SRR13927737#TTGCGGGTCAGGTCTA-1 SRR13927737#GGTTGCGCAAAGCTTC-1
[1,]                     .                              .         
[2,]                     .                              .         
[3,]                     0.09113499                     0.08396973
[4,]                     0.02096304                     0.01727603
[5,]                     .                              .         
     SRR13927737#CTTAATCTCGCGCTGA-1
[1,]                     .         
[2,]                     .         
[3,]                     0.08327649
[4,]                     0.01368917
[5,]                     .         


2023-11-14 08:19:31.799916 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927737.arrow, 0.363 mins elapsed.

2023-11-14 08:19:32.668782 : se SRR13927737, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 4439 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(4439): SRR13927737#GTCACGGGTCAACAGG-1
  SRR13927737#TGCCTGTTCTTACTCA-1 ... SRR13927737#TTGCTATAGATGTTCC-1
  SRR13927737#GCTCACTAGTTAGCGG-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-14 08:19:32.68269 : Completed GeneIntegrationMatrix : SRR13927737(1 of 4), 0.447 mins elapsed.
2023-11-14 08:19:32.683653 : Reading GeneIntegrationMatrix : SRR13927738(2 of 4), 0.447 mins elapsed.

2023-11-14 08:19:32.709244 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927738 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927738.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3736
[1] "SRR13927738#TAGCATGAGGCTAAAT-1" "SRR13927738#TGTTAGGTCGGTCAGC-1"
[3] "SRR13927738#CTCTCAGCATGTCCCT-1" "SRR13927738#ACTATTCGTCCAACCG-1"
[5] "SRR13927738#TGGTCCTTCTACTGCC-1" "SRR13927738#ATTGTCTAGACTAGGC-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-14 08:19:33.022253 : featureDF SRR13927738, Class = DFrame
DataFrame with 23525 rows and 6 columns
      seqnames     start       end    strand        name       idx
         <Rle> <integer> <integer> <integer> <character> <integer>
1         chr1     69091     70008         1       OR4F5         1
2         chr1    817371    819837         1      FAM87B         2
3         chr1    825138    859446         1   LINC01128         3
4         chr1    827522    826206         2   LINC00115         4
5         chr1    876903    868071         2      FAM41C         5
...        ...       ...       ...       ...         ...       ...
23521     chrX 155264589 155258241         2      RAB39B       923
23522     chrX 155334657 155276211         2       CLIC2       924
23523     chrX 155380787 155381134         1      H2AFB1       925
23524     chrX 155382115 155383230         1        F8A2       926
23525     chrX 155669944 155490115         2       TMLHE       927

2023-11-14 08:19:33.046858 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927738.arrow, 0.006 mins elapsed.

2023-11-14 08:19:52.6369 : mat SRR13927738, Class = dgCMatrix
mat SRR13927738: nRows = 23525, nCols = 3736
mat SRR13927738: NonZeroEntries = 41898897, EntryRange = [ 9.73975135001421e-05 , 758.763031303812 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927738#TAGCATGAGGCTAAAT-1 SRR13927738#TGTTAGGTCGGTCAGC-1
[1,]                    .                                .        
[2,]                    .                                .        
[3,]                    0.002346352                      0.0084575
[4,]                    0.005836413                      .        
[5,]                    .                                .        
     SRR13927738#CTCTCAGCATGTCCCT-1 SRR13927738#ACTATTCGTCCAACCG-1
[1,]                     .                             .          
[2,]                     .                             .          
[3,]                     .                             0.004081077
[4,]                     0.02001583                    .          
[5,]                     .                             .          
     SRR13927738#TGGTCCTTCTACTGCC-1
[1,]                    .          
[2,]                    .          
[3,]                    0.007640877
[4,]                    0.015048796
[5,]                    .          


2023-11-14 08:19:52.970758 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927738.arrow, 0.338 mins elapsed.

2023-11-14 08:19:53.227532 : se SRR13927738, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 3736 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(3736): SRR13927738#TAGCATGAGGCTAAAT-1
  SRR13927738#TGTTAGGTCGGTCAGC-1 ... SRR13927738#GCTCCTATCCGTGCGA-1
  SRR13927738#GGCACGTAGTCGATAA-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-14 08:19:53.23263 : Completed GeneIntegrationMatrix : SRR13927738(2 of 4), 0.79 mins elapsed.
2023-11-14 08:19:53.233045 : Reading GeneIntegrationMatrix : SRR13927735(3 of 4), 0.79 mins elapsed.

2023-11-14 08:19:53.245229 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927735 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927735.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3813
[1] "SRR13927735#TTATGTCTCCAGGTAT-1" "SRR13927735#TATTGCTCATCAGAAA-1"
[3] "SRR13927735#TTCGATTGTAGGGTTG-1" "SRR13927735#CATTCATTCGGATGTT-1"
[5] "SRR13927735#ACGTTAGGTCAACTGT-1" "SRR13927735#AAATGCCCAGCAATGG-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-14 08:19:53.430932 : featureDF SRR13927735, Class = DFrame
DataFrame with 23525 rows and 6 columns
      seqnames     start       end    strand        name       idx
         <Rle> <integer> <integer> <integer> <character> <integer>
1         chr1     69091     70008         1       OR4F5         1
2         chr1    817371    819837         1      FAM87B         2
3         chr1    825138    859446         1   LINC01128         3
4         chr1    827522    826206         2   LINC00115         4
5         chr1    876903    868071         2      FAM41C         5
...        ...       ...       ...       ...         ...       ...
23521     chrX 155264589 155258241         2      RAB39B       923
23522     chrX 155334657 155276211         2       CLIC2       924
23523     chrX 155380787 155381134         1      H2AFB1       925
23524     chrX 155382115 155383230         1        F8A2       926
23525     chrX 155669944 155490115         2       TMLHE       927

2023-11-14 08:19:53.44502 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927735.arrow, 0.003 mins elapsed.

2023-11-14 08:20:13.49085 : mat SRR13927735, Class = dgCMatrix
mat SRR13927735: nRows = 23525, nCols = 3813
mat SRR13927735: NonZeroEntries = 45535634, EntryRange = [ 7.20715653308254e-05 , 673.990433360078 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927735#TTATGTCTCCAGGTAT-1 SRR13927735#TATTGCTCATCAGAAA-1
[1,]                    .                              .          
[2,]                    .                              .          
[3,]                    0.070892732                    0.057009268
[4,]                    0.008800039                    0.002978408
[5,]                    .                              .          
     SRR13927735#TTCGATTGTAGGGTTG-1 SRR13927735#CATTCATTCGGATGTT-1
[1,]                     .                               .        
[2,]                     .                               .        
[3,]                     0.05523011                      0.0248301
[4,]                     .                               .        
[5,]                     .                               .        
     SRR13927735#ACGTTAGGTCAACTGT-1
[1,]                     .         
[2,]                     .         
[3,]                     0.05384571
[4,]                     0.01228688
[5,]                     .         


2023-11-14 08:20:13.795937 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927735.arrow, 0.343 mins elapsed.

2023-11-14 08:20:14.215464 : se SRR13927735, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 3813 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(3813): SRR13927735#TTATGTCTCCAGGTAT-1
  SRR13927735#TATTGCTCATCAGAAA-1 ... SRR13927735#GATCATGGTGTTGTTG-1
  SRR13927735#GTCCATCCATTTCACT-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-14 08:20:14.221046 : Completed GeneIntegrationMatrix : SRR13927735(3 of 4), 1.14 mins elapsed.
2023-11-14 08:20:14.221506 : Reading GeneIntegrationMatrix : SRR13927736(4 of 4), 1.14 mins elapsed.

2023-11-14 08:20:14.234202 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927736 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927736.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3534
[1] "SRR13927736#GCTGTTCTCCTTGACC-1" "SRR13927736#AGTTACGAGTAGGTCG-1"
[3] "SRR13927736#AGTTTGGAGCATTCCA-1" "SRR13927736#GCACGGTTCGGCAATT-1"
[5] "SRR13927736#TCAGCTCTCTCTATTG-1" "SRR13927736#CCTTAATTCTACTTTG-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-14 08:20:14.448937 : featureDF SRR13927736, Class = DFrame
DataFrame with 23525 rows and 6 columns
      seqnames     start       end    strand        name       idx
         <Rle> <integer> <integer> <integer> <character> <integer>
1         chr1     69091     70008         1       OR4F5         1
2         chr1    817371    819837         1      FAM87B         2
3         chr1    825138    859446         1   LINC01128         3
4         chr1    827522    826206         2   LINC00115         4
5         chr1    876903    868071         2      FAM41C         5
...        ...       ...       ...       ...         ...       ...
23521     chrX 155264589 155258241         2      RAB39B       923
23522     chrX 155334657 155276211         2       CLIC2       924
23523     chrX 155380787 155381134         1      H2AFB1       925
23524     chrX 155382115 155383230         1        F8A2       926
23525     chrX 155669944 155490115         2       TMLHE       927

2023-11-14 08:20:14.467012 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927736.arrow, 0.004 mins elapsed.

2023-11-14 08:20:31.891275 : mat SRR13927736, Class = dgCMatrix
mat SRR13927736: nRows = 23525, nCols = 3534
mat SRR13927736: NonZeroEntries = 39003616, EntryRange = [ 8.77627031078293e-05 , 770.387276519725 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927736#GCTGTTCTCCTTGACC-1 SRR13927736#AGTTACGAGTAGGTCG-1
[1,]                     .                             .          
[2,]                     .                             .          
[3,]                     0.07845381                    0.028344539
[4,]                     0.01395289                    0.008499189
[5,]                     .                             .          
     SRR13927736#AGTTTGGAGCATTCCA-1 SRR13927736#GCACGGTTCGGCAATT-1
[1,]                     .                              .         
[2,]                     .                              .         
[3,]                     0.06535130                     0.06576222
[4,]                     0.01810177                     0.01315728
[5,]                     .                              .         
     SRR13927736#TCAGCTCTCTCTATTG-1
[1,]                    .          
[2,]                    .          
[3,]                    0.016229526
[4,]                    0.003862592
[5,]                    .          


2023-11-14 08:20:32.208107 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927736.arrow, 0.3 mins elapsed.

2023-11-14 08:20:32.478863 : se SRR13927736, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 3534 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(3534): SRR13927736#GCTGTTCTCCTTGACC-1
  SRR13927736#AGTTACGAGTAGGTCG-1 ... SRR13927736#TAAACCGCAAGTCTGT-1
  SRR13927736#TGCATTTGTTCTACCC-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-14 08:20:32.483757 : Completed GeneIntegrationMatrix : SRR13927736(4 of 4), 1.444 mins elapsed.
2023-11-14 08:20:32.484817 : Organizing colData, 1.444 mins elapsed.
2023-11-14 08:20:32.62989 : Organizing rowData, 1.447 mins elapsed.
2023-11-14 08:20:32.633138 : Organizing rowRanges, 1.447 mins elapsed.
2023-11-14 08:20:32.638663 : Organizing Assays (1 of 1), 1.447 mins elapsed.
2023-11-14 08:20:41.927447 : Constructing SummarizedExperiment, 1.602 mins elapsed.
2023-11-14 08:20:44.655792 : Finished Matrix Creation, 1.647 mins elapsed.

2023-11-14 08:20:44.668539 : getMatrixFromProject Input-Parameters, Class = list

getMatrixFromProject Input-Parameters$ArchRProj: length = 1

getMatrixFromProject Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromProject Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromProject Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromProject Input-Parameters$verbose: length = 1
[1] TRUE


getMatrixFromProject Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromProject Input-Parameters$threads: length = 1
[1] 1


getMatrixFromProject Input-Parameters$logFile: length = 1
[1] "x"


2023-11-14 08:20:46.811227 : Reading PeakMatrix : SRR13927737(1 of 4), 0.036 mins elapsed.

2023-11-14 08:20:46.839305 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927737 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927737.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 4439
[1] "SRR13927737#GTCACGGGTCAACAGG-1" "SRR13927737#TGCCTGTTCTTACTCA-1"
[3] "SRR13927737#TTGCGGGTCAGGTCTA-1" "SRR13927737#GGTTGCGCAAAGCTTC-1"
[5] "SRR13927737#CTTAATCTCGCGCTGA-1" "SRR13927737#TGTACGATCGGAAAGC-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-14 08:20:47.370255 : featureDF SRR13927737, Class = DFrame
DataFrame with 210721 rows and 4 columns
       seqnames       idx     start       end
          <Rle> <integer> <integer> <integer>
1          chr1         1    804677    805177
2          chr1         2    811866    812366
3          chr1         3    817090    817590
4          chr1         4    818831    819331
5          chr1         5    826543    827043
...         ...       ...       ...       ...
210717     chrX      5777 155898525 155899025
210718     chrX      5778 155899804 155900304
210719     chrX      5779 155956608 155957108
210720     chrX      5780 155966827 155967327
210721     chrX      5781 155971625 155972125

2023-11-14 08:20:47.442056 : Getting PeakMatrix from ArrowFile : SRR13927737.arrow, 0.01 mins elapsed.

2023-11-14 08:21:00.364377 : mat SRR13927737, Class = dgCMatrix
mat SRR13927737: nRows = 210721, nCols = 4439
mat SRR13927737: NonZeroEntries = 29963541, EntryRange = [ 1 , 4 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927737#GTCACGGGTCAACAGG-1 SRR13927737#TGCCTGTTCTTACTCA-1
[1,]                              .                              .
[2,]                              .                              2
[3,]                              .                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927737#TTGCGGGTCAGGTCTA-1 SRR13927737#GGTTGCGCAAAGCTTC-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              2                              .
[4,]                              .                              .
[5,]                              4                              .
     SRR13927737#CTTAATCTCGCGCTGA-1
[1,]                              .
[2,]                              .
[3,]                              .
[4,]                              .
[5,]                              .


2023-11-14 08:21:00.605284 : Organizing SE PeakMatrix from ArrowFile : SRR13927737.arrow, 0.229 mins elapsed.

2023-11-14 08:21:00.905508 : se SRR13927737, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 4439 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(4439): SRR13927737#GTCACGGGTCAACAGG-1
  SRR13927737#TGCCTGTTCTTACTCA-1 ... SRR13927737#TTGCTATAGATGTTCC-1
  SRR13927737#GCTCACTAGTTAGCGG-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-14 08:21:00.938643 : Completed PeakMatrix : SRR13927737(1 of 4), 0.271 mins elapsed.
2023-11-14 08:21:00.939086 : Reading PeakMatrix : SRR13927738(2 of 4), 0.271 mins elapsed.

2023-11-14 08:21:00.95168 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927738 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927738.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3736
[1] "SRR13927738#TAGCATGAGGCTAAAT-1" "SRR13927738#TGTTAGGTCGGTCAGC-1"
[3] "SRR13927738#CTCTCAGCATGTCCCT-1" "SRR13927738#ACTATTCGTCCAACCG-1"
[5] "SRR13927738#TGGTCCTTCTACTGCC-1" "SRR13927738#ATTGTCTAGACTAGGC-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-14 08:21:01.278755 : featureDF SRR13927738, Class = DFrame
DataFrame with 210721 rows and 4 columns
       seqnames       idx     start       end
          <Rle> <integer> <integer> <integer>
1          chr1         1    804677    805177
2          chr1         2    811866    812366
3          chr1         3    817090    817590
4          chr1         4    818831    819331
5          chr1         5    826543    827043
...         ...       ...       ...       ...
210717     chrX      5777 155898525 155899025
210718     chrX      5778 155899804 155900304
210719     chrX      5779 155956608 155957108
210720     chrX      5780 155966827 155967327
210721     chrX      5781 155971625 155972125

2023-11-14 08:21:01.319089 : Getting PeakMatrix from ArrowFile : SRR13927738.arrow, 0.006 mins elapsed.

2023-11-14 08:21:25.02217 : mat SRR13927738, Class = dgCMatrix
mat SRR13927738: nRows = 210721, nCols = 3736
mat SRR13927738: NonZeroEntries = 30165921, EntryRange = [ 1 , 4 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927738#TAGCATGAGGCTAAAT-1 SRR13927738#TGTTAGGTCGGTCAGC-1
[1,]                              2                              .
[2,]                              .                              .
[3,]                              .                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927738#CTCTCAGCATGTCCCT-1 SRR13927738#ACTATTCGTCCAACCG-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              .                              .
[4,]                              .                              .
[5,]                              2                              .
     SRR13927738#TGGTCCTTCTACTGCC-1
[1,]                              .
[2,]                              .
[3,]                              2
[4,]                              .
[5,]                              .


2023-11-14 08:21:25.267001 : Organizing SE PeakMatrix from ArrowFile : SRR13927738.arrow, 0.405 mins elapsed.

2023-11-14 08:21:25.555339 : se SRR13927738, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 3736 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(3736): SRR13927738#TAGCATGAGGCTAAAT-1
  SRR13927738#TGTTAGGTCGGTCAGC-1 ... SRR13927738#GCTCCTATCCGTGCGA-1
  SRR13927738#GGCACGTAGTCGATAA-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-14 08:21:25.562651 : Completed PeakMatrix : SRR13927738(2 of 4), 0.682 mins elapsed.
2023-11-14 08:21:25.563078 : Reading PeakMatrix : SRR13927735(3 of 4), 0.682 mins elapsed.

2023-11-14 08:21:25.576284 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927735 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927735.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3813
[1] "SRR13927735#TTATGTCTCCAGGTAT-1" "SRR13927735#TATTGCTCATCAGAAA-1"
[3] "SRR13927735#TTCGATTGTAGGGTTG-1" "SRR13927735#CATTCATTCGGATGTT-1"
[5] "SRR13927735#ACGTTAGGTCAACTGT-1" "SRR13927735#AAATGCCCAGCAATGG-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-14 08:21:25.840976 : featureDF SRR13927735, Class = DFrame
DataFrame with 210721 rows and 4 columns
       seqnames       idx     start       end
          <Rle> <integer> <integer> <integer>
1          chr1         1    804677    805177
2          chr1         2    811866    812366
3          chr1         3    817090    817590
4          chr1         4    818831    819331
5          chr1         5    826543    827043
...         ...       ...       ...       ...
210717     chrX      5777 155898525 155899025
210718     chrX      5778 155899804 155900304
210719     chrX      5779 155956608 155957108
210720     chrX      5780 155966827 155967327
210721     chrX      5781 155971625 155972125

2023-11-14 08:21:25.884959 : Getting PeakMatrix from ArrowFile : SRR13927735.arrow, 0.005 mins elapsed.

2023-11-14 08:21:37.088665 : mat SRR13927735, Class = dgCMatrix
mat SRR13927735: nRows = 210721, nCols = 3813
mat SRR13927735: NonZeroEntries = 26038084, EntryRange = [ 1 , 4 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927735#TTATGTCTCCAGGTAT-1 SRR13927735#TATTGCTCATCAGAAA-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              1                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927735#TTCGATTGTAGGGTTG-1 SRR13927735#CATTCATTCGGATGTT-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              .                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927735#ACGTTAGGTCAACTGT-1
[1,]                              2
[2,]                              .
[3,]                              .
[4,]                              .
[5,]                              .


2023-11-14 08:21:37.264688 : Organizing SE PeakMatrix from ArrowFile : SRR13927735.arrow, 0.195 mins elapsed.

2023-11-14 08:21:37.553277 : se SRR13927735, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 3813 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(3813): SRR13927735#TTATGTCTCCAGGTAT-1
  SRR13927735#TATTGCTCATCAGAAA-1 ... SRR13927735#GATCATGGTGTTGTTG-1
  SRR13927735#GTCCATCCATTTCACT-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-14 08:21:37.5606 : Completed PeakMatrix : SRR13927735(3 of 4), 0.882 mins elapsed.
2023-11-14 08:21:37.561035 : Reading PeakMatrix : SRR13927736(4 of 4), 0.882 mins elapsed.

2023-11-14 08:21:37.573837 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927736 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927736.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3534
[1] "SRR13927736#GCTGTTCTCCTTGACC-1" "SRR13927736#AGTTACGAGTAGGTCG-1"
[3] "SRR13927736#AGTTTGGAGCATTCCA-1" "SRR13927736#GCACGGTTCGGCAATT-1"
[5] "SRR13927736#TCAGCTCTCTCTATTG-1" "SRR13927736#CCTTAATTCTACTTTG-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-14 08:21:37.854557 : featureDF SRR13927736, Class = DFrame
DataFrame with 210721 rows and 4 columns
       seqnames       idx     start       end
          <Rle> <integer> <integer> <integer>
1          chr1         1    804677    805177
2          chr1         2    811866    812366
3          chr1         3    817090    817590
4          chr1         4    818831    819331
5          chr1         5    826543    827043
...         ...       ...       ...       ...
210717     chrX      5777 155898525 155899025
210718     chrX      5778 155899804 155900304
210719     chrX      5779 155956608 155957108
210720     chrX      5780 155966827 155967327
210721     chrX      5781 155971625 155972125

2023-11-14 08:21:37.897677 : Getting PeakMatrix from ArrowFile : SRR13927736.arrow, 0.005 mins elapsed.

2023-11-14 08:21:48.07122 : mat SRR13927736, Class = dgCMatrix
mat SRR13927736: nRows = 210721, nCols = 3534
mat SRR13927736: NonZeroEntries = 25583439, EntryRange = [ 1 , 4 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927736#GCTGTTCTCCTTGACC-1 SRR13927736#AGTTACGAGTAGGTCG-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              .                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927736#AGTTTGGAGCATTCCA-1 SRR13927736#GCACGGTTCGGCAATT-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              1                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927736#TCAGCTCTCTCTATTG-1
[1,]                              .
[2,]                              .
[3,]                              .
[4,]                              .
[5,]                              .


2023-11-14 08:21:48.274805 : Organizing SE PeakMatrix from ArrowFile : SRR13927736.arrow, 0.178 mins elapsed.

2023-11-14 08:21:48.548519 : se SRR13927736, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 3534 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(3534): SRR13927736#GCTGTTCTCCTTGACC-1
  SRR13927736#AGTTACGAGTAGGTCG-1 ... SRR13927736#TAAACCGCAAGTCTGT-1
  SRR13927736#TGCATTTGTTCTACCC-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-14 08:21:48.555337 : Completed PeakMatrix : SRR13927736(4 of 4), 1.065 mins elapsed.
2023-11-14 08:21:48.556258 : Organizing colData, 1.065 mins elapsed.
2023-11-14 08:21:48.700895 : Organizing rowData, 1.067 mins elapsed.
2023-11-14 08:21:48.709363 : Organizing rowRanges, 1.067 mins elapsed.
2023-11-14 08:21:48.720807 : Organizing Assays (1 of 1), 1.068 mins elapsed.
2023-11-14 08:21:52.051728 : Constructing SummarizedExperiment, 1.123 mins elapsed.
2023-11-14 08:22:08.945223 : Finished Matrix Creation, 1.405 mins elapsed.

2023-11-15 09:21:44.520218 : getMatrixFromProject Input-Parameters, Class = list

getMatrixFromProject Input-Parameters$ArchRProj: length = 1

getMatrixFromProject Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromProject Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromProject Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromProject Input-Parameters$verbose: length = 1
[1] TRUE


getMatrixFromProject Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromProject Input-Parameters$threads: length = 1
[1] 1


getMatrixFromProject Input-Parameters$logFile: length = 1
[1] "x"


2023-11-15 09:21:45.459273 : Reading GeneIntegrationMatrix : SRR13927737(1 of 4), 0.016 mins elapsed.

2023-11-15 09:21:45.516526 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927737 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927737.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 4439
[1] "SRR13927737#GTCACGGGTCAACAGG-1" "SRR13927737#TGCCTGTTCTTACTCA-1" "SRR13927737#TTGCGGGTCAGGTCTA-1"
[4] "SRR13927737#GGTTGCGCAAAGCTTC-1" "SRR13927737#CTTAATCTCGCGCTGA-1" "SRR13927737#TGTACGATCGGAAAGC-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-15 09:21:48.455108 : featureDF SRR13927737, Class = DFrame
DataFrame with 23525 rows and 6 columns
      seqnames     start       end    strand        name       idx
         <Rle> <integer> <integer> <integer> <character> <integer>
1         chr1     69091     70008         1       OR4F5         1
2         chr1    817371    819837         1      FAM87B         2
3         chr1    825138    859446         1   LINC01128         3
4         chr1    827522    826206         2   LINC00115         4
5         chr1    876903    868071         2      FAM41C         5
...        ...       ...       ...       ...         ...       ...
23521     chrX 155264589 155258241         2      RAB39B       923
23522     chrX 155334657 155276211         2       CLIC2       924
23523     chrX 155380787 155381134         1      H2AFB1       925
23524     chrX 155382115 155383230         1        F8A2       926
23525     chrX 155669944 155490115         2       TMLHE       927

2023-11-15 09:21:48.490438 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927737.arrow, 0.05 mins elapsed.

2023-11-15 09:22:23.562233 : mat SRR13927737, Class = dgCMatrix
mat SRR13927737: nRows = 23525, nCols = 4439
mat SRR13927737: NonZeroEntries = 50277774, EntryRange = [ 0.000148130054904683 , 376.345870564644 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927737#GTCACGGGTCAACAGG-1 SRR13927737#TGCCTGTTCTTACTCA-1 SRR13927737#TTGCGGGTCAGGTCTA-1
[1,]                      .                             .                              .         
[2,]                      .                             .                              .         
[3,]                      0.1057576                     0.09630696                     0.09113499
[4,]                      .                             0.01806508                     0.02096304
[5,]                      .                             .                              .         
     SRR13927737#GGTTGCGCAAAGCTTC-1 SRR13927737#CTTAATCTCGCGCTGA-1
[1,]                     .                              .         
[2,]                     .                              .         
[3,]                     0.08396973                     0.08327649
[4,]                     0.01727603                     0.01368917
[5,]                     .                              .         


2023-11-15 09:22:24.430589 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927737.arrow, 0.649 mins elapsed.

2023-11-15 09:22:25.297926 : se SRR13927737, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 4439 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(4439): SRR13927737#GTCACGGGTCAACAGG-1 SRR13927737#TGCCTGTTCTTACTCA-1 ...
  SRR13927737#TTGCTATAGATGTTCC-1 SRR13927737#GCTCACTAGTTAGCGG-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-15 09:22:25.309865 : Completed GeneIntegrationMatrix : SRR13927737(1 of 4), 0.68 mins elapsed.
2023-11-15 09:22:25.310719 : Reading GeneIntegrationMatrix : SRR13927738(2 of 4), 0.68 mins elapsed.

2023-11-15 09:22:25.333894 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927738 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927738.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3736
[1] "SRR13927738#TAGCATGAGGCTAAAT-1" "SRR13927738#TGTTAGGTCGGTCAGC-1" "SRR13927738#CTCTCAGCATGTCCCT-1"
[4] "SRR13927738#ACTATTCGTCCAACCG-1" "SRR13927738#TGGTCCTTCTACTGCC-1" "SRR13927738#ATTGTCTAGACTAGGC-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-15 09:22:25.838992 : featureDF SRR13927738, Class = DFrame
DataFrame with 23525 rows and 6 columns
      seqnames     start       end    strand        name       idx
         <Rle> <integer> <integer> <integer> <character> <integer>
1         chr1     69091     70008         1       OR4F5         1
2         chr1    817371    819837         1      FAM87B         2
3         chr1    825138    859446         1   LINC01128         3
4         chr1    827522    826206         2   LINC00115         4
5         chr1    876903    868071         2      FAM41C         5
...        ...       ...       ...       ...         ...       ...
23521     chrX 155264589 155258241         2      RAB39B       923
23522     chrX 155334657 155276211         2       CLIC2       924
23523     chrX 155380787 155381134         1      H2AFB1       925
23524     chrX 155382115 155383230         1        F8A2       926
23525     chrX 155669944 155490115         2       TMLHE       927

2023-11-15 09:22:25.875649 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927738.arrow, 0.009 mins elapsed.

2023-11-15 09:22:44.334784 : mat SRR13927738, Class = dgCMatrix
mat SRR13927738: nRows = 23525, nCols = 3736
mat SRR13927738: NonZeroEntries = 41898897, EntryRange = [ 9.73975135001421e-05 , 758.763031303812 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927738#TAGCATGAGGCTAAAT-1 SRR13927738#TGTTAGGTCGGTCAGC-1 SRR13927738#CTCTCAGCATGTCCCT-1
[1,]                    .                                .                             .         
[2,]                    .                                .                             .         
[3,]                    0.002346352                      0.0084575                     .         
[4,]                    0.005836413                      .                             0.02001583
[5,]                    .                                .                             .         
     SRR13927738#ACTATTCGTCCAACCG-1 SRR13927738#TGGTCCTTCTACTGCC-1
[1,]                    .                              .          
[2,]                    .                              .          
[3,]                    0.004081077                    0.007640877
[4,]                    .                              0.015048796
[5,]                    .                              .          


2023-11-15 09:22:45.049611 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927738.arrow, 0.329 mins elapsed.

2023-11-15 09:22:45.972306 : se SRR13927738, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 3736 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(3736): SRR13927738#TAGCATGAGGCTAAAT-1 SRR13927738#TGTTAGGTCGGTCAGC-1 ...
  SRR13927738#GCTCCTATCCGTGCGA-1 SRR13927738#GGCACGTAGTCGATAA-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-15 09:22:45.981096 : Completed GeneIntegrationMatrix : SRR13927738(2 of 4), 1.024 mins elapsed.
2023-11-15 09:22:45.981858 : Reading GeneIntegrationMatrix : SRR13927735(3 of 4), 1.024 mins elapsed.

2023-11-15 09:22:46.001994 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927735 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927735.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3813
[1] "SRR13927735#TTATGTCTCCAGGTAT-1" "SRR13927735#TATTGCTCATCAGAAA-1" "SRR13927735#TTCGATTGTAGGGTTG-1"
[4] "SRR13927735#CATTCATTCGGATGTT-1" "SRR13927735#ACGTTAGGTCAACTGT-1" "SRR13927735#AAATGCCCAGCAATGG-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-15 09:22:46.462347 : featureDF SRR13927735, Class = DFrame
DataFrame with 23525 rows and 6 columns
      seqnames     start       end    strand        name       idx
         <Rle> <integer> <integer> <integer> <character> <integer>
1         chr1     69091     70008         1       OR4F5         1
2         chr1    817371    819837         1      FAM87B         2
3         chr1    825138    859446         1   LINC01128         3
4         chr1    827522    826206         2   LINC00115         4
5         chr1    876903    868071         2      FAM41C         5
...        ...       ...       ...       ...         ...       ...
23521     chrX 155264589 155258241         2      RAB39B       923
23522     chrX 155334657 155276211         2       CLIC2       924
23523     chrX 155380787 155381134         1      H2AFB1       925
23524     chrX 155382115 155383230         1        F8A2       926
23525     chrX 155669944 155490115         2       TMLHE       927

2023-11-15 09:22:46.499866 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927735.arrow, 0.008 mins elapsed.

2023-11-15 09:23:05.128243 : mat SRR13927735, Class = dgCMatrix
mat SRR13927735: nRows = 23525, nCols = 3813
mat SRR13927735: NonZeroEntries = 45535634, EntryRange = [ 7.20715653308254e-05 , 673.990433360078 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927735#TTATGTCTCCAGGTAT-1 SRR13927735#TATTGCTCATCAGAAA-1 SRR13927735#TTCGATTGTAGGGTTG-1
[1,]                    .                              .                               .         
[2,]                    .                              .                               .         
[3,]                    0.070892732                    0.057009268                     0.05523011
[4,]                    0.008800039                    0.002978408                     .         
[5,]                    .                              .                               .         
     SRR13927735#CATTCATTCGGATGTT-1 SRR13927735#ACGTTAGGTCAACTGT-1
[1,]                      .                             .         
[2,]                      .                             .         
[3,]                      0.0248301                     0.05384571
[4,]                      .                             0.01228688
[5,]                      .                             .         


2023-11-15 09:23:05.887511 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927735.arrow, 0.331 mins elapsed.

2023-11-15 09:23:06.679346 : se SRR13927735, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 3813 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(3813): SRR13927735#TTATGTCTCCAGGTAT-1 SRR13927735#TATTGCTCATCAGAAA-1 ...
  SRR13927735#GATCATGGTGTTGTTG-1 SRR13927735#GTCCATCCATTTCACT-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-15 09:23:06.689478 : Completed GeneIntegrationMatrix : SRR13927735(3 of 4), 1.37 mins elapsed.
2023-11-15 09:23:06.690413 : Reading GeneIntegrationMatrix : SRR13927736(4 of 4), 1.37 mins elapsed.

2023-11-15 09:23:06.713163 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927736 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927736.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3534
[1] "SRR13927736#GCTGTTCTCCTTGACC-1" "SRR13927736#AGTTACGAGTAGGTCG-1" "SRR13927736#AGTTTGGAGCATTCCA-1"
[4] "SRR13927736#GCACGGTTCGGCAATT-1" "SRR13927736#TCAGCTCTCTCTATTG-1" "SRR13927736#CCTTAATTCTACTTTG-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-15 09:23:07.281434 : featureDF SRR13927736, Class = DFrame
DataFrame with 23525 rows and 6 columns
      seqnames     start       end    strand        name       idx
         <Rle> <integer> <integer> <integer> <character> <integer>
1         chr1     69091     70008         1       OR4F5         1
2         chr1    817371    819837         1      FAM87B         2
3         chr1    825138    859446         1   LINC01128         3
4         chr1    827522    826206         2   LINC00115         4
5         chr1    876903    868071         2      FAM41C         5
...        ...       ...       ...       ...         ...       ...
23521     chrX 155264589 155258241         2      RAB39B       923
23522     chrX 155334657 155276211         2       CLIC2       924
23523     chrX 155380787 155381134         1      H2AFB1       925
23524     chrX 155382115 155383230         1        F8A2       926
23525     chrX 155669944 155490115         2       TMLHE       927

2023-11-15 09:23:07.31785 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927736.arrow, 0.01 mins elapsed.

2023-11-15 09:23:23.372725 : mat SRR13927736, Class = dgCMatrix
mat SRR13927736: nRows = 23525, nCols = 3534
mat SRR13927736: NonZeroEntries = 39003616, EntryRange = [ 8.77627031078293e-05 , 770.387276519725 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927736#GCTGTTCTCCTTGACC-1 SRR13927736#AGTTACGAGTAGGTCG-1 SRR13927736#AGTTTGGAGCATTCCA-1
[1,]                     .                             .                               .         
[2,]                     .                             .                               .         
[3,]                     0.07845381                    0.028344539                     0.06535130
[4,]                     0.01395289                    0.008499189                     0.01810177
[5,]                     .                             .                               .         
     SRR13927736#GCACGGTTCGGCAATT-1 SRR13927736#TCAGCTCTCTCTATTG-1
[1,]                     .                             .          
[2,]                     .                             .          
[3,]                     0.06576222                    0.016229526
[4,]                     0.01315728                    0.003862592
[5,]                     .                             .          


2023-11-15 09:23:24.065885 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927736.arrow, 0.289 mins elapsed.

2023-11-15 09:23:25.1068 : se SRR13927736, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 3534 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(3534): SRR13927736#GCTGTTCTCCTTGACC-1 SRR13927736#AGTTACGAGTAGGTCG-1 ...
  SRR13927736#TAAACCGCAAGTCTGT-1 SRR13927736#TGCATTTGTTCTACCC-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-15 09:23:25.315057 : Completed GeneIntegrationMatrix : SRR13927736(4 of 4), 1.68 mins elapsed.
2023-11-15 09:23:25.324095 : Organizing colData, 1.68 mins elapsed.
2023-11-15 09:23:25.920098 : Organizing rowData, 1.69 mins elapsed.
2023-11-15 09:23:25.927802 : Organizing rowRanges, 1.69 mins elapsed.
2023-11-15 09:23:25.940156 : Organizing Assays (1 of 1), 1.69 mins elapsed.
2023-11-15 09:23:34.084817 : Constructing SummarizedExperiment, 1.826 mins elapsed.
2023-11-15 09:23:35.796669 : Finished Matrix Creation, 1.853 mins elapsed.

2023-11-15 09:23:35.802217 : getMatrixFromProject Input-Parameters, Class = list

getMatrixFromProject Input-Parameters$ArchRProj: length = 1

getMatrixFromProject Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromProject Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromProject Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromProject Input-Parameters$verbose: length = 1
[1] TRUE


getMatrixFromProject Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromProject Input-Parameters$threads: length = 1
[1] 1


getMatrixFromProject Input-Parameters$logFile: length = 1
[1] "x"


2023-11-15 09:23:36.470466 : Reading PeakMatrix : SRR13927737(1 of 4), 0.011 mins elapsed.

2023-11-15 09:23:36.502595 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927737 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927737.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 4439
[1] "SRR13927737#GTCACGGGTCAACAGG-1" "SRR13927737#TGCCTGTTCTTACTCA-1" "SRR13927737#TTGCGGGTCAGGTCTA-1"
[4] "SRR13927737#GGTTGCGCAAAGCTTC-1" "SRR13927737#CTTAATCTCGCGCTGA-1" "SRR13927737#TGTACGATCGGAAAGC-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-15 09:23:37.108704 : featureDF SRR13927737, Class = DFrame
DataFrame with 210721 rows and 4 columns
       seqnames       idx     start       end
          <Rle> <integer> <integer> <integer>
1          chr1         1    804677    805177
2          chr1         2    811866    812366
3          chr1         3    817090    817590
4          chr1         4    818831    819331
5          chr1         5    826543    827043
...         ...       ...       ...       ...
210717     chrX      5777 155898525 155899025
210718     chrX      5778 155899804 155900304
210719     chrX      5779 155956608 155957108
210720     chrX      5780 155966827 155967327
210721     chrX      5781 155971625 155972125

2023-11-15 09:23:37.354035 : Getting PeakMatrix from ArrowFile : SRR13927737.arrow, 0.014 mins elapsed.

2023-11-15 09:23:49.090359 : mat SRR13927737, Class = dgCMatrix
mat SRR13927737: nRows = 210721, nCols = 4439
mat SRR13927737: NonZeroEntries = 29963541, EntryRange = [ 1 , 4 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927737#GTCACGGGTCAACAGG-1 SRR13927737#TGCCTGTTCTTACTCA-1 SRR13927737#TTGCGGGTCAGGTCTA-1
[1,]                              .                              .                              .
[2,]                              .                              2                              .
[3,]                              .                              .                              2
[4,]                              .                              .                              .
[5,]                              .                              .                              4
     SRR13927737#GGTTGCGCAAAGCTTC-1 SRR13927737#CTTAATCTCGCGCTGA-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              .                              .
[4,]                              .                              .
[5,]                              .                              .


2023-11-15 09:23:49.824112 : Organizing SE PeakMatrix from ArrowFile : SRR13927737.arrow, 0.222 mins elapsed.

2023-11-15 09:23:50.791615 : se SRR13927737, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 4439 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(4439): SRR13927737#GTCACGGGTCAACAGG-1 SRR13927737#TGCCTGTTCTTACTCA-1 ...
  SRR13927737#TTGCTATAGATGTTCC-1 SRR13927737#GCTCACTAGTTAGCGG-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-15 09:23:50.81196 : Completed PeakMatrix : SRR13927737(1 of 4), 0.25 mins elapsed.
2023-11-15 09:23:50.812994 : Reading PeakMatrix : SRR13927738(2 of 4), 0.25 mins elapsed.

2023-11-15 09:23:50.839487 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927738 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927738.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3736
[1] "SRR13927738#TAGCATGAGGCTAAAT-1" "SRR13927738#TGTTAGGTCGGTCAGC-1" "SRR13927738#CTCTCAGCATGTCCCT-1"
[4] "SRR13927738#ACTATTCGTCCAACCG-1" "SRR13927738#TGGTCCTTCTACTGCC-1" "SRR13927738#ATTGTCTAGACTAGGC-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-15 09:23:51.601959 : featureDF SRR13927738, Class = DFrame
DataFrame with 210721 rows and 4 columns
       seqnames       idx     start       end
          <Rle> <integer> <integer> <integer>
1          chr1         1    804677    805177
2          chr1         2    811866    812366
3          chr1         3    817090    817590
4          chr1         4    818831    819331
5          chr1         5    826543    827043
...         ...       ...       ...       ...
210717     chrX      5777 155898525 155899025
210718     chrX      5778 155899804 155900304
210719     chrX      5779 155956608 155957108
210720     chrX      5780 155966827 155967327
210721     chrX      5781 155971625 155972125

2023-11-15 09:23:51.715241 : Getting PeakMatrix from ArrowFile : SRR13927738.arrow, 0.015 mins elapsed.

2023-11-15 09:24:03.460555 : mat SRR13927738, Class = dgCMatrix
mat SRR13927738: nRows = 210721, nCols = 3736
mat SRR13927738: NonZeroEntries = 30165921, EntryRange = [ 1 , 4 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927738#TAGCATGAGGCTAAAT-1 SRR13927738#TGTTAGGTCGGTCAGC-1 SRR13927738#CTCTCAGCATGTCCCT-1
[1,]                              2                              .                              .
[2,]                              .                              .                              .
[3,]                              .                              .                              .
[4,]                              .                              .                              .
[5,]                              .                              .                              2
     SRR13927738#ACTATTCGTCCAACCG-1 SRR13927738#TGGTCCTTCTACTGCC-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              .                              2
[4,]                              .                              .
[5,]                              .                              .


2023-11-15 09:24:04.176016 : Organizing SE PeakMatrix from ArrowFile : SRR13927738.arrow, 0.222 mins elapsed.

2023-11-15 09:24:05.145676 : se SRR13927738, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 3736 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(3736): SRR13927738#TAGCATGAGGCTAAAT-1 SRR13927738#TGTTAGGTCGGTCAGC-1 ...
  SRR13927738#GCTCCTATCCGTGCGA-1 SRR13927738#GGCACGTAGTCGATAA-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-15 09:24:05.163945 : Completed PeakMatrix : SRR13927738(2 of 4), 0.489 mins elapsed.
2023-11-15 09:24:05.1652 : Reading PeakMatrix : SRR13927735(3 of 4), 0.489 mins elapsed.

2023-11-15 09:24:05.196896 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927735 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927735.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3813
[1] "SRR13927735#TTATGTCTCCAGGTAT-1" "SRR13927735#TATTGCTCATCAGAAA-1" "SRR13927735#TTCGATTGTAGGGTTG-1"
[4] "SRR13927735#CATTCATTCGGATGTT-1" "SRR13927735#ACGTTAGGTCAACTGT-1" "SRR13927735#AAATGCCCAGCAATGG-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-15 09:24:05.88683 : featureDF SRR13927735, Class = DFrame
DataFrame with 210721 rows and 4 columns
       seqnames       idx     start       end
          <Rle> <integer> <integer> <integer>
1          chr1         1    804677    805177
2          chr1         2    811866    812366
3          chr1         3    817090    817590
4          chr1         4    818831    819331
5          chr1         5    826543    827043
...         ...       ...       ...       ...
210717     chrX      5777 155898525 155899025
210718     chrX      5778 155899804 155900304
210719     chrX      5779 155956608 155957108
210720     chrX      5780 155966827 155967327
210721     chrX      5781 155971625 155972125

2023-11-15 09:24:05.952454 : Getting PeakMatrix from ArrowFile : SRR13927735.arrow, 0.013 mins elapsed.

2023-11-15 09:24:16.355906 : mat SRR13927735, Class = dgCMatrix
mat SRR13927735: nRows = 210721, nCols = 3813
mat SRR13927735: NonZeroEntries = 26038084, EntryRange = [ 1 , 4 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927735#TTATGTCTCCAGGTAT-1 SRR13927735#TATTGCTCATCAGAAA-1 SRR13927735#TTCGATTGTAGGGTTG-1
[1,]                              .                              .                              .
[2,]                              .                              .                              .
[3,]                              1                              .                              .
[4,]                              .                              .                              .
[5,]                              .                              .                              .
     SRR13927735#CATTCATTCGGATGTT-1 SRR13927735#ACGTTAGGTCAACTGT-1
[1,]                              .                              2
[2,]                              .                              .
[3,]                              .                              .
[4,]                              .                              .
[5,]                              .                              .


2023-11-15 09:24:16.519599 : Organizing SE PeakMatrix from ArrowFile : SRR13927735.arrow, 0.189 mins elapsed.

2023-11-15 09:24:17.298327 : se SRR13927735, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 3813 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(3813): SRR13927735#TTATGTCTCCAGGTAT-1 SRR13927735#TATTGCTCATCAGAAA-1 ...
  SRR13927735#GATCATGGTGTTGTTG-1 SRR13927735#GTCCATCCATTTCACT-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-15 09:24:17.317067 : Completed PeakMatrix : SRR13927735(3 of 4), 0.692 mins elapsed.
2023-11-15 09:24:17.318301 : Reading PeakMatrix : SRR13927736(4 of 4), 0.692 mins elapsed.

2023-11-15 09:24:17.350024 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927736 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927736.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3534
[1] "SRR13927736#GCTGTTCTCCTTGACC-1" "SRR13927736#AGTTACGAGTAGGTCG-1" "SRR13927736#AGTTTGGAGCATTCCA-1"
[4] "SRR13927736#GCACGGTTCGGCAATT-1" "SRR13927736#TCAGCTCTCTCTATTG-1" "SRR13927736#CCTTAATTCTACTTTG-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-15 09:24:18.133683 : featureDF SRR13927736, Class = DFrame
DataFrame with 210721 rows and 4 columns
       seqnames       idx     start       end
          <Rle> <integer> <integer> <integer>
1          chr1         1    804677    805177
2          chr1         2    811866    812366
3          chr1         3    817090    817590
4          chr1         4    818831    819331
5          chr1         5    826543    827043
...         ...       ...       ...       ...
210717     chrX      5777 155898525 155899025
210718     chrX      5778 155899804 155900304
210719     chrX      5779 155956608 155957108
210720     chrX      5780 155966827 155967327
210721     chrX      5781 155971625 155972125

2023-11-15 09:24:18.24724 : Getting PeakMatrix from ArrowFile : SRR13927736.arrow, 0.015 mins elapsed.

2023-11-15 09:24:27.80719 : mat SRR13927736, Class = dgCMatrix
mat SRR13927736: nRows = 210721, nCols = 3534
mat SRR13927736: NonZeroEntries = 25583439, EntryRange = [ 1 , 4 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927736#GCTGTTCTCCTTGACC-1 SRR13927736#AGTTACGAGTAGGTCG-1 SRR13927736#AGTTTGGAGCATTCCA-1
[1,]                              .                              .                              .
[2,]                              .                              .                              .
[3,]                              .                              .                              1
[4,]                              .                              .                              .
[5,]                              .                              .                              .
     SRR13927736#GCACGGTTCGGCAATT-1 SRR13927736#TCAGCTCTCTCTATTG-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              .                              .
[4,]                              .                              .
[5,]                              .                              .


2023-11-15 09:24:28.008235 : Organizing SE PeakMatrix from ArrowFile : SRR13927736.arrow, 0.178 mins elapsed.

2023-11-15 09:24:28.794451 : se SRR13927736, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 3534 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(3534): SRR13927736#GCTGTTCTCCTTGACC-1 SRR13927736#AGTTACGAGTAGGTCG-1 ...
  SRR13927736#TAAACCGCAAGTCTGT-1 SRR13927736#TGCATTTGTTCTACCC-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-15 09:24:28.813183 : Completed PeakMatrix : SRR13927736(4 of 4), 0.884 mins elapsed.
2023-11-15 09:24:28.815712 : Organizing colData, 0.884 mins elapsed.
2023-11-15 09:24:29.360011 : Organizing rowData, 0.893 mins elapsed.
2023-11-15 09:24:29.382219 : Organizing rowRanges, 0.893 mins elapsed.
2023-11-15 09:24:29.410318 : Organizing Assays (1 of 1), 0.893 mins elapsed.
2023-11-15 09:24:31.800544 : Constructing SummarizedExperiment, 0.933 mins elapsed.
2023-11-15 09:24:48.627058 : Finished Matrix Creation, 1.214 mins elapsed.

2023-11-16 03:47:39.593726 : getMatrixFromProject Input-Parameters, Class = list

getMatrixFromProject Input-Parameters$ArchRProj: length = 1

getMatrixFromProject Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromProject Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromProject Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromProject Input-Parameters$verbose: length = 1
[1] TRUE


getMatrixFromProject Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromProject Input-Parameters$threads: length = 1
[1] 1


getMatrixFromProject Input-Parameters$logFile: length = 1
[1] "x"


2023-11-16 03:47:44.355768 : Reading GeneIntegrationMatrix : SRR13927737(1 of 4), 0.079 mins elapsed.

2023-11-16 03:47:44.403947 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927737 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927737.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 4439
[1] "SRR13927737#GTCACGGGTCAACAGG-1" "SRR13927737#TGCCTGTTCTTACTCA-1"
[3] "SRR13927737#TTGCGGGTCAGGTCTA-1" "SRR13927737#GGTTGCGCAAAGCTTC-1"
[5] "SRR13927737#CTTAATCTCGCGCTGA-1" "SRR13927737#TGTACGATCGGAAAGC-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-16 03:47:45.495101 : featureDF SRR13927737, Class = DFrame
DataFrame with 23525 rows and 6 columns
      seqnames     start       end    strand        name       idx
         <Rle> <integer> <integer> <integer> <character> <integer>
1         chr1     69091     70008         1       OR4F5         1
2         chr1    817371    819837         1      FAM87B         2
3         chr1    825138    859446         1   LINC01128         3
4         chr1    827522    826206         2   LINC00115         4
5         chr1    876903    868071         2      FAM41C         5
...        ...       ...       ...       ...         ...       ...
23521     chrX 155264589 155258241         2      RAB39B       923
23522     chrX 155334657 155276211         2       CLIC2       924
23523     chrX 155380787 155381134         1      H2AFB1       925
23524     chrX 155382115 155383230         1        F8A2       926
23525     chrX 155669944 155490115         2       TMLHE       927

2023-11-16 03:47:45.526563 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927737.arrow, 0.019 mins elapsed.

2023-11-16 03:48:07.643353 : mat SRR13927737, Class = dgCMatrix
mat SRR13927737: nRows = 23525, nCols = 4439
mat SRR13927737: NonZeroEntries = 50277774, EntryRange = [ 0.000148130054904683 , 376.345870564644 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927737#GTCACGGGTCAACAGG-1 SRR13927737#TGCCTGTTCTTACTCA-1
[1,]                      .                             .         
[2,]                      .                             .         
[3,]                      0.1057576                     0.09630696
[4,]                      .                             0.01806508
[5,]                      .                             .         
     SRR13927737#TTGCGGGTCAGGTCTA-1 SRR13927737#GGTTGCGCAAAGCTTC-1
[1,]                     .                              .         
[2,]                     .                              .         
[3,]                     0.09113499                     0.08396973
[4,]                     0.02096304                     0.01727603
[5,]                     .                              .         
     SRR13927737#CTTAATCTCGCGCTGA-1
[1,]                     .         
[2,]                     .         
[3,]                     0.08327649
[4,]                     0.01368917
[5,]                     .         


2023-11-16 03:48:08.115118 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927737.arrow, 0.395 mins elapsed.

2023-11-16 03:48:08.450639 : se SRR13927737, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 4439 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(4439): SRR13927737#GTCACGGGTCAACAGG-1
  SRR13927737#TGCCTGTTCTTACTCA-1 ... SRR13927737#TTGCTATAGATGTTCC-1
  SRR13927737#GCTCACTAGTTAGCGG-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-16 03:48:08.456956 : Completed GeneIntegrationMatrix : SRR13927737(1 of 4), 0.481 mins elapsed.
2023-11-16 03:48:08.457374 : Reading GeneIntegrationMatrix : SRR13927738(2 of 4), 0.481 mins elapsed.

2023-11-16 03:48:08.481437 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927738 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927738.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3736
[1] "SRR13927738#TAGCATGAGGCTAAAT-1" "SRR13927738#TGTTAGGTCGGTCAGC-1"
[3] "SRR13927738#CTCTCAGCATGTCCCT-1" "SRR13927738#ACTATTCGTCCAACCG-1"
[5] "SRR13927738#TGGTCCTTCTACTGCC-1" "SRR13927738#ATTGTCTAGACTAGGC-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-16 03:48:09.495113 : featureDF SRR13927738, Class = DFrame
DataFrame with 23525 rows and 6 columns
      seqnames     start       end    strand        name       idx
         <Rle> <integer> <integer> <integer> <character> <integer>
1         chr1     69091     70008         1       OR4F5         1
2         chr1    817371    819837         1      FAM87B         2
3         chr1    825138    859446         1   LINC01128         3
4         chr1    827522    826206         2   LINC00115         4
5         chr1    876903    868071         2      FAM41C         5
...        ...       ...       ...       ...         ...       ...
23521     chrX 155264589 155258241         2      RAB39B       923
23522     chrX 155334657 155276211         2       CLIC2       924
23523     chrX 155380787 155381134         1      H2AFB1       925
23524     chrX 155382115 155383230         1        F8A2       926
23525     chrX 155669944 155490115         2       TMLHE       927

2023-11-16 03:48:09.531124 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927738.arrow, 0.017 mins elapsed.

2023-11-16 03:48:27.665817 : mat SRR13927738, Class = dgCMatrix
mat SRR13927738: nRows = 23525, nCols = 3736
mat SRR13927738: NonZeroEntries = 41898897, EntryRange = [ 9.73975135001421e-05 , 758.763031303812 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927738#TAGCATGAGGCTAAAT-1 SRR13927738#TGTTAGGTCGGTCAGC-1
[1,]                    .                                .        
[2,]                    .                                .        
[3,]                    0.002346352                      0.0084575
[4,]                    0.005836413                      .        
[5,]                    .                                .        
     SRR13927738#CTCTCAGCATGTCCCT-1 SRR13927738#ACTATTCGTCCAACCG-1
[1,]                     .                             .          
[2,]                     .                             .          
[3,]                     .                             0.004081077
[4,]                     0.02001583                    .          
[5,]                     .                             .          
     SRR13927738#TGGTCCTTCTACTGCC-1
[1,]                    .          
[2,]                    .          
[3,]                    0.007640877
[4,]                    0.015048796
[5,]                    .          


2023-11-16 03:48:28.005979 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927738.arrow, 0.325 mins elapsed.

2023-11-16 03:48:28.309291 : se SRR13927738, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 3736 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(3736): SRR13927738#TAGCATGAGGCTAAAT-1
  SRR13927738#TGTTAGGTCGGTCAGC-1 ... SRR13927738#GCTCCTATCCGTGCGA-1
  SRR13927738#GGCACGTAGTCGATAA-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-16 03:48:28.314307 : Completed GeneIntegrationMatrix : SRR13927738(2 of 4), 0.812 mins elapsed.
2023-11-16 03:48:28.314714 : Reading GeneIntegrationMatrix : SRR13927735(3 of 4), 0.812 mins elapsed.

2023-11-16 03:48:28.337926 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927735 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927735.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3813
[1] "SRR13927735#TTATGTCTCCAGGTAT-1" "SRR13927735#TATTGCTCATCAGAAA-1"
[3] "SRR13927735#TTCGATTGTAGGGTTG-1" "SRR13927735#CATTCATTCGGATGTT-1"
[5] "SRR13927735#ACGTTAGGTCAACTGT-1" "SRR13927735#AAATGCCCAGCAATGG-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-16 03:48:29.316366 : featureDF SRR13927735, Class = DFrame
DataFrame with 23525 rows and 6 columns
      seqnames     start       end    strand        name       idx
         <Rle> <integer> <integer> <integer> <character> <integer>
1         chr1     69091     70008         1       OR4F5         1
2         chr1    817371    819837         1      FAM87B         2
3         chr1    825138    859446         1   LINC01128         3
4         chr1    827522    826206         2   LINC00115         4
5         chr1    876903    868071         2      FAM41C         5
...        ...       ...       ...       ...         ...       ...
23521     chrX 155264589 155258241         2      RAB39B       923
23522     chrX 155334657 155276211         2       CLIC2       924
23523     chrX 155380787 155381134         1      H2AFB1       925
23524     chrX 155382115 155383230         1        F8A2       926
23525     chrX 155669944 155490115         2       TMLHE       927

2023-11-16 03:48:29.349478 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927735.arrow, 0.017 mins elapsed.

2023-11-16 03:48:49.043076 : mat SRR13927735, Class = dgCMatrix
mat SRR13927735: nRows = 23525, nCols = 3813
mat SRR13927735: NonZeroEntries = 45535634, EntryRange = [ 7.20715653308254e-05 , 673.990433360078 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927735#TTATGTCTCCAGGTAT-1 SRR13927735#TATTGCTCATCAGAAA-1
[1,]                    .                              .          
[2,]                    .                              .          
[3,]                    0.070892732                    0.057009268
[4,]                    0.008800039                    0.002978408
[5,]                    .                              .          
     SRR13927735#TTCGATTGTAGGGTTG-1 SRR13927735#CATTCATTCGGATGTT-1
[1,]                     .                               .        
[2,]                     .                               .        
[3,]                     0.05523011                      0.0248301
[4,]                     .                               .        
[5,]                     .                               .        
     SRR13927735#ACGTTAGGTCAACTGT-1
[1,]                     .         
[2,]                     .         
[3,]                     0.05384571
[4,]                     0.01228688
[5,]                     .         


2023-11-16 03:48:49.409805 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927735.arrow, 0.351 mins elapsed.

2023-11-16 03:48:49.733125 : se SRR13927735, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 3813 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(3813): SRR13927735#TTATGTCTCCAGGTAT-1
  SRR13927735#TATTGCTCATCAGAAA-1 ... SRR13927735#GATCATGGTGTTGTTG-1
  SRR13927735#GTCCATCCATTTCACT-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-16 03:48:49.738704 : Completed GeneIntegrationMatrix : SRR13927735(3 of 4), 1.169 mins elapsed.
2023-11-16 03:48:49.739149 : Reading GeneIntegrationMatrix : SRR13927736(4 of 4), 1.169 mins elapsed.

2023-11-16 03:48:49.762331 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927736 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927736.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3534
[1] "SRR13927736#GCTGTTCTCCTTGACC-1" "SRR13927736#AGTTACGAGTAGGTCG-1"
[3] "SRR13927736#AGTTTGGAGCATTCCA-1" "SRR13927736#GCACGGTTCGGCAATT-1"
[5] "SRR13927736#TCAGCTCTCTCTATTG-1" "SRR13927736#CCTTAATTCTACTTTG-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-16 03:48:50.799886 : featureDF SRR13927736, Class = DFrame
DataFrame with 23525 rows and 6 columns
      seqnames     start       end    strand        name       idx
         <Rle> <integer> <integer> <integer> <character> <integer>
1         chr1     69091     70008         1       OR4F5         1
2         chr1    817371    819837         1      FAM87B         2
3         chr1    825138    859446         1   LINC01128         3
4         chr1    827522    826206         2   LINC00115         4
5         chr1    876903    868071         2      FAM41C         5
...        ...       ...       ...       ...         ...       ...
23521     chrX 155264589 155258241         2      RAB39B       923
23522     chrX 155334657 155276211         2       CLIC2       924
23523     chrX 155380787 155381134         1      H2AFB1       925
23524     chrX 155382115 155383230         1        F8A2       926
23525     chrX 155669944 155490115         2       TMLHE       927

2023-11-16 03:48:50.83644 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927736.arrow, 0.018 mins elapsed.

2023-11-16 03:49:07.262917 : mat SRR13927736, Class = dgCMatrix
mat SRR13927736: nRows = 23525, nCols = 3534
mat SRR13927736: NonZeroEntries = 39003616, EntryRange = [ 8.77627031078293e-05 , 770.387276519725 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927736#GCTGTTCTCCTTGACC-1 SRR13927736#AGTTACGAGTAGGTCG-1
[1,]                     .                             .          
[2,]                     .                             .          
[3,]                     0.07845381                    0.028344539
[4,]                     0.01395289                    0.008499189
[5,]                     .                             .          
     SRR13927736#AGTTTGGAGCATTCCA-1 SRR13927736#GCACGGTTCGGCAATT-1
[1,]                     .                              .         
[2,]                     .                              .         
[3,]                     0.06535130                     0.06576222
[4,]                     0.01810177                     0.01315728
[5,]                     .                              .         
     SRR13927736#TCAGCTCTCTCTATTG-1
[1,]                    .          
[2,]                    .          
[3,]                    0.016229526
[4,]                    0.003862592
[5,]                    .          


2023-11-16 03:49:07.574581 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927736.arrow, 0.297 mins elapsed.

2023-11-16 03:49:07.83126 : se SRR13927736, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 3534 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(3534): SRR13927736#GCTGTTCTCCTTGACC-1
  SRR13927736#AGTTACGAGTAGGTCG-1 ... SRR13927736#TAAACCGCAAGTCTGT-1
  SRR13927736#TGCATTTGTTCTACCC-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-16 03:49:07.836251 : Completed GeneIntegrationMatrix : SRR13927736(4 of 4), 1.471 mins elapsed.
2023-11-16 03:49:07.837314 : Organizing colData, 1.471 mins elapsed.
2023-11-16 03:49:07.986126 : Organizing rowData, 1.473 mins elapsed.
2023-11-16 03:49:07.989334 : Organizing rowRanges, 1.473 mins elapsed.
2023-11-16 03:49:07.995036 : Organizing Assays (1 of 1), 1.473 mins elapsed.
2023-11-16 03:49:17.016996 : Constructing SummarizedExperiment, 1.624 mins elapsed.
2023-11-16 03:49:19.770077 : Finished Matrix Creation, 1.67 mins elapsed.

2023-11-16 03:49:19.772529 : getMatrixFromProject Input-Parameters, Class = list

getMatrixFromProject Input-Parameters$ArchRProj: length = 1

getMatrixFromProject Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromProject Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromProject Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromProject Input-Parameters$verbose: length = 1
[1] TRUE


getMatrixFromProject Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromProject Input-Parameters$threads: length = 1
[1] 1


getMatrixFromProject Input-Parameters$logFile: length = 1
[1] "x"


2023-11-16 03:49:20.911326 : Reading PeakMatrix : SRR13927737(1 of 4), 0.019 mins elapsed.

2023-11-16 03:49:20.942228 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927737 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927737.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 4439
[1] "SRR13927737#GTCACGGGTCAACAGG-1" "SRR13927737#TGCCTGTTCTTACTCA-1"
[3] "SRR13927737#TTGCGGGTCAGGTCTA-1" "SRR13927737#GGTTGCGCAAAGCTTC-1"
[5] "SRR13927737#CTTAATCTCGCGCTGA-1" "SRR13927737#TGTACGATCGGAAAGC-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-16 03:49:21.539495 : featureDF SRR13927737, Class = DFrame
DataFrame with 210721 rows and 4 columns
       seqnames       idx     start       end
          <Rle> <integer> <integer> <integer>
1          chr1         1    804677    805177
2          chr1         2    811866    812366
3          chr1         3    817090    817590
4          chr1         4    818831    819331
5          chr1         5    826543    827043
...         ...       ...       ...       ...
210717     chrX      5777 155898525 155899025
210718     chrX      5778 155899804 155900304
210719     chrX      5779 155956608 155957108
210720     chrX      5780 155966827 155967327
210721     chrX      5781 155971625 155972125

2023-11-16 03:49:21.602484 : Getting PeakMatrix from ArrowFile : SRR13927737.arrow, 0.011 mins elapsed.

2023-11-16 03:49:34.394531 : mat SRR13927737, Class = dgCMatrix
mat SRR13927737: nRows = 210721, nCols = 4439
mat SRR13927737: NonZeroEntries = 29963541, EntryRange = [ 1 , 4 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927737#GTCACGGGTCAACAGG-1 SRR13927737#TGCCTGTTCTTACTCA-1
[1,]                              .                              .
[2,]                              .                              2
[3,]                              .                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927737#TTGCGGGTCAGGTCTA-1 SRR13927737#GGTTGCGCAAAGCTTC-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              2                              .
[4,]                              .                              .
[5,]                              4                              .
     SRR13927737#CTTAATCTCGCGCTGA-1
[1,]                              .
[2,]                              .
[3,]                              .
[4,]                              .
[5,]                              .


2023-11-16 03:49:34.586517 : Organizing SE PeakMatrix from ArrowFile : SRR13927737.arrow, 0.227 mins elapsed.

2023-11-16 03:49:34.904989 : se SRR13927737, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 4439 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(4439): SRR13927737#GTCACGGGTCAACAGG-1
  SRR13927737#TGCCTGTTCTTACTCA-1 ... SRR13927737#TTGCTATAGATGTTCC-1
  SRR13927737#GCTCACTAGTTAGCGG-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-16 03:49:34.914519 : Completed PeakMatrix : SRR13927737(1 of 4), 0.252 mins elapsed.
2023-11-16 03:49:34.914945 : Reading PeakMatrix : SRR13927738(2 of 4), 0.252 mins elapsed.

2023-11-16 03:49:34.930557 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927738 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927738.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3736
[1] "SRR13927738#TAGCATGAGGCTAAAT-1" "SRR13927738#TGTTAGGTCGGTCAGC-1"
[3] "SRR13927738#CTCTCAGCATGTCCCT-1" "SRR13927738#ACTATTCGTCCAACCG-1"
[5] "SRR13927738#TGGTCCTTCTACTGCC-1" "SRR13927738#ATTGTCTAGACTAGGC-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-16 03:49:35.232424 : featureDF SRR13927738, Class = DFrame
DataFrame with 210721 rows and 4 columns
       seqnames       idx     start       end
          <Rle> <integer> <integer> <integer>
1          chr1         1    804677    805177
2          chr1         2    811866    812366
3          chr1         3    817090    817590
4          chr1         4    818831    819331
5          chr1         5    826543    827043
...         ...       ...       ...       ...
210717     chrX      5777 155898525 155899025
210718     chrX      5778 155899804 155900304
210719     chrX      5779 155956608 155957108
210720     chrX      5780 155966827 155967327
210721     chrX      5781 155971625 155972125

2023-11-16 03:49:35.274956 : Getting PeakMatrix from ArrowFile : SRR13927738.arrow, 0.006 mins elapsed.

2023-11-16 03:49:47.374978 : mat SRR13927738, Class = dgCMatrix
mat SRR13927738: nRows = 210721, nCols = 3736
mat SRR13927738: NonZeroEntries = 30165921, EntryRange = [ 1 , 4 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927738#TAGCATGAGGCTAAAT-1 SRR13927738#TGTTAGGTCGGTCAGC-1
[1,]                              2                              .
[2,]                              .                              .
[3,]                              .                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927738#CTCTCAGCATGTCCCT-1 SRR13927738#ACTATTCGTCCAACCG-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              .                              .
[4,]                              .                              .
[5,]                              2                              .
     SRR13927738#TGGTCCTTCTACTGCC-1
[1,]                              .
[2,]                              .
[3,]                              2
[4,]                              .
[5,]                              .


2023-11-16 03:49:47.616376 : Organizing SE PeakMatrix from ArrowFile : SRR13927738.arrow, 0.211 mins elapsed.

2023-11-16 03:49:47.911508 : se SRR13927738, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 3736 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(3736): SRR13927738#TAGCATGAGGCTAAAT-1
  SRR13927738#TGTTAGGTCGGTCAGC-1 ... SRR13927738#GCTCCTATCCGTGCGA-1
  SRR13927738#GGCACGTAGTCGATAA-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-16 03:49:47.918817 : Completed PeakMatrix : SRR13927738(2 of 4), 0.469 mins elapsed.
2023-11-16 03:49:47.91925 : Reading PeakMatrix : SRR13927735(3 of 4), 0.469 mins elapsed.

2023-11-16 03:49:47.93203 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927735 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927735.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3813
[1] "SRR13927735#TTATGTCTCCAGGTAT-1" "SRR13927735#TATTGCTCATCAGAAA-1"
[3] "SRR13927735#TTCGATTGTAGGGTTG-1" "SRR13927735#CATTCATTCGGATGTT-1"
[5] "SRR13927735#ACGTTAGGTCAACTGT-1" "SRR13927735#AAATGCCCAGCAATGG-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-16 03:49:48.178533 : featureDF SRR13927735, Class = DFrame
DataFrame with 210721 rows and 4 columns
       seqnames       idx     start       end
          <Rle> <integer> <integer> <integer>
1          chr1         1    804677    805177
2          chr1         2    811866    812366
3          chr1         3    817090    817590
4          chr1         4    818831    819331
5          chr1         5    826543    827043
...         ...       ...       ...       ...
210717     chrX      5777 155898525 155899025
210718     chrX      5778 155899804 155900304
210719     chrX      5779 155956608 155957108
210720     chrX      5780 155966827 155967327
210721     chrX      5781 155971625 155972125

2023-11-16 03:49:48.218843 : Getting PeakMatrix from ArrowFile : SRR13927735.arrow, 0.005 mins elapsed.

2023-11-16 03:49:59.082146 : mat SRR13927735, Class = dgCMatrix
mat SRR13927735: nRows = 210721, nCols = 3813
mat SRR13927735: NonZeroEntries = 26038084, EntryRange = [ 1 , 4 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927735#TTATGTCTCCAGGTAT-1 SRR13927735#TATTGCTCATCAGAAA-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              1                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927735#TTCGATTGTAGGGTTG-1 SRR13927735#CATTCATTCGGATGTT-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              .                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927735#ACGTTAGGTCAACTGT-1
[1,]                              2
[2,]                              .
[3,]                              .
[4,]                              .
[5,]                              .


2023-11-16 03:49:59.292242 : Organizing SE PeakMatrix from ArrowFile : SRR13927735.arrow, 0.189 mins elapsed.

2023-11-16 03:49:59.578011 : se SRR13927735, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 3813 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(3813): SRR13927735#TTATGTCTCCAGGTAT-1
  SRR13927735#TATTGCTCATCAGAAA-1 ... SRR13927735#GATCATGGTGTTGTTG-1
  SRR13927735#GTCCATCCATTTCACT-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-16 03:49:59.585167 : Completed PeakMatrix : SRR13927735(3 of 4), 0.664 mins elapsed.
2023-11-16 03:49:59.585611 : Reading PeakMatrix : SRR13927736(4 of 4), 0.664 mins elapsed.

2023-11-16 03:49:59.598188 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927736 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927736.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3534
[1] "SRR13927736#GCTGTTCTCCTTGACC-1" "SRR13927736#AGTTACGAGTAGGTCG-1"
[3] "SRR13927736#AGTTTGGAGCATTCCA-1" "SRR13927736#GCACGGTTCGGCAATT-1"
[5] "SRR13927736#TCAGCTCTCTCTATTG-1" "SRR13927736#CCTTAATTCTACTTTG-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-16 03:49:59.830064 : featureDF SRR13927736, Class = DFrame
DataFrame with 210721 rows and 4 columns
       seqnames       idx     start       end
          <Rle> <integer> <integer> <integer>
1          chr1         1    804677    805177
2          chr1         2    811866    812366
3          chr1         3    817090    817590
4          chr1         4    818831    819331
5          chr1         5    826543    827043
...         ...       ...       ...       ...
210717     chrX      5777 155898525 155899025
210718     chrX      5778 155899804 155900304
210719     chrX      5779 155956608 155957108
210720     chrX      5780 155966827 155967327
210721     chrX      5781 155971625 155972125

2023-11-16 03:49:59.87186 : Getting PeakMatrix from ArrowFile : SRR13927736.arrow, 0.005 mins elapsed.

2023-11-16 03:50:10.612028 : mat SRR13927736, Class = dgCMatrix
mat SRR13927736: nRows = 210721, nCols = 3534
mat SRR13927736: NonZeroEntries = 25583439, EntryRange = [ 1 , 4 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927736#GCTGTTCTCCTTGACC-1 SRR13927736#AGTTACGAGTAGGTCG-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              .                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927736#AGTTTGGAGCATTCCA-1 SRR13927736#GCACGGTTCGGCAATT-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              1                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927736#TCAGCTCTCTCTATTG-1
[1,]                              .
[2,]                              .
[3,]                              .
[4,]                              .
[5,]                              .


2023-11-16 03:50:10.821592 : Organizing SE PeakMatrix from ArrowFile : SRR13927736.arrow, 0.187 mins elapsed.

2023-11-16 03:50:11.112301 : se SRR13927736, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 3534 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(3534): SRR13927736#GCTGTTCTCCTTGACC-1
  SRR13927736#AGTTACGAGTAGGTCG-1 ... SRR13927736#TAAACCGCAAGTCTGT-1
  SRR13927736#TGCATTTGTTCTACCC-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-16 03:50:11.119145 : Completed PeakMatrix : SRR13927736(4 of 4), 0.856 mins elapsed.
2023-11-16 03:50:11.120033 : Organizing colData, 0.856 mins elapsed.
2023-11-16 03:50:11.257588 : Organizing rowData, 0.858 mins elapsed.
2023-11-16 03:50:11.266003 : Organizing rowRanges, 0.858 mins elapsed.
2023-11-16 03:50:11.277679 : Organizing Assays (1 of 1), 0.858 mins elapsed.
2023-11-16 03:50:14.947278 : Constructing SummarizedExperiment, 0.92 mins elapsed.
2023-11-16 03:50:32.137995 : Finished Matrix Creation, 1.206 mins elapsed.

2023-11-16 04:56:51.01661 : getMatrixFromProject Input-Parameters, Class = list

getMatrixFromProject Input-Parameters$ArchRProj: length = 1

getMatrixFromProject Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromProject Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromProject Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromProject Input-Parameters$verbose: length = 1
[1] TRUE


getMatrixFromProject Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromProject Input-Parameters$threads: length = 1
[1] 1


getMatrixFromProject Input-Parameters$logFile: length = 1
[1] "x"


2023-11-16 04:56:55.964795 : Reading GeneIntegrationMatrix : SRR13927737(1 of 4), 0.083 mins elapsed.

2023-11-16 04:56:56.072272 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927737 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927737.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 4439
[1] "SRR13927737#GTCACGGGTCAACAGG-1" "SRR13927737#TGCCTGTTCTTACTCA-1" "SRR13927737#TTGCGGGTCAGGTCTA-1"
[4] "SRR13927737#GGTTGCGCAAAGCTTC-1" "SRR13927737#CTTAATCTCGCGCTGA-1" "SRR13927737#TGTACGATCGGAAAGC-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-16 04:56:57.933668 : featureDF SRR13927737, Class = DFrame
DataFrame with 23525 rows and 6 columns
      seqnames     start       end    strand        name       idx
         <Rle> <integer> <integer> <integer> <character> <integer>
1         chr1     69091     70008         1       OR4F5         1
2         chr1    817371    819837         1      FAM87B         2
3         chr1    825138    859446         1   LINC01128         3
4         chr1    827522    826206         2   LINC00115         4
5         chr1    876903    868071         2      FAM41C         5
...        ...       ...       ...       ...         ...       ...
23521     chrX 155264589 155258241         2      RAB39B       923
23522     chrX 155334657 155276211         2       CLIC2       924
23523     chrX 155380787 155381134         1      H2AFB1       925
23524     chrX 155382115 155383230         1        F8A2       926
23525     chrX 155669944 155490115         2       TMLHE       927

2023-11-16 04:56:57.976414 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927737.arrow, 0.032 mins elapsed.

2023-11-16 04:57:30.448378 : mat SRR13927737, Class = dgCMatrix
mat SRR13927737: nRows = 23525, nCols = 4439
mat SRR13927737: NonZeroEntries = 50277774, EntryRange = [ 0.000148130054904683 , 376.345870564644 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927737#GTCACGGGTCAACAGG-1 SRR13927737#TGCCTGTTCTTACTCA-1 SRR13927737#TTGCGGGTCAGGTCTA-1
[1,]                      .                             .                              .         
[2,]                      .                             .                              .         
[3,]                      0.1057576                     0.09630696                     0.09113499
[4,]                      .                             0.01806508                     0.02096304
[5,]                      .                             .                              .         
     SRR13927737#GGTTGCGCAAAGCTTC-1 SRR13927737#CTTAATCTCGCGCTGA-1
[1,]                     .                              .         
[2,]                     .                              .         
[3,]                     0.08396973                     0.08327649
[4,]                     0.01727603                     0.01368917
[5,]                     .                              .         


2023-11-16 04:57:32.614036 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927737.arrow, 0.609 mins elapsed.

2023-11-16 04:57:34.697214 : se SRR13927737, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 4439 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(4439): SRR13927737#GTCACGGGTCAACAGG-1 SRR13927737#TGCCTGTTCTTACTCA-1 ...
  SRR13927737#TTGCTATAGATGTTCC-1 SRR13927737#GCTCACTAGTTAGCGG-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-16 04:57:34.712385 : Completed GeneIntegrationMatrix : SRR13927737(1 of 4), 0.728 mins elapsed.
2023-11-16 04:57:34.713542 : Reading GeneIntegrationMatrix : SRR13927738(2 of 4), 0.728 mins elapsed.

2023-11-16 04:57:34.743128 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927738 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927738.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3736
[1] "SRR13927738#TAGCATGAGGCTAAAT-1" "SRR13927738#TGTTAGGTCGGTCAGC-1" "SRR13927738#CTCTCAGCATGTCCCT-1"
[4] "SRR13927738#ACTATTCGTCCAACCG-1" "SRR13927738#TGGTCCTTCTACTGCC-1" "SRR13927738#ATTGTCTAGACTAGGC-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-16 04:57:36.78884 : featureDF SRR13927738, Class = DFrame
DataFrame with 23525 rows and 6 columns
      seqnames     start       end    strand        name       idx
         <Rle> <integer> <integer> <integer> <character> <integer>
1         chr1     69091     70008         1       OR4F5         1
2         chr1    817371    819837         1      FAM87B         2
3         chr1    825138    859446         1   LINC01128         3
4         chr1    827522    826206         2   LINC00115         4
5         chr1    876903    868071         2      FAM41C         5
...        ...       ...       ...       ...         ...       ...
23521     chrX 155264589 155258241         2      RAB39B       923
23522     chrX 155334657 155276211         2       CLIC2       924
23523     chrX 155380787 155381134         1      H2AFB1       925
23524     chrX 155382115 155383230         1        F8A2       926
23525     chrX 155669944 155490115         2       TMLHE       927

2023-11-16 04:57:38.524873 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927738.arrow, 0.063 mins elapsed.

2023-11-16 04:58:04.055872 : mat SRR13927738, Class = dgCMatrix
mat SRR13927738: nRows = 23525, nCols = 3736
mat SRR13927738: NonZeroEntries = 41898897, EntryRange = [ 9.73975135001421e-05 , 758.763031303812 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927738#TAGCATGAGGCTAAAT-1 SRR13927738#TGTTAGGTCGGTCAGC-1 SRR13927738#CTCTCAGCATGTCCCT-1
[1,]                    .                                .                             .         
[2,]                    .                                .                             .         
[3,]                    0.002346352                      0.0084575                     .         
[4,]                    0.005836413                      .                             0.02001583
[5,]                    .                                .                             .         
     SRR13927738#ACTATTCGTCCAACCG-1 SRR13927738#TGGTCCTTCTACTGCC-1
[1,]                    .                              .          
[2,]                    .                              .          
[3,]                    0.004081077                    0.007640877
[4,]                    .                              0.015048796
[5,]                    .                              .          


2023-11-16 04:58:06.431276 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927738.arrow, 0.528 mins elapsed.

2023-11-16 04:58:07.958327 : se SRR13927738, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 3736 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(3736): SRR13927738#TAGCATGAGGCTAAAT-1 SRR13927738#TGTTAGGTCGGTCAGC-1 ...
  SRR13927738#GCTCCTATCCGTGCGA-1 SRR13927738#GGCACGTAGTCGATAA-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-16 04:58:07.96853 : Completed GeneIntegrationMatrix : SRR13927738(2 of 4), 1.283 mins elapsed.
2023-11-16 04:58:07.969355 : Reading GeneIntegrationMatrix : SRR13927735(3 of 4), 1.283 mins elapsed.

2023-11-16 04:58:08.051576 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927735 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927735.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3813
[1] "SRR13927735#TTATGTCTCCAGGTAT-1" "SRR13927735#TATTGCTCATCAGAAA-1" "SRR13927735#TTCGATTGTAGGGTTG-1"
[4] "SRR13927735#CATTCATTCGGATGTT-1" "SRR13927735#ACGTTAGGTCAACTGT-1" "SRR13927735#AAATGCCCAGCAATGG-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-16 04:58:10.538479 : featureDF SRR13927735, Class = DFrame
DataFrame with 23525 rows and 6 columns
      seqnames     start       end    strand        name       idx
         <Rle> <integer> <integer> <integer> <character> <integer>
1         chr1     69091     70008         1       OR4F5         1
2         chr1    817371    819837         1      FAM87B         2
3         chr1    825138    859446         1   LINC01128         3
4         chr1    827522    826206         2   LINC00115         4
5         chr1    876903    868071         2      FAM41C         5
...        ...       ...       ...       ...         ...       ...
23521     chrX 155264589 155258241         2      RAB39B       923
23522     chrX 155334657 155276211         2       CLIC2       924
23523     chrX 155380787 155381134         1      H2AFB1       925
23524     chrX 155382115 155383230         1        F8A2       926
23525     chrX 155669944 155490115         2       TMLHE       927

2023-11-16 04:58:10.575021 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927735.arrow, 0.042 mins elapsed.

2023-11-16 04:58:37.654355 : mat SRR13927735, Class = dgCMatrix
mat SRR13927735: nRows = 23525, nCols = 3813
mat SRR13927735: NonZeroEntries = 45535634, EntryRange = [ 7.20715653308254e-05 , 673.990433360078 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927735#TTATGTCTCCAGGTAT-1 SRR13927735#TATTGCTCATCAGAAA-1 SRR13927735#TTCGATTGTAGGGTTG-1
[1,]                    .                              .                               .         
[2,]                    .                              .                               .         
[3,]                    0.070892732                    0.057009268                     0.05523011
[4,]                    0.008800039                    0.002978408                     .         
[5,]                    .                              .                               .         
     SRR13927735#CATTCATTCGGATGTT-1 SRR13927735#ACGTTAGGTCAACTGT-1
[1,]                      .                             .         
[2,]                      .                             .         
[3,]                      0.0248301                     0.05384571
[4,]                      .                             0.01228688
[5,]                      .                             .         


2023-11-16 04:58:39.615747 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927735.arrow, 0.526 mins elapsed.

2023-11-16 04:58:41.431431 : se SRR13927735, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 3813 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(3813): SRR13927735#TTATGTCTCCAGGTAT-1 SRR13927735#TATTGCTCATCAGAAA-1 ...
  SRR13927735#GATCATGGTGTTGTTG-1 SRR13927735#GTCCATCCATTTCACT-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-16 04:58:41.444293 : Completed GeneIntegrationMatrix : SRR13927735(3 of 4), 1.841 mins elapsed.
2023-11-16 04:58:41.445424 : Reading GeneIntegrationMatrix : SRR13927736(4 of 4), 1.841 mins elapsed.

2023-11-16 04:58:42.261657 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927736 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927736.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3534
[1] "SRR13927736#GCTGTTCTCCTTGACC-1" "SRR13927736#AGTTACGAGTAGGTCG-1" "SRR13927736#AGTTTGGAGCATTCCA-1"
[4] "SRR13927736#GCACGGTTCGGCAATT-1" "SRR13927736#TCAGCTCTCTCTATTG-1" "SRR13927736#CCTTAATTCTACTTTG-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-16 04:58:45.169224 : featureDF SRR13927736, Class = DFrame
DataFrame with 23525 rows and 6 columns
      seqnames     start       end    strand        name       idx
         <Rle> <integer> <integer> <integer> <character> <integer>
1         chr1     69091     70008         1       OR4F5         1
2         chr1    817371    819837         1      FAM87B         2
3         chr1    825138    859446         1   LINC01128         3
4         chr1    827522    826206         2   LINC00115         4
5         chr1    876903    868071         2      FAM41C         5
...        ...       ...       ...       ...         ...       ...
23521     chrX 155264589 155258241         2      RAB39B       923
23522     chrX 155334657 155276211         2       CLIC2       924
23523     chrX 155380787 155381134         1      H2AFB1       925
23524     chrX 155382115 155383230         1        F8A2       926
23525     chrX 155669944 155490115         2       TMLHE       927

2023-11-16 04:58:45.206305 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927736.arrow, 0.049 mins elapsed.

2023-11-16 04:59:08.622289 : mat SRR13927736, Class = dgCMatrix
mat SRR13927736: nRows = 23525, nCols = 3534
mat SRR13927736: NonZeroEntries = 39003616, EntryRange = [ 8.77627031078293e-05 , 770.387276519725 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927736#GCTGTTCTCCTTGACC-1 SRR13927736#AGTTACGAGTAGGTCG-1 SRR13927736#AGTTTGGAGCATTCCA-1
[1,]                     .                             .                               .         
[2,]                     .                             .                               .         
[3,]                     0.07845381                    0.028344539                     0.06535130
[4,]                     0.01395289                    0.008499189                     0.01810177
[5,]                     .                             .                               .         
     SRR13927736#GCACGGTTCGGCAATT-1 SRR13927736#TCAGCTCTCTCTATTG-1
[1,]                     .                             .          
[2,]                     .                             .          
[3,]                     0.06576222                    0.016229526
[4,]                     0.01315728                    0.003862592
[5,]                     .                             .          


2023-11-16 04:59:10.721316 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927736.arrow, 0.474 mins elapsed.

2023-11-16 04:59:12.304729 : se SRR13927736, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 3534 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(3534): SRR13927736#GCTGTTCTCCTTGACC-1 SRR13927736#AGTTACGAGTAGGTCG-1 ...
  SRR13927736#TAAACCGCAAGTCTGT-1 SRR13927736#TGCATTTGTTCTACCC-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-16 04:59:12.317865 : Completed GeneIntegrationMatrix : SRR13927736(4 of 4), 2.355 mins elapsed.
2023-11-16 04:59:12.325937 : Organizing colData, 2.355 mins elapsed.
2023-11-16 04:59:13.667262 : Organizing rowData, 2.378 mins elapsed.
2023-11-16 04:59:13.675348 : Organizing rowRanges, 2.378 mins elapsed.
2023-11-16 04:59:13.689287 : Organizing Assays (1 of 1), 2.378 mins elapsed.
2023-11-16 04:59:24.952601 : Constructing SummarizedExperiment, 2.566 mins elapsed.
2023-11-16 04:59:28.55604 : Finished Matrix Creation, 2.616 mins elapsed.

2023-11-16 04:59:28.565143 : getMatrixFromProject Input-Parameters, Class = list

getMatrixFromProject Input-Parameters$ArchRProj: length = 1

getMatrixFromProject Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromProject Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromProject Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromProject Input-Parameters$verbose: length = 1
[1] TRUE


getMatrixFromProject Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromProject Input-Parameters$threads: length = 1
[1] 1


getMatrixFromProject Input-Parameters$logFile: length = 1
[1] "x"


2023-11-16 04:59:31.793331 : Reading PeakMatrix : SRR13927737(1 of 4), 0.054 mins elapsed.

2023-11-16 04:59:31.899011 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927737 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927737.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 4439
[1] "SRR13927737#GTCACGGGTCAACAGG-1" "SRR13927737#TGCCTGTTCTTACTCA-1" "SRR13927737#TTGCGGGTCAGGTCTA-1"
[4] "SRR13927737#GGTTGCGCAAAGCTTC-1" "SRR13927737#CTTAATCTCGCGCTGA-1" "SRR13927737#TGTACGATCGGAAAGC-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-16 04:59:33.46228 : featureDF SRR13927737, Class = DFrame
DataFrame with 210721 rows and 4 columns
       seqnames       idx     start       end
          <Rle> <integer> <integer> <integer>
1          chr1         1    804677    805177
2          chr1         2    811866    812366
3          chr1         3    817090    817590
4          chr1         4    818831    819331
5          chr1         5    826543    827043
...         ...       ...       ...       ...
210717     chrX      5777 155898525 155899025
210718     chrX      5778 155899804 155900304
210719     chrX      5779 155956608 155957108
210720     chrX      5780 155966827 155967327
210721     chrX      5781 155971625 155972125

2023-11-16 04:59:33.543992 : Getting PeakMatrix from ArrowFile : SRR13927737.arrow, 0.027 mins elapsed.

2023-11-16 04:59:52.722138 : mat SRR13927737, Class = dgCMatrix
mat SRR13927737: nRows = 210721, nCols = 4439
mat SRR13927737: NonZeroEntries = 29963541, EntryRange = [ 1 , 4 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927737#GTCACGGGTCAACAGG-1 SRR13927737#TGCCTGTTCTTACTCA-1 SRR13927737#TTGCGGGTCAGGTCTA-1
[1,]                              .                              .                              .
[2,]                              .                              2                              .
[3,]                              .                              .                              2
[4,]                              .                              .                              .
[5,]                              .                              .                              4
     SRR13927737#GGTTGCGCAAAGCTTC-1 SRR13927737#CTTAATCTCGCGCTGA-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              .                              .
[4,]                              .                              .
[5,]                              .                              .


2023-11-16 04:59:52.930115 : Organizing SE PeakMatrix from ArrowFile : SRR13927737.arrow, 0.351 mins elapsed.

2023-11-16 04:59:56.53894 : se SRR13927737, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 4439 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(4439): SRR13927737#GTCACGGGTCAACAGG-1 SRR13927737#TGCCTGTTCTTACTCA-1 ...
  SRR13927737#TTGCTATAGATGTTCC-1 SRR13927737#GCTCACTAGTTAGCGG-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-16 04:59:56.555345 : Completed PeakMatrix : SRR13927737(1 of 4), 0.467 mins elapsed.
2023-11-16 04:59:56.556111 : Reading PeakMatrix : SRR13927738(2 of 4), 0.467 mins elapsed.

2023-11-16 04:59:56.576986 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927738 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927738.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3736
[1] "SRR13927738#TAGCATGAGGCTAAAT-1" "SRR13927738#TGTTAGGTCGGTCAGC-1" "SRR13927738#CTCTCAGCATGTCCCT-1"
[4] "SRR13927738#ACTATTCGTCCAACCG-1" "SRR13927738#TGGTCCTTCTACTGCC-1" "SRR13927738#ATTGTCTAGACTAGGC-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-16 04:59:59.100333 : featureDF SRR13927738, Class = DFrame
DataFrame with 210721 rows and 4 columns
       seqnames       idx     start       end
          <Rle> <integer> <integer> <integer>
1          chr1         1    804677    805177
2          chr1         2    811866    812366
3          chr1         3    817090    817590
4          chr1         4    818831    819331
5          chr1         5    826543    827043
...         ...       ...       ...       ...
210717     chrX      5777 155898525 155899025
210718     chrX      5778 155899804 155900304
210719     chrX      5779 155956608 155957108
210720     chrX      5780 155966827 155967327
210721     chrX      5781 155971625 155972125

2023-11-16 04:59:59.212683 : Getting PeakMatrix from ArrowFile : SRR13927738.arrow, 0.044 mins elapsed.

2023-11-16 05:00:15.75531 : mat SRR13927738, Class = dgCMatrix
mat SRR13927738: nRows = 210721, nCols = 3736
mat SRR13927738: NonZeroEntries = 30165921, EntryRange = [ 1 , 4 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927738#TAGCATGAGGCTAAAT-1 SRR13927738#TGTTAGGTCGGTCAGC-1 SRR13927738#CTCTCAGCATGTCCCT-1
[1,]                              2                              .                              .
[2,]                              .                              .                              .
[3,]                              .                              .                              .
[4,]                              .                              .                              .
[5,]                              .                              .                              2
     SRR13927738#ACTATTCGTCCAACCG-1 SRR13927738#TGGTCCTTCTACTGCC-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              .                              2
[4,]                              .                              .
[5,]                              .                              .


2023-11-16 05:00:18.7787 : Organizing SE PeakMatrix from ArrowFile : SRR13927738.arrow, 0.37 mins elapsed.

2023-11-16 05:00:22.732098 : se SRR13927738, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 3736 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(3736): SRR13927738#TAGCATGAGGCTAAAT-1 SRR13927738#TGTTAGGTCGGTCAGC-1 ...
  SRR13927738#GCTCCTATCCGTGCGA-1 SRR13927738#GGCACGTAGTCGATAA-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-16 05:00:22.749622 : Completed PeakMatrix : SRR13927738(2 of 4), 0.903 mins elapsed.
2023-11-16 05:00:22.750735 : Reading PeakMatrix : SRR13927735(3 of 4), 0.903 mins elapsed.

2023-11-16 05:00:22.779887 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927735 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927735.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3813
[1] "SRR13927735#TTATGTCTCCAGGTAT-1" "SRR13927735#TATTGCTCATCAGAAA-1" "SRR13927735#TTCGATTGTAGGGTTG-1"
[4] "SRR13927735#CATTCATTCGGATGTT-1" "SRR13927735#ACGTTAGGTCAACTGT-1" "SRR13927735#AAATGCCCAGCAATGG-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-16 05:00:26.84875 : featureDF SRR13927735, Class = DFrame
DataFrame with 210721 rows and 4 columns
       seqnames       idx     start       end
          <Rle> <integer> <integer> <integer>
1          chr1         1    804677    805177
2          chr1         2    811866    812366
3          chr1         3    817090    817590
4          chr1         4    818831    819331
5          chr1         5    826543    827043
...         ...       ...       ...       ...
210717     chrX      5777 155898525 155899025
210718     chrX      5778 155899804 155900304
210719     chrX      5779 155956608 155957108
210720     chrX      5780 155966827 155967327
210721     chrX      5781 155971625 155972125

2023-11-16 05:00:28.400246 : Getting PeakMatrix from ArrowFile : SRR13927735.arrow, 0.094 mins elapsed.

2023-11-16 05:00:43.761042 : mat SRR13927735, Class = dgCMatrix
mat SRR13927735: nRows = 210721, nCols = 3813
mat SRR13927735: NonZeroEntries = 26038084, EntryRange = [ 1 , 4 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927735#TTATGTCTCCAGGTAT-1 SRR13927735#TATTGCTCATCAGAAA-1 SRR13927735#TTCGATTGTAGGGTTG-1
[1,]                              .                              .                              .
[2,]                              .                              .                              .
[3,]                              1                              .                              .
[4,]                              .                              .                              .
[5,]                              .                              .                              .
     SRR13927735#CATTCATTCGGATGTT-1 SRR13927735#ACGTTAGGTCAACTGT-1
[1,]                              .                              2
[2,]                              .                              .
[3,]                              .                              .
[4,]                              .                              .
[5,]                              .                              .


2023-11-16 05:00:43.976018 : Organizing SE PeakMatrix from ArrowFile : SRR13927735.arrow, 0.353 mins elapsed.

2023-11-16 05:00:47.203347 : se SRR13927735, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 3813 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(3813): SRR13927735#TTATGTCTCCAGGTAT-1 SRR13927735#TATTGCTCATCAGAAA-1 ...
  SRR13927735#GATCATGGTGTTGTTG-1 SRR13927735#GTCCATCCATTTCACT-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-16 05:00:47.221219 : Completed PeakMatrix : SRR13927735(3 of 4), 1.311 mins elapsed.
2023-11-16 05:00:47.222293 : Reading PeakMatrix : SRR13927736(4 of 4), 1.311 mins elapsed.

2023-11-16 05:00:47.25162 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927736 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927736.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3534
[1] "SRR13927736#GCTGTTCTCCTTGACC-1" "SRR13927736#AGTTACGAGTAGGTCG-1" "SRR13927736#AGTTTGGAGCATTCCA-1"
[4] "SRR13927736#GCACGGTTCGGCAATT-1" "SRR13927736#TCAGCTCTCTCTATTG-1" "SRR13927736#CCTTAATTCTACTTTG-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-16 05:00:49.95911 : featureDF SRR13927736, Class = DFrame
DataFrame with 210721 rows and 4 columns
       seqnames       idx     start       end
          <Rle> <integer> <integer> <integer>
1          chr1         1    804677    805177
2          chr1         2    811866    812366
3          chr1         3    817090    817590
4          chr1         4    818831    819331
5          chr1         5    826543    827043
...         ...       ...       ...       ...
210717     chrX      5777 155898525 155899025
210718     chrX      5778 155899804 155900304
210719     chrX      5779 155956608 155957108
210720     chrX      5780 155966827 155967327
210721     chrX      5781 155971625 155972125

2023-11-16 05:00:50.038154 : Getting PeakMatrix from ArrowFile : SRR13927736.arrow, 0.046 mins elapsed.

2023-11-16 05:01:05.021386 : mat SRR13927736, Class = dgCMatrix
mat SRR13927736: nRows = 210721, nCols = 3534
mat SRR13927736: NonZeroEntries = 25583439, EntryRange = [ 1 , 4 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927736#GCTGTTCTCCTTGACC-1 SRR13927736#AGTTACGAGTAGGTCG-1 SRR13927736#AGTTTGGAGCATTCCA-1
[1,]                              .                              .                              .
[2,]                              .                              .                              .
[3,]                              .                              .                              1
[4,]                              .                              .                              .
[5,]                              .                              .                              .
     SRR13927736#GCACGGTTCGGCAATT-1 SRR13927736#TCAGCTCTCTCTATTG-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              .                              .
[4,]                              .                              .
[5,]                              .                              .


2023-11-16 05:01:08.450836 : Organizing SE PeakMatrix from ArrowFile : SRR13927736.arrow, 0.353 mins elapsed.

2023-11-16 05:01:10.798464 : se SRR13927736, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 3534 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(3534): SRR13927736#GCTGTTCTCCTTGACC-1 SRR13927736#AGTTACGAGTAGGTCG-1 ...
  SRR13927736#TAAACCGCAAGTCTGT-1 SRR13927736#TGCATTTGTTCTACCC-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2023-11-16 05:01:10.815889 : Completed PeakMatrix : SRR13927736(4 of 4), 1.704 mins elapsed.
2023-11-16 05:01:10.818319 : Organizing colData, 1.704 mins elapsed.
2023-11-16 05:01:17.737973 : Organizing rowData, 1.82 mins elapsed.
2023-11-16 05:01:17.760749 : Organizing rowRanges, 1.82 mins elapsed.
2023-11-16 05:01:17.790576 : Organizing Assays (1 of 1), 1.82 mins elapsed.
2023-11-16 05:01:25.965777 : Constructing SummarizedExperiment, 1.957 mins elapsed.
2023-11-16 05:01:47.176593 : Finished Matrix Creation, 2.31 mins elapsed.

2023-11-20 02:50:56.935308 : getMatrixFromProject Input-Parameters, Class = list

getMatrixFromProject Input-Parameters$ArchRProj: length = 1

getMatrixFromProject Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromProject Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromProject Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromProject Input-Parameters$verbose: length = 1
[1] TRUE


getMatrixFromProject Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromProject Input-Parameters$threads: length = 1
[1] 1


getMatrixFromProject Input-Parameters$logFile: length = 1
[1] "x"


2023-11-20 02:50:57.981784 : Reading GeneIntegrationMatrix : SRR13927737(1 of 4), 0.018 mins elapsed.

2023-11-20 02:50:58.016622 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927737 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927737.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 4439
[1] "SRR13927737#GTCACGGGTCAACAGG-1" "SRR13927737#TGCCTGTTCTTACTCA-1"
[3] "SRR13927737#TTGCGGGTCAGGTCTA-1" "SRR13927737#GGTTGCGCAAAGCTTC-1"
[5] "SRR13927737#CTTAATCTCGCGCTGA-1" "SRR13927737#TGTACGATCGGAAAGC-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-20 02:50:58.394957 : featureDF SRR13927737, Class = DFrame
DataFrame with 23525 rows and 6 columns
      seqnames     start       end    strand        name       idx
         <Rle> <integer> <integer> <integer> <character> <integer>
1         chr1     69091     70008         1       OR4F5         1
2         chr1    817371    819837         1      FAM87B         2
3         chr1    825138    859446         1   LINC01128         3
4         chr1    827522    826206         2   LINC00115         4
5         chr1    876903    868071         2      FAM41C         5
...        ...       ...       ...       ...         ...       ...
23521     chrX 155264589 155258241         2      RAB39B       923
23522     chrX 155334657 155276211         2       CLIC2       924
23523     chrX 155380787 155381134         1      H2AFB1       925
23524     chrX 155382115 155383230         1        F8A2       926
23525     chrX 155669944 155490115         2       TMLHE       927

2023-11-20 02:50:58.423724 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927737.arrow, 0.007 mins elapsed.

2023-11-20 02:51:22.580848 : mat SRR13927737, Class = dgCMatrix
mat SRR13927737: nRows = 23525, nCols = 4439
mat SRR13927737: NonZeroEntries = 50277774, EntryRange = [ 0.000148130054904683 , 376.345870564644 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927737#GTCACGGGTCAACAGG-1 SRR13927737#TGCCTGTTCTTACTCA-1
[1,]                      .                             .         
[2,]                      .                             .         
[3,]                      0.1057576                     0.09630696
[4,]                      .                             0.01806508
[5,]                      .                             .         
     SRR13927737#TTGCGGGTCAGGTCTA-1 SRR13927737#GGTTGCGCAAAGCTTC-1
[1,]                     .                              .         
[2,]                     .                              .         
[3,]                     0.09113499                     0.08396973
[4,]                     0.02096304                     0.01727603
[5,]                     .                              .         
     SRR13927737#CTTAATCTCGCGCTGA-1
[1,]                     .         
[2,]                     .         
[3,]                     0.08327649
[4,]                     0.01368917
[5,]                     .         


2023-11-20 02:51:23.151749 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927737.arrow, 0.419 mins elapsed.

2023-11-20 02:51:23.465941 : se SRR13927737, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 4439 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(4439): SRR13927737#GTCACGGGTCAACAGG-1
  SRR13927737#TGCCTGTTCTTACTCA-1 ... SRR13927737#TTGCTATAGATGTTCC-1
  SRR13927737#GCTCACTAGTTAGCGG-1
colData names(21): BlacklistRatio nDiFrags ... ReadsInPeaks FRIP

2023-11-20 02:51:23.472637 : Completed GeneIntegrationMatrix : SRR13927737(1 of 4), 0.442 mins elapsed.
2023-11-20 02:51:23.473205 : Reading GeneIntegrationMatrix : SRR13927738(2 of 4), 0.442 mins elapsed.

2023-11-20 02:51:23.487179 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927738 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927738.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3736
[1] "SRR13927738#TAGCATGAGGCTAAAT-1" "SRR13927738#TGTTAGGTCGGTCAGC-1"
[3] "SRR13927738#CTCTCAGCATGTCCCT-1" "SRR13927738#ACTATTCGTCCAACCG-1"
[5] "SRR13927738#TGGTCCTTCTACTGCC-1" "SRR13927738#ATTGTCTAGACTAGGC-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-20 02:51:23.682309 : featureDF SRR13927738, Class = DFrame
DataFrame with 23525 rows and 6 columns
      seqnames     start       end    strand        name       idx
         <Rle> <integer> <integer> <integer> <character> <integer>
1         chr1     69091     70008         1       OR4F5         1
2         chr1    817371    819837         1      FAM87B         2
3         chr1    825138    859446         1   LINC01128         3
4         chr1    827522    826206         2   LINC00115         4
5         chr1    876903    868071         2      FAM41C         5
...        ...       ...       ...       ...         ...       ...
23521     chrX 155264589 155258241         2      RAB39B       923
23522     chrX 155334657 155276211         2       CLIC2       924
23523     chrX 155380787 155381134         1      H2AFB1       925
23524     chrX 155382115 155383230         1        F8A2       926
23525     chrX 155669944 155490115         2       TMLHE       927

2023-11-20 02:51:23.699161 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927738.arrow, 0.004 mins elapsed.

2023-11-20 02:51:46.444423 : mat SRR13927738, Class = dgCMatrix
mat SRR13927738: nRows = 23525, nCols = 3736
mat SRR13927738: NonZeroEntries = 41898897, EntryRange = [ 9.73975135001421e-05 , 758.763031303812 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927738#TAGCATGAGGCTAAAT-1 SRR13927738#TGTTAGGTCGGTCAGC-1
[1,]                    .                                .        
[2,]                    .                                .        
[3,]                    0.002346352                      0.0084575
[4,]                    0.005836413                      .        
[5,]                    .                                .        
     SRR13927738#CTCTCAGCATGTCCCT-1 SRR13927738#ACTATTCGTCCAACCG-1
[1,]                     .                             .          
[2,]                     .                             .          
[3,]                     .                             0.004081077
[4,]                     0.02001583                    .          
[5,]                     .                             .          
     SRR13927738#TGGTCCTTCTACTGCC-1
[1,]                    .          
[2,]                    .          
[3,]                    0.007640877
[4,]                    0.015048796
[5,]                    .          


2023-11-20 02:51:46.855319 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927738.arrow, 0.389 mins elapsed.

2023-11-20 02:51:47.161915 : se SRR13927738, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 3736 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(3736): SRR13927738#TAGCATGAGGCTAAAT-1
  SRR13927738#TGTTAGGTCGGTCAGC-1 ... SRR13927738#GCTCCTATCCGTGCGA-1
  SRR13927738#GGCACGTAGTCGATAA-1
colData names(21): BlacklistRatio nDiFrags ... ReadsInPeaks FRIP

2023-11-20 02:51:47.16772 : Completed GeneIntegrationMatrix : SRR13927738(2 of 4), 0.837 mins elapsed.
2023-11-20 02:51:47.16829 : Reading GeneIntegrationMatrix : SRR13927735(3 of 4), 0.837 mins elapsed.

2023-11-20 02:51:47.182037 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927735 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927735.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3813
[1] "SRR13927735#TTATGTCTCCAGGTAT-1" "SRR13927735#TATTGCTCATCAGAAA-1"
[3] "SRR13927735#TTCGATTGTAGGGTTG-1" "SRR13927735#CATTCATTCGGATGTT-1"
[5] "SRR13927735#ACGTTAGGTCAACTGT-1" "SRR13927735#AAATGCCCAGCAATGG-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-20 02:51:47.377951 : featureDF SRR13927735, Class = DFrame
DataFrame with 23525 rows and 6 columns
      seqnames     start       end    strand        name       idx
         <Rle> <integer> <integer> <integer> <character> <integer>
1         chr1     69091     70008         1       OR4F5         1
2         chr1    817371    819837         1      FAM87B         2
3         chr1    825138    859446         1   LINC01128         3
4         chr1    827522    826206         2   LINC00115         4
5         chr1    876903    868071         2      FAM41C         5
...        ...       ...       ...       ...         ...       ...
23521     chrX 155264589 155258241         2      RAB39B       923
23522     chrX 155334657 155276211         2       CLIC2       924
23523     chrX 155380787 155381134         1      H2AFB1       925
23524     chrX 155382115 155383230         1        F8A2       926
23525     chrX 155669944 155490115         2       TMLHE       927

2023-11-20 02:51:47.393748 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927735.arrow, 0.004 mins elapsed.

2023-11-20 02:52:08.084442 : mat SRR13927735, Class = dgCMatrix
mat SRR13927735: nRows = 23525, nCols = 3813
mat SRR13927735: NonZeroEntries = 45535634, EntryRange = [ 7.20715653308254e-05 , 673.990433360078 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927735#TTATGTCTCCAGGTAT-1 SRR13927735#TATTGCTCATCAGAAA-1
[1,]                    .                              .          
[2,]                    .                              .          
[3,]                    0.070892732                    0.057009268
[4,]                    0.008800039                    0.002978408
[5,]                    .                              .          
     SRR13927735#TTCGATTGTAGGGTTG-1 SRR13927735#CATTCATTCGGATGTT-1
[1,]                     .                               .        
[2,]                     .                               .        
[3,]                     0.05523011                      0.0248301
[4,]                     .                               .        
[5,]                     .                               .        
     SRR13927735#ACGTTAGGTCAACTGT-1
[1,]                     .         
[2,]                     .         
[3,]                     0.05384571
[4,]                     0.01228688
[5,]                     .         


2023-11-20 02:52:08.550973 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927735.arrow, 0.356 mins elapsed.

2023-11-20 02:52:08.827218 : se SRR13927735, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 3813 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(3813): SRR13927735#TTATGTCTCCAGGTAT-1
  SRR13927735#TATTGCTCATCAGAAA-1 ... SRR13927735#GATCATGGTGTTGTTG-1
  SRR13927735#GTCCATCCATTTCACT-1
colData names(21): BlacklistRatio nDiFrags ... ReadsInPeaks FRIP

2023-11-20 02:52:08.832862 : Completed GeneIntegrationMatrix : SRR13927735(3 of 4), 1.198 mins elapsed.
2023-11-20 02:52:08.833427 : Reading GeneIntegrationMatrix : SRR13927736(4 of 4), 1.198 mins elapsed.

2023-11-20 02:52:08.846299 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927736 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927736.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3534
[1] "SRR13927736#GCTGTTCTCCTTGACC-1" "SRR13927736#AGTTACGAGTAGGTCG-1"
[3] "SRR13927736#AGTTTGGAGCATTCCA-1" "SRR13927736#GCACGGTTCGGCAATT-1"
[5] "SRR13927736#TCAGCTCTCTCTATTG-1" "SRR13927736#CCTTAATTCTACTTTG-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-20 02:52:09.028031 : featureDF SRR13927736, Class = DFrame
DataFrame with 23525 rows and 6 columns
      seqnames     start       end    strand        name       idx
         <Rle> <integer> <integer> <integer> <character> <integer>
1         chr1     69091     70008         1       OR4F5         1
2         chr1    817371    819837         1      FAM87B         2
3         chr1    825138    859446         1   LINC01128         3
4         chr1    827522    826206         2   LINC00115         4
5         chr1    876903    868071         2      FAM41C         5
...        ...       ...       ...       ...         ...       ...
23521     chrX 155264589 155258241         2      RAB39B       923
23522     chrX 155334657 155276211         2       CLIC2       924
23523     chrX 155380787 155381134         1      H2AFB1       925
23524     chrX 155382115 155383230         1        F8A2       926
23525     chrX 155669944 155490115         2       TMLHE       927

2023-11-20 02:52:09.043906 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927736.arrow, 0.003 mins elapsed.

2023-11-20 02:52:27.50027 : mat SRR13927736, Class = dgCMatrix
mat SRR13927736: nRows = 23525, nCols = 3534
mat SRR13927736: NonZeroEntries = 39003616, EntryRange = [ 8.77627031078293e-05 , 770.387276519725 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927736#GCTGTTCTCCTTGACC-1 SRR13927736#AGTTACGAGTAGGTCG-1
[1,]                     .                             .          
[2,]                     .                             .          
[3,]                     0.07845381                    0.028344539
[4,]                     0.01395289                    0.008499189
[5,]                     .                             .          
     SRR13927736#AGTTTGGAGCATTCCA-1 SRR13927736#GCACGGTTCGGCAATT-1
[1,]                     .                              .         
[2,]                     .                              .         
[3,]                     0.06535130                     0.06576222
[4,]                     0.01810177                     0.01315728
[5,]                     .                              .         
     SRR13927736#TCAGCTCTCTCTATTG-1
[1,]                    .          
[2,]                    .          
[3,]                    0.016229526
[4,]                    0.003862592
[5,]                    .          


2023-11-20 02:52:27.919782 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927736.arrow, 0.318 mins elapsed.

2023-11-20 02:52:28.210406 : se SRR13927736, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 3534 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(3534): SRR13927736#GCTGTTCTCCTTGACC-1
  SRR13927736#AGTTACGAGTAGGTCG-1 ... SRR13927736#TAAACCGCAAGTCTGT-1
  SRR13927736#TGCATTTGTTCTACCC-1
colData names(21): BlacklistRatio nDiFrags ... ReadsInPeaks FRIP

2023-11-20 02:52:28.216552 : Completed GeneIntegrationMatrix : SRR13927736(4 of 4), 1.521 mins elapsed.
2023-11-20 02:52:28.224223 : Organizing colData, 1.521 mins elapsed.
2023-11-20 02:52:28.392253 : Organizing rowData, 1.524 mins elapsed.
2023-11-20 02:52:28.396249 : Organizing rowRanges, 1.524 mins elapsed.
2023-11-20 02:52:28.403021 : Organizing Assays (1 of 1), 1.525 mins elapsed.
2023-11-20 02:52:37.530073 : Constructing SummarizedExperiment, 1.677 mins elapsed.
2023-11-20 02:52:40.771331 : Finished Matrix Creation, 1.731 mins elapsed.

2023-11-20 02:52:40.774164 : getMatrixFromProject Input-Parameters, Class = list

getMatrixFromProject Input-Parameters$ArchRProj: length = 1

getMatrixFromProject Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromProject Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromProject Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromProject Input-Parameters$verbose: length = 1
[1] TRUE


getMatrixFromProject Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromProject Input-Parameters$threads: length = 1
[1] 1


getMatrixFromProject Input-Parameters$logFile: length = 1
[1] "x"


2023-11-20 02:52:41.723747 : Reading PeakMatrix : SRR13927737(1 of 4), 0.016 mins elapsed.

2023-11-20 02:52:41.745449 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927737 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927737.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 4439
[1] "SRR13927737#GTCACGGGTCAACAGG-1" "SRR13927737#TGCCTGTTCTTACTCA-1"
[3] "SRR13927737#TTGCGGGTCAGGTCTA-1" "SRR13927737#GGTTGCGCAAAGCTTC-1"
[5] "SRR13927737#CTTAATCTCGCGCTGA-1" "SRR13927737#TGTACGATCGGAAAGC-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-20 02:52:42.224352 : featureDF SRR13927737, Class = DFrame
DataFrame with 210721 rows and 4 columns
       seqnames       idx     start       end
          <Rle> <integer> <integer> <integer>
1          chr1         1    804677    805177
2          chr1         2    811866    812366
3          chr1         3    817090    817590
4          chr1         4    818831    819331
5          chr1         5    826543    827043
...         ...       ...       ...       ...
210717     chrX      5777 155898525 155899025
210718     chrX      5778 155899804 155900304
210719     chrX      5779 155956608 155957108
210720     chrX      5780 155966827 155967327
210721     chrX      5781 155971625 155972125

2023-11-20 02:52:42.295717 : Getting PeakMatrix from ArrowFile : SRR13927737.arrow, 0.009 mins elapsed.

2023-11-20 02:52:55.343688 : mat SRR13927737, Class = dgCMatrix
mat SRR13927737: nRows = 210721, nCols = 4439
mat SRR13927737: NonZeroEntries = 29963541, EntryRange = [ 1 , 4 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927737#GTCACGGGTCAACAGG-1 SRR13927737#TGCCTGTTCTTACTCA-1
[1,]                              .                              .
[2,]                              .                              2
[3,]                              .                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927737#TTGCGGGTCAGGTCTA-1 SRR13927737#GGTTGCGCAAAGCTTC-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              2                              .
[4,]                              .                              .
[5,]                              4                              .
     SRR13927737#CTTAATCTCGCGCTGA-1
[1,]                              .
[2,]                              .
[3,]                              .
[4,]                              .
[5,]                              .


2023-11-20 02:52:55.586143 : Organizing SE PeakMatrix from ArrowFile : SRR13927737.arrow, 0.231 mins elapsed.

2023-11-20 02:52:55.896003 : se SRR13927737, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 4439 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(4439): SRR13927737#GTCACGGGTCAACAGG-1
  SRR13927737#TGCCTGTTCTTACTCA-1 ... SRR13927737#TTGCTATAGATGTTCC-1
  SRR13927737#GCTCACTAGTTAGCGG-1
colData names(21): BlacklistRatio nDiFrags ... ReadsInPeaks FRIP

2023-11-20 02:52:55.906015 : Completed PeakMatrix : SRR13927737(1 of 4), 0.252 mins elapsed.
2023-11-20 02:52:55.906669 : Reading PeakMatrix : SRR13927738(2 of 4), 0.252 mins elapsed.

2023-11-20 02:52:55.92014 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927738 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927738.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3736
[1] "SRR13927738#TAGCATGAGGCTAAAT-1" "SRR13927738#TGTTAGGTCGGTCAGC-1"
[3] "SRR13927738#CTCTCAGCATGTCCCT-1" "SRR13927738#ACTATTCGTCCAACCG-1"
[5] "SRR13927738#TGGTCCTTCTACTGCC-1" "SRR13927738#ATTGTCTAGACTAGGC-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-20 02:52:56.243222 : featureDF SRR13927738, Class = DFrame
DataFrame with 210721 rows and 4 columns
       seqnames       idx     start       end
          <Rle> <integer> <integer> <integer>
1          chr1         1    804677    805177
2          chr1         2    811866    812366
3          chr1         3    817090    817590
4          chr1         4    818831    819331
5          chr1         5    826543    827043
...         ...       ...       ...       ...
210717     chrX      5777 155898525 155899025
210718     chrX      5778 155899804 155900304
210719     chrX      5779 155956608 155957108
210720     chrX      5780 155966827 155967327
210721     chrX      5781 155971625 155972125

2023-11-20 02:52:56.289641 : Getting PeakMatrix from ArrowFile : SRR13927738.arrow, 0.006 mins elapsed.

2023-11-20 02:53:09.24545 : mat SRR13927738, Class = dgCMatrix
mat SRR13927738: nRows = 210721, nCols = 3736
mat SRR13927738: NonZeroEntries = 30165921, EntryRange = [ 1 , 4 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927738#TAGCATGAGGCTAAAT-1 SRR13927738#TGTTAGGTCGGTCAGC-1
[1,]                              2                              .
[2,]                              .                              .
[3,]                              .                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927738#CTCTCAGCATGTCCCT-1 SRR13927738#ACTATTCGTCCAACCG-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              .                              .
[4,]                              .                              .
[5,]                              2                              .
     SRR13927738#TGGTCCTTCTACTGCC-1
[1,]                              .
[2,]                              .
[3,]                              2
[4,]                              .
[5,]                              .


2023-11-20 02:53:09.485727 : Organizing SE PeakMatrix from ArrowFile : SRR13927738.arrow, 0.226 mins elapsed.

2023-11-20 02:53:09.778052 : se SRR13927738, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 3736 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(3736): SRR13927738#TAGCATGAGGCTAAAT-1
  SRR13927738#TGTTAGGTCGGTCAGC-1 ... SRR13927738#GCTCCTATCCGTGCGA-1
  SRR13927738#GGCACGTAGTCGATAA-1
colData names(21): BlacklistRatio nDiFrags ... ReadsInPeaks FRIP

2023-11-20 02:53:09.785746 : Completed PeakMatrix : SRR13927738(2 of 4), 0.484 mins elapsed.
2023-11-20 02:53:09.786399 : Reading PeakMatrix : SRR13927735(3 of 4), 0.484 mins elapsed.

2023-11-20 02:53:09.800049 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927735 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927735.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3813
[1] "SRR13927735#TTATGTCTCCAGGTAT-1" "SRR13927735#TATTGCTCATCAGAAA-1"
[3] "SRR13927735#TTCGATTGTAGGGTTG-1" "SRR13927735#CATTCATTCGGATGTT-1"
[5] "SRR13927735#ACGTTAGGTCAACTGT-1" "SRR13927735#AAATGCCCAGCAATGG-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-20 02:53:10.478585 : featureDF SRR13927735, Class = DFrame
DataFrame with 210721 rows and 4 columns
       seqnames       idx     start       end
          <Rle> <integer> <integer> <integer>
1          chr1         1    804677    805177
2          chr1         2    811866    812366
3          chr1         3    817090    817590
4          chr1         4    818831    819331
5          chr1         5    826543    827043
...         ...       ...       ...       ...
210717     chrX      5777 155898525 155899025
210718     chrX      5778 155899804 155900304
210719     chrX      5779 155956608 155957108
210720     chrX      5780 155966827 155967327
210721     chrX      5781 155971625 155972125

2023-11-20 02:53:10.565153 : Getting PeakMatrix from ArrowFile : SRR13927735.arrow, 0.013 mins elapsed.

2023-11-20 02:53:22.914511 : mat SRR13927735, Class = dgCMatrix
mat SRR13927735: nRows = 210721, nCols = 3813
mat SRR13927735: NonZeroEntries = 26038084, EntryRange = [ 1 , 4 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927735#TTATGTCTCCAGGTAT-1 SRR13927735#TATTGCTCATCAGAAA-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              1                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927735#TTCGATTGTAGGGTTG-1 SRR13927735#CATTCATTCGGATGTT-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              .                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927735#ACGTTAGGTCAACTGT-1
[1,]                              2
[2,]                              .
[3,]                              .
[4,]                              .
[5,]                              .


2023-11-20 02:53:23.103083 : Organizing SE PeakMatrix from ArrowFile : SRR13927735.arrow, 0.222 mins elapsed.

2023-11-20 02:53:23.429681 : se SRR13927735, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 3813 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(3813): SRR13927735#TTATGTCTCCAGGTAT-1
  SRR13927735#TATTGCTCATCAGAAA-1 ... SRR13927735#GATCATGGTGTTGTTG-1
  SRR13927735#GTCCATCCATTTCACT-1
colData names(21): BlacklistRatio nDiFrags ... ReadsInPeaks FRIP

2023-11-20 02:53:23.437839 : Completed PeakMatrix : SRR13927735(3 of 4), 0.711 mins elapsed.
2023-11-20 02:53:23.438472 : Reading PeakMatrix : SRR13927736(4 of 4), 0.711 mins elapsed.

2023-11-20 02:53:23.452515 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927736 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927736.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3534
[1] "SRR13927736#GCTGTTCTCCTTGACC-1" "SRR13927736#AGTTACGAGTAGGTCG-1"
[3] "SRR13927736#AGTTTGGAGCATTCCA-1" "SRR13927736#GCACGGTTCGGCAATT-1"
[5] "SRR13927736#TCAGCTCTCTCTATTG-1" "SRR13927736#CCTTAATTCTACTTTG-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-20 02:53:23.977116 : featureDF SRR13927736, Class = DFrame
DataFrame with 210721 rows and 4 columns
       seqnames       idx     start       end
          <Rle> <integer> <integer> <integer>
1          chr1         1    804677    805177
2          chr1         2    811866    812366
3          chr1         3    817090    817590
4          chr1         4    818831    819331
5          chr1         5    826543    827043
...         ...       ...       ...       ...
210717     chrX      5777 155898525 155899025
210718     chrX      5778 155899804 155900304
210719     chrX      5779 155956608 155957108
210720     chrX      5780 155966827 155967327
210721     chrX      5781 155971625 155972125

2023-11-20 02:53:24.065688 : Getting PeakMatrix from ArrowFile : SRR13927736.arrow, 0.01 mins elapsed.

2023-11-20 02:53:36.635261 : mat SRR13927736, Class = dgCMatrix
mat SRR13927736: nRows = 210721, nCols = 3534
mat SRR13927736: NonZeroEntries = 25583439, EntryRange = [ 1 , 4 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927736#GCTGTTCTCCTTGACC-1 SRR13927736#AGTTACGAGTAGGTCG-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              .                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927736#AGTTTGGAGCATTCCA-1 SRR13927736#GCACGGTTCGGCAATT-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              1                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927736#TCAGCTCTCTCTATTG-1
[1,]                              .
[2,]                              .
[3,]                              .
[4,]                              .
[5,]                              .


2023-11-20 02:53:36.886649 : Organizing SE PeakMatrix from ArrowFile : SRR13927736.arrow, 0.224 mins elapsed.

2023-11-20 02:53:37.223824 : se SRR13927736, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 3534 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(3534): SRR13927736#GCTGTTCTCCTTGACC-1
  SRR13927736#AGTTACGAGTAGGTCG-1 ... SRR13927736#TAAACCGCAAGTCTGT-1
  SRR13927736#TGCATTTGTTCTACCC-1
colData names(21): BlacklistRatio nDiFrags ... ReadsInPeaks FRIP

2023-11-20 02:53:37.233066 : Completed PeakMatrix : SRR13927736(4 of 4), 0.941 mins elapsed.
2023-11-20 02:53:37.234599 : Organizing colData, 0.941 mins elapsed.
2023-11-20 02:53:37.425753 : Organizing rowData, 0.944 mins elapsed.
2023-11-20 02:53:37.437838 : Organizing rowRanges, 0.944 mins elapsed.
2023-11-20 02:53:37.453424 : Organizing Assays (1 of 1), 0.945 mins elapsed.
2023-11-20 02:53:41.002239 : Constructing SummarizedExperiment, 1.004 mins elapsed.
2023-11-20 02:54:02.612444 : Finished Matrix Creation, 1.364 mins elapsed.

2023-11-20 04:28:20.904416 : getMatrixFromProject Input-Parameters, Class = list

getMatrixFromProject Input-Parameters$ArchRProj: length = 1

getMatrixFromProject Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromProject Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromProject Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromProject Input-Parameters$verbose: length = 1
[1] TRUE


getMatrixFromProject Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromProject Input-Parameters$threads: length = 1
[1] 1


getMatrixFromProject Input-Parameters$logFile: length = 1
[1] "x"


2023-11-20 04:28:22.088192 : Reading GeneIntegrationMatrix : SRR13927737(1 of 4), 0.02 mins elapsed.

2023-11-20 04:28:22.123441 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927737 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927737.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 4439
[1] "SRR13927737#GTCACGGGTCAACAGG-1" "SRR13927737#TGCCTGTTCTTACTCA-1"
[3] "SRR13927737#TTGCGGGTCAGGTCTA-1" "SRR13927737#GGTTGCGCAAAGCTTC-1"
[5] "SRR13927737#CTTAATCTCGCGCTGA-1" "SRR13927737#TGTACGATCGGAAAGC-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-20 04:28:22.516602 : featureDF SRR13927737, Class = DFrame
DataFrame with 23525 rows and 6 columns
      seqnames     start       end    strand        name       idx
         <Rle> <integer> <integer> <integer> <character> <integer>
1         chr1     69091     70008         1       OR4F5         1
2         chr1    817371    819837         1      FAM87B         2
3         chr1    825138    859446         1   LINC01128         3
4         chr1    827522    826206         2   LINC00115         4
5         chr1    876903    868071         2      FAM41C         5
...        ...       ...       ...       ...         ...       ...
23521     chrX 155264589 155258241         2      RAB39B       923
23522     chrX 155334657 155276211         2       CLIC2       924
23523     chrX 155380787 155381134         1      H2AFB1       925
23524     chrX 155382115 155383230         1        F8A2       926
23525     chrX 155669944 155490115         2       TMLHE       927

2023-11-20 04:28:22.543623 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927737.arrow, 0.007 mins elapsed.

2023-11-20 04:28:48.132433 : mat SRR13927737, Class = dgCMatrix
mat SRR13927737: nRows = 23525, nCols = 4439
mat SRR13927737: NonZeroEntries = 50277774, EntryRange = [ 0.000148130054904683 , 376.345870564644 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927737#GTCACGGGTCAACAGG-1 SRR13927737#TGCCTGTTCTTACTCA-1
[1,]                      .                             .         
[2,]                      .                             .         
[3,]                      0.1057576                     0.09630696
[4,]                      .                             0.01806508
[5,]                      .                             .         
     SRR13927737#TTGCGGGTCAGGTCTA-1 SRR13927737#GGTTGCGCAAAGCTTC-1
[1,]                     .                              .         
[2,]                     .                              .         
[3,]                     0.09113499                     0.08396973
[4,]                     0.02096304                     0.01727603
[5,]                     .                              .         
     SRR13927737#CTTAATCTCGCGCTGA-1
[1,]                     .         
[2,]                     .         
[3,]                     0.08327649
[4,]                     0.01368917
[5,]                     .         


2023-11-20 04:28:48.700859 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927737.arrow, 0.443 mins elapsed.

2023-11-20 04:28:49.067316 : se SRR13927737, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 4439 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(4439): SRR13927737#GTCACGGGTCAACAGG-1
  SRR13927737#TGCCTGTTCTTACTCA-1 ... SRR13927737#TTGCTATAGATGTTCC-1
  SRR13927737#GCTCACTAGTTAGCGG-1
colData names(21): BlacklistRatio nDiFrags ... ReadsInPeaks FRIP

2023-11-20 04:28:49.07477 : Completed GeneIntegrationMatrix : SRR13927737(1 of 4), 0.47 mins elapsed.
2023-11-20 04:28:49.07546 : Reading GeneIntegrationMatrix : SRR13927738(2 of 4), 0.47 mins elapsed.

2023-11-20 04:28:49.090478 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927738 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927738.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3736
[1] "SRR13927738#TAGCATGAGGCTAAAT-1" "SRR13927738#TGTTAGGTCGGTCAGC-1"
[3] "SRR13927738#CTCTCAGCATGTCCCT-1" "SRR13927738#ACTATTCGTCCAACCG-1"
[5] "SRR13927738#TGGTCCTTCTACTGCC-1" "SRR13927738#ATTGTCTAGACTAGGC-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-20 04:28:49.303923 : featureDF SRR13927738, Class = DFrame
DataFrame with 23525 rows and 6 columns
      seqnames     start       end    strand        name       idx
         <Rle> <integer> <integer> <integer> <character> <integer>
1         chr1     69091     70008         1       OR4F5         1
2         chr1    817371    819837         1      FAM87B         2
3         chr1    825138    859446         1   LINC01128         3
4         chr1    827522    826206         2   LINC00115         4
5         chr1    876903    868071         2      FAM41C         5
...        ...       ...       ...       ...         ...       ...
23521     chrX 155264589 155258241         2      RAB39B       923
23522     chrX 155334657 155276211         2       CLIC2       924
23523     chrX 155380787 155381134         1      H2AFB1       925
23524     chrX 155382115 155383230         1        F8A2       926
23525     chrX 155669944 155490115         2       TMLHE       927

2023-11-20 04:28:49.32258 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927738.arrow, 0.004 mins elapsed.

2023-11-20 04:29:11.725497 : mat SRR13927738, Class = dgCMatrix
mat SRR13927738: nRows = 23525, nCols = 3736
mat SRR13927738: NonZeroEntries = 41898897, EntryRange = [ 9.73975135001421e-05 , 758.763031303812 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927738#TAGCATGAGGCTAAAT-1 SRR13927738#TGTTAGGTCGGTCAGC-1
[1,]                    .                                .        
[2,]                    .                                .        
[3,]                    0.002346352                      0.0084575
[4,]                    0.005836413                      .        
[5,]                    .                                .        
     SRR13927738#CTCTCAGCATGTCCCT-1 SRR13927738#ACTATTCGTCCAACCG-1
[1,]                     .                             .          
[2,]                     .                             .          
[3,]                     .                             0.004081077
[4,]                     0.02001583                    .          
[5,]                     .                             .          
     SRR13927738#TGGTCCTTCTACTGCC-1
[1,]                    .          
[2,]                    .          
[3,]                    0.007640877
[4,]                    0.015048796
[5,]                    .          


2023-11-20 04:29:12.134392 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927738.arrow, 0.384 mins elapsed.

2023-11-20 04:29:12.447331 : se SRR13927738, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 3736 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(3736): SRR13927738#TAGCATGAGGCTAAAT-1
  SRR13927738#TGTTAGGTCGGTCAGC-1 ... SRR13927738#GCTCCTATCCGTGCGA-1
  SRR13927738#GGCACGTAGTCGATAA-1
colData names(21): BlacklistRatio nDiFrags ... ReadsInPeaks FRIP

2023-11-20 04:29:12.453533 : Completed GeneIntegrationMatrix : SRR13927738(2 of 4), 0.859 mins elapsed.
2023-11-20 04:29:12.454229 : Reading GeneIntegrationMatrix : SRR13927735(3 of 4), 0.859 mins elapsed.

2023-11-20 04:29:12.468629 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927735 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927735.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3813
[1] "SRR13927735#TTATGTCTCCAGGTAT-1" "SRR13927735#TATTGCTCATCAGAAA-1"
[3] "SRR13927735#TTCGATTGTAGGGTTG-1" "SRR13927735#CATTCATTCGGATGTT-1"
[5] "SRR13927735#ACGTTAGGTCAACTGT-1" "SRR13927735#AAATGCCCAGCAATGG-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-20 04:29:12.675363 : featureDF SRR13927735, Class = DFrame
DataFrame with 23525 rows and 6 columns
      seqnames     start       end    strand        name       idx
         <Rle> <integer> <integer> <integer> <character> <integer>
1         chr1     69091     70008         1       OR4F5         1
2         chr1    817371    819837         1      FAM87B         2
3         chr1    825138    859446         1   LINC01128         3
4         chr1    827522    826206         2   LINC00115         4
5         chr1    876903    868071         2      FAM41C         5
...        ...       ...       ...       ...         ...       ...
23521     chrX 155264589 155258241         2      RAB39B       923
23522     chrX 155334657 155276211         2       CLIC2       924
23523     chrX 155380787 155381134         1      H2AFB1       925
23524     chrX 155382115 155383230         1        F8A2       926
23525     chrX 155669944 155490115         2       TMLHE       927

2023-11-20 04:29:12.693356 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927735.arrow, 0.004 mins elapsed.

2023-11-20 04:29:32.526881 : mat SRR13927735, Class = dgCMatrix
mat SRR13927735: nRows = 23525, nCols = 3813
mat SRR13927735: NonZeroEntries = 45535634, EntryRange = [ 7.20715653308254e-05 , 673.990433360078 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927735#TTATGTCTCCAGGTAT-1 SRR13927735#TATTGCTCATCAGAAA-1
[1,]                    .                              .          
[2,]                    .                              .          
[3,]                    0.070892732                    0.057009268
[4,]                    0.008800039                    0.002978408
[5,]                    .                              .          
     SRR13927735#TTCGATTGTAGGGTTG-1 SRR13927735#CATTCATTCGGATGTT-1
[1,]                     .                               .        
[2,]                     .                               .        
[3,]                     0.05523011                      0.0248301
[4,]                     .                               .        
[5,]                     .                               .        
     SRR13927735#ACGTTAGGTCAACTGT-1
[1,]                     .         
[2,]                     .         
[3,]                     0.05384571
[4,]                     0.01228688
[5,]                     .         


2023-11-20 04:29:32.978095 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927735.arrow, 0.342 mins elapsed.

2023-11-20 04:29:33.262895 : se SRR13927735, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 3813 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(3813): SRR13927735#TTATGTCTCCAGGTAT-1
  SRR13927735#TATTGCTCATCAGAAA-1 ... SRR13927735#GATCATGGTGTTGTTG-1
  SRR13927735#GTCCATCCATTTCACT-1
colData names(21): BlacklistRatio nDiFrags ... ReadsInPeaks FRIP

2023-11-20 04:29:33.268652 : Completed GeneIntegrationMatrix : SRR13927735(3 of 4), 1.206 mins elapsed.
2023-11-20 04:29:33.269336 : Reading GeneIntegrationMatrix : SRR13927736(4 of 4), 1.206 mins elapsed.

2023-11-20 04:29:33.283543 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927736 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927736.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3534
[1] "SRR13927736#GCTGTTCTCCTTGACC-1" "SRR13927736#AGTTACGAGTAGGTCG-1"
[3] "SRR13927736#AGTTTGGAGCATTCCA-1" "SRR13927736#GCACGGTTCGGCAATT-1"
[5] "SRR13927736#TCAGCTCTCTCTATTG-1" "SRR13927736#CCTTAATTCTACTTTG-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-20 04:29:33.477752 : featureDF SRR13927736, Class = DFrame
DataFrame with 23525 rows and 6 columns
      seqnames     start       end    strand        name       idx
         <Rle> <integer> <integer> <integer> <character> <integer>
1         chr1     69091     70008         1       OR4F5         1
2         chr1    817371    819837         1      FAM87B         2
3         chr1    825138    859446         1   LINC01128         3
4         chr1    827522    826206         2   LINC00115         4
5         chr1    876903    868071         2      FAM41C         5
...        ...       ...       ...       ...         ...       ...
23521     chrX 155264589 155258241         2      RAB39B       923
23522     chrX 155334657 155276211         2       CLIC2       924
23523     chrX 155380787 155381134         1      H2AFB1       925
23524     chrX 155382115 155383230         1        F8A2       926
23525     chrX 155669944 155490115         2       TMLHE       927

2023-11-20 04:29:33.494557 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927736.arrow, 0.004 mins elapsed.

2023-11-20 04:29:49.569222 : mat SRR13927736, Class = dgCMatrix
mat SRR13927736: nRows = 23525, nCols = 3534
mat SRR13927736: NonZeroEntries = 39003616, EntryRange = [ 8.77627031078293e-05 , 770.387276519725 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927736#GCTGTTCTCCTTGACC-1 SRR13927736#AGTTACGAGTAGGTCG-1
[1,]                     .                             .          
[2,]                     .                             .          
[3,]                     0.07845381                    0.028344539
[4,]                     0.01395289                    0.008499189
[5,]                     .                             .          
     SRR13927736#AGTTTGGAGCATTCCA-1 SRR13927736#GCACGGTTCGGCAATT-1
[1,]                     .                              .         
[2,]                     .                              .         
[3,]                     0.06535130                     0.06576222
[4,]                     0.01810177                     0.01315728
[5,]                     .                              .         
     SRR13927736#TCAGCTCTCTCTATTG-1
[1,]                    .          
[2,]                    .          
[3,]                    0.016229526
[4,]                    0.003862592
[5,]                    .          


2023-11-20 04:29:49.956011 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927736.arrow, 0.278 mins elapsed.

2023-11-20 04:29:50.230959 : se SRR13927736, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 3534 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(3534): SRR13927736#GCTGTTCTCCTTGACC-1
  SRR13927736#AGTTACGAGTAGGTCG-1 ... SRR13927736#TAAACCGCAAGTCTGT-1
  SRR13927736#TGCATTTGTTCTACCC-1
colData names(21): BlacklistRatio nDiFrags ... ReadsInPeaks FRIP

2023-11-20 04:29:50.236664 : Completed GeneIntegrationMatrix : SRR13927736(4 of 4), 1.489 mins elapsed.
2023-11-20 04:29:50.243573 : Organizing colData, 1.489 mins elapsed.
2023-11-20 04:29:50.399392 : Organizing rowData, 1.492 mins elapsed.
2023-11-20 04:29:50.403297 : Organizing rowRanges, 1.492 mins elapsed.
2023-11-20 04:29:50.409203 : Organizing Assays (1 of 1), 1.492 mins elapsed.
2023-11-20 04:29:58.212804 : Constructing SummarizedExperiment, 1.622 mins elapsed.
2023-11-20 04:30:01.008341 : Finished Matrix Creation, 1.668 mins elapsed.

2023-11-20 04:30:01.011143 : getMatrixFromProject Input-Parameters, Class = list

getMatrixFromProject Input-Parameters$ArchRProj: length = 1

getMatrixFromProject Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromProject Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromProject Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromProject Input-Parameters$verbose: length = 1
[1] TRUE


getMatrixFromProject Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromProject Input-Parameters$threads: length = 1
[1] 1


getMatrixFromProject Input-Parameters$logFile: length = 1
[1] "x"


2023-11-20 04:30:01.694738 : Reading PeakMatrix : SRR13927737(1 of 4), 0.011 mins elapsed.

2023-11-20 04:30:01.711034 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927737 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927737.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 4439
[1] "SRR13927737#GTCACGGGTCAACAGG-1" "SRR13927737#TGCCTGTTCTTACTCA-1"
[3] "SRR13927737#TTGCGGGTCAGGTCTA-1" "SRR13927737#GGTTGCGCAAAGCTTC-1"
[5] "SRR13927737#CTTAATCTCGCGCTGA-1" "SRR13927737#TGTACGATCGGAAAGC-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-20 04:30:02.024963 : featureDF SRR13927737, Class = DFrame
DataFrame with 210721 rows and 4 columns
       seqnames       idx     start       end
          <Rle> <integer> <integer> <integer>
1          chr1         1    804677    805177
2          chr1         2    811866    812366
3          chr1         3    817090    817590
4          chr1         4    818831    819331
5          chr1         5    826543    827043
...         ...       ...       ...       ...
210717     chrX      5777 155898525 155899025
210718     chrX      5778 155899804 155900304
210719     chrX      5779 155956608 155957108
210720     chrX      5780 155966827 155967327
210721     chrX      5781 155971625 155972125

2023-11-20 04:30:02.082472 : Getting PeakMatrix from ArrowFile : SRR13927737.arrow, 0.006 mins elapsed.

2023-11-20 04:30:14.780796 : mat SRR13927737, Class = dgCMatrix
mat SRR13927737: nRows = 210721, nCols = 4439
mat SRR13927737: NonZeroEntries = 29963541, EntryRange = [ 1 , 4 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927737#GTCACGGGTCAACAGG-1 SRR13927737#TGCCTGTTCTTACTCA-1
[1,]                              .                              .
[2,]                              .                              2
[3,]                              .                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927737#TTGCGGGTCAGGTCTA-1 SRR13927737#GGTTGCGCAAAGCTTC-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              2                              .
[4,]                              .                              .
[5,]                              4                              .
     SRR13927737#CTTAATCTCGCGCTGA-1
[1,]                              .
[2,]                              .
[3,]                              .
[4,]                              .
[5,]                              .


2023-11-20 04:30:15.008364 : Organizing SE PeakMatrix from ArrowFile : SRR13927737.arrow, 0.222 mins elapsed.

2023-11-20 04:30:15.327249 : se SRR13927737, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 4439 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(4439): SRR13927737#GTCACGGGTCAACAGG-1
  SRR13927737#TGCCTGTTCTTACTCA-1 ... SRR13927737#TTGCTATAGATGTTCC-1
  SRR13927737#GCTCACTAGTTAGCGG-1
colData names(21): BlacklistRatio nDiFrags ... ReadsInPeaks FRIP

2023-11-20 04:30:15.337894 : Completed PeakMatrix : SRR13927737(1 of 4), 0.239 mins elapsed.
2023-11-20 04:30:15.338667 : Reading PeakMatrix : SRR13927738(2 of 4), 0.239 mins elapsed.

2023-11-20 04:30:15.355561 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927738 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927738.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3736
[1] "SRR13927738#TAGCATGAGGCTAAAT-1" "SRR13927738#TGTTAGGTCGGTCAGC-1"
[3] "SRR13927738#CTCTCAGCATGTCCCT-1" "SRR13927738#ACTATTCGTCCAACCG-1"
[5] "SRR13927738#TGGTCCTTCTACTGCC-1" "SRR13927738#ATTGTCTAGACTAGGC-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-20 04:30:15.621952 : featureDF SRR13927738, Class = DFrame
DataFrame with 210721 rows and 4 columns
       seqnames       idx     start       end
          <Rle> <integer> <integer> <integer>
1          chr1         1    804677    805177
2          chr1         2    811866    812366
3          chr1         3    817090    817590
4          chr1         4    818831    819331
5          chr1         5    826543    827043
...         ...       ...       ...       ...
210717     chrX      5777 155898525 155899025
210718     chrX      5778 155899804 155900304
210719     chrX      5779 155956608 155957108
210720     chrX      5780 155966827 155967327
210721     chrX      5781 155971625 155972125

2023-11-20 04:30:15.668766 : Getting PeakMatrix from ArrowFile : SRR13927738.arrow, 0.005 mins elapsed.

2023-11-20 04:30:28.708271 : mat SRR13927738, Class = dgCMatrix
mat SRR13927738: nRows = 210721, nCols = 3736
mat SRR13927738: NonZeroEntries = 30165921, EntryRange = [ 1 , 4 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927738#TAGCATGAGGCTAAAT-1 SRR13927738#TGTTAGGTCGGTCAGC-1
[1,]                              2                              .
[2,]                              .                              .
[3,]                              .                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927738#CTCTCAGCATGTCCCT-1 SRR13927738#ACTATTCGTCCAACCG-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              .                              .
[4,]                              .                              .
[5,]                              2                              .
     SRR13927738#TGGTCCTTCTACTGCC-1
[1,]                              .
[2,]                              .
[3,]                              2
[4,]                              .
[5,]                              .


2023-11-20 04:30:28.93938 : Organizing SE PeakMatrix from ArrowFile : SRR13927738.arrow, 0.226 mins elapsed.

2023-11-20 04:30:29.290008 : se SRR13927738, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 3736 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(3736): SRR13927738#TAGCATGAGGCTAAAT-1
  SRR13927738#TGTTAGGTCGGTCAGC-1 ... SRR13927738#GCTCCTATCCGTGCGA-1
  SRR13927738#GGCACGTAGTCGATAA-1
colData names(21): BlacklistRatio nDiFrags ... ReadsInPeaks FRIP

2023-11-20 04:30:29.298723 : Completed PeakMatrix : SRR13927738(2 of 4), 0.471 mins elapsed.
2023-11-20 04:30:29.299445 : Reading PeakMatrix : SRR13927735(3 of 4), 0.471 mins elapsed.

2023-11-20 04:30:29.314633 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927735 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927735.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3813
[1] "SRR13927735#TTATGTCTCCAGGTAT-1" "SRR13927735#TATTGCTCATCAGAAA-1"
[3] "SRR13927735#TTCGATTGTAGGGTTG-1" "SRR13927735#CATTCATTCGGATGTT-1"
[5] "SRR13927735#ACGTTAGGTCAACTGT-1" "SRR13927735#AAATGCCCAGCAATGG-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-20 04:30:29.617368 : featureDF SRR13927735, Class = DFrame
DataFrame with 210721 rows and 4 columns
       seqnames       idx     start       end
          <Rle> <integer> <integer> <integer>
1          chr1         1    804677    805177
2          chr1         2    811866    812366
3          chr1         3    817090    817590
4          chr1         4    818831    819331
5          chr1         5    826543    827043
...         ...       ...       ...       ...
210717     chrX      5777 155898525 155899025
210718     chrX      5778 155899804 155900304
210719     chrX      5779 155956608 155957108
210720     chrX      5780 155966827 155967327
210721     chrX      5781 155971625 155972125

2023-11-20 04:30:29.667375 : Getting PeakMatrix from ArrowFile : SRR13927735.arrow, 0.006 mins elapsed.

2023-11-20 04:30:41.863375 : mat SRR13927735, Class = dgCMatrix
mat SRR13927735: nRows = 210721, nCols = 3813
mat SRR13927735: NonZeroEntries = 26038084, EntryRange = [ 1 , 4 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927735#TTATGTCTCCAGGTAT-1 SRR13927735#TATTGCTCATCAGAAA-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              1                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927735#TTCGATTGTAGGGTTG-1 SRR13927735#CATTCATTCGGATGTT-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              .                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927735#ACGTTAGGTCAACTGT-1
[1,]                              2
[2,]                              .
[3,]                              .
[4,]                              .
[5,]                              .


2023-11-20 04:30:42.041122 : Organizing SE PeakMatrix from ArrowFile : SRR13927735.arrow, 0.212 mins elapsed.

2023-11-20 04:30:42.356211 : se SRR13927735, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 3813 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(3813): SRR13927735#TTATGTCTCCAGGTAT-1
  SRR13927735#TATTGCTCATCAGAAA-1 ... SRR13927735#GATCATGGTGTTGTTG-1
  SRR13927735#GTCCATCCATTTCACT-1
colData names(21): BlacklistRatio nDiFrags ... ReadsInPeaks FRIP

2023-11-20 04:30:42.364067 : Completed PeakMatrix : SRR13927735(3 of 4), 0.689 mins elapsed.
2023-11-20 04:30:42.364762 : Reading PeakMatrix : SRR13927736(4 of 4), 0.689 mins elapsed.

2023-11-20 04:30:42.378834 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927736 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927736.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3534
[1] "SRR13927736#GCTGTTCTCCTTGACC-1" "SRR13927736#AGTTACGAGTAGGTCG-1"
[3] "SRR13927736#AGTTTGGAGCATTCCA-1" "SRR13927736#GCACGGTTCGGCAATT-1"
[5] "SRR13927736#TCAGCTCTCTCTATTG-1" "SRR13927736#CCTTAATTCTACTTTG-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-20 04:30:42.664991 : featureDF SRR13927736, Class = DFrame
DataFrame with 210721 rows and 4 columns
       seqnames       idx     start       end
          <Rle> <integer> <integer> <integer>
1          chr1         1    804677    805177
2          chr1         2    811866    812366
3          chr1         3    817090    817590
4          chr1         4    818831    819331
5          chr1         5    826543    827043
...         ...       ...       ...       ...
210717     chrX      5777 155898525 155899025
210718     chrX      5778 155899804 155900304
210719     chrX      5779 155956608 155957108
210720     chrX      5780 155966827 155967327
210721     chrX      5781 155971625 155972125

2023-11-20 04:30:42.711766 : Getting PeakMatrix from ArrowFile : SRR13927736.arrow, 0.006 mins elapsed.

2023-11-20 04:30:53.878183 : mat SRR13927736, Class = dgCMatrix
mat SRR13927736: nRows = 210721, nCols = 3534
mat SRR13927736: NonZeroEntries = 25583439, EntryRange = [ 1 , 4 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927736#GCTGTTCTCCTTGACC-1 SRR13927736#AGTTACGAGTAGGTCG-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              .                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927736#AGTTTGGAGCATTCCA-1 SRR13927736#GCACGGTTCGGCAATT-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              1                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927736#TCAGCTCTCTCTATTG-1
[1,]                              .
[2,]                              .
[3,]                              .
[4,]                              .
[5,]                              .


2023-11-20 04:30:54.131762 : Organizing SE PeakMatrix from ArrowFile : SRR13927736.arrow, 0.196 mins elapsed.

2023-11-20 04:30:54.438867 : se SRR13927736, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 3534 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(3534): SRR13927736#GCTGTTCTCCTTGACC-1
  SRR13927736#AGTTACGAGTAGGTCG-1 ... SRR13927736#TAAACCGCAAGTCTGT-1
  SRR13927736#TGCATTTGTTCTACCC-1
colData names(21): BlacklistRatio nDiFrags ... ReadsInPeaks FRIP

2023-11-20 04:30:54.446644 : Completed PeakMatrix : SRR13927736(4 of 4), 0.891 mins elapsed.
2023-11-20 04:30:54.448056 : Organizing colData, 0.891 mins elapsed.
2023-11-20 04:30:54.600953 : Organizing rowData, 0.893 mins elapsed.
2023-11-20 04:30:54.610404 : Organizing rowRanges, 0.893 mins elapsed.
2023-11-20 04:30:54.623472 : Organizing Assays (1 of 1), 0.894 mins elapsed.
2023-11-20 04:30:57.773187 : Constructing SummarizedExperiment, 0.946 mins elapsed.
2023-11-20 04:31:24.000274 : Finished Matrix Creation, 1.383 mins elapsed.

2023-11-20 06:03:05.869139 : getMatrixFromProject Input-Parameters, Class = list

getMatrixFromProject Input-Parameters$ArchRProj: length = 1

getMatrixFromProject Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromProject Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromProject Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromProject Input-Parameters$verbose: length = 1
[1] TRUE


getMatrixFromProject Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromProject Input-Parameters$threads: length = 1
[1] 1


getMatrixFromProject Input-Parameters$logFile: length = 1
[1] "x"


2023-11-20 06:03:06.622311 : Reading GeneIntegrationMatrix : SRR13927737(1 of 4), 0.013 mins elapsed.

2023-11-20 06:03:06.646712 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927737 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927737.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 4439
[1] "SRR13927737#GTCACGGGTCAACAGG-1" "SRR13927737#TGCCTGTTCTTACTCA-1"
[3] "SRR13927737#TTGCGGGTCAGGTCTA-1" "SRR13927737#GGTTGCGCAAAGCTTC-1"
[5] "SRR13927737#CTTAATCTCGCGCTGA-1" "SRR13927737#TGTACGATCGGAAAGC-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-20 06:03:07.016081 : featureDF SRR13927737, Class = DFrame
DataFrame with 23525 rows and 6 columns
      seqnames     start       end    strand        name       idx
         <Rle> <integer> <integer> <integer> <character> <integer>
1         chr1     69091     70008         1       OR4F5         1
2         chr1    817371    819837         1      FAM87B         2
3         chr1    825138    859446         1   LINC01128         3
4         chr1    827522    826206         2   LINC00115         4
5         chr1    876903    868071         2      FAM41C         5
...        ...       ...       ...       ...         ...       ...
23521     chrX 155264589 155258241         2      RAB39B       923
23522     chrX 155334657 155276211         2       CLIC2       924
23523     chrX 155380787 155381134         1      H2AFB1       925
23524     chrX 155382115 155383230         1        F8A2       926
23525     chrX 155669944 155490115         2       TMLHE       927

2023-11-20 06:03:07.044062 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927737.arrow, 0.007 mins elapsed.

2023-11-20 06:03:31.779162 : mat SRR13927737, Class = dgCMatrix
mat SRR13927737: nRows = 23525, nCols = 4439
mat SRR13927737: NonZeroEntries = 50277774, EntryRange = [ 0.000148130054904683 , 376.345870564644 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927737#GTCACGGGTCAACAGG-1 SRR13927737#TGCCTGTTCTTACTCA-1
[1,]                      .                             .         
[2,]                      .                             .         
[3,]                      0.1057576                     0.09630696
[4,]                      .                             0.01806508
[5,]                      .                             .         
     SRR13927737#TTGCGGGTCAGGTCTA-1 SRR13927737#GGTTGCGCAAAGCTTC-1
[1,]                     .                              .         
[2,]                     .                              .         
[3,]                     0.09113499                     0.08396973
[4,]                     0.02096304                     0.01727603
[5,]                     .                              .         
     SRR13927737#CTTAATCTCGCGCTGA-1
[1,]                     .         
[2,]                     .         
[3,]                     0.08327649
[4,]                     0.01368917
[5,]                     .         


2023-11-20 06:03:32.348403 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927737.arrow, 0.428 mins elapsed.

2023-11-20 06:03:32.705898 : se SRR13927737, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 4439 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(4439): SRR13927737#GTCACGGGTCAACAGG-1
  SRR13927737#TGCCTGTTCTTACTCA-1 ... SRR13927737#TTGCTATAGATGTTCC-1
  SRR13927737#GCTCACTAGTTAGCGG-1
colData names(21): BlacklistRatio nDiFrags ... ReadsInPeaks FRIP

2023-11-20 06:03:32.713342 : Completed GeneIntegrationMatrix : SRR13927737(1 of 4), 0.447 mins elapsed.
2023-11-20 06:03:32.713995 : Reading GeneIntegrationMatrix : SRR13927738(2 of 4), 0.448 mins elapsed.

2023-11-20 06:03:32.727944 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927738 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927738.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3736
[1] "SRR13927738#TAGCATGAGGCTAAAT-1" "SRR13927738#TGTTAGGTCGGTCAGC-1"
[3] "SRR13927738#CTCTCAGCATGTCCCT-1" "SRR13927738#ACTATTCGTCCAACCG-1"
[5] "SRR13927738#TGGTCCTTCTACTGCC-1" "SRR13927738#ATTGTCTAGACTAGGC-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-20 06:03:32.936903 : featureDF SRR13927738, Class = DFrame
DataFrame with 23525 rows and 6 columns
      seqnames     start       end    strand        name       idx
         <Rle> <integer> <integer> <integer> <character> <integer>
1         chr1     69091     70008         1       OR4F5         1
2         chr1    817371    819837         1      FAM87B         2
3         chr1    825138    859446         1   LINC01128         3
4         chr1    827522    826206         2   LINC00115         4
5         chr1    876903    868071         2      FAM41C         5
...        ...       ...       ...       ...         ...       ...
23521     chrX 155264589 155258241         2      RAB39B       923
23522     chrX 155334657 155276211         2       CLIC2       924
23523     chrX 155380787 155381134         1      H2AFB1       925
23524     chrX 155382115 155383230         1        F8A2       926
23525     chrX 155669944 155490115         2       TMLHE       927

2023-11-20 06:03:32.955416 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927738.arrow, 0.004 mins elapsed.

2023-11-20 06:03:53.779294 : mat SRR13927738, Class = dgCMatrix
mat SRR13927738: nRows = 23525, nCols = 3736
mat SRR13927738: NonZeroEntries = 41898897, EntryRange = [ 9.73975135001421e-05 , 758.763031303812 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927738#TAGCATGAGGCTAAAT-1 SRR13927738#TGTTAGGTCGGTCAGC-1
[1,]                    .                                .        
[2,]                    .                                .        
[3,]                    0.002346352                      0.0084575
[4,]                    0.005836413                      .        
[5,]                    .                                .        
     SRR13927738#CTCTCAGCATGTCCCT-1 SRR13927738#ACTATTCGTCCAACCG-1
[1,]                     .                             .          
[2,]                     .                             .          
[3,]                     .                             0.004081077
[4,]                     0.02001583                    .          
[5,]                     .                             .          
     SRR13927738#TGGTCCTTCTACTGCC-1
[1,]                    .          
[2,]                    .          
[3,]                    0.007640877
[4,]                    0.015048796
[5,]                    .          


2023-11-20 06:03:54.15635 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927738.arrow, 0.357 mins elapsed.

2023-11-20 06:03:54.44609 : se SRR13927738, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 3736 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(3736): SRR13927738#TAGCATGAGGCTAAAT-1
  SRR13927738#TGTTAGGTCGGTCAGC-1 ... SRR13927738#GCTCCTATCCGTGCGA-1
  SRR13927738#GGCACGTAGTCGATAA-1
colData names(21): BlacklistRatio nDiFrags ... ReadsInPeaks FRIP

2023-11-20 06:03:54.460993 : Completed GeneIntegrationMatrix : SRR13927738(2 of 4), 0.81 mins elapsed.
2023-11-20 06:03:54.461686 : Reading GeneIntegrationMatrix : SRR13927735(3 of 4), 0.81 mins elapsed.

2023-11-20 06:03:54.474933 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927735 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927735.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3813
[1] "SRR13927735#TTATGTCTCCAGGTAT-1" "SRR13927735#TATTGCTCATCAGAAA-1"
[3] "SRR13927735#TTCGATTGTAGGGTTG-1" "SRR13927735#CATTCATTCGGATGTT-1"
[5] "SRR13927735#ACGTTAGGTCAACTGT-1" "SRR13927735#AAATGCCCAGCAATGG-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-20 06:03:54.66554 : featureDF SRR13927735, Class = DFrame
DataFrame with 23525 rows and 6 columns
      seqnames     start       end    strand        name       idx
         <Rle> <integer> <integer> <integer> <character> <integer>
1         chr1     69091     70008         1       OR4F5         1
2         chr1    817371    819837         1      FAM87B         2
3         chr1    825138    859446         1   LINC01128         3
4         chr1    827522    826206         2   LINC00115         4
5         chr1    876903    868071         2      FAM41C         5
...        ...       ...       ...       ...         ...       ...
23521     chrX 155264589 155258241         2      RAB39B       923
23522     chrX 155334657 155276211         2       CLIC2       924
23523     chrX 155380787 155381134         1      H2AFB1       925
23524     chrX 155382115 155383230         1        F8A2       926
23525     chrX 155669944 155490115         2       TMLHE       927

2023-11-20 06:03:54.681568 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927735.arrow, 0.003 mins elapsed.

2023-11-20 06:04:12.728576 : mat SRR13927735, Class = dgCMatrix
mat SRR13927735: nRows = 23525, nCols = 3813
mat SRR13927735: NonZeroEntries = 45535634, EntryRange = [ 7.20715653308254e-05 , 673.990433360078 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927735#TTATGTCTCCAGGTAT-1 SRR13927735#TATTGCTCATCAGAAA-1
[1,]                    .                              .          
[2,]                    .                              .          
[3,]                    0.070892732                    0.057009268
[4,]                    0.008800039                    0.002978408
[5,]                    .                              .          
     SRR13927735#TTCGATTGTAGGGTTG-1 SRR13927735#CATTCATTCGGATGTT-1
[1,]                     .                               .        
[2,]                     .                               .        
[3,]                     0.05523011                      0.0248301
[4,]                     .                               .        
[5,]                     .                               .        
     SRR13927735#ACGTTAGGTCAACTGT-1
[1,]                     .         
[2,]                     .         
[3,]                     0.05384571
[4,]                     0.01228688
[5,]                     .         


2023-11-20 06:04:13.182362 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927735.arrow, 0.312 mins elapsed.

2023-11-20 06:04:13.471398 : se SRR13927735, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 3813 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(3813): SRR13927735#TTATGTCTCCAGGTAT-1
  SRR13927735#TATTGCTCATCAGAAA-1 ... SRR13927735#GATCATGGTGTTGTTG-1
  SRR13927735#GTCCATCCATTTCACT-1
colData names(21): BlacklistRatio nDiFrags ... ReadsInPeaks FRIP

2023-11-20 06:04:13.477752 : Completed GeneIntegrationMatrix : SRR13927735(3 of 4), 1.127 mins elapsed.
2023-11-20 06:04:13.478457 : Reading GeneIntegrationMatrix : SRR13927736(4 of 4), 1.127 mins elapsed.

2023-11-20 06:04:13.491682 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927736 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927736.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3534
[1] "SRR13927736#GCTGTTCTCCTTGACC-1" "SRR13927736#AGTTACGAGTAGGTCG-1"
[3] "SRR13927736#AGTTTGGAGCATTCCA-1" "SRR13927736#GCACGGTTCGGCAATT-1"
[5] "SRR13927736#TCAGCTCTCTCTATTG-1" "SRR13927736#CCTTAATTCTACTTTG-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-20 06:04:13.688476 : featureDF SRR13927736, Class = DFrame
DataFrame with 23525 rows and 6 columns
      seqnames     start       end    strand        name       idx
         <Rle> <integer> <integer> <integer> <character> <integer>
1         chr1     69091     70008         1       OR4F5         1
2         chr1    817371    819837         1      FAM87B         2
3         chr1    825138    859446         1   LINC01128         3
4         chr1    827522    826206         2   LINC00115         4
5         chr1    876903    868071         2      FAM41C         5
...        ...       ...       ...       ...         ...       ...
23521     chrX 155264589 155258241         2      RAB39B       923
23522     chrX 155334657 155276211         2       CLIC2       924
23523     chrX 155380787 155381134         1      H2AFB1       925
23524     chrX 155382115 155383230         1        F8A2       926
23525     chrX 155669944 155490115         2       TMLHE       927

2023-11-20 06:04:13.705046 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927736.arrow, 0.004 mins elapsed.

2023-11-20 06:04:28.871834 : mat SRR13927736, Class = dgCMatrix
mat SRR13927736: nRows = 23525, nCols = 3534
mat SRR13927736: NonZeroEntries = 39003616, EntryRange = [ 8.77627031078293e-05 , 770.387276519725 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927736#GCTGTTCTCCTTGACC-1 SRR13927736#AGTTACGAGTAGGTCG-1
[1,]                     .                             .          
[2,]                     .                             .          
[3,]                     0.07845381                    0.028344539
[4,]                     0.01395289                    0.008499189
[5,]                     .                             .          
     SRR13927736#AGTTTGGAGCATTCCA-1 SRR13927736#GCACGGTTCGGCAATT-1
[1,]                     .                              .         
[2,]                     .                              .         
[3,]                     0.06535130                     0.06576222
[4,]                     0.01810177                     0.01315728
[5,]                     .                              .         
     SRR13927736#TCAGCTCTCTCTATTG-1
[1,]                    .          
[2,]                    .          
[3,]                    0.016229526
[4,]                    0.003862592
[5,]                    .          


2023-11-20 06:04:29.183567 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927736.arrow, 0.262 mins elapsed.

2023-11-20 06:04:29.451146 : se SRR13927736, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 3534 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(3534): SRR13927736#GCTGTTCTCCTTGACC-1
  SRR13927736#AGTTACGAGTAGGTCG-1 ... SRR13927736#TAAACCGCAAGTCTGT-1
  SRR13927736#TGCATTTGTTCTACCC-1
colData names(21): BlacklistRatio nDiFrags ... ReadsInPeaks FRIP

2023-11-20 06:04:29.457202 : Completed GeneIntegrationMatrix : SRR13927736(4 of 4), 1.393 mins elapsed.
2023-11-20 06:04:29.463722 : Organizing colData, 1.393 mins elapsed.
2023-11-20 06:04:29.617879 : Organizing rowData, 1.396 mins elapsed.
2023-11-20 06:04:29.621622 : Organizing rowRanges, 1.396 mins elapsed.
2023-11-20 06:04:29.627183 : Organizing Assays (1 of 1), 1.396 mins elapsed.
2023-11-20 06:04:36.792517 : Constructing SummarizedExperiment, 1.515 mins elapsed.
2023-11-20 06:04:39.386923 : Finished Matrix Creation, 1.559 mins elapsed.

2023-11-20 06:04:39.389701 : getMatrixFromProject Input-Parameters, Class = list

getMatrixFromProject Input-Parameters$ArchRProj: length = 1

getMatrixFromProject Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromProject Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromProject Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromProject Input-Parameters$verbose: length = 1
[1] TRUE


getMatrixFromProject Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromProject Input-Parameters$threads: length = 1
[1] 1


getMatrixFromProject Input-Parameters$logFile: length = 1
[1] "x"


2023-11-20 06:04:40.057927 : Reading PeakMatrix : SRR13927737(1 of 4), 0.011 mins elapsed.

2023-11-20 06:04:40.073146 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927737 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927737.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 4439
[1] "SRR13927737#GTCACGGGTCAACAGG-1" "SRR13927737#TGCCTGTTCTTACTCA-1"
[3] "SRR13927737#TTGCGGGTCAGGTCTA-1" "SRR13927737#GGTTGCGCAAAGCTTC-1"
[5] "SRR13927737#CTTAATCTCGCGCTGA-1" "SRR13927737#TGTACGATCGGAAAGC-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-20 06:04:40.375232 : featureDF SRR13927737, Class = DFrame
DataFrame with 210721 rows and 4 columns
       seqnames       idx     start       end
          <Rle> <integer> <integer> <integer>
1          chr1         1    804677    805177
2          chr1         2    811866    812366
3          chr1         3    817090    817590
4          chr1         4    818831    819331
5          chr1         5    826543    827043
...         ...       ...       ...       ...
210717     chrX      5777 155898525 155899025
210718     chrX      5778 155899804 155900304
210719     chrX      5779 155956608 155957108
210720     chrX      5780 155966827 155967327
210721     chrX      5781 155971625 155972125

2023-11-20 06:04:40.431509 : Getting PeakMatrix from ArrowFile : SRR13927737.arrow, 0.006 mins elapsed.

2023-11-20 06:04:52.49571 : mat SRR13927737, Class = dgCMatrix
mat SRR13927737: nRows = 210721, nCols = 4439
mat SRR13927737: NonZeroEntries = 29963541, EntryRange = [ 1 , 4 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927737#GTCACGGGTCAACAGG-1 SRR13927737#TGCCTGTTCTTACTCA-1
[1,]                              .                              .
[2,]                              .                              2
[3,]                              .                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927737#TTGCGGGTCAGGTCTA-1 SRR13927737#GGTTGCGCAAAGCTTC-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              2                              .
[4,]                              .                              .
[5,]                              4                              .
     SRR13927737#CTTAATCTCGCGCTGA-1
[1,]                              .
[2,]                              .
[3,]                              .
[4,]                              .
[5,]                              .


2023-11-20 06:04:52.729185 : Organizing SE PeakMatrix from ArrowFile : SRR13927737.arrow, 0.211 mins elapsed.

2023-11-20 06:04:53.030306 : se SRR13927737, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 4439 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(4439): SRR13927737#GTCACGGGTCAACAGG-1
  SRR13927737#TGCCTGTTCTTACTCA-1 ... SRR13927737#TTGCTATAGATGTTCC-1
  SRR13927737#GCTCACTAGTTAGCGG-1
colData names(21): BlacklistRatio nDiFrags ... ReadsInPeaks FRIP

2023-11-20 06:04:53.041356 : Completed PeakMatrix : SRR13927737(1 of 4), 0.228 mins elapsed.
2023-11-20 06:04:53.042007 : Reading PeakMatrix : SRR13927738(2 of 4), 0.228 mins elapsed.

2023-11-20 06:04:53.057036 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927738 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927738.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3736
[1] "SRR13927738#TAGCATGAGGCTAAAT-1" "SRR13927738#TGTTAGGTCGGTCAGC-1"
[3] "SRR13927738#CTCTCAGCATGTCCCT-1" "SRR13927738#ACTATTCGTCCAACCG-1"
[5] "SRR13927738#TGGTCCTTCTACTGCC-1" "SRR13927738#ATTGTCTAGACTAGGC-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-20 06:04:53.351287 : featureDF SRR13927738, Class = DFrame
DataFrame with 210721 rows and 4 columns
       seqnames       idx     start       end
          <Rle> <integer> <integer> <integer>
1          chr1         1    804677    805177
2          chr1         2    811866    812366
3          chr1         3    817090    817590
4          chr1         4    818831    819331
5          chr1         5    826543    827043
...         ...       ...       ...       ...
210717     chrX      5777 155898525 155899025
210718     chrX      5778 155899804 155900304
210719     chrX      5779 155956608 155957108
210720     chrX      5780 155966827 155967327
210721     chrX      5781 155971625 155972125

2023-11-20 06:04:53.399143 : Getting PeakMatrix from ArrowFile : SRR13927738.arrow, 0.006 mins elapsed.

2023-11-20 06:05:05.367527 : mat SRR13927738, Class = dgCMatrix
mat SRR13927738: nRows = 210721, nCols = 3736
mat SRR13927738: NonZeroEntries = 30165921, EntryRange = [ 1 , 4 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927738#TAGCATGAGGCTAAAT-1 SRR13927738#TGTTAGGTCGGTCAGC-1
[1,]                              2                              .
[2,]                              .                              .
[3,]                              .                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927738#CTCTCAGCATGTCCCT-1 SRR13927738#ACTATTCGTCCAACCG-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              .                              .
[4,]                              .                              .
[5,]                              2                              .
     SRR13927738#TGGTCCTTCTACTGCC-1
[1,]                              .
[2,]                              .
[3,]                              2
[4,]                              .
[5,]                              .


2023-11-20 06:05:05.603421 : Organizing SE PeakMatrix from ArrowFile : SRR13927738.arrow, 0.209 mins elapsed.

2023-11-20 06:05:05.9255 : se SRR13927738, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 3736 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(3736): SRR13927738#TAGCATGAGGCTAAAT-1
  SRR13927738#TGTTAGGTCGGTCAGC-1 ... SRR13927738#GCTCCTATCCGTGCGA-1
  SRR13927738#GGCACGTAGTCGATAA-1
colData names(21): BlacklistRatio nDiFrags ... ReadsInPeaks FRIP

2023-11-20 06:05:05.933712 : Completed PeakMatrix : SRR13927738(2 of 4), 0.442 mins elapsed.
2023-11-20 06:05:05.934405 : Reading PeakMatrix : SRR13927735(3 of 4), 0.442 mins elapsed.

2023-11-20 06:05:05.953171 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927735 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927735.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3813
[1] "SRR13927735#TTATGTCTCCAGGTAT-1" "SRR13927735#TATTGCTCATCAGAAA-1"
[3] "SRR13927735#TTCGATTGTAGGGTTG-1" "SRR13927735#CATTCATTCGGATGTT-1"
[5] "SRR13927735#ACGTTAGGTCAACTGT-1" "SRR13927735#AAATGCCCAGCAATGG-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-20 06:05:06.250616 : featureDF SRR13927735, Class = DFrame
DataFrame with 210721 rows and 4 columns
       seqnames       idx     start       end
          <Rle> <integer> <integer> <integer>
1          chr1         1    804677    805177
2          chr1         2    811866    812366
3          chr1         3    817090    817590
4          chr1         4    818831    819331
5          chr1         5    826543    827043
...         ...       ...       ...       ...
210717     chrX      5777 155898525 155899025
210718     chrX      5778 155899804 155900304
210719     chrX      5779 155956608 155957108
210720     chrX      5780 155966827 155967327
210721     chrX      5781 155971625 155972125

2023-11-20 06:05:06.298104 : Getting PeakMatrix from ArrowFile : SRR13927735.arrow, 0.006 mins elapsed.

2023-11-20 06:05:17.056483 : mat SRR13927735, Class = dgCMatrix
mat SRR13927735: nRows = 210721, nCols = 3813
mat SRR13927735: NonZeroEntries = 26038084, EntryRange = [ 1 , 4 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927735#TTATGTCTCCAGGTAT-1 SRR13927735#TATTGCTCATCAGAAA-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              1                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927735#TTCGATTGTAGGGTTG-1 SRR13927735#CATTCATTCGGATGTT-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              .                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927735#ACGTTAGGTCAACTGT-1
[1,]                              2
[2,]                              .
[3,]                              .
[4,]                              .
[5,]                              .


2023-11-20 06:05:17.230603 : Organizing SE PeakMatrix from ArrowFile : SRR13927735.arrow, 0.188 mins elapsed.

2023-11-20 06:05:17.544396 : se SRR13927735, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 3813 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(3813): SRR13927735#TTATGTCTCCAGGTAT-1
  SRR13927735#TATTGCTCATCAGAAA-1 ... SRR13927735#GATCATGGTGTTGTTG-1
  SRR13927735#GTCCATCCATTTCACT-1
colData names(21): BlacklistRatio nDiFrags ... ReadsInPeaks FRIP

2023-11-20 06:05:17.552302 : Completed PeakMatrix : SRR13927735(3 of 4), 0.636 mins elapsed.
2023-11-20 06:05:17.55298 : Reading PeakMatrix : SRR13927736(4 of 4), 0.636 mins elapsed.

2023-11-20 06:05:17.566419 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927736 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927736.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3534
[1] "SRR13927736#GCTGTTCTCCTTGACC-1" "SRR13927736#AGTTACGAGTAGGTCG-1"
[3] "SRR13927736#AGTTTGGAGCATTCCA-1" "SRR13927736#GCACGGTTCGGCAATT-1"
[5] "SRR13927736#TCAGCTCTCTCTATTG-1" "SRR13927736#CCTTAATTCTACTTTG-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-20 06:05:17.853637 : featureDF SRR13927736, Class = DFrame
DataFrame with 210721 rows and 4 columns
       seqnames       idx     start       end
          <Rle> <integer> <integer> <integer>
1          chr1         1    804677    805177
2          chr1         2    811866    812366
3          chr1         3    817090    817590
4          chr1         4    818831    819331
5          chr1         5    826543    827043
...         ...       ...       ...       ...
210717     chrX      5777 155898525 155899025
210718     chrX      5778 155899804 155900304
210719     chrX      5779 155956608 155957108
210720     chrX      5780 155966827 155967327
210721     chrX      5781 155971625 155972125

2023-11-20 06:05:17.901064 : Getting PeakMatrix from ArrowFile : SRR13927736.arrow, 0.006 mins elapsed.

2023-11-20 06:05:28.409052 : mat SRR13927736, Class = dgCMatrix
mat SRR13927736: nRows = 210721, nCols = 3534
mat SRR13927736: NonZeroEntries = 25583439, EntryRange = [ 1 , 4 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927736#GCTGTTCTCCTTGACC-1 SRR13927736#AGTTACGAGTAGGTCG-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              .                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927736#AGTTTGGAGCATTCCA-1 SRR13927736#GCACGGTTCGGCAATT-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              1                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927736#TCAGCTCTCTCTATTG-1
[1,]                              .
[2,]                              .
[3,]                              .
[4,]                              .
[5,]                              .


2023-11-20 06:05:28.655134 : Organizing SE PeakMatrix from ArrowFile : SRR13927736.arrow, 0.185 mins elapsed.

2023-11-20 06:05:28.968614 : se SRR13927736, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 3534 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(3534): SRR13927736#GCTGTTCTCCTTGACC-1
  SRR13927736#AGTTACGAGTAGGTCG-1 ... SRR13927736#TAAACCGCAAGTCTGT-1
  SRR13927736#TGCATTTGTTCTACCC-1
colData names(21): BlacklistRatio nDiFrags ... ReadsInPeaks FRIP

2023-11-20 06:05:28.976268 : Completed PeakMatrix : SRR13927736(4 of 4), 0.826 mins elapsed.
2023-11-20 06:05:28.977615 : Organizing colData, 0.826 mins elapsed.
2023-11-20 06:05:29.124874 : Organizing rowData, 0.829 mins elapsed.
2023-11-20 06:05:29.133808 : Organizing rowRanges, 0.829 mins elapsed.
2023-11-20 06:05:29.145468 : Organizing Assays (1 of 1), 0.829 mins elapsed.
2023-11-20 06:05:31.819311 : Constructing SummarizedExperiment, 0.874 mins elapsed.
2023-11-20 06:05:53.437045 : Finished Matrix Creation, 1.234 mins elapsed.

2023-11-21 03:25:14.375371 : getMatrixFromProject Input-Parameters, Class = list

getMatrixFromProject Input-Parameters$ArchRProj: length = 1

getMatrixFromProject Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromProject Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromProject Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromProject Input-Parameters$verbose: length = 1
[1] TRUE


getMatrixFromProject Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromProject Input-Parameters$threads: length = 1
[1] 1


getMatrixFromProject Input-Parameters$logFile: length = 1
[1] "x"


2023-11-21 03:25:15.034003 : Reading GeneIntegrationMatrix : SRR13927737(1 of 4), 0.011 mins elapsed.

2023-11-21 03:25:15.058499 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927737 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927737.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 4439
[1] "SRR13927737#GTCACGGGTCAACAGG-1" "SRR13927737#TGCCTGTTCTTACTCA-1"
[3] "SRR13927737#TTGCGGGTCAGGTCTA-1" "SRR13927737#GGTTGCGCAAAGCTTC-1"
[5] "SRR13927737#CTTAATCTCGCGCTGA-1" "SRR13927737#TGTACGATCGGAAAGC-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-21 03:25:15.387388 : featureDF SRR13927737, Class = DFrame
DataFrame with 23525 rows and 6 columns
      seqnames     start       end    strand        name       idx
         <Rle> <integer> <integer> <integer> <character> <integer>
1         chr1     69091     70008         1       OR4F5         1
2         chr1    817371    819837         1      FAM87B         2
3         chr1    825138    859446         1   LINC01128         3
4         chr1    827522    826206         2   LINC00115         4
5         chr1    876903    868071         2      FAM41C         5
...        ...       ...       ...       ...         ...       ...
23521     chrX 155264589 155258241         2      RAB39B       923
23522     chrX 155334657 155276211         2       CLIC2       924
23523     chrX 155380787 155381134         1      H2AFB1       925
23524     chrX 155382115 155383230         1        F8A2       926
23525     chrX 155669944 155490115         2       TMLHE       927

2023-11-21 03:25:15.411147 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927737.arrow, 0.006 mins elapsed.

2023-11-21 03:25:43.41652 : mat SRR13927737, Class = dgCMatrix
mat SRR13927737: nRows = 23525, nCols = 4439
mat SRR13927737: NonZeroEntries = 50277774, EntryRange = [ 0.000148130054904683 , 376.345870564644 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927737#GTCACGGGTCAACAGG-1 SRR13927737#TGCCTGTTCTTACTCA-1
[1,]                      .                             .         
[2,]                      .                             .         
[3,]                      0.1057576                     0.09630696
[4,]                      .                             0.01806508
[5,]                      .                             .         
     SRR13927737#TTGCGGGTCAGGTCTA-1 SRR13927737#GGTTGCGCAAAGCTTC-1
[1,]                     .                              .         
[2,]                     .                              .         
[3,]                     0.09113499                     0.08396973
[4,]                     0.02096304                     0.01727603
[5,]                     .                              .         
     SRR13927737#CTTAATCTCGCGCTGA-1
[1,]                     .         
[2,]                     .         
[3,]                     0.08327649
[4,]                     0.01368917
[5,]                     .         


2023-11-21 03:25:43.964872 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927737.arrow, 0.482 mins elapsed.

2023-11-21 03:25:44.346436 : se SRR13927737, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 4439 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(4439): SRR13927737#GTCACGGGTCAACAGG-1
  SRR13927737#TGCCTGTTCTTACTCA-1 ... SRR13927737#TTGCTATAGATGTTCC-1
  SRR13927737#GCTCACTAGTTAGCGG-1
colData names(21): BlacklistRatio nDiFrags ... ReadsInPeaks FRIP

2023-11-21 03:25:44.353371 : Completed GeneIntegrationMatrix : SRR13927737(1 of 4), 0.5 mins elapsed.
2023-11-21 03:25:44.354074 : Reading GeneIntegrationMatrix : SRR13927738(2 of 4), 0.5 mins elapsed.

2023-11-21 03:25:44.367218 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927738 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927738.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3736
[1] "SRR13927738#TAGCATGAGGCTAAAT-1" "SRR13927738#TGTTAGGTCGGTCAGC-1"
[3] "SRR13927738#CTCTCAGCATGTCCCT-1" "SRR13927738#ACTATTCGTCCAACCG-1"
[5] "SRR13927738#TGGTCCTTCTACTGCC-1" "SRR13927738#ATTGTCTAGACTAGGC-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-21 03:25:44.55618 : featureDF SRR13927738, Class = DFrame
DataFrame with 23525 rows and 6 columns
      seqnames     start       end    strand        name       idx
         <Rle> <integer> <integer> <integer> <character> <integer>
1         chr1     69091     70008         1       OR4F5         1
2         chr1    817371    819837         1      FAM87B         2
3         chr1    825138    859446         1   LINC01128         3
4         chr1    827522    826206         2   LINC00115         4
5         chr1    876903    868071         2      FAM41C         5
...        ...       ...       ...       ...         ...       ...
23521     chrX 155264589 155258241         2      RAB39B       923
23522     chrX 155334657 155276211         2       CLIC2       924
23523     chrX 155380787 155381134         1      H2AFB1       925
23524     chrX 155382115 155383230         1        F8A2       926
23525     chrX 155669944 155490115         2       TMLHE       927

2023-11-21 03:25:44.572918 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927738.arrow, 0.003 mins elapsed.

2023-11-21 03:26:05.792476 : mat SRR13927738, Class = dgCMatrix
mat SRR13927738: nRows = 23525, nCols = 3736
mat SRR13927738: NonZeroEntries = 41898897, EntryRange = [ 9.73975135001421e-05 , 758.763031303812 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927738#TAGCATGAGGCTAAAT-1 SRR13927738#TGTTAGGTCGGTCAGC-1
[1,]                    .                                .        
[2,]                    .                                .        
[3,]                    0.002346352                      0.0084575
[4,]                    0.005836413                      .        
[5,]                    .                                .        
     SRR13927738#CTCTCAGCATGTCCCT-1 SRR13927738#ACTATTCGTCCAACCG-1
[1,]                     .                             .          
[2,]                     .                             .          
[3,]                     .                             0.004081077
[4,]                     0.02001583                    .          
[5,]                     .                             .          
     SRR13927738#TGGTCCTTCTACTGCC-1
[1,]                    .          
[2,]                    .          
[3,]                    0.007640877
[4,]                    0.015048796
[5,]                    .          


2023-11-21 03:26:06.790455 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927738.arrow, 0.374 mins elapsed.

2023-11-21 03:26:07.067214 : se SRR13927738, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 3736 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(3736): SRR13927738#TAGCATGAGGCTAAAT-1
  SRR13927738#TGTTAGGTCGGTCAGC-1 ... SRR13927738#GCTCCTATCCGTGCGA-1
  SRR13927738#GGCACGTAGTCGATAA-1
colData names(21): BlacklistRatio nDiFrags ... ReadsInPeaks FRIP

2023-11-21 03:26:07.072843 : Completed GeneIntegrationMatrix : SRR13927738(2 of 4), 0.879 mins elapsed.
2023-11-21 03:26:07.073501 : Reading GeneIntegrationMatrix : SRR13927735(3 of 4), 0.879 mins elapsed.

2023-11-21 03:26:07.086785 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927735 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927735.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3813
[1] "SRR13927735#TTATGTCTCCAGGTAT-1" "SRR13927735#TATTGCTCATCAGAAA-1"
[3] "SRR13927735#TTCGATTGTAGGGTTG-1" "SRR13927735#CATTCATTCGGATGTT-1"
[5] "SRR13927735#ACGTTAGGTCAACTGT-1" "SRR13927735#AAATGCCCAGCAATGG-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-21 03:26:07.279948 : featureDF SRR13927735, Class = DFrame
DataFrame with 23525 rows and 6 columns
      seqnames     start       end    strand        name       idx
         <Rle> <integer> <integer> <integer> <character> <integer>
1         chr1     69091     70008         1       OR4F5         1
2         chr1    817371    819837         1      FAM87B         2
3         chr1    825138    859446         1   LINC01128         3
4         chr1    827522    826206         2   LINC00115         4
5         chr1    876903    868071         2      FAM41C         5
...        ...       ...       ...       ...         ...       ...
23521     chrX 155264589 155258241         2      RAB39B       923
23522     chrX 155334657 155276211         2       CLIC2       924
23523     chrX 155380787 155381134         1      H2AFB1       925
23524     chrX 155382115 155383230         1        F8A2       926
23525     chrX 155669944 155490115         2       TMLHE       927

2023-11-21 03:26:07.297044 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927735.arrow, 0.004 mins elapsed.

2023-11-21 03:26:30.25996 : mat SRR13927735, Class = dgCMatrix
mat SRR13927735: nRows = 23525, nCols = 3813
mat SRR13927735: NonZeroEntries = 45535634, EntryRange = [ 7.20715653308254e-05 , 673.990433360078 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927735#TTATGTCTCCAGGTAT-1 SRR13927735#TATTGCTCATCAGAAA-1
[1,]                    .                              .          
[2,]                    .                              .          
[3,]                    0.070892732                    0.057009268
[4,]                    0.008800039                    0.002978408
[5,]                    .                              .          
     SRR13927735#TTCGATTGTAGGGTTG-1 SRR13927735#CATTCATTCGGATGTT-1
[1,]                     .                               .        
[2,]                     .                               .        
[3,]                     0.05523011                      0.0248301
[4,]                     .                               .        
[5,]                     .                               .        
     SRR13927735#ACGTTAGGTCAACTGT-1
[1,]                     .         
[2,]                     .         
[3,]                     0.05384571
[4,]                     0.01228688
[5,]                     .         


2023-11-21 03:26:30.711076 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927735.arrow, 0.394 mins elapsed.

2023-11-21 03:26:30.977214 : se SRR13927735, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 3813 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(3813): SRR13927735#TTATGTCTCCAGGTAT-1
  SRR13927735#TATTGCTCATCAGAAA-1 ... SRR13927735#GATCATGGTGTTGTTG-1
  SRR13927735#GTCCATCCATTTCACT-1
colData names(21): BlacklistRatio nDiFrags ... ReadsInPeaks FRIP

2023-11-21 03:26:30.982827 : Completed GeneIntegrationMatrix : SRR13927735(3 of 4), 1.277 mins elapsed.
2023-11-21 03:26:30.983464 : Reading GeneIntegrationMatrix : SRR13927736(4 of 4), 1.277 mins elapsed.

2023-11-21 03:26:30.995954 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927736 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927736.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3534
[1] "SRR13927736#GCTGTTCTCCTTGACC-1" "SRR13927736#AGTTACGAGTAGGTCG-1"
[3] "SRR13927736#AGTTTGGAGCATTCCA-1" "SRR13927736#GCACGGTTCGGCAATT-1"
[5] "SRR13927736#TCAGCTCTCTCTATTG-1" "SRR13927736#CCTTAATTCTACTTTG-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-21 03:26:31.179098 : featureDF SRR13927736, Class = DFrame
DataFrame with 23525 rows and 6 columns
      seqnames     start       end    strand        name       idx
         <Rle> <integer> <integer> <integer> <character> <integer>
1         chr1     69091     70008         1       OR4F5         1
2         chr1    817371    819837         1      FAM87B         2
3         chr1    825138    859446         1   LINC01128         3
4         chr1    827522    826206         2   LINC00115         4
5         chr1    876903    868071         2      FAM41C         5
...        ...       ...       ...       ...         ...       ...
23521     chrX 155264589 155258241         2      RAB39B       923
23522     chrX 155334657 155276211         2       CLIC2       924
23523     chrX 155380787 155381134         1      H2AFB1       925
23524     chrX 155382115 155383230         1        F8A2       926
23525     chrX 155669944 155490115         2       TMLHE       927

2023-11-21 03:26:31.195157 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927736.arrow, 0.003 mins elapsed.

2023-11-21 03:26:50.309563 : mat SRR13927736, Class = dgCMatrix
mat SRR13927736: nRows = 23525, nCols = 3534
mat SRR13927736: NonZeroEntries = 39003616, EntryRange = [ 8.77627031078293e-05 , 770.387276519725 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927736#GCTGTTCTCCTTGACC-1 SRR13927736#AGTTACGAGTAGGTCG-1
[1,]                     .                             .          
[2,]                     .                             .          
[3,]                     0.07845381                    0.028344539
[4,]                     0.01395289                    0.008499189
[5,]                     .                             .          
     SRR13927736#AGTTTGGAGCATTCCA-1 SRR13927736#GCACGGTTCGGCAATT-1
[1,]                     .                              .         
[2,]                     .                              .         
[3,]                     0.06535130                     0.06576222
[4,]                     0.01810177                     0.01315728
[5,]                     .                              .         
     SRR13927736#TCAGCTCTCTCTATTG-1
[1,]                    .          
[2,]                    .          
[3,]                    0.016229526
[4,]                    0.003862592
[5,]                    .          


2023-11-21 03:26:50.691432 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927736.arrow, 0.328 mins elapsed.

2023-11-21 03:26:50.950267 : se SRR13927736, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 3534 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(3534): SRR13927736#GCTGTTCTCCTTGACC-1
  SRR13927736#AGTTACGAGTAGGTCG-1 ... SRR13927736#TAAACCGCAAGTCTGT-1
  SRR13927736#TGCATTTGTTCTACCC-1
colData names(21): BlacklistRatio nDiFrags ... ReadsInPeaks FRIP

2023-11-21 03:26:50.955755 : Completed GeneIntegrationMatrix : SRR13927736(4 of 4), 1.61 mins elapsed.
2023-11-21 03:26:50.96148 : Organizing colData, 1.61 mins elapsed.
2023-11-21 03:26:51.108646 : Organizing rowData, 1.613 mins elapsed.
2023-11-21 03:26:51.112044 : Organizing rowRanges, 1.613 mins elapsed.
2023-11-21 03:26:51.117186 : Organizing Assays (1 of 1), 1.613 mins elapsed.
2023-11-21 03:26:58.954271 : Constructing SummarizedExperiment, 1.743 mins elapsed.
2023-11-21 03:27:01.681513 : Finished Matrix Creation, 1.789 mins elapsed.

2023-11-21 03:27:01.684143 : getMatrixFromProject Input-Parameters, Class = list

getMatrixFromProject Input-Parameters$ArchRProj: length = 1

getMatrixFromProject Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromProject Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromProject Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromProject Input-Parameters$verbose: length = 1
[1] TRUE


getMatrixFromProject Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromProject Input-Parameters$threads: length = 1
[1] 1


getMatrixFromProject Input-Parameters$logFile: length = 1
[1] "x"


2023-11-21 03:27:02.256907 : Reading PeakMatrix : SRR13927737(1 of 4), 0.01 mins elapsed.

2023-11-21 03:27:02.271931 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927737 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927737.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 4439
[1] "SRR13927737#GTCACGGGTCAACAGG-1" "SRR13927737#TGCCTGTTCTTACTCA-1"
[3] "SRR13927737#TTGCGGGTCAGGTCTA-1" "SRR13927737#GGTTGCGCAAAGCTTC-1"
[5] "SRR13927737#CTTAATCTCGCGCTGA-1" "SRR13927737#TGTACGATCGGAAAGC-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-21 03:27:02.579271 : featureDF SRR13927737, Class = DFrame
DataFrame with 210721 rows and 4 columns
       seqnames       idx     start       end
          <Rle> <integer> <integer> <integer>
1          chr1         1    804677    805177
2          chr1         2    811866    812366
3          chr1         3    817090    817590
4          chr1         4    818831    819331
5          chr1         5    826543    827043
...         ...       ...       ...       ...
210717     chrX      5777 155898525 155899025
210718     chrX      5778 155899804 155900304
210719     chrX      5779 155956608 155957108
210720     chrX      5780 155966827 155967327
210721     chrX      5781 155971625 155972125

2023-11-21 03:27:02.634441 : Getting PeakMatrix from ArrowFile : SRR13927737.arrow, 0.006 mins elapsed.

2023-11-21 03:27:16.525204 : mat SRR13927737, Class = dgCMatrix
mat SRR13927737: nRows = 210721, nCols = 4439
mat SRR13927737: NonZeroEntries = 29963541, EntryRange = [ 1 , 4 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927737#GTCACGGGTCAACAGG-1 SRR13927737#TGCCTGTTCTTACTCA-1
[1,]                              .                              .
[2,]                              .                              2
[3,]                              .                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927737#TTGCGGGTCAGGTCTA-1 SRR13927737#GGTTGCGCAAAGCTTC-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              2                              .
[4,]                              .                              .
[5,]                              4                              .
     SRR13927737#CTTAATCTCGCGCTGA-1
[1,]                              .
[2,]                              .
[3,]                              .
[4,]                              .
[5,]                              .


2023-11-21 03:27:16.822167 : Organizing SE PeakMatrix from ArrowFile : SRR13927737.arrow, 0.243 mins elapsed.

2023-11-21 03:27:17.120456 : se SRR13927737, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 4439 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(4439): SRR13927737#GTCACGGGTCAACAGG-1
  SRR13927737#TGCCTGTTCTTACTCA-1 ... SRR13927737#TTGCTATAGATGTTCC-1
  SRR13927737#GCTCACTAGTTAGCGG-1
colData names(21): BlacklistRatio nDiFrags ... ReadsInPeaks FRIP

2023-11-21 03:27:17.130763 : Completed PeakMatrix : SRR13927737(1 of 4), 0.257 mins elapsed.
2023-11-21 03:27:17.131412 : Reading PeakMatrix : SRR13927738(2 of 4), 0.257 mins elapsed.

2023-11-21 03:27:17.144051 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927738 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927738.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3736
[1] "SRR13927738#TAGCATGAGGCTAAAT-1" "SRR13927738#TGTTAGGTCGGTCAGC-1"
[3] "SRR13927738#CTCTCAGCATGTCCCT-1" "SRR13927738#ACTATTCGTCCAACCG-1"
[5] "SRR13927738#TGGTCCTTCTACTGCC-1" "SRR13927738#ATTGTCTAGACTAGGC-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-21 03:27:17.418992 : featureDF SRR13927738, Class = DFrame
DataFrame with 210721 rows and 4 columns
       seqnames       idx     start       end
          <Rle> <integer> <integer> <integer>
1          chr1         1    804677    805177
2          chr1         2    811866    812366
3          chr1         3    817090    817590
4          chr1         4    818831    819331
5          chr1         5    826543    827043
...         ...       ...       ...       ...
210717     chrX      5777 155898525 155899025
210718     chrX      5778 155899804 155900304
210719     chrX      5779 155956608 155957108
210720     chrX      5780 155966827 155967327
210721     chrX      5781 155971625 155972125

2023-11-21 03:27:17.465513 : Getting PeakMatrix from ArrowFile : SRR13927738.arrow, 0.005 mins elapsed.

2023-11-21 03:27:31.525353 : mat SRR13927738, Class = dgCMatrix
mat SRR13927738: nRows = 210721, nCols = 3736
mat SRR13927738: NonZeroEntries = 30165921, EntryRange = [ 1 , 4 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927738#TAGCATGAGGCTAAAT-1 SRR13927738#TGTTAGGTCGGTCAGC-1
[1,]                              2                              .
[2,]                              .                              .
[3,]                              .                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927738#CTCTCAGCATGTCCCT-1 SRR13927738#ACTATTCGTCCAACCG-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              .                              .
[4,]                              .                              .
[5,]                              2                              .
     SRR13927738#TGGTCCTTCTACTGCC-1
[1,]                              .
[2,]                              .
[3,]                              2
[4,]                              .
[5,]                              .


2023-11-21 03:27:31.824317 : Organizing SE PeakMatrix from ArrowFile : SRR13927738.arrow, 0.245 mins elapsed.

2023-11-21 03:27:32.128091 : se SRR13927738, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 3736 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(3736): SRR13927738#TAGCATGAGGCTAAAT-1
  SRR13927738#TGTTAGGTCGGTCAGC-1 ... SRR13927738#GCTCCTATCCGTGCGA-1
  SRR13927738#GGCACGTAGTCGATAA-1
colData names(21): BlacklistRatio nDiFrags ... ReadsInPeaks FRIP

2023-11-21 03:27:32.135894 : Completed PeakMatrix : SRR13927738(2 of 4), 0.508 mins elapsed.
2023-11-21 03:27:32.136565 : Reading PeakMatrix : SRR13927735(3 of 4), 0.508 mins elapsed.

2023-11-21 03:27:32.149765 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927735 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927735.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3813
[1] "SRR13927735#TTATGTCTCCAGGTAT-1" "SRR13927735#TATTGCTCATCAGAAA-1"
[3] "SRR13927735#TTCGATTGTAGGGTTG-1" "SRR13927735#CATTCATTCGGATGTT-1"
[5] "SRR13927735#ACGTTAGGTCAACTGT-1" "SRR13927735#AAATGCCCAGCAATGG-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-21 03:27:32.429347 : featureDF SRR13927735, Class = DFrame
DataFrame with 210721 rows and 4 columns
       seqnames       idx     start       end
          <Rle> <integer> <integer> <integer>
1          chr1         1    804677    805177
2          chr1         2    811866    812366
3          chr1         3    817090    817590
4          chr1         4    818831    819331
5          chr1         5    826543    827043
...         ...       ...       ...       ...
210717     chrX      5777 155898525 155899025
210718     chrX      5778 155899804 155900304
210719     chrX      5779 155956608 155957108
210720     chrX      5780 155966827 155967327
210721     chrX      5781 155971625 155972125

2023-11-21 03:27:32.476429 : Getting PeakMatrix from ArrowFile : SRR13927735.arrow, 0.005 mins elapsed.

2023-11-21 03:27:46.074732 : mat SRR13927735, Class = dgCMatrix
mat SRR13927735: nRows = 210721, nCols = 3813
mat SRR13927735: NonZeroEntries = 26038084, EntryRange = [ 1 , 4 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927735#TTATGTCTCCAGGTAT-1 SRR13927735#TATTGCTCATCAGAAA-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              1                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927735#TTCGATTGTAGGGTTG-1 SRR13927735#CATTCATTCGGATGTT-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              .                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927735#ACGTTAGGTCAACTGT-1
[1,]                              2
[2,]                              .
[3,]                              .
[4,]                              .
[5,]                              .


2023-11-21 03:27:46.249349 : Organizing SE PeakMatrix from ArrowFile : SRR13927735.arrow, 0.235 mins elapsed.

2023-11-21 03:27:46.574306 : se SRR13927735, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 3813 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(3813): SRR13927735#TTATGTCTCCAGGTAT-1
  SRR13927735#TATTGCTCATCAGAAA-1 ... SRR13927735#GATCATGGTGTTGTTG-1
  SRR13927735#GTCCATCCATTTCACT-1
colData names(21): BlacklistRatio nDiFrags ... ReadsInPeaks FRIP

2023-11-21 03:27:46.582173 : Completed PeakMatrix : SRR13927735(3 of 4), 0.748 mins elapsed.
2023-11-21 03:27:46.582868 : Reading PeakMatrix : SRR13927736(4 of 4), 0.748 mins elapsed.

2023-11-21 03:27:46.595871 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927736 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927736.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3534
[1] "SRR13927736#GCTGTTCTCCTTGACC-1" "SRR13927736#AGTTACGAGTAGGTCG-1"
[3] "SRR13927736#AGTTTGGAGCATTCCA-1" "SRR13927736#GCACGGTTCGGCAATT-1"
[5] "SRR13927736#TCAGCTCTCTCTATTG-1" "SRR13927736#CCTTAATTCTACTTTG-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-21 03:27:46.905879 : featureDF SRR13927736, Class = DFrame
DataFrame with 210721 rows and 4 columns
       seqnames       idx     start       end
          <Rle> <integer> <integer> <integer>
1          chr1         1    804677    805177
2          chr1         2    811866    812366
3          chr1         3    817090    817590
4          chr1         4    818831    819331
5          chr1         5    826543    827043
...         ...       ...       ...       ...
210717     chrX      5777 155898525 155899025
210718     chrX      5778 155899804 155900304
210719     chrX      5779 155956608 155957108
210720     chrX      5780 155966827 155967327
210721     chrX      5781 155971625 155972125

2023-11-21 03:27:46.952953 : Getting PeakMatrix from ArrowFile : SRR13927736.arrow, 0.006 mins elapsed.

2023-11-21 03:27:58.044644 : mat SRR13927736, Class = dgCMatrix
mat SRR13927736: nRows = 210721, nCols = 3534
mat SRR13927736: NonZeroEntries = 25583439, EntryRange = [ 1 , 4 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927736#GCTGTTCTCCTTGACC-1 SRR13927736#AGTTACGAGTAGGTCG-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              .                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927736#AGTTTGGAGCATTCCA-1 SRR13927736#GCACGGTTCGGCAATT-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              1                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927736#TCAGCTCTCTCTATTG-1
[1,]                              .
[2,]                              .
[3,]                              .
[4,]                              .
[5,]                              .


2023-11-21 03:27:58.302411 : Organizing SE PeakMatrix from ArrowFile : SRR13927736.arrow, 0.195 mins elapsed.

2023-11-21 03:27:58.598793 : se SRR13927736, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 3534 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(3534): SRR13927736#GCTGTTCTCCTTGACC-1
  SRR13927736#AGTTACGAGTAGGTCG-1 ... SRR13927736#TAAACCGCAAGTCTGT-1
  SRR13927736#TGCATTTGTTCTACCC-1
colData names(21): BlacklistRatio nDiFrags ... ReadsInPeaks FRIP

2023-11-21 03:27:58.606504 : Completed PeakMatrix : SRR13927736(4 of 4), 0.949 mins elapsed.
2023-11-21 03:27:58.607841 : Organizing colData, 0.949 mins elapsed.
2023-11-21 03:27:58.755675 : Organizing rowData, 0.951 mins elapsed.
2023-11-21 03:27:58.76477 : Organizing rowRanges, 0.951 mins elapsed.
2023-11-21 03:27:58.776977 : Organizing Assays (1 of 1), 0.952 mins elapsed.
2023-11-21 03:28:01.851628 : Constructing SummarizedExperiment, 1.003 mins elapsed.
2023-11-21 03:28:27.021569 : Finished Matrix Creation, 1.422 mins elapsed.

2023-11-21 05:16:30.697845 : getMatrixFromProject Input-Parameters, Class = list

getMatrixFromProject Input-Parameters$ArchRProj: length = 1

getMatrixFromProject Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromProject Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromProject Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromProject Input-Parameters$verbose: length = 1
[1] TRUE


getMatrixFromProject Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromProject Input-Parameters$threads: length = 1
[1] 1


getMatrixFromProject Input-Parameters$logFile: length = 1
[1] "x"


2023-11-21 05:16:32.098408 : Reading GeneIntegrationMatrix : SRR13927737(1 of 4), 0.023 mins elapsed.

2023-11-21 05:16:32.132706 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927737 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927737.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 4439
[1] "SRR13927737#GTCACGGGTCAACAGG-1" "SRR13927737#TGCCTGTTCTTACTCA-1"
[3] "SRR13927737#TTGCGGGTCAGGTCTA-1" "SRR13927737#GGTTGCGCAAAGCTTC-1"
[5] "SRR13927737#CTTAATCTCGCGCTGA-1" "SRR13927737#TGTACGATCGGAAAGC-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-21 05:16:32.660944 : featureDF SRR13927737, Class = DFrame
DataFrame with 23525 rows and 6 columns
      seqnames     start       end    strand        name       idx
         <Rle> <integer> <integer> <integer> <character> <integer>
1         chr1     69091     70008         1       OR4F5         1
2         chr1    817371    819837         1      FAM87B         2
3         chr1    825138    859446         1   LINC01128         3
4         chr1    827522    826206         2   LINC00115         4
5         chr1    876903    868071         2      FAM41C         5
...        ...       ...       ...       ...         ...       ...
23521     chrX 155264589 155258241         2      RAB39B       923
23522     chrX 155334657 155276211         2       CLIC2       924
23523     chrX 155380787 155381134         1      H2AFB1       925
23524     chrX 155382115 155383230         1        F8A2       926
23525     chrX 155669944 155490115         2       TMLHE       927

2023-11-21 05:16:32.700668 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927737.arrow, 0.009 mins elapsed.

2023-11-21 05:16:57.556199 : mat SRR13927737, Class = dgCMatrix
mat SRR13927737: nRows = 23525, nCols = 4439
mat SRR13927737: NonZeroEntries = 50277774, EntryRange = [ 0.000148130054904683 , 376.345870564644 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927737#GTCACGGGTCAACAGG-1 SRR13927737#TGCCTGTTCTTACTCA-1
[1,]                      .                             .         
[2,]                      .                             .         
[3,]                      0.1057576                     0.09630696
[4,]                      .                             0.01806508
[5,]                      .                             .         
     SRR13927737#TTGCGGGTCAGGTCTA-1 SRR13927737#GGTTGCGCAAAGCTTC-1
[1,]                     .                              .         
[2,]                     .                              .         
[3,]                     0.09113499                     0.08396973
[4,]                     0.02096304                     0.01727603
[5,]                     .                              .         
     SRR13927737#CTTAATCTCGCGCTGA-1
[1,]                     .         
[2,]                     .         
[3,]                     0.08327649
[4,]                     0.01368917
[5,]                     .         


2023-11-21 05:16:57.968987 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927737.arrow, 0.431 mins elapsed.

2023-11-21 05:16:58.29443 : se SRR13927737, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 4439 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(4439): SRR13927737#GTCACGGGTCAACAGG-1
  SRR13927737#TGCCTGTTCTTACTCA-1 ... SRR13927737#TTGCTATAGATGTTCC-1
  SRR13927737#GCTCACTAGTTAGCGG-1
colData names(21): BlacklistRatio nDiFrags ... ReadsInPeaks FRIP

2023-11-21 05:16:58.300516 : Completed GeneIntegrationMatrix : SRR13927737(1 of 4), 0.46 mins elapsed.
2023-11-21 05:16:58.301136 : Reading GeneIntegrationMatrix : SRR13927738(2 of 4), 0.46 mins elapsed.

2023-11-21 05:16:58.312776 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927738 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927738.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3736
[1] "SRR13927738#TAGCATGAGGCTAAAT-1" "SRR13927738#TGTTAGGTCGGTCAGC-1"
[3] "SRR13927738#CTCTCAGCATGTCCCT-1" "SRR13927738#ACTATTCGTCCAACCG-1"
[5] "SRR13927738#TGGTCCTTCTACTGCC-1" "SRR13927738#ATTGTCTAGACTAGGC-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-21 05:16:58.47728 : featureDF SRR13927738, Class = DFrame
DataFrame with 23525 rows and 6 columns
      seqnames     start       end    strand        name       idx
         <Rle> <integer> <integer> <integer> <character> <integer>
1         chr1     69091     70008         1       OR4F5         1
2         chr1    817371    819837         1      FAM87B         2
3         chr1    825138    859446         1   LINC01128         3
4         chr1    827522    826206         2   LINC00115         4
5         chr1    876903    868071         2      FAM41C         5
...        ...       ...       ...       ...         ...       ...
23521     chrX 155264589 155258241         2      RAB39B       923
23522     chrX 155334657 155276211         2       CLIC2       924
23523     chrX 155380787 155381134         1      H2AFB1       925
23524     chrX 155382115 155383230         1        F8A2       926
23525     chrX 155669944 155490115         2       TMLHE       927

2023-11-21 05:16:58.490962 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927738.arrow, 0.003 mins elapsed.

2023-11-21 05:17:15.994451 : mat SRR13927738, Class = dgCMatrix
mat SRR13927738: nRows = 23525, nCols = 3736
mat SRR13927738: NonZeroEntries = 41898897, EntryRange = [ 9.73975135001421e-05 , 758.763031303812 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927738#TAGCATGAGGCTAAAT-1 SRR13927738#TGTTAGGTCGGTCAGC-1
[1,]                    .                                .        
[2,]                    .                                .        
[3,]                    0.002346352                      0.0084575
[4,]                    0.005836413                      .        
[5,]                    .                                .        
     SRR13927738#CTCTCAGCATGTCCCT-1 SRR13927738#ACTATTCGTCCAACCG-1
[1,]                     .                             .          
[2,]                     .                             .          
[3,]                     .                             0.004081077
[4,]                     0.02001583                    .          
[5,]                     .                             .          
     SRR13927738#TGGTCCTTCTACTGCC-1
[1,]                    .          
[2,]                    .          
[3,]                    0.007640877
[4,]                    0.015048796
[5,]                    .          


2023-11-21 05:17:16.775753 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927738.arrow, 0.308 mins elapsed.

2023-11-21 05:17:17.028579 : se SRR13927738, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 3736 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(3736): SRR13927738#TAGCATGAGGCTAAAT-1
  SRR13927738#TGTTAGGTCGGTCAGC-1 ... SRR13927738#GCTCCTATCCGTGCGA-1
  SRR13927738#GGCACGTAGTCGATAA-1
colData names(21): BlacklistRatio nDiFrags ... ReadsInPeaks FRIP

2023-11-21 05:17:17.033851 : Completed GeneIntegrationMatrix : SRR13927738(2 of 4), 0.772 mins elapsed.
2023-11-21 05:17:17.034462 : Reading GeneIntegrationMatrix : SRR13927735(3 of 4), 0.772 mins elapsed.

2023-11-21 05:17:17.047098 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927735 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927735.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3813
[1] "SRR13927735#TTATGTCTCCAGGTAT-1" "SRR13927735#TATTGCTCATCAGAAA-1"
[3] "SRR13927735#TTCGATTGTAGGGTTG-1" "SRR13927735#CATTCATTCGGATGTT-1"
[5] "SRR13927735#ACGTTAGGTCAACTGT-1" "SRR13927735#AAATGCCCAGCAATGG-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-21 05:17:17.21752 : featureDF SRR13927735, Class = DFrame
DataFrame with 23525 rows and 6 columns
      seqnames     start       end    strand        name       idx
         <Rle> <integer> <integer> <integer> <character> <integer>
1         chr1     69091     70008         1       OR4F5         1
2         chr1    817371    819837         1      FAM87B         2
3         chr1    825138    859446         1   LINC01128         3
4         chr1    827522    826206         2   LINC00115         4
5         chr1    876903    868071         2      FAM41C         5
...        ...       ...       ...       ...         ...       ...
23521     chrX 155264589 155258241         2      RAB39B       923
23522     chrX 155334657 155276211         2       CLIC2       924
23523     chrX 155380787 155381134         1      H2AFB1       925
23524     chrX 155382115 155383230         1        F8A2       926
23525     chrX 155669944 155490115         2       TMLHE       927

2023-11-21 05:17:17.231944 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927735.arrow, 0.003 mins elapsed.

2023-11-21 05:17:36.728151 : mat SRR13927735, Class = dgCMatrix
mat SRR13927735: nRows = 23525, nCols = 3813
mat SRR13927735: NonZeroEntries = 45535634, EntryRange = [ 7.20715653308254e-05 , 673.990433360078 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927735#TTATGTCTCCAGGTAT-1 SRR13927735#TATTGCTCATCAGAAA-1
[1,]                    .                              .          
[2,]                    .                              .          
[3,]                    0.070892732                    0.057009268
[4,]                    0.008800039                    0.002978408
[5,]                    .                              .          
     SRR13927735#TTCGATTGTAGGGTTG-1 SRR13927735#CATTCATTCGGATGTT-1
[1,]                     .                               .        
[2,]                     .                               .        
[3,]                     0.05523011                      0.0248301
[4,]                     .                               .        
[5,]                     .                               .        
     SRR13927735#ACGTTAGGTCAACTGT-1
[1,]                     .         
[2,]                     .         
[3,]                     0.05384571
[4,]                     0.01228688
[5,]                     .         


2023-11-21 05:17:37.06962 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927735.arrow, 0.334 mins elapsed.

2023-11-21 05:17:37.313357 : se SRR13927735, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 3813 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(3813): SRR13927735#TTATGTCTCCAGGTAT-1
  SRR13927735#TATTGCTCATCAGAAA-1 ... SRR13927735#GATCATGGTGTTGTTG-1
  SRR13927735#GTCCATCCATTTCACT-1
colData names(21): BlacklistRatio nDiFrags ... ReadsInPeaks FRIP

2023-11-21 05:17:37.318566 : Completed GeneIntegrationMatrix : SRR13927735(3 of 4), 1.11 mins elapsed.
2023-11-21 05:17:37.31918 : Reading GeneIntegrationMatrix : SRR13927736(4 of 4), 1.11 mins elapsed.

2023-11-21 05:17:37.330523 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927736 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927736.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "GeneIntegrationMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3534
[1] "SRR13927736#GCTGTTCTCCTTGACC-1" "SRR13927736#AGTTACGAGTAGGTCG-1"
[3] "SRR13927736#AGTTTGGAGCATTCCA-1" "SRR13927736#GCACGGTTCGGCAATT-1"
[5] "SRR13927736#TCAGCTCTCTCTATTG-1" "SRR13927736#CCTTAATTCTACTTTG-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-21 05:17:37.495765 : featureDF SRR13927736, Class = DFrame
DataFrame with 23525 rows and 6 columns
      seqnames     start       end    strand        name       idx
         <Rle> <integer> <integer> <integer> <character> <integer>
1         chr1     69091     70008         1       OR4F5         1
2         chr1    817371    819837         1      FAM87B         2
3         chr1    825138    859446         1   LINC01128         3
4         chr1    827522    826206         2   LINC00115         4
5         chr1    876903    868071         2      FAM41C         5
...        ...       ...       ...       ...         ...       ...
23521     chrX 155264589 155258241         2      RAB39B       923
23522     chrX 155334657 155276211         2       CLIC2       924
23523     chrX 155380787 155381134         1      H2AFB1       925
23524     chrX 155382115 155383230         1        F8A2       926
23525     chrX 155669944 155490115         2       TMLHE       927

2023-11-21 05:17:37.50942 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927736.arrow, 0.003 mins elapsed.

2023-11-21 05:17:53.675331 : mat SRR13927736, Class = dgCMatrix
mat SRR13927736: nRows = 23525, nCols = 3534
mat SRR13927736: NonZeroEntries = 39003616, EntryRange = [ 8.77627031078293e-05 , 770.387276519725 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927736#GCTGTTCTCCTTGACC-1 SRR13927736#AGTTACGAGTAGGTCG-1
[1,]                     .                             .          
[2,]                     .                             .          
[3,]                     0.07845381                    0.028344539
[4,]                     0.01395289                    0.008499189
[5,]                     .                             .          
     SRR13927736#AGTTTGGAGCATTCCA-1 SRR13927736#GCACGGTTCGGCAATT-1
[1,]                     .                              .         
[2,]                     .                              .         
[3,]                     0.06535130                     0.06576222
[4,]                     0.01810177                     0.01315728
[5,]                     .                              .         
     SRR13927736#TCAGCTCTCTCTATTG-1
[1,]                    .          
[2,]                    .          
[3,]                    0.016229526
[4,]                    0.003862592
[5,]                    .          


2023-11-21 05:17:53.962265 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927736.arrow, 0.277 mins elapsed.

2023-11-21 05:17:54.195783 : se SRR13927736, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 3534 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(3534): SRR13927736#GCTGTTCTCCTTGACC-1
  SRR13927736#AGTTACGAGTAGGTCG-1 ... SRR13927736#TAAACCGCAAGTCTGT-1
  SRR13927736#TGCATTTGTTCTACCC-1
colData names(21): BlacklistRatio nDiFrags ... ReadsInPeaks FRIP

2023-11-21 05:17:54.200783 : Completed GeneIntegrationMatrix : SRR13927736(4 of 4), 1.392 mins elapsed.
2023-11-21 05:17:54.20489 : Organizing colData, 1.392 mins elapsed.
2023-11-21 05:17:54.341611 : Organizing rowData, 1.394 mins elapsed.
2023-11-21 05:17:54.344947 : Organizing rowRanges, 1.394 mins elapsed.
2023-11-21 05:17:54.349657 : Organizing Assays (1 of 1), 1.394 mins elapsed.
2023-11-21 05:17:59.429618 : Constructing SummarizedExperiment, 1.479 mins elapsed.
2023-11-21 05:18:01.836366 : Finished Matrix Creation, 1.519 mins elapsed.

2023-11-21 05:18:01.839461 : getMatrixFromProject Input-Parameters, Class = list

getMatrixFromProject Input-Parameters$ArchRProj: length = 1

getMatrixFromProject Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromProject Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromProject Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromProject Input-Parameters$verbose: length = 1
[1] TRUE


getMatrixFromProject Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromProject Input-Parameters$threads: length = 1
[1] 1


getMatrixFromProject Input-Parameters$logFile: length = 1
[1] "x"


2023-11-21 05:18:02.555487 : Reading PeakMatrix : SRR13927737(1 of 4), 0.012 mins elapsed.

2023-11-21 05:18:02.57233 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927737 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927737.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 4439
[1] "SRR13927737#GTCACGGGTCAACAGG-1" "SRR13927737#TGCCTGTTCTTACTCA-1"
[3] "SRR13927737#TTGCGGGTCAGGTCTA-1" "SRR13927737#GGTTGCGCAAAGCTTC-1"
[5] "SRR13927737#CTTAATCTCGCGCTGA-1" "SRR13927737#TGTACGATCGGAAAGC-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-21 05:18:02.903448 : featureDF SRR13927737, Class = DFrame
DataFrame with 210721 rows and 4 columns
       seqnames       idx     start       end
          <Rle> <integer> <integer> <integer>
1          chr1         1    804677    805177
2          chr1         2    811866    812366
3          chr1         3    817090    817590
4          chr1         4    818831    819331
5          chr1         5    826543    827043
...         ...       ...       ...       ...
210717     chrX      5777 155898525 155899025
210718     chrX      5778 155899804 155900304
210719     chrX      5779 155956608 155957108
210720     chrX      5780 155966827 155967327
210721     chrX      5781 155971625 155972125

2023-11-21 05:18:02.956607 : Getting PeakMatrix from ArrowFile : SRR13927737.arrow, 0.006 mins elapsed.

2023-11-21 05:18:13.804756 : mat SRR13927737, Class = dgCMatrix
mat SRR13927737: nRows = 210721, nCols = 4439
mat SRR13927737: NonZeroEntries = 29963541, EntryRange = [ 1 , 4 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927737#GTCACGGGTCAACAGG-1 SRR13927737#TGCCTGTTCTTACTCA-1
[1,]                              .                              .
[2,]                              .                              2
[3,]                              .                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927737#TTGCGGGTCAGGTCTA-1 SRR13927737#GGTTGCGCAAAGCTTC-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              2                              .
[4,]                              .                              .
[5,]                              4                              .
     SRR13927737#CTTAATCTCGCGCTGA-1
[1,]                              .
[2,]                              .
[3,]                              .
[4,]                              .
[5,]                              .


2023-11-21 05:18:14.02683 : Organizing SE PeakMatrix from ArrowFile : SRR13927737.arrow, 0.191 mins elapsed.

2023-11-21 05:18:14.298368 : se SRR13927737, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 4439 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(4439): SRR13927737#GTCACGGGTCAACAGG-1
  SRR13927737#TGCCTGTTCTTACTCA-1 ... SRR13927737#TTGCTATAGATGTTCC-1
  SRR13927737#GCTCACTAGTTAGCGG-1
colData names(21): BlacklistRatio nDiFrags ... ReadsInPeaks FRIP

2023-11-21 05:18:14.307543 : Completed PeakMatrix : SRR13927737(1 of 4), 0.208 mins elapsed.
2023-11-21 05:18:14.308154 : Reading PeakMatrix : SRR13927738(2 of 4), 0.208 mins elapsed.

2023-11-21 05:18:14.319844 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927738 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927738.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3736
[1] "SRR13927738#TAGCATGAGGCTAAAT-1" "SRR13927738#TGTTAGGTCGGTCAGC-1"
[3] "SRR13927738#CTCTCAGCATGTCCCT-1" "SRR13927738#ACTATTCGTCCAACCG-1"
[5] "SRR13927738#TGGTCCTTCTACTGCC-1" "SRR13927738#ATTGTCTAGACTAGGC-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-21 05:18:14.547801 : featureDF SRR13927738, Class = DFrame
DataFrame with 210721 rows and 4 columns
       seqnames       idx     start       end
          <Rle> <integer> <integer> <integer>
1          chr1         1    804677    805177
2          chr1         2    811866    812366
3          chr1         3    817090    817590
4          chr1         4    818831    819331
5          chr1         5    826543    827043
...         ...       ...       ...       ...
210717     chrX      5777 155898525 155899025
210718     chrX      5778 155899804 155900304
210719     chrX      5779 155956608 155957108
210720     chrX      5780 155966827 155967327
210721     chrX      5781 155971625 155972125

2023-11-21 05:18:14.587556 : Getting PeakMatrix from ArrowFile : SRR13927738.arrow, 0.004 mins elapsed.

2023-11-21 05:18:25.919451 : mat SRR13927738, Class = dgCMatrix
mat SRR13927738: nRows = 210721, nCols = 3736
mat SRR13927738: NonZeroEntries = 30165921, EntryRange = [ 1 , 4 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927738#TAGCATGAGGCTAAAT-1 SRR13927738#TGTTAGGTCGGTCAGC-1
[1,]                              2                              .
[2,]                              .                              .
[3,]                              .                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927738#CTCTCAGCATGTCCCT-1 SRR13927738#ACTATTCGTCCAACCG-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              .                              .
[4,]                              .                              .
[5,]                              2                              .
     SRR13927738#TGGTCCTTCTACTGCC-1
[1,]                              .
[2,]                              .
[3,]                              2
[4,]                              .
[5,]                              .


2023-11-21 05:18:26.151361 : Organizing SE PeakMatrix from ArrowFile : SRR13927738.arrow, 0.197 mins elapsed.

2023-11-21 05:18:26.436878 : se SRR13927738, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 3736 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(3736): SRR13927738#TAGCATGAGGCTAAAT-1
  SRR13927738#TGTTAGGTCGGTCAGC-1 ... SRR13927738#GCTCCTATCCGTGCGA-1
  SRR13927738#GGCACGTAGTCGATAA-1
colData names(21): BlacklistRatio nDiFrags ... ReadsInPeaks FRIP

2023-11-21 05:18:26.444698 : Completed PeakMatrix : SRR13927738(2 of 4), 0.41 mins elapsed.
2023-11-21 05:18:26.4453 : Reading PeakMatrix : SRR13927735(3 of 4), 0.41 mins elapsed.

2023-11-21 05:18:26.457078 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927735 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927735.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3813
[1] "SRR13927735#TTATGTCTCCAGGTAT-1" "SRR13927735#TATTGCTCATCAGAAA-1"
[3] "SRR13927735#TTCGATTGTAGGGTTG-1" "SRR13927735#CATTCATTCGGATGTT-1"
[5] "SRR13927735#ACGTTAGGTCAACTGT-1" "SRR13927735#AAATGCCCAGCAATGG-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-21 05:18:26.830597 : featureDF SRR13927735, Class = DFrame
DataFrame with 210721 rows and 4 columns
       seqnames       idx     start       end
          <Rle> <integer> <integer> <integer>
1          chr1         1    804677    805177
2          chr1         2    811866    812366
3          chr1         3    817090    817590
4          chr1         4    818831    819331
5          chr1         5    826543    827043
...         ...       ...       ...       ...
210717     chrX      5777 155898525 155899025
210718     chrX      5778 155899804 155900304
210719     chrX      5779 155956608 155957108
210720     chrX      5780 155966827 155967327
210721     chrX      5781 155971625 155972125

2023-11-21 05:18:26.870041 : Getting PeakMatrix from ArrowFile : SRR13927735.arrow, 0.007 mins elapsed.

2023-11-21 05:18:38.931096 : mat SRR13927735, Class = dgCMatrix
mat SRR13927735: nRows = 210721, nCols = 3813
mat SRR13927735: NonZeroEntries = 26038084, EntryRange = [ 1 , 4 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927735#TTATGTCTCCAGGTAT-1 SRR13927735#TATTGCTCATCAGAAA-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              1                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927735#TTCGATTGTAGGGTTG-1 SRR13927735#CATTCATTCGGATGTT-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              .                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927735#ACGTTAGGTCAACTGT-1
[1,]                              2
[2,]                              .
[3,]                              .
[4,]                              .
[5,]                              .


2023-11-21 05:18:39.097337 : Organizing SE PeakMatrix from ArrowFile : SRR13927735.arrow, 0.211 mins elapsed.

2023-11-21 05:18:39.3669 : se SRR13927735, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 3813 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(3813): SRR13927735#TTATGTCTCCAGGTAT-1
  SRR13927735#TATTGCTCATCAGAAA-1 ... SRR13927735#GATCATGGTGTTGTTG-1
  SRR13927735#GTCCATCCATTTCACT-1
colData names(21): BlacklistRatio nDiFrags ... ReadsInPeaks FRIP

2023-11-21 05:18:39.374169 : Completed PeakMatrix : SRR13927735(3 of 4), 0.626 mins elapsed.
2023-11-21 05:18:39.374777 : Reading PeakMatrix : SRR13927736(4 of 4), 0.626 mins elapsed.

2023-11-21 05:18:39.387093 : getMatrixFromArrow Input-Parameters, Class = list

getMatrixFromArrow Input-Parameters$ArrowFile: length = 1
                                                                              SRR13927736 
"/gstore/project/lineage/prostate/GSE168667/OUTPUT/multiome/ArrowFiles/SRR13927736.arrow" 


getMatrixFromArrow Input-Parameters$useMatrix: length = 1
[1] "PeakMatrix"


getMatrixFromArrow Input-Parameters$useSeqnames: length = 0
NULL


getMatrixFromArrow Input-Parameters$excludeChr: length = 0
NULL


getMatrixFromArrow Input-Parameters$cellNames: length = 3534
[1] "SRR13927736#GCTGTTCTCCTTGACC-1" "SRR13927736#AGTTACGAGTAGGTCG-1"
[3] "SRR13927736#AGTTTGGAGCATTCCA-1" "SRR13927736#GCACGGTTCGGCAATT-1"
[5] "SRR13927736#TCAGCTCTCTCTATTG-1" "SRR13927736#CCTTAATTCTACTTTG-1"


getMatrixFromArrow Input-Parameters$ArchRProj: length = 1

getMatrixFromArrow Input-Parameters$verbose: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$binarize: length = 1
[1] FALSE


getMatrixFromArrow Input-Parameters$logFile: length = 1
[1] "x"



2023-11-21 05:18:39.62753 : featureDF SRR13927736, Class = DFrame
DataFrame with 210721 rows and 4 columns
       seqnames       idx     start       end
          <Rle> <integer> <integer> <integer>
1          chr1         1    804677    805177
2          chr1         2    811866    812366
3          chr1         3    817090    817590
4          chr1         4    818831    819331
5          chr1         5    826543    827043
...         ...       ...       ...       ...
210717     chrX      5777 155898525 155899025
210718     chrX      5778 155899804 155900304
210719     chrX      5779 155956608 155957108
210720     chrX      5780 155966827 155967327
210721     chrX      5781 155971625 155972125

2023-11-21 05:18:39.667704 : Getting PeakMatrix from ArrowFile : SRR13927736.arrow, 0.005 mins elapsed.

2023-11-21 05:18:49.018342 : mat SRR13927736, Class = dgCMatrix
mat SRR13927736: nRows = 210721, nCols = 3534
mat SRR13927736: NonZeroEntries = 25583439, EntryRange = [ 1 , 4 ]
5 x 5 sparse Matrix of class "dgCMatrix"
     SRR13927736#GCTGTTCTCCTTGACC-1 SRR13927736#AGTTACGAGTAGGTCG-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              .                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927736#AGTTTGGAGCATTCCA-1 SRR13927736#GCACGGTTCGGCAATT-1
[1,]                              .                              .
[2,]                              .                              .
[3,]                              1                              .
[4,]                              .                              .
[5,]                              .                              .
     SRR13927736#TCAGCTCTCTCTATTG-1
[1,]                              .
[2,]                              .
[3,]                              .
[4,]                              .
[5,]                              .


2023-11-21 05:18:49.207265 : Organizing SE PeakMatrix from ArrowFile : SRR13927736.arrow, 0.164 mins elapsed.

2023-11-21 05:18:49.465388 : se SRR13927736, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 3534 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(3534): SRR13927736#GCTGTTCTCCTTGACC-1
  SRR13927736#AGTTACGAGTAGGTCG-1 ... SRR13927736#TAAACCGCAAGTCTGT-1
  SRR13927736#TGCATTTGTTCTACCC-1
colData names(21): BlacklistRatio nDiFrags ... ReadsInPeaks FRIP

2023-11-21 05:18:49.472285 : Completed PeakMatrix : SRR13927736(4 of 4), 0.794 mins elapsed.
2023-11-21 05:18:49.473431 : Organizing colData, 0.794 mins elapsed.
2023-11-21 05:18:49.606192 : Organizing rowData, 0.796 mins elapsed.
2023-11-21 05:18:49.614302 : Organizing rowRanges, 0.796 mins elapsed.
2023-11-21 05:18:49.624927 : Organizing Assays (1 of 1), 0.796 mins elapsed.
2023-11-21 05:18:51.401084 : Constructing SummarizedExperiment, 0.826 mins elapsed.
2023-11-21 05:19:05.165389 : Finished Matrix Creation, 1.055 mins elapsed.
