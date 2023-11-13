
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
