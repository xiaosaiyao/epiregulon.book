
           ___      .______        ______  __    __  .______      
          /   \     |   _  \      /      ||  |  |  | |   _  \     
         /  ^  \    |  |_)  |    |  ,----'|  |__|  | |  |_)  |    
        /  /_\  \   |      /     |  |     |   __   | |      /     
       /  _____  \  |  |\  \\___ |  `----.|  |  |  | |  |\  \\___.
      /__/     \__\ | _| `._____| \______||__|  |__| | _| `._____|
    
Logging With ArchR!

Start Time : 2025-01-13 23:35:50.356463

------- ArchR Info

ArchRThreads = 14
ArchRGenome = Hg19test2

------- System Info

Computer OS = unix
Total Cores = 28

------- Session Info

R Under development (unstable) (2023-12-04 r85659)
Platform: x86_64-pc-linux-gnu
Running under: Ubuntu 18.04.6 LTS

Matrix products: default
BLAS:   /usr/local/lib/R/lib/libRblas.so 
LAPACK: /usr/local/lib/R/lib/libRlapack.so;  LAPACK version 3.11.0

Random number generation:
 RNG:     L'Ecuyer-CMRG 
 Normal:  Inversion 
 Sample:  Rejection 
 
locale:
 [1] LC_CTYPE=en_US.UTF-8 LC_NUMERIC=C         LC_TIME=C            LC_COLLATE=C         LC_MONETARY=C        LC_MESSAGES=C       
 [7] LC_PAPER=C           LC_NAME=C            LC_ADDRESS=C         LC_TELEPHONE=C       LC_MEASUREMENT=C     LC_IDENTIFICATION=C 

time zone: Etc/UTC
tzcode source: system (glibc)

attached base packages:
 [1] parallel  grid      stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] igraph_2.1.2                      msigdbr_7.5.1                     epiregulon.extra_1.1.4           
 [4] BiocStyle_2.33.1                  DropletUtils_1.25.2               zellkonverter_1.15.4             
 [7] ensembldb_2.29.1                  AnnotationFilter_1.29.0           GenomicFeatures_1.57.1           
[10] AnnotationDbi_1.67.0              pbmcMultiome.SeuratData_0.1.4     SeuratData_0.2.2.9001            
[13] Signac_1.14.0                     Seurat_5.1.0                      SeuratObject_5.0.2               
[16] sp_2.1-4                          epiregulon.archr_0.99.5           rhdf5_2.48.0                     
[19] RcppArmadillo_14.2.2-1            Rcpp_1.0.13-1                     sparseMatrixStats_1.16.0         
[22] data.table_1.16.4                 stringr_1.5.1                     plyr_1.8.9                       
[25] magrittr_2.0.3                    ggplot2_3.5.1                     gtable_0.3.6                     
[28] gtools_3.9.5                      gridExtra_2.3                     devtools_2.4.5                   
[31] usethis_3.1.0                     ArchR_1.0.3                       Matrix_1.6-4                     
[34] BSgenome.Hsapiens.UCSC.hg38_1.4.5 BSgenome_1.72.0                   rtracklayer_1.64.0               
[37] BiocIO_1.14.0                     Biostrings_2.72.1                 XVector_0.44.0                   
[40] scMultiome_1.4.2                  MultiAssayExperiment_1.30.3       ExperimentHub_2.12.0             
[43] AnnotationHub_3.12.0              BiocFileCache_2.12.0              dbplyr_2.5.0                     
[46] epiregulon_1.0.1                  SingleCellExperiment_1.26.0       SummarizedExperiment_1.35.5      
[49] Biobase_2.64.0                    GenomicRanges_1.56.2              GenomeInfoDb_1.40.1              
[52] IRanges_2.38.1                    S4Vectors_0.42.1                  BiocGenerics_0.51.3              
[55] MatrixGenerics_1.16.0             matrixStats_1.4.1                

loaded via a namespace (and not attached):
  [1] R.methodsS3_1.8.2                 GSEABase_1.66.0                   urlchecker_1.0.1                 
  [4] poweRlaw_0.80.0                   goftest_1.2-3                     HDF5Array_1.33.8                 
  [7] vctrs_0.6.5                       ggtangle_0.0.3                    spatstat.random_3.3-2            
 [10] digest_0.6.37                     png_0.1-8                         ggrepel_0.9.6                    
 [13] deldir_2.0-4                      parallelly_1.41.0                 MASS_7.3-61                      
 [16] reshape2_1.4.4                    qvalue_2.37.0                     httpuv_1.6.15                    
 [19] withr_3.0.2                       ggrastr_1.0.2                     ggfun_0.1.7                      
 [22] xfun_0.49                         ellipsis_0.3.2                    survival_3.7-0                   
 [25] memoise_2.0.1                     ggbeeswarm_0.7.2                  gson_0.1.0                       
 [28] clusterProfiler_4.13.4            profvis_0.4.0                     tidytree_0.4.6                   
 [31] KEGGgraph_1.65.0                  zoo_1.8-12                        pbapply_1.7-2                    
 [34] R.oo_1.27.0                       KEGGREST_1.44.1                   promises_1.3.2                   
 [37] httr_1.4.7                        restfulr_0.0.15                   globals_0.16.3                   
 [40] fitdistrplus_1.2-1                rhdf5filters_1.16.0               rstudioapi_0.17.1                
 [43] DOSE_3.99.1                       UCSC.utils_1.0.0                  miniUI_0.1.1.1                   
 [46] generics_0.1.3                    dir.expiry_1.13.1                 babelgene_22.9                   
 [49] curl_6.0.1                        zlibbioc_1.50.0                   ggraph_2.2.1                     
 [52] ScaledMatrix_1.12.0               polyclip_1.10-7                   GenomeInfoDbData_1.2.12          
 [55] SparseArray_1.5.45                xtable_1.8-4                      pracma_2.4.4                     
 [58] evaluate_1.0.1                    S4Arrays_1.5.11                   hms_1.1.3                        
 [61] bookdown_0.41                     irlba_2.3.5.1                     colorspace_2.1-1                 
 [64] filelock_1.0.3                    ROCR_1.0-11                       reticulate_1.39.0                
 [67] spatstat.data_3.1-2               Rgraphviz_2.49.1                  lmtest_0.9-40                    
 [70] readr_2.1.5                       ggtree_3.13.2                     viridis_0.6.5                    
 [73] later_1.4.1                       lattice_0.22-6                    spatstat.geom_3.3-3              
 [76] future.apply_1.11.3               scattermore_1.2                   XML_3.99-0.17                    
 [79] scuttle_1.14.0                    cowplot_1.1.3                     RcppAnnoy_0.0.22                 
 [82] pillar_1.9.0                      nlme_3.1-166                      pwalign_1.0.0                    
 [85] caTools_1.18.3                    compiler_4.4.0                    beachmat_2.20.0                  
 [88] RSpectra_0.16-2                   stringi_1.8.4                     tensor_1.5                       
 [91] GenomicAlignments_1.40.0          crayon_1.5.3                      abind_1.4-8                      
 [94] scater_1.33.4                     gridGraphics_0.5-1                locfit_1.5-9.10                  
 [97] graphlayouts_1.2.0                bit_4.5.0.1                       dplyr_1.1.4                      
[100] fastmatch_1.1-4                   codetools_0.2-20                  BiocSingular_1.20.0              
[103] bslib_0.8.0                       plotly_4.10.4                     mime_0.12                        
[106] splines_4.4.0                     fastDummies_1.7.4                 basilisk_1.17.2                  
[109] EnrichmentBrowser_2.35.1          knitr_1.49                        blob_1.2.4                       
[112] utf8_1.2.4                        BiocVersion_3.19.1                seqLogo_1.70.0                   
[115] fs_1.6.5                          listenv_0.9.1                     checkmate_2.3.2                  
[118] DelayedMatrixStats_1.26.0         pkgbuild_1.4.5                    ggplotify_0.1.2                  
[121] tibble_3.2.1                      statmod_1.5.0                     tzdb_0.4.0                       
[124] tweenr_2.0.3                      pkgconfig_2.0.3                   BSgenome.Hsapiens.UCSC.hg19_1.4.3
[127] tools_4.4.0                       cachem_1.1.0                      RSQLite_2.3.9                    
[130] viridisLite_0.4.2                 DBI_1.2.3                         fastmap_1.2.0                    
[133] rmarkdown_2.29                    scales_1.3.0                      ica_1.0-3                        
[136] Rsamtools_2.20.0                  sass_0.4.9                        patchwork_1.3.0                  
[139] BiocManager_1.30.25               dotCall64_1.2                     graph_1.82.0                     
[142] RANN_2.6.2                        farver_2.1.2                      tidygraph_1.3.1                  
[145] yaml_2.3.10                       cli_3.6.3                         purrr_1.0.2                      
[148] motifmatchr_1.26.0                leiden_0.4.3.1                    lifecycle_1.0.4                  
[151] uwot_0.2.2                        bluster_1.14.0                    sessioninfo_1.2.2                
[154] backports_1.5.0                   BiocParallel_1.38.0               annotate_1.82.0                  
[157] rjson_0.2.23                      ggridges_0.5.6                    progressr_0.15.1                 
[160] ape_5.8                           limma_3.60.6                      jsonlite_1.8.9                   
[163] edgeR_4.2.2                       RcppHNSW_0.6.0                    TFBSTools_1.31.2                 
[166] bitops_1.0-9                      bit64_4.5.2                       Rtsne_0.17                       
[169] yulab.utils_0.1.7                 BiocNeighbors_1.22.0              spatstat.utils_3.1-0             
[172] CNEr_1.40.0                       metapod_1.12.0                    jquerylib_0.1.4                  
[175] highr_0.11                        GOSemSim_2.31.2                   dqrng_0.4.1                      
[178] spatstat.univar_3.0-1             R.utils_2.12.3                    lazyeval_0.2.2                   
[181] shiny_1.9.1                       enrichplot_1.25.5                 htmltools_0.5.8.1                
[184] GO.db_3.19.1                      sctransform_0.4.1                 rappdirs_0.3.3                   
[187] basilisk.utils_1.17.3             glue_1.8.0                        TFMPvalue_0.0.9                  
[190] spam_2.11-0                       RCurl_1.98-1.16                   treeio_1.29.2                    
[193] scran_1.32.0                      R6_2.5.1                          tidyr_1.3.1                      
[196] labeling_0.4.3                    RcppRoll_0.3.1                    cluster_2.1.6                    
[199] pkgload_1.4.0                     Rhdf5lib_1.26.0                   aplot_0.2.3                      
[202] vipor_0.4.7                       DelayedArray_0.31.14              tidyselect_1.2.1                 
[205] ProtGenerics_1.37.1               ggforce_0.4.2                     future_1.34.0                    
[208] rsvd_1.0.5                        munsell_0.5.1                     KernSmooth_2.23-22               
[211] fgsea_1.31.6                      htmlwidgets_1.6.4                 RColorBrewer_1.1-3               
[214] rlang_1.1.4                       spatstat.sparse_3.1-0             spatstat.explore_3.3-3           
[217] remotes_2.5.0                     fansi_1.0.6                       Cairo_1.6-2                      
[220] beeswarm_0.4.0                   


------- Log Info


2025-01-13 23:35:53.648163 : getMatrixFromProject Input-Parameters, Class = list

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


2025-01-13 23:35:57.542904 : Reading GeneIntegrationMatrix : SRR13927737(1 of 4), 0.12 mins elapsed.

2025-01-13 23:35:58.937917 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-13 23:36:01.692159 : featureDF SRR13927737, Class = DFrame
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

2025-01-13 23:36:02.673401 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927737.arrow, 0.062 mins elapsed.

2025-01-13 23:36:28.507799 : mat SRR13927737, Class = dgCMatrix
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


2025-01-13 23:36:30.962281 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927737.arrow, 0.534 mins elapsed.

2025-01-13 23:36:32.734032 : se SRR13927737, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 4439 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(4439): SRR13927737#GTCACGGGTCAACAGG-1 SRR13927737#TGCCTGTTCTTACTCA-1 ... SRR13927737#TTGCTATAGATGTTCC-1
  SRR13927737#GCTCACTAGTTAGCGG-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2025-01-13 23:36:32.744521 : Completed GeneIntegrationMatrix : SRR13927737(1 of 4), 0.707 mins elapsed.
2025-01-13 23:36:32.745417 : Reading GeneIntegrationMatrix : SRR13927738(2 of 4), 0.707 mins elapsed.

2025-01-13 23:36:32.785628 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-13 23:36:35.004172 : featureDF SRR13927738, Class = DFrame
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

2025-01-13 23:36:35.047241 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927738.arrow, 0.038 mins elapsed.

2025-01-13 23:36:58.841415 : mat SRR13927738, Class = dgCMatrix
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


2025-01-13 23:37:01.028763 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927738.arrow, 0.471 mins elapsed.

2025-01-13 23:37:02.721933 : se SRR13927738, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 3736 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(3736): SRR13927738#TAGCATGAGGCTAAAT-1 SRR13927738#TGTTAGGTCGGTCAGC-1 ... SRR13927738#GCTCCTATCCGTGCGA-1
  SRR13927738#GGCACGTAGTCGATAA-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2025-01-13 23:37:02.730872 : Completed GeneIntegrationMatrix : SRR13927738(2 of 4), 1.207 mins elapsed.
2025-01-13 23:37:02.73173 : Reading GeneIntegrationMatrix : SRR13927735(3 of 4), 1.207 mins elapsed.

2025-01-13 23:37:02.770508 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-13 23:37:05.091062 : featureDF SRR13927735, Class = DFrame
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

2025-01-13 23:37:05.134717 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927735.arrow, 0.039 mins elapsed.

2025-01-13 23:37:29.362289 : mat SRR13927735, Class = dgCMatrix
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


2025-01-13 23:37:31.371714 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927735.arrow, 0.477 mins elapsed.

2025-01-13 23:37:33.035023 : se SRR13927735, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 3813 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(3813): SRR13927735#TTATGTCTCCAGGTAT-1 SRR13927735#TATTGCTCATCAGAAA-1 ... SRR13927735#GATCATGGTGTTGTTG-1
  SRR13927735#GTCCATCCATTTCACT-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2025-01-13 23:37:33.04832 : Completed GeneIntegrationMatrix : SRR13927735(3 of 4), 1.712 mins elapsed.
2025-01-13 23:37:33.049582 : Reading GeneIntegrationMatrix : SRR13927736(4 of 4), 1.712 mins elapsed.

2025-01-13 23:37:33.09933 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-13 23:37:35.176793 : featureDF SRR13927736, Class = DFrame
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

2025-01-13 23:37:35.203375 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927736.arrow, 0.035 mins elapsed.

2025-01-13 23:37:51.802981 : mat SRR13927736, Class = dgCMatrix
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


2025-01-13 23:37:54.125041 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927736.arrow, 0.35 mins elapsed.

2025-01-13 23:37:55.637863 : se SRR13927736, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 3534 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(3534): SRR13927736#GCTGTTCTCCTTGACC-1 SRR13927736#AGTTACGAGTAGGTCG-1 ... SRR13927736#TAAACCGCAAGTCTGT-1
  SRR13927736#TGCATTTGTTCTACCC-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2025-01-13 23:37:55.644674 : Completed GeneIntegrationMatrix : SRR13927736(4 of 4), 2.089 mins elapsed.
2025-01-13 23:37:55.646678 : Organizing colData, 2.089 mins elapsed.
2025-01-13 23:37:57.098636 : Organizing rowData, 2.113 mins elapsed.
2025-01-13 23:37:57.107224 : Organizing rowRanges, 2.113 mins elapsed.
2025-01-13 23:37:57.121918 : Organizing Assays (1 of 1), 2.113 mins elapsed.
2025-01-13 23:38:07.621079 : Constructing SummarizedExperiment, 2.288 mins elapsed.
2025-01-13 23:38:11.535641 : Finished Matrix Creation, 2.353 mins elapsed.

2025-01-13 23:38:11.540077 : getMatrixFromProject Input-Parameters, Class = list

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


2025-01-13 23:38:14.322589 : Reading PeakMatrix : SRR13927737(1 of 4), 0.046 mins elapsed.

2025-01-13 23:38:14.442307 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-13 23:38:16.581611 : featureDF SRR13927737, Class = DFrame
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

2025-01-13 23:38:16.663918 : Getting PeakMatrix from ArrowFile : SRR13927737.arrow, 0.037 mins elapsed.

2025-01-13 23:38:30.855815 : mat SRR13927737, Class = dgCMatrix
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


2025-01-13 23:38:32.854327 : Organizing SE PeakMatrix from ArrowFile : SRR13927737.arrow, 0.307 mins elapsed.

2025-01-13 23:38:34.313659 : se SRR13927737, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 4439 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(4439): SRR13927737#GTCACGGGTCAACAGG-1 SRR13927737#TGCCTGTTCTTACTCA-1 ... SRR13927737#TTGCTATAGATGTTCC-1
  SRR13927737#GCTCACTAGTTAGCGG-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2025-01-13 23:38:34.325836 : Completed PeakMatrix : SRR13927737(1 of 4), 0.38 mins elapsed.
2025-01-13 23:38:34.326523 : Reading PeakMatrix : SRR13927738(2 of 4), 0.38 mins elapsed.

2025-01-13 23:38:35.080224 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-13 23:38:37.428038 : featureDF SRR13927738, Class = DFrame
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

2025-01-13 23:38:37.517093 : Getting PeakMatrix from ArrowFile : SRR13927738.arrow, 0.041 mins elapsed.

2025-01-13 23:38:51.289794 : mat SRR13927738, Class = dgCMatrix
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


2025-01-13 23:38:52.549667 : Organizing SE PeakMatrix from ArrowFile : SRR13927738.arrow, 0.291 mins elapsed.

2025-01-13 23:38:54.61032 : se SRR13927738, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 3736 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(3736): SRR13927738#TAGCATGAGGCTAAAT-1 SRR13927738#TGTTAGGTCGGTCAGC-1 ... SRR13927738#GCTCCTATCCGTGCGA-1
  SRR13927738#GGCACGTAGTCGATAA-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2025-01-13 23:38:54.628472 : Completed PeakMatrix : SRR13927738(2 of 4), 0.718 mins elapsed.
2025-01-13 23:38:54.62977 : Reading PeakMatrix : SRR13927735(3 of 4), 0.718 mins elapsed.

2025-01-13 23:38:54.684032 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-13 23:38:57.018433 : featureDF SRR13927735, Class = DFrame
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

2025-01-13 23:38:57.105129 : Getting PeakMatrix from ArrowFile : SRR13927735.arrow, 0.04 mins elapsed.

2025-01-13 23:39:09.593212 : mat SRR13927735, Class = dgCMatrix
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


2025-01-13 23:39:11.717335 : Organizing SE PeakMatrix from ArrowFile : SRR13927735.arrow, 0.284 mins elapsed.

2025-01-13 23:39:13.610332 : se SRR13927735, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 3813 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(3813): SRR13927735#TTATGTCTCCAGGTAT-1 SRR13927735#TATTGCTCATCAGAAA-1 ... SRR13927735#GATCATGGTGTTGTTG-1
  SRR13927735#GTCCATCCATTTCACT-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2025-01-13 23:39:13.622503 : Completed PeakMatrix : SRR13927735(3 of 4), 1.035 mins elapsed.
2025-01-13 23:39:13.623372 : Reading PeakMatrix : SRR13927736(4 of 4), 1.035 mins elapsed.

2025-01-13 23:39:13.664442 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-13 23:39:15.957479 : featureDF SRR13927736, Class = DFrame
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

2025-01-13 23:39:16.04032 : Getting PeakMatrix from ArrowFile : SRR13927736.arrow, 0.04 mins elapsed.

2025-01-13 23:39:28.652755 : mat SRR13927736, Class = dgCMatrix
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


2025-01-13 23:39:28.821428 : Organizing SE PeakMatrix from ArrowFile : SRR13927736.arrow, 0.253 mins elapsed.

2025-01-13 23:39:30.317678 : se SRR13927736, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 3534 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(3534): SRR13927736#GCTGTTCTCCTTGACC-1 SRR13927736#AGTTACGAGTAGGTCG-1 ... SRR13927736#TAAACCGCAAGTCTGT-1
  SRR13927736#TGCATTTGTTCTACCC-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2025-01-13 23:39:30.326887 : Completed PeakMatrix : SRR13927736(4 of 4), 1.313 mins elapsed.
2025-01-13 23:39:30.328993 : Organizing colData, 1.313 mins elapsed.
2025-01-13 23:39:31.675667 : Organizing rowData, 1.336 mins elapsed.
2025-01-13 23:39:31.696454 : Organizing rowRanges, 1.336 mins elapsed.
2025-01-13 23:39:33.827352 : Organizing Assays (1 of 1), 1.371 mins elapsed.
2025-01-13 23:39:38.508786 : Constructing SummarizedExperiment, 1.449 mins elapsed.
2025-01-13 23:40:01.310688 : Finished Matrix Creation, 1.829 mins elapsed.

2025-01-14 20:07:13.035798 : getMatrixFromProject Input-Parameters, Class = list

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


2025-01-14 20:07:20.874547 : Reading GeneIntegrationMatrix : SRR13927737(1 of 4), 0.131 mins elapsed.

2025-01-14 20:07:21.295392 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-14 20:07:23.206053 : featureDF SRR13927737, Class = DFrame
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

2025-01-14 20:07:23.264658 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927737.arrow, 0.033 mins elapsed.

2025-01-14 20:07:58.667581 : mat SRR13927737, Class = dgCMatrix
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


2025-01-14 20:07:59.717292 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927737.arrow, 0.64 mins elapsed.

2025-01-14 20:08:00.687486 : se SRR13927737, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 4439 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(4439): SRR13927737#GTCACGGGTCAACAGG-1 SRR13927737#TGCCTGTTCTTACTCA-1 ... SRR13927737#TTGCTATAGATGTTCC-1
  SRR13927737#GCTCACTAGTTAGCGG-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2025-01-14 20:08:00.697446 : Completed GeneIntegrationMatrix : SRR13927737(1 of 4), 0.794 mins elapsed.
2025-01-14 20:08:00.698179 : Reading GeneIntegrationMatrix : SRR13927738(2 of 4), 0.794 mins elapsed.

2025-01-14 20:08:00.939736 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-14 20:08:02.670048 : featureDF SRR13927738, Class = DFrame
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

2025-01-14 20:08:02.694539 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927738.arrow, 0.029 mins elapsed.

2025-01-14 20:08:21.382918 : mat SRR13927738, Class = dgCMatrix
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


2025-01-14 20:08:23.238219 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927738.arrow, 0.372 mins elapsed.

2025-01-14 20:08:23.998352 : se SRR13927738, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 3736 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(3736): SRR13927738#TAGCATGAGGCTAAAT-1 SRR13927738#TGTTAGGTCGGTCAGC-1 ... SRR13927738#GCTCCTATCCGTGCGA-1
  SRR13927738#GGCACGTAGTCGATAA-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2025-01-14 20:08:24.005853 : Completed GeneIntegrationMatrix : SRR13927738(2 of 4), 1.183 mins elapsed.
2025-01-14 20:08:24.006537 : Reading GeneIntegrationMatrix : SRR13927735(3 of 4), 1.183 mins elapsed.

2025-01-14 20:08:24.042767 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-14 20:08:25.501925 : featureDF SRR13927735, Class = DFrame
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

2025-01-14 20:08:25.538582 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927735.arrow, 0.025 mins elapsed.

2025-01-14 20:08:43.820917 : mat SRR13927735, Class = dgCMatrix
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


2025-01-14 20:08:44.881876 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927735.arrow, 0.347 mins elapsed.

2025-01-14 20:08:45.848256 : se SRR13927735, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 3813 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(3813): SRR13927735#TTATGTCTCCAGGTAT-1 SRR13927735#TATTGCTCATCAGAAA-1 ... SRR13927735#GATCATGGTGTTGTTG-1
  SRR13927735#GTCCATCCATTTCACT-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2025-01-14 20:08:45.858253 : Completed GeneIntegrationMatrix : SRR13927735(3 of 4), 1.547 mins elapsed.
2025-01-14 20:08:45.859147 : Reading GeneIntegrationMatrix : SRR13927736(4 of 4), 1.547 mins elapsed.

2025-01-14 20:08:45.895812 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-14 20:08:48.138025 : featureDF SRR13927736, Class = DFrame
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

2025-01-14 20:08:48.164647 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927736.arrow, 0.038 mins elapsed.

2025-01-14 20:09:04.43288 : mat SRR13927736, Class = dgCMatrix
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


2025-01-14 20:09:05.441078 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927736.arrow, 0.326 mins elapsed.

2025-01-14 20:09:06.315999 : se SRR13927736, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 3534 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(3534): SRR13927736#GCTGTTCTCCTTGACC-1 SRR13927736#AGTTACGAGTAGGTCG-1 ... SRR13927736#TAAACCGCAAGTCTGT-1
  SRR13927736#TGCATTTGTTCTACCC-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2025-01-14 20:09:06.328203 : Completed GeneIntegrationMatrix : SRR13927736(4 of 4), 1.888 mins elapsed.
2025-01-14 20:09:06.331526 : Organizing colData, 1.888 mins elapsed.
2025-01-14 20:09:07.072175 : Organizing rowData, 1.901 mins elapsed.
2025-01-14 20:09:07.080179 : Organizing rowRanges, 1.901 mins elapsed.
2025-01-14 20:09:07.094821 : Organizing Assays (1 of 1), 1.901 mins elapsed.
2025-01-14 20:09:14.474174 : Constructing SummarizedExperiment, 2.024 mins elapsed.
2025-01-14 20:09:15.972996 : Finished Matrix Creation, 2.046 mins elapsed.

2025-01-14 20:09:15.9953 : getMatrixFromProject Input-Parameters, Class = list

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


2025-01-14 20:09:17.653626 : Reading PeakMatrix : SRR13927737(1 of 4), 0.028 mins elapsed.

2025-01-14 20:09:17.723287 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-14 20:09:18.92223 : featureDF SRR13927737, Class = DFrame
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

2025-01-14 20:09:19.556652 : Getting PeakMatrix from ArrowFile : SRR13927737.arrow, 0.031 mins elapsed.

2025-01-14 20:09:31.073322 : mat SRR13927737, Class = dgCMatrix
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


2025-01-14 20:09:32.234351 : Organizing SE PeakMatrix from ArrowFile : SRR13927737.arrow, 0.242 mins elapsed.

2025-01-14 20:09:33.33241 : se SRR13927737, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 4439 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(4439): SRR13927737#GTCACGGGTCAACAGG-1 SRR13927737#TGCCTGTTCTTACTCA-1 ... SRR13927737#TTGCTATAGATGTTCC-1
  SRR13927737#GCTCACTAGTTAGCGG-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2025-01-14 20:09:33.356122 : Completed PeakMatrix : SRR13927737(1 of 4), 0.289 mins elapsed.
2025-01-14 20:09:33.357242 : Reading PeakMatrix : SRR13927738(2 of 4), 0.289 mins elapsed.

2025-01-14 20:09:33.394624 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-14 20:09:34.947654 : featureDF SRR13927738, Class = DFrame
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

2025-01-14 20:09:35.025526 : Getting PeakMatrix from ArrowFile : SRR13927738.arrow, 0.027 mins elapsed.

2025-01-14 20:09:46.764805 : mat SRR13927738, Class = dgCMatrix
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


2025-01-14 20:09:47.725021 : Organizing SE PeakMatrix from ArrowFile : SRR13927738.arrow, 0.239 mins elapsed.

2025-01-14 20:09:48.902606 : se SRR13927738, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 3736 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(3736): SRR13927738#TAGCATGAGGCTAAAT-1 SRR13927738#TGTTAGGTCGGTCAGC-1 ... SRR13927738#GCTCCTATCCGTGCGA-1
  SRR13927738#GGCACGTAGTCGATAA-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2025-01-14 20:09:48.918988 : Completed PeakMatrix : SRR13927738(2 of 4), 0.549 mins elapsed.
2025-01-14 20:09:48.920131 : Reading PeakMatrix : SRR13927735(3 of 4), 0.549 mins elapsed.

2025-01-14 20:09:48.963756 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-14 20:09:50.777003 : featureDF SRR13927735, Class = DFrame
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

2025-01-14 20:09:50.85492 : Getting PeakMatrix from ArrowFile : SRR13927735.arrow, 0.032 mins elapsed.

2025-01-14 20:10:01.066351 : mat SRR13927735, Class = dgCMatrix
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


2025-01-14 20:10:02.075826 : Organizing SE PeakMatrix from ArrowFile : SRR13927735.arrow, 0.219 mins elapsed.

2025-01-14 20:10:03.670799 : se SRR13927735, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 3813 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(3813): SRR13927735#TTATGTCTCCAGGTAT-1 SRR13927735#TATTGCTCATCAGAAA-1 ... SRR13927735#GATCATGGTGTTGTTG-1
  SRR13927735#GTCCATCCATTTCACT-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2025-01-14 20:10:03.682132 : Completed PeakMatrix : SRR13927735(3 of 4), 0.795 mins elapsed.
2025-01-14 20:10:03.682894 : Reading PeakMatrix : SRR13927736(4 of 4), 0.795 mins elapsed.

2025-01-14 20:10:03.71717 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-14 20:10:05.953436 : featureDF SRR13927736, Class = DFrame
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

2025-01-14 20:10:06.068831 : Getting PeakMatrix from ArrowFile : SRR13927736.arrow, 0.039 mins elapsed.

2025-01-14 20:10:17.574597 : mat SRR13927736, Class = dgCMatrix
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


2025-01-14 20:10:18.648875 : Organizing SE PeakMatrix from ArrowFile : SRR13927736.arrow, 0.249 mins elapsed.

2025-01-14 20:10:19.879697 : se SRR13927736, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 3534 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(3534): SRR13927736#GCTGTTCTCCTTGACC-1 SRR13927736#AGTTACGAGTAGGTCG-1 ... SRR13927736#TAAACCGCAAGTCTGT-1
  SRR13927736#TGCATTTGTTCTACCC-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2025-01-14 20:10:19.896054 : Completed PeakMatrix : SRR13927736(4 of 4), 1.065 mins elapsed.
2025-01-14 20:10:19.899168 : Organizing colData, 1.065 mins elapsed.
2025-01-14 20:10:20.787257 : Organizing rowData, 1.08 mins elapsed.
2025-01-14 20:10:20.806713 : Organizing rowRanges, 1.08 mins elapsed.
2025-01-14 20:10:20.832217 : Organizing Assays (1 of 1), 1.081 mins elapsed.
2025-01-14 20:10:25.616054 : Constructing SummarizedExperiment, 1.16 mins elapsed.
2025-01-14 20:10:43.911105 : Finished Matrix Creation, 1.465 mins elapsed.

2025-01-14 20:34:11.301916 : getMatrixFromProject Input-Parameters, Class = list

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


2025-01-14 20:34:13.336435 : Reading GeneIntegrationMatrix : SRR13927737(1 of 4), 0.034 mins elapsed.

2025-01-14 20:34:13.758144 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-14 20:34:15.373418 : featureDF SRR13927737, Class = DFrame
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

2025-01-14 20:34:15.433107 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927737.arrow, 0.028 mins elapsed.

2025-01-14 20:34:50.629205 : mat SRR13927737, Class = dgCMatrix
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


2025-01-14 20:34:52.55746 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927737.arrow, 0.647 mins elapsed.

2025-01-14 20:34:53.43373 : se SRR13927737, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 4439 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(4439): SRR13927737#GTCACGGGTCAACAGG-1 SRR13927737#TGCCTGTTCTTACTCA-1 ... SRR13927737#TTGCTATAGATGTTCC-1
  SRR13927737#GCTCACTAGTTAGCGG-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2025-01-14 20:34:53.443342 : Completed GeneIntegrationMatrix : SRR13927737(1 of 4), 0.702 mins elapsed.
2025-01-14 20:34:53.444049 : Reading GeneIntegrationMatrix : SRR13927738(2 of 4), 0.702 mins elapsed.

2025-01-14 20:34:53.474466 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-14 20:34:54.655238 : featureDF SRR13927738, Class = DFrame
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

2025-01-14 20:34:54.683558 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927738.arrow, 0.02 mins elapsed.

2025-01-14 20:35:11.211797 : mat SRR13927738, Class = dgCMatrix
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


2025-01-14 20:35:12.040621 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927738.arrow, 0.309 mins elapsed.

2025-01-14 20:35:12.775224 : se SRR13927738, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 3736 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(3736): SRR13927738#TAGCATGAGGCTAAAT-1 SRR13927738#TGTTAGGTCGGTCAGC-1 ... SRR13927738#GCTCCTATCCGTGCGA-1
  SRR13927738#GGCACGTAGTCGATAA-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2025-01-14 20:35:12.782349 : Completed GeneIntegrationMatrix : SRR13927738(2 of 4), 1.025 mins elapsed.
2025-01-14 20:35:12.782976 : Reading GeneIntegrationMatrix : SRR13927735(3 of 4), 1.025 mins elapsed.

2025-01-14 20:35:12.808828 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-14 20:35:14.018905 : featureDF SRR13927735, Class = DFrame
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

2025-01-14 20:35:14.055355 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927735.arrow, 0.021 mins elapsed.

2025-01-14 20:35:31.869528 : mat SRR13927735, Class = dgCMatrix
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


2025-01-14 20:35:32.834243 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927735.arrow, 0.334 mins elapsed.

2025-01-14 20:35:33.772421 : se SRR13927735, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 3813 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(3813): SRR13927735#TTATGTCTCCAGGTAT-1 SRR13927735#TATTGCTCATCAGAAA-1 ... SRR13927735#GATCATGGTGTTGTTG-1
  SRR13927735#GTCCATCCATTTCACT-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2025-01-14 20:35:33.782015 : Completed GeneIntegrationMatrix : SRR13927735(3 of 4), 1.375 mins elapsed.
2025-01-14 20:35:33.782887 : Reading GeneIntegrationMatrix : SRR13927736(4 of 4), 1.375 mins elapsed.

2025-01-14 20:35:33.816864 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-14 20:35:35.327145 : featureDF SRR13927736, Class = DFrame
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

2025-01-14 20:35:35.377897 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927736.arrow, 0.026 mins elapsed.

2025-01-14 20:35:50.554471 : mat SRR13927736, Class = dgCMatrix
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


2025-01-14 20:35:52.553724 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927736.arrow, 0.312 mins elapsed.

2025-01-14 20:35:53.338934 : se SRR13927736, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 3534 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(3534): SRR13927736#GCTGTTCTCCTTGACC-1 SRR13927736#AGTTACGAGTAGGTCG-1 ... SRR13927736#TAAACCGCAAGTCTGT-1
  SRR13927736#TGCATTTGTTCTACCC-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2025-01-14 20:35:53.347047 : Completed GeneIntegrationMatrix : SRR13927736(4 of 4), 1.701 mins elapsed.
2025-01-14 20:35:53.34921 : Organizing colData, 1.701 mins elapsed.
2025-01-14 20:35:53.832946 : Organizing rowData, 1.709 mins elapsed.
2025-01-14 20:35:53.840908 : Organizing rowRanges, 1.709 mins elapsed.
2025-01-14 20:35:53.857009 : Organizing Assays (1 of 1), 1.709 mins elapsed.
2025-01-14 20:36:00.285365 : Constructing SummarizedExperiment, 1.816 mins elapsed.
2025-01-14 20:36:01.717017 : Finished Matrix Creation, 1.838 mins elapsed.

2025-01-14 20:36:01.886709 : getMatrixFromProject Input-Parameters, Class = list

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


2025-01-14 20:36:03.321941 : Reading PeakMatrix : SRR13927737(1 of 4), 0.024 mins elapsed.

2025-01-14 20:36:03.391763 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-14 20:36:04.236233 : featureDF SRR13927737, Class = DFrame
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

2025-01-14 20:36:04.47121 : Getting PeakMatrix from ArrowFile : SRR13927737.arrow, 0.018 mins elapsed.

2025-01-14 20:36:17.174231 : mat SRR13927737, Class = dgCMatrix
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


2025-01-14 20:36:17.799058 : Organizing SE PeakMatrix from ArrowFile : SRR13927737.arrow, 0.24 mins elapsed.

2025-01-14 20:36:18.955333 : se SRR13927737, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 4439 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(4439): SRR13927737#GTCACGGGTCAACAGG-1 SRR13927737#TGCCTGTTCTTACTCA-1 ... SRR13927737#TTGCTATAGATGTTCC-1
  SRR13927737#GCTCACTAGTTAGCGG-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2025-01-14 20:36:18.972303 : Completed PeakMatrix : SRR13927737(1 of 4), 0.285 mins elapsed.
2025-01-14 20:36:18.973178 : Reading PeakMatrix : SRR13927738(2 of 4), 0.285 mins elapsed.

2025-01-14 20:36:19.00835 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-14 20:36:20.276983 : featureDF SRR13927738, Class = DFrame
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

2025-01-14 20:36:20.350539 : Getting PeakMatrix from ArrowFile : SRR13927738.arrow, 0.022 mins elapsed.

2025-01-14 20:36:32.100372 : mat SRR13927738, Class = dgCMatrix
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


2025-01-14 20:36:33.079599 : Organizing SE PeakMatrix from ArrowFile : SRR13927738.arrow, 0.235 mins elapsed.

2025-01-14 20:36:34.126273 : se SRR13927738, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 3736 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(3736): SRR13927738#TAGCATGAGGCTAAAT-1 SRR13927738#TGTTAGGTCGGTCAGC-1 ... SRR13927738#GCTCCTATCCGTGCGA-1
  SRR13927738#GGCACGTAGTCGATAA-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2025-01-14 20:36:34.135851 : Completed PeakMatrix : SRR13927738(2 of 4), 0.537 mins elapsed.
2025-01-14 20:36:34.136568 : Reading PeakMatrix : SRR13927735(3 of 4), 0.538 mins elapsed.

2025-01-14 20:36:34.156478 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-14 20:36:35.83752 : featureDF SRR13927735, Class = DFrame
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

2025-01-14 20:36:35.906529 : Getting PeakMatrix from ArrowFile : SRR13927735.arrow, 0.029 mins elapsed.

2025-01-14 20:36:46.547737 : mat SRR13927735, Class = dgCMatrix
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


2025-01-14 20:36:47.593592 : Organizing SE PeakMatrix from ArrowFile : SRR13927735.arrow, 0.224 mins elapsed.

2025-01-14 20:36:48.896896 : se SRR13927735, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 3813 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(3813): SRR13927735#TTATGTCTCCAGGTAT-1 SRR13927735#TATTGCTCATCAGAAA-1 ... SRR13927735#GATCATGGTGTTGTTG-1
  SRR13927735#GTCCATCCATTTCACT-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2025-01-14 20:36:48.906291 : Completed PeakMatrix : SRR13927735(3 of 4), 0.784 mins elapsed.
2025-01-14 20:36:48.906989 : Reading PeakMatrix : SRR13927736(4 of 4), 0.784 mins elapsed.

2025-01-14 20:36:49.263128 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-14 20:36:50.961748 : featureDF SRR13927736, Class = DFrame
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

2025-01-14 20:36:51.038388 : Getting PeakMatrix from ArrowFile : SRR13927736.arrow, 0.03 mins elapsed.

2025-01-14 20:37:01.031487 : mat SRR13927736, Class = dgCMatrix
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


2025-01-14 20:37:02.095684 : Organizing SE PeakMatrix from ArrowFile : SRR13927736.arrow, 0.214 mins elapsed.

2025-01-14 20:37:03.046221 : se SRR13927736, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 3534 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(3534): SRR13927736#GCTGTTCTCCTTGACC-1 SRR13927736#AGTTACGAGTAGGTCG-1 ... SRR13927736#TAAACCGCAAGTCTGT-1
  SRR13927736#TGCATTTGTTCTACCC-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2025-01-14 20:37:03.063607 : Completed PeakMatrix : SRR13927736(4 of 4), 1.02 mins elapsed.
2025-01-14 20:37:03.066911 : Organizing colData, 1.02 mins elapsed.
2025-01-14 20:37:03.770622 : Organizing rowData, 1.031 mins elapsed.
2025-01-14 20:37:03.790703 : Organizing rowRanges, 1.032 mins elapsed.
2025-01-14 20:37:03.816741 : Organizing Assays (1 of 1), 1.032 mins elapsed.
2025-01-14 20:37:07.242708 : Constructing SummarizedExperiment, 1.089 mins elapsed.
2025-01-14 20:37:27.958435 : Finished Matrix Creation, 1.435 mins elapsed.

2025-01-15 13:08:52.073469 : getMatrixFromProject Input-Parameters, Class = list

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


2025-01-15 13:08:54.506767 : Reading GeneIntegrationMatrix : SRR13927737(1 of 4), 0.041 mins elapsed.

2025-01-15 13:08:54.601988 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-15 13:08:58.003195 : featureDF SRR13927737, Class = DFrame
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

2025-01-15 13:08:58.044121 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927737.arrow, 0.057 mins elapsed.

2025-01-15 13:09:21.208524 : mat SRR13927737, Class = dgCMatrix
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


2025-01-15 13:09:23.079231 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927737.arrow, 0.475 mins elapsed.

2025-01-15 13:09:24.508709 : se SRR13927737, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 4439 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(4439): SRR13927737#GTCACGGGTCAACAGG-1 SRR13927737#TGCCTGTTCTTACTCA-1 ...
  SRR13927737#TTGCTATAGATGTTCC-1 SRR13927737#GCTCACTAGTTAGCGG-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2025-01-15 13:09:24.52487 : Completed GeneIntegrationMatrix : SRR13927737(1 of 4), 0.541 mins elapsed.
2025-01-15 13:09:24.526221 : Reading GeneIntegrationMatrix : SRR13927738(2 of 4), 0.541 mins elapsed.

2025-01-15 13:09:24.563021 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-15 13:09:26.260519 : featureDF SRR13927738, Class = DFrame
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

2025-01-15 13:09:26.30193 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927738.arrow, 0.029 mins elapsed.

2025-01-15 13:09:45.006574 : mat SRR13927738, Class = dgCMatrix
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


2025-01-15 13:09:46.749738 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927738.arrow, 0.37 mins elapsed.

2025-01-15 13:09:48.313793 : se SRR13927738, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 3736 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(3736): SRR13927738#TAGCATGAGGCTAAAT-1 SRR13927738#TGTTAGGTCGGTCAGC-1 ...
  SRR13927738#GCTCCTATCCGTGCGA-1 SRR13927738#GGCACGTAGTCGATAA-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2025-01-15 13:09:48.327331 : Completed GeneIntegrationMatrix : SRR13927738(2 of 4), 0.938 mins elapsed.
2025-01-15 13:09:48.328705 : Reading GeneIntegrationMatrix : SRR13927735(3 of 4), 0.938 mins elapsed.

2025-01-15 13:09:48.376842 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-15 13:09:50.460793 : featureDF SRR13927735, Class = DFrame
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

2025-01-15 13:09:50.50265 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927735.arrow, 0.035 mins elapsed.

2025-01-15 13:10:16.598433 : mat SRR13927735, Class = dgCMatrix
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


2025-01-15 13:10:18.107507 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927735.arrow, 0.496 mins elapsed.

2025-01-15 13:10:20.02038 : se SRR13927735, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 3813 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(3813): SRR13927735#TTATGTCTCCAGGTAT-1 SRR13927735#TATTGCTCATCAGAAA-1 ...
  SRR13927735#GATCATGGTGTTGTTG-1 SRR13927735#GTCCATCCATTTCACT-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2025-01-15 13:10:20.033488 : Completed GeneIntegrationMatrix : SRR13927735(3 of 4), 1.466 mins elapsed.
2025-01-15 13:10:20.064071 : Reading GeneIntegrationMatrix : SRR13927736(4 of 4), 1.467 mins elapsed.

2025-01-15 13:10:20.099953 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-15 13:10:21.97544 : featureDF SRR13927736, Class = DFrame
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

2025-01-15 13:10:22.015865 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927736.arrow, 0.032 mins elapsed.

2025-01-15 13:10:38.830887 : mat SRR13927736, Class = dgCMatrix
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


2025-01-15 13:10:40.694101 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927736.arrow, 0.343 mins elapsed.

2025-01-15 13:10:42.046066 : se SRR13927736, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 3534 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(3534): SRR13927736#GCTGTTCTCCTTGACC-1 SRR13927736#AGTTACGAGTAGGTCG-1 ...
  SRR13927736#TAAACCGCAAGTCTGT-1 SRR13927736#TGCATTTGTTCTACCC-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2025-01-15 13:10:42.055953 : Completed GeneIntegrationMatrix : SRR13927736(4 of 4), 1.833 mins elapsed.
2025-01-15 13:10:42.058766 : Organizing colData, 1.833 mins elapsed.
2025-01-15 13:10:42.89501 : Organizing rowData, 1.847 mins elapsed.
2025-01-15 13:10:42.904487 : Organizing rowRanges, 1.848 mins elapsed.
2025-01-15 13:10:42.919944 : Organizing Assays (1 of 1), 1.848 mins elapsed.
2025-01-15 13:10:52.93399 : Constructing SummarizedExperiment, 2.015 mins elapsed.
2025-01-15 13:10:56.588112 : Finished Matrix Creation, 2.076 mins elapsed.

2025-01-15 13:10:56.592596 : getMatrixFromProject Input-Parameters, Class = list

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


2025-01-15 13:11:14.452027 : Reading PeakMatrix : SRR13927737(1 of 4), 0.298 mins elapsed.

2025-01-15 13:11:14.490062 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-15 13:11:16.159262 : featureDF SRR13927737, Class = DFrame
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

2025-01-15 13:11:16.265506 : Getting PeakMatrix from ArrowFile : SRR13927737.arrow, 0.03 mins elapsed.

2025-01-15 13:11:31.281547 : mat SRR13927737, Class = dgCMatrix
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


2025-01-15 13:11:32.838624 : Organizing SE PeakMatrix from ArrowFile : SRR13927737.arrow, 0.306 mins elapsed.

2025-01-15 13:11:34.16229 : se SRR13927737, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 4439 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(4439): SRR13927737#GTCACGGGTCAACAGG-1 SRR13927737#TGCCTGTTCTTACTCA-1 ...
  SRR13927737#TTGCTATAGATGTTCC-1 SRR13927737#GCTCACTAGTTAGCGG-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2025-01-15 13:11:34.179621 : Completed PeakMatrix : SRR13927737(1 of 4), 0.626 mins elapsed.
2025-01-15 13:11:34.180669 : Reading PeakMatrix : SRR13927738(2 of 4), 0.626 mins elapsed.

2025-01-15 13:11:34.211897 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-15 13:11:36.216133 : featureDF SRR13927738, Class = DFrame
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

2025-01-15 13:11:37.081265 : Getting PeakMatrix from ArrowFile : SRR13927738.arrow, 0.048 mins elapsed.

2025-01-15 13:11:52.101394 : mat SRR13927738, Class = dgCMatrix
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


2025-01-15 13:11:53.700849 : Organizing SE PeakMatrix from ArrowFile : SRR13927738.arrow, 0.325 mins elapsed.

2025-01-15 13:11:55.782731 : se SRR13927738, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 3736 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(3736): SRR13927738#TAGCATGAGGCTAAAT-1 SRR13927738#TGTTAGGTCGGTCAGC-1 ...
  SRR13927738#GCTCCTATCCGTGCGA-1 SRR13927738#GGCACGTAGTCGATAA-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2025-01-15 13:11:55.799991 : Completed PeakMatrix : SRR13927738(2 of 4), 0.987 mins elapsed.
2025-01-15 13:11:55.801302 : Reading PeakMatrix : SRR13927735(3 of 4), 0.987 mins elapsed.

2025-01-15 13:11:55.856596 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-15 13:11:57.865708 : featureDF SRR13927735, Class = DFrame
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

2025-01-15 13:11:57.980444 : Getting PeakMatrix from ArrowFile : SRR13927735.arrow, 0.035 mins elapsed.

2025-01-15 13:12:10.665905 : mat SRR13927735, Class = dgCMatrix
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


2025-01-15 13:12:12.349972 : Organizing SE PeakMatrix from ArrowFile : SRR13927735.arrow, 0.275 mins elapsed.

2025-01-15 13:12:14.572149 : se SRR13927735, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 3813 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(3813): SRR13927735#TTATGTCTCCAGGTAT-1 SRR13927735#TATTGCTCATCAGAAA-1 ...
  SRR13927735#GATCATGGTGTTGTTG-1 SRR13927735#GTCCATCCATTTCACT-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2025-01-15 13:12:14.581542 : Completed PeakMatrix : SRR13927735(3 of 4), 1.3 mins elapsed.
2025-01-15 13:12:14.582232 : Reading PeakMatrix : SRR13927736(4 of 4), 1.3 mins elapsed.

2025-01-15 13:12:14.616158 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-15 13:12:16.696691 : featureDF SRR13927736, Class = DFrame
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

2025-01-15 13:12:16.755594 : Getting PeakMatrix from ArrowFile : SRR13927736.arrow, 0.036 mins elapsed.

2025-01-15 13:12:28.665842 : mat SRR13927736, Class = dgCMatrix
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


2025-01-15 13:12:30.222614 : Organizing SE PeakMatrix from ArrowFile : SRR13927736.arrow, 0.26 mins elapsed.

2025-01-15 13:12:31.802909 : se SRR13927736, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 3534 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(3534): SRR13927736#GCTGTTCTCCTTGACC-1 SRR13927736#AGTTACGAGTAGGTCG-1 ...
  SRR13927736#TAAACCGCAAGTCTGT-1 SRR13927736#TGCATTTGTTCTACCC-1
colData names(21): BlacklistRatio NucleosomeRatio ... ReadsInPeaks FRIP

2025-01-15 13:12:31.81567 : Completed PeakMatrix : SRR13927736(4 of 4), 1.587 mins elapsed.
2025-01-15 13:12:31.818898 : Organizing colData, 1.587 mins elapsed.
2025-01-15 13:12:32.780385 : Organizing rowData, 1.603 mins elapsed.
2025-01-15 13:12:32.804039 : Organizing rowRanges, 1.603 mins elapsed.
2025-01-15 13:12:32.83654 : Organizing Assays (1 of 1), 1.604 mins elapsed.
2025-01-15 13:12:37.315318 : Constructing SummarizedExperiment, 1.679 mins elapsed.
2025-01-15 13:13:00.373452 : Finished Matrix Creation, 2.063 mins elapsed.

2025-01-15 13:50:06.447681 : getMatrixFromProject Input-Parameters, Class = list

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


2025-01-15 13:50:11.854682 : Reading GeneIntegrationMatrix : SRR13927737(1 of 4), 0.09 mins elapsed.

2025-01-15 13:50:11.926521 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-15 13:50:16.543164 : featureDF SRR13927737, Class = DFrame
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

2025-01-15 13:50:16.559761 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927737.arrow, 0.077 mins elapsed.

2025-01-15 13:50:33.302807 : mat SRR13927737, Class = dgCMatrix
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


2025-01-15 13:50:34.303396 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927737.arrow, 0.373 mins elapsed.

2025-01-15 13:50:35.614869 : se SRR13927737, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 4439 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(4439): SRR13927737#GTCACGGGTCAACAGG-1 SRR13927737#TGCCTGTTCTTACTCA-1 ... SRR13927737#TTGCTATAGATGTTCC-1
  SRR13927737#GCTCACTAGTTAGCGG-1
colData names(21): BlacklistRatio nDiFrags ... ReadsInPeaks FRIP

2025-01-15 13:50:35.620651 : Completed GeneIntegrationMatrix : SRR13927737(1 of 4), 0.486 mins elapsed.
2025-01-15 13:50:35.621379 : Reading GeneIntegrationMatrix : SRR13927738(2 of 4), 0.486 mins elapsed.

2025-01-15 13:50:35.682119 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-15 13:50:40.15843 : featureDF SRR13927738, Class = DFrame
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

2025-01-15 13:50:40.175604 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927738.arrow, 0.075 mins elapsed.

2025-01-15 13:50:54.942866 : mat SRR13927738, Class = dgCMatrix
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


2025-01-15 13:50:55.783097 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927738.arrow, 0.335 mins elapsed.

2025-01-15 13:50:58.143372 : se SRR13927738, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 3736 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(3736): SRR13927738#TAGCATGAGGCTAAAT-1 SRR13927738#TGTTAGGTCGGTCAGC-1 ... SRR13927738#GCTCCTATCCGTGCGA-1
  SRR13927738#GGCACGTAGTCGATAA-1
colData names(21): BlacklistRatio nDiFrags ... ReadsInPeaks FRIP

2025-01-15 13:50:58.149585 : Completed GeneIntegrationMatrix : SRR13927738(2 of 4), 0.862 mins elapsed.
2025-01-15 13:50:58.150493 : Reading GeneIntegrationMatrix : SRR13927735(3 of 4), 0.862 mins elapsed.

2025-01-15 13:50:58.271513 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-15 13:51:01.290965 : featureDF SRR13927735, Class = DFrame
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

2025-01-15 13:51:01.307231 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927735.arrow, 0.051 mins elapsed.

2025-01-15 13:51:16.695758 : mat SRR13927735, Class = dgCMatrix
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


2025-01-15 13:51:17.611996 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927735.arrow, 0.322 mins elapsed.

2025-01-15 13:51:19.271425 : se SRR13927735, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 3813 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(3813): SRR13927735#TTATGTCTCCAGGTAT-1 SRR13927735#TATTGCTCATCAGAAA-1 ... SRR13927735#GATCATGGTGTTGTTG-1
  SRR13927735#GTCCATCCATTTCACT-1
colData names(21): BlacklistRatio nDiFrags ... ReadsInPeaks FRIP

2025-01-15 13:51:19.27779 : Completed GeneIntegrationMatrix : SRR13927735(3 of 4), 1.214 mins elapsed.
2025-01-15 13:51:19.278686 : Reading GeneIntegrationMatrix : SRR13927736(4 of 4), 1.214 mins elapsed.

2025-01-15 13:51:19.772795 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-15 13:51:24.446342 : featureDF SRR13927736, Class = DFrame
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

2025-01-15 13:51:24.461356 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927736.arrow, 0.078 mins elapsed.

2025-01-15 13:51:42.100786 : mat SRR13927736, Class = dgCMatrix
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


2025-01-15 13:51:42.935054 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927736.arrow, 0.386 mins elapsed.

2025-01-15 13:51:45.976214 : se SRR13927736, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 3534 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(3534): SRR13927736#GCTGTTCTCCTTGACC-1 SRR13927736#AGTTACGAGTAGGTCG-1 ... SRR13927736#TAAACCGCAAGTCTGT-1
  SRR13927736#TGCATTTGTTCTACCC-1
colData names(21): BlacklistRatio nDiFrags ... ReadsInPeaks FRIP

2025-01-15 13:51:45.981067 : Completed GeneIntegrationMatrix : SRR13927736(4 of 4), 1.659 mins elapsed.
2025-01-15 13:51:45.982721 : Organizing colData, 1.659 mins elapsed.
2025-01-15 13:51:46.676579 : Organizing rowData, 1.671 mins elapsed.
2025-01-15 13:51:46.680061 : Organizing rowRanges, 1.671 mins elapsed.
2025-01-15 13:51:46.685169 : Organizing Assays (1 of 1), 1.671 mins elapsed.
2025-01-15 13:51:51.917075 : Constructing SummarizedExperiment, 1.758 mins elapsed.
2025-01-15 13:51:54.751247 : Finished Matrix Creation, 1.801 mins elapsed.

2025-01-15 13:51:54.769736 : getMatrixFromProject Input-Parameters, Class = list

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


2025-01-15 13:51:59.513337 : Reading PeakMatrix : SRR13927737(1 of 4), 0.079 mins elapsed.

2025-01-15 13:51:59.597148 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-15 13:52:04.508653 : featureDF SRR13927737, Class = DFrame
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

2025-01-15 13:52:04.554778 : Getting PeakMatrix from ArrowFile : SRR13927737.arrow, 0.083 mins elapsed.

2025-01-15 13:52:15.266669 : mat SRR13927737, Class = dgCMatrix
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


2025-01-15 13:52:15.482879 : Organizing SE PeakMatrix from ArrowFile : SRR13927737.arrow, 0.265 mins elapsed.

2025-01-15 13:52:16.493877 : se SRR13927737, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 4439 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(4439): SRR13927737#GTCACGGGTCAACAGG-1 SRR13927737#TGCCTGTTCTTACTCA-1 ... SRR13927737#TTGCTATAGATGTTCC-1
  SRR13927737#GCTCACTAGTTAGCGG-1
colData names(21): BlacklistRatio nDiFrags ... ReadsInPeaks FRIP

2025-01-15 13:52:16.502035 : Completed PeakMatrix : SRR13927737(1 of 4), 0.362 mins elapsed.
2025-01-15 13:52:16.502791 : Reading PeakMatrix : SRR13927738(2 of 4), 0.362 mins elapsed.

2025-01-15 13:52:17.210998 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-15 13:52:21.782705 : featureDF SRR13927738, Class = DFrame
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

2025-01-15 13:52:21.842691 : Getting PeakMatrix from ArrowFile : SRR13927738.arrow, 0.077 mins elapsed.

2025-01-15 13:52:33.529044 : mat SRR13927738, Class = dgCMatrix
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


2025-01-15 13:52:33.746312 : Organizing SE PeakMatrix from ArrowFile : SRR13927738.arrow, 0.276 mins elapsed.

2025-01-15 13:52:35.223044 : se SRR13927738, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 3736 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(3736): SRR13927738#TAGCATGAGGCTAAAT-1 SRR13927738#TGTTAGGTCGGTCAGC-1 ... SRR13927738#GCTCCTATCCGTGCGA-1
  SRR13927738#GGCACGTAGTCGATAA-1
colData names(21): BlacklistRatio nDiFrags ... ReadsInPeaks FRIP

2025-01-15 13:52:35.229708 : Completed PeakMatrix : SRR13927738(2 of 4), 0.674 mins elapsed.
2025-01-15 13:52:35.230435 : Reading PeakMatrix : SRR13927735(3 of 4), 0.674 mins elapsed.

2025-01-15 13:52:35.949733 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-15 13:52:41.279432 : featureDF SRR13927735, Class = DFrame
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

2025-01-15 13:52:41.323362 : Getting PeakMatrix from ArrowFile : SRR13927735.arrow, 0.09 mins elapsed.

2025-01-15 13:52:51.30818 : mat SRR13927735, Class = dgCMatrix
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


2025-01-15 13:52:52.048168 : Organizing SE PeakMatrix from ArrowFile : SRR13927735.arrow, 0.268 mins elapsed.

2025-01-15 13:52:53.325868 : se SRR13927735, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 3813 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(3813): SRR13927735#TTATGTCTCCAGGTAT-1 SRR13927735#TATTGCTCATCAGAAA-1 ... SRR13927735#GATCATGGTGTTGTTG-1
  SRR13927735#GTCCATCCATTTCACT-1
colData names(21): BlacklistRatio nDiFrags ... ReadsInPeaks FRIP

2025-01-15 13:52:53.3324 : Completed PeakMatrix : SRR13927735(3 of 4), 0.976 mins elapsed.
2025-01-15 13:52:53.333161 : Reading PeakMatrix : SRR13927736(4 of 4), 0.976 mins elapsed.

2025-01-15 13:52:53.90035 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-15 13:52:58.662097 : featureDF SRR13927736, Class = DFrame
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

2025-01-15 13:52:58.707628 : Getting PeakMatrix from ArrowFile : SRR13927736.arrow, 0.08 mins elapsed.

2025-01-15 13:53:08.301543 : mat SRR13927736, Class = dgCMatrix
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


2025-01-15 13:53:08.494918 : Organizing SE PeakMatrix from ArrowFile : SRR13927736.arrow, 0.243 mins elapsed.

2025-01-15 13:53:09.159264 : se SRR13927736, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 3534 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(3534): SRR13927736#GCTGTTCTCCTTGACC-1 SRR13927736#AGTTACGAGTAGGTCG-1 ... SRR13927736#TAAACCGCAAGTCTGT-1
  SRR13927736#TGCATTTGTTCTACCC-1
colData names(21): BlacklistRatio nDiFrags ... ReadsInPeaks FRIP

2025-01-15 13:53:09.193579 : Completed PeakMatrix : SRR13927736(4 of 4), 1.24 mins elapsed.
2025-01-15 13:53:09.195686 : Organizing colData, 1.24 mins elapsed.
2025-01-15 13:53:09.60848 : Organizing rowData, 1.247 mins elapsed.
2025-01-15 13:53:09.617087 : Organizing rowRanges, 1.247 mins elapsed.
2025-01-15 13:53:09.628994 : Organizing Assays (1 of 1), 1.248 mins elapsed.
2025-01-15 13:53:12.356164 : Constructing SummarizedExperiment, 1.293 mins elapsed.
2025-01-15 13:53:22.27318 : Finished Matrix Creation, 1.458 mins elapsed.

2025-01-15 14:29:13.890402 : getMatrixFromProject Input-Parameters, Class = list

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


2025-01-15 14:29:20.136499 : Reading GeneIntegrationMatrix : SRR13927737(1 of 4), 0.104 mins elapsed.

2025-01-15 14:29:20.233609 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-15 14:29:24.620471 : featureDF SRR13927737, Class = DFrame
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

2025-01-15 14:29:24.65243 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927737.arrow, 0.074 mins elapsed.

2025-01-15 14:29:55.02686 : mat SRR13927737, Class = dgCMatrix
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


2025-01-15 14:29:55.680845 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927737.arrow, 0.591 mins elapsed.

2025-01-15 14:29:58.846697 : se SRR13927737, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 4439 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(4439): SRR13927737#GTCACGGGTCAACAGG-1 SRR13927737#TGCCTGTTCTTACTCA-1 ... SRR13927737#TTGCTATAGATGTTCC-1
  SRR13927737#GCTCACTAGTTAGCGG-1
colData names(21): BlacklistRatio nDiFrags ... ReadsInPeaks FRIP

2025-01-15 14:29:58.852227 : Completed GeneIntegrationMatrix : SRR13927737(1 of 4), 0.749 mins elapsed.
2025-01-15 14:29:58.852769 : Reading GeneIntegrationMatrix : SRR13927738(2 of 4), 0.749 mins elapsed.

2025-01-15 14:29:58.906046 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-15 14:30:03.490195 : featureDF SRR13927738, Class = DFrame
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

2025-01-15 14:30:03.508149 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927738.arrow, 0.077 mins elapsed.

2025-01-15 14:30:19.294295 : mat SRR13927738, Class = dgCMatrix
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


2025-01-15 14:30:19.87516 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927738.arrow, 0.349 mins elapsed.

2025-01-15 14:30:20.932859 : se SRR13927738, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 3736 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(3736): SRR13927738#TAGCATGAGGCTAAAT-1 SRR13927738#TGTTAGGTCGGTCAGC-1 ... SRR13927738#GCTCCTATCCGTGCGA-1
  SRR13927738#GGCACGTAGTCGATAA-1
colData names(21): BlacklistRatio nDiFrags ... ReadsInPeaks FRIP

2025-01-15 14:30:20.938056 : Completed GeneIntegrationMatrix : SRR13927738(2 of 4), 1.117 mins elapsed.
2025-01-15 14:30:20.93864 : Reading GeneIntegrationMatrix : SRR13927735(3 of 4), 1.118 mins elapsed.

2025-01-15 14:30:21.046259 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-15 14:30:24.642595 : featureDF SRR13927735, Class = DFrame
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

2025-01-15 14:30:24.759096 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927735.arrow, 0.062 mins elapsed.

2025-01-15 14:30:39.815155 : mat SRR13927735, Class = dgCMatrix
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


2025-01-15 14:30:40.342207 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927735.arrow, 0.322 mins elapsed.

2025-01-15 14:30:42.330884 : se SRR13927735, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 3813 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(3813): SRR13927735#TTATGTCTCCAGGTAT-1 SRR13927735#TATTGCTCATCAGAAA-1 ... SRR13927735#GATCATGGTGTTGTTG-1
  SRR13927735#GTCCATCCATTTCACT-1
colData names(21): BlacklistRatio nDiFrags ... ReadsInPeaks FRIP

2025-01-15 14:30:42.335284 : Completed GeneIntegrationMatrix : SRR13927735(3 of 4), 1.474 mins elapsed.
2025-01-15 14:30:42.335897 : Reading GeneIntegrationMatrix : SRR13927736(4 of 4), 1.474 mins elapsed.

2025-01-15 14:30:42.400312 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-15 14:30:47.406859 : featureDF SRR13927736, Class = DFrame
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

2025-01-15 14:30:47.42118 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927736.arrow, 0.084 mins elapsed.

2025-01-15 14:30:58.589366 : mat SRR13927736, Class = dgCMatrix
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


2025-01-15 14:30:59.898842 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927736.arrow, 0.292 mins elapsed.

2025-01-15 14:31:00.507728 : se SRR13927736, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 3534 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(3534): SRR13927736#GCTGTTCTCCTTGACC-1 SRR13927736#AGTTACGAGTAGGTCG-1 ... SRR13927736#TAAACCGCAAGTCTGT-1
  SRR13927736#TGCATTTGTTCTACCC-1
colData names(21): BlacklistRatio nDiFrags ... ReadsInPeaks FRIP

2025-01-15 14:31:00.512745 : Completed GeneIntegrationMatrix : SRR13927736(4 of 4), 1.777 mins elapsed.
2025-01-15 14:31:00.514318 : Organizing colData, 1.777 mins elapsed.
2025-01-15 14:31:00.638406 : Organizing rowData, 1.779 mins elapsed.
2025-01-15 14:31:00.641506 : Organizing rowRanges, 1.779 mins elapsed.
2025-01-15 14:31:00.646749 : Organizing Assays (1 of 1), 1.779 mins elapsed.
2025-01-15 14:31:05.712767 : Constructing SummarizedExperiment, 1.864 mins elapsed.
2025-01-15 14:31:07.374287 : Finished Matrix Creation, 1.89 mins elapsed.

2025-01-15 14:31:07.390563 : getMatrixFromProject Input-Parameters, Class = list

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


2025-01-15 14:31:14.310321 : Reading PeakMatrix : SRR13927737(1 of 4), 0.115 mins elapsed.

2025-01-15 14:31:14.335425 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-15 14:31:16.021257 : featureDF SRR13927737, Class = DFrame
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

2025-01-15 14:31:16.077767 : Getting PeakMatrix from ArrowFile : SRR13927737.arrow, 0.029 mins elapsed.

2025-01-15 14:31:26.791203 : mat SRR13927737, Class = dgCMatrix
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


2025-01-15 14:31:27.482953 : Organizing SE PeakMatrix from ArrowFile : SRR13927737.arrow, 0.219 mins elapsed.

2025-01-15 14:31:28.524453 : se SRR13927737, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 4439 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(4439): SRR13927737#GTCACGGGTCAACAGG-1 SRR13927737#TGCCTGTTCTTACTCA-1 ... SRR13927737#TTGCTATAGATGTTCC-1
  SRR13927737#GCTCACTAGTTAGCGG-1
colData names(21): BlacklistRatio nDiFrags ... ReadsInPeaks FRIP

2025-01-15 14:31:28.534234 : Completed PeakMatrix : SRR13927737(1 of 4), 0.352 mins elapsed.
2025-01-15 14:31:28.534858 : Reading PeakMatrix : SRR13927738(2 of 4), 0.352 mins elapsed.

2025-01-15 14:31:28.602911 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-15 14:31:31.472711 : featureDF SRR13927738, Class = DFrame
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

2025-01-15 14:31:31.653487 : Getting PeakMatrix from ArrowFile : SRR13927738.arrow, 0.051 mins elapsed.

2025-01-15 14:31:40.575832 : mat SRR13927738, Class = dgCMatrix
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


2025-01-15 14:31:40.937273 : Organizing SE PeakMatrix from ArrowFile : SRR13927738.arrow, 0.206 mins elapsed.

2025-01-15 14:31:41.43126 : se SRR13927738, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 3736 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(3736): SRR13927738#TAGCATGAGGCTAAAT-1 SRR13927738#TGTTAGGTCGGTCAGC-1 ... SRR13927738#GCTCCTATCCGTGCGA-1
  SRR13927738#GGCACGTAGTCGATAA-1
colData names(21): BlacklistRatio nDiFrags ... ReadsInPeaks FRIP

2025-01-15 14:31:41.438072 : Completed PeakMatrix : SRR13927738(2 of 4), 0.567 mins elapsed.
2025-01-15 14:31:41.438689 : Reading PeakMatrix : SRR13927735(3 of 4), 0.567 mins elapsed.

2025-01-15 14:31:41.494656 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-15 14:31:46.455504 : featureDF SRR13927735, Class = DFrame
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

2025-01-15 14:31:46.498116 : Getting PeakMatrix from ArrowFile : SRR13927735.arrow, 0.083 mins elapsed.

2025-01-15 14:31:54.336324 : mat SRR13927735, Class = dgCMatrix
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


2025-01-15 14:31:54.494797 : Organizing SE PeakMatrix from ArrowFile : SRR13927735.arrow, 0.217 mins elapsed.

2025-01-15 14:31:55.035551 : se SRR13927735, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 3813 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(3813): SRR13927735#TTATGTCTCCAGGTAT-1 SRR13927735#TATTGCTCATCAGAAA-1 ... SRR13927735#GATCATGGTGTTGTTG-1
  SRR13927735#GTCCATCCATTTCACT-1
colData names(21): BlacklistRatio nDiFrags ... ReadsInPeaks FRIP

2025-01-15 14:31:55.041387 : Completed PeakMatrix : SRR13927735(3 of 4), 0.794 mins elapsed.
2025-01-15 14:31:55.041958 : Reading PeakMatrix : SRR13927736(4 of 4), 0.794 mins elapsed.

2025-01-15 14:31:55.808048 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-15 14:31:59.771028 : featureDF SRR13927736, Class = DFrame
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

2025-01-15 14:31:59.817923 : Getting PeakMatrix from ArrowFile : SRR13927736.arrow, 0.067 mins elapsed.

2025-01-15 14:32:07.202194 : mat SRR13927736, Class = dgCMatrix
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


2025-01-15 14:32:07.53902 : Organizing SE PeakMatrix from ArrowFile : SRR13927736.arrow, 0.196 mins elapsed.

2025-01-15 14:32:08.060415 : se SRR13927736, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 3534 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(3534): SRR13927736#GCTGTTCTCCTTGACC-1 SRR13927736#AGTTACGAGTAGGTCG-1 ... SRR13927736#TAAACCGCAAGTCTGT-1
  SRR13927736#TGCATTTGTTCTACCC-1
colData names(21): BlacklistRatio nDiFrags ... ReadsInPeaks FRIP

2025-01-15 14:32:08.066626 : Completed PeakMatrix : SRR13927736(4 of 4), 1.011 mins elapsed.
2025-01-15 14:32:08.068305 : Organizing colData, 1.011 mins elapsed.
2025-01-15 14:32:08.201024 : Organizing rowData, 1.013 mins elapsed.
2025-01-15 14:32:08.207996 : Organizing rowRanges, 1.014 mins elapsed.
2025-01-15 14:32:08.217479 : Organizing Assays (1 of 1), 1.014 mins elapsed.
2025-01-15 14:32:09.62476 : Constructing SummarizedExperiment, 1.037 mins elapsed.
2025-01-15 14:32:19.383537 : Finished Matrix Creation, 1.2 mins elapsed.

2025-01-15 15:58:24.138047 : getMatrixFromProject Input-Parameters, Class = list

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


2025-01-15 15:58:29.441488 : Reading GeneIntegrationMatrix : SRR13927737(1 of 4), 0.088 mins elapsed.

2025-01-15 15:58:30.188501 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-15 15:58:35.123216 : featureDF SRR13927737, Class = DFrame
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

2025-01-15 15:58:35.139672 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927737.arrow, 0.083 mins elapsed.

2025-01-15 15:58:51.789695 : mat SRR13927737, Class = dgCMatrix
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


2025-01-15 15:58:52.980108 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927737.arrow, 0.38 mins elapsed.

2025-01-15 15:58:54.463151 : se SRR13927737, Class = SummarizedExperiment
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

2025-01-15 15:58:54.468902 : Completed GeneIntegrationMatrix : SRR13927737(1 of 4), 0.506 mins elapsed.
2025-01-15 15:58:54.469602 : Reading GeneIntegrationMatrix : SRR13927738(2 of 4), 0.506 mins elapsed.

2025-01-15 15:58:54.510936 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-15 15:58:59.434281 : featureDF SRR13927738, Class = DFrame
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

2025-01-15 15:58:59.45125 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927738.arrow, 0.082 mins elapsed.

2025-01-15 15:59:13.63129 : mat SRR13927738, Class = dgCMatrix
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


2025-01-15 15:59:15.169901 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927738.arrow, 0.344 mins elapsed.

2025-01-15 15:59:17.925105 : se SRR13927738, Class = SummarizedExperiment
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

2025-01-15 15:59:17.930209 : Completed GeneIntegrationMatrix : SRR13927738(2 of 4), 0.897 mins elapsed.
2025-01-15 15:59:17.930924 : Reading GeneIntegrationMatrix : SRR13927735(3 of 4), 0.897 mins elapsed.

2025-01-15 15:59:18.076242 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-15 15:59:22.416329 : featureDF SRR13927735, Class = DFrame
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

2025-01-15 15:59:22.433826 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927735.arrow, 0.073 mins elapsed.

2025-01-15 15:59:39.394294 : mat SRR13927735, Class = dgCMatrix
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


2025-01-15 15:59:40.533908 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927735.arrow, 0.374 mins elapsed.

2025-01-15 15:59:44.474904 : se SRR13927735, Class = SummarizedExperiment
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

2025-01-15 15:59:44.480752 : Completed GeneIntegrationMatrix : SRR13927735(3 of 4), 1.339 mins elapsed.
2025-01-15 15:59:44.482625 : Reading GeneIntegrationMatrix : SRR13927736(4 of 4), 1.339 mins elapsed.

2025-01-15 15:59:44.568064 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-15 15:59:49.556014 : featureDF SRR13927736, Class = DFrame
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

2025-01-15 15:59:49.578213 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927736.arrow, 0.084 mins elapsed.

2025-01-15 16:00:05.373183 : mat SRR13927736, Class = dgCMatrix
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


2025-01-15 16:00:06.463752 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927736.arrow, 0.365 mins elapsed.

2025-01-15 16:00:07.780287 : se SRR13927736, Class = SummarizedExperiment
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

2025-01-15 16:00:07.786269 : Completed GeneIntegrationMatrix : SRR13927736(4 of 4), 1.727 mins elapsed.
2025-01-15 16:00:07.788707 : Organizing colData, 1.728 mins elapsed.
2025-01-15 16:00:08.329959 : Organizing rowData, 1.737 mins elapsed.
2025-01-15 16:00:08.334394 : Organizing rowRanges, 1.737 mins elapsed.
2025-01-15 16:00:08.340796 : Organizing Assays (1 of 1), 1.737 mins elapsed.
2025-01-15 16:00:14.957623 : Constructing SummarizedExperiment, 1.847 mins elapsed.
2025-01-15 16:00:18.155966 : Finished Matrix Creation, 1.9 mins elapsed.

2025-01-15 16:00:18.15831 : getMatrixFromProject Input-Parameters, Class = list

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


2025-01-15 16:00:24.739678 : Reading PeakMatrix : SRR13927737(1 of 4), 0.11 mins elapsed.

2025-01-15 16:00:24.80841 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-15 16:00:27.749207 : featureDF SRR13927737, Class = DFrame
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

2025-01-15 16:00:27.797919 : Getting PeakMatrix from ArrowFile : SRR13927737.arrow, 0.05 mins elapsed.

2025-01-15 16:00:38.932401 : mat SRR13927737, Class = dgCMatrix
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


2025-01-15 16:00:39.744876 : Organizing SE PeakMatrix from ArrowFile : SRR13927737.arrow, 0.249 mins elapsed.

2025-01-15 16:00:40.474612 : se SRR13927737, Class = RangedSummarizedExperiment
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

2025-01-15 16:00:40.483062 : Completed PeakMatrix : SRR13927737(1 of 4), 0.372 mins elapsed.
2025-01-15 16:00:40.483768 : Reading PeakMatrix : SRR13927738(2 of 4), 0.372 mins elapsed.

2025-01-15 16:00:40.58234 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-15 16:00:44.390681 : featureDF SRR13927738, Class = DFrame
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

2025-01-15 16:00:44.436444 : Getting PeakMatrix from ArrowFile : SRR13927738.arrow, 0.064 mins elapsed.

2025-01-15 16:00:55.290354 : mat SRR13927738, Class = dgCMatrix
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


2025-01-15 16:00:56.077422 : Organizing SE PeakMatrix from ArrowFile : SRR13927738.arrow, 0.258 mins elapsed.

2025-01-15 16:00:56.978274 : se SRR13927738, Class = RangedSummarizedExperiment
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

2025-01-15 16:00:56.984734 : Completed PeakMatrix : SRR13927738(2 of 4), 0.647 mins elapsed.
2025-01-15 16:00:56.985478 : Reading PeakMatrix : SRR13927735(3 of 4), 0.647 mins elapsed.

2025-01-15 16:00:58.059389 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-15 16:01:01.218904 : featureDF SRR13927735, Class = DFrame
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

2025-01-15 16:01:01.262895 : Getting PeakMatrix from ArrowFile : SRR13927735.arrow, 0.053 mins elapsed.

2025-01-15 16:01:13.253499 : mat SRR13927735, Class = dgCMatrix
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


2025-01-15 16:01:13.446213 : Organizing SE PeakMatrix from ArrowFile : SRR13927735.arrow, 0.256 mins elapsed.

2025-01-15 16:01:14.896731 : se SRR13927735, Class = RangedSummarizedExperiment
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

2025-01-15 16:01:14.904023 : Completed PeakMatrix : SRR13927735(3 of 4), 0.946 mins elapsed.
2025-01-15 16:01:14.904792 : Reading PeakMatrix : SRR13927736(4 of 4), 0.946 mins elapsed.

2025-01-15 16:01:14.982094 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-15 16:01:19.642665 : featureDF SRR13927736, Class = DFrame
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

2025-01-15 16:01:19.686303 : Getting PeakMatrix from ArrowFile : SRR13927736.arrow, 0.078 mins elapsed.

2025-01-15 16:01:31.276731 : mat SRR13927736, Class = dgCMatrix
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


2025-01-15 16:01:32.149876 : Organizing SE PeakMatrix from ArrowFile : SRR13927736.arrow, 0.286 mins elapsed.

2025-01-15 16:01:47.711194 : se SRR13927736, Class = RangedSummarizedExperiment
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

2025-01-15 16:01:47.71875 : Completed PeakMatrix : SRR13927736(4 of 4), 1.493 mins elapsed.
2025-01-15 16:01:47.720766 : Organizing colData, 1.493 mins elapsed.
2025-01-15 16:01:47.862855 : Organizing rowData, 1.495 mins elapsed.
2025-01-15 16:01:47.870394 : Organizing rowRanges, 1.495 mins elapsed.
2025-01-15 16:01:47.881086 : Organizing Assays (1 of 1), 1.495 mins elapsed.
2025-01-15 16:01:53.141679 : Constructing SummarizedExperiment, 1.583 mins elapsed.
2025-01-15 16:02:02.863703 : Finished Matrix Creation, 1.745 mins elapsed.

2025-01-15 18:53:55.838521 : getMatrixFromProject Input-Parameters, Class = list

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


2025-01-15 18:54:09.425138 : Reading GeneIntegrationMatrix : SRR13927737(1 of 4), 0.226 mins elapsed.

2025-01-15 18:54:09.494374 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-15 18:54:12.169848 : featureDF SRR13927737, Class = DFrame
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

2025-01-15 18:54:12.19198 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927737.arrow, 0.045 mins elapsed.

2025-01-15 18:54:31.563575 : mat SRR13927737, Class = dgCMatrix
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


2025-01-15 18:54:31.922139 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927737.arrow, 0.374 mins elapsed.

2025-01-15 18:54:32.486392 : se SRR13927737, Class = SummarizedExperiment
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

2025-01-15 18:54:32.491001 : Completed GeneIntegrationMatrix : SRR13927737(1 of 4), 0.611 mins elapsed.
2025-01-15 18:54:32.491516 : Reading GeneIntegrationMatrix : SRR13927738(2 of 4), 0.611 mins elapsed.

2025-01-15 18:54:32.559007 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-15 18:54:35.379161 : featureDF SRR13927738, Class = DFrame
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

2025-01-15 18:54:35.391131 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927738.arrow, 0.047 mins elapsed.

2025-01-15 18:54:46.037519 : mat SRR13927738, Class = dgCMatrix
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


2025-01-15 18:54:46.333771 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927738.arrow, 0.23 mins elapsed.

2025-01-15 18:54:46.815906 : se SRR13927738, Class = SummarizedExperiment
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

2025-01-15 18:54:46.819417 : Completed GeneIntegrationMatrix : SRR13927738(2 of 4), 0.85 mins elapsed.
2025-01-15 18:54:46.819919 : Reading GeneIntegrationMatrix : SRR13927735(3 of 4), 0.85 mins elapsed.

2025-01-15 18:54:46.949928 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-15 18:54:49.58146 : featureDF SRR13927735, Class = DFrame
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

2025-01-15 18:54:49.592982 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927735.arrow, 0.044 mins elapsed.

2025-01-15 18:55:00.510371 : mat SRR13927735, Class = dgCMatrix
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


2025-01-15 18:55:00.825464 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927735.arrow, 0.231 mins elapsed.

2025-01-15 18:55:01.230707 : se SRR13927735, Class = SummarizedExperiment
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

2025-01-15 18:55:01.234194 : Completed GeneIntegrationMatrix : SRR13927735(3 of 4), 1.09 mins elapsed.
2025-01-15 18:55:01.234691 : Reading GeneIntegrationMatrix : SRR13927736(4 of 4), 1.09 mins elapsed.

2025-01-15 18:55:01.295041 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-15 18:55:04.163533 : featureDF SRR13927736, Class = DFrame
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

2025-01-15 18:55:04.174103 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927736.arrow, 0.048 mins elapsed.

2025-01-15 18:55:13.405236 : mat SRR13927736, Class = dgCMatrix
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


2025-01-15 18:55:13.677128 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927736.arrow, 0.206 mins elapsed.

2025-01-15 18:55:13.945042 : se SRR13927736, Class = SummarizedExperiment
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

2025-01-15 18:55:13.948746 : Completed GeneIntegrationMatrix : SRR13927736(4 of 4), 1.302 mins elapsed.
2025-01-15 18:55:13.949964 : Organizing colData, 1.302 mins elapsed.
2025-01-15 18:55:14.043864 : Organizing rowData, 1.303 mins elapsed.
2025-01-15 18:55:14.046527 : Organizing rowRanges, 1.304 mins elapsed.
2025-01-15 18:55:14.050771 : Organizing Assays (1 of 1), 1.304 mins elapsed.
2025-01-15 18:55:16.070254 : Constructing SummarizedExperiment, 1.337 mins elapsed.
2025-01-15 18:55:16.473266 : Finished Matrix Creation, 1.344 mins elapsed.

2025-01-15 18:55:16.474529 : getMatrixFromProject Input-Parameters, Class = list

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


2025-01-15 18:55:29.400108 : Reading PeakMatrix : SRR13927737(1 of 4), 0.215 mins elapsed.

2025-01-15 18:55:29.517923 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-15 18:55:32.556952 : featureDF SRR13927737, Class = DFrame
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

2025-01-15 18:55:32.586051 : Getting PeakMatrix from ArrowFile : SRR13927737.arrow, 0.051 mins elapsed.

2025-01-15 18:55:40.338528 : mat SRR13927737, Class = dgCMatrix
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


2025-01-15 18:55:40.538311 : Organizing SE PeakMatrix from ArrowFile : SRR13927737.arrow, 0.184 mins elapsed.

2025-01-15 18:55:41.456243 : se SRR13927737, Class = RangedSummarizedExperiment
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

2025-01-15 18:55:41.464232 : Completed PeakMatrix : SRR13927737(1 of 4), 0.416 mins elapsed.
2025-01-15 18:55:41.464752 : Reading PeakMatrix : SRR13927738(2 of 4), 0.417 mins elapsed.

2025-01-15 18:55:41.544453 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-15 18:55:44.098079 : featureDF SRR13927738, Class = DFrame
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

2025-01-15 18:55:44.131566 : Getting PeakMatrix from ArrowFile : SRR13927738.arrow, 0.043 mins elapsed.

2025-01-15 18:55:51.03116 : mat SRR13927738, Class = dgCMatrix
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


2025-01-15 18:55:51.24681 : Organizing SE PeakMatrix from ArrowFile : SRR13927738.arrow, 0.162 mins elapsed.

2025-01-15 18:55:51.522104 : se SRR13927738, Class = RangedSummarizedExperiment
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

2025-01-15 18:55:51.527005 : Completed PeakMatrix : SRR13927738(2 of 4), 0.584 mins elapsed.
2025-01-15 18:55:51.527496 : Reading PeakMatrix : SRR13927735(3 of 4), 0.584 mins elapsed.

2025-01-15 18:55:51.61679 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-15 18:55:54.586575 : featureDF SRR13927735, Class = DFrame
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

2025-01-15 18:55:54.616296 : Getting PeakMatrix from ArrowFile : SRR13927735.arrow, 0.05 mins elapsed.

2025-01-15 18:56:00.717666 : mat SRR13927735, Class = dgCMatrix
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


2025-01-15 18:56:00.890053 : Organizing SE PeakMatrix from ArrowFile : SRR13927735.arrow, 0.155 mins elapsed.

2025-01-15 18:56:01.159515 : se SRR13927735, Class = RangedSummarizedExperiment
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

2025-01-15 18:56:01.164272 : Completed PeakMatrix : SRR13927735(3 of 4), 0.745 mins elapsed.
2025-01-15 18:56:01.164764 : Reading PeakMatrix : SRR13927736(4 of 4), 0.745 mins elapsed.

2025-01-15 18:56:01.240911 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-15 18:56:04.910011 : featureDF SRR13927736, Class = DFrame
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

2025-01-15 18:56:04.93759 : Getting PeakMatrix from ArrowFile : SRR13927736.arrow, 0.062 mins elapsed.

2025-01-15 18:56:10.554295 : mat SRR13927736, Class = dgCMatrix
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


2025-01-15 18:56:10.726918 : Organizing SE PeakMatrix from ArrowFile : SRR13927736.arrow, 0.158 mins elapsed.

2025-01-15 18:56:11.053748 : se SRR13927736, Class = RangedSummarizedExperiment
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

2025-01-15 18:56:11.058354 : Completed PeakMatrix : SRR13927736(4 of 4), 0.91 mins elapsed.
2025-01-15 18:56:11.059645 : Organizing colData, 0.91 mins elapsed.
2025-01-15 18:56:11.146832 : Organizing rowData, 0.911 mins elapsed.
2025-01-15 18:56:11.152109 : Organizing rowRanges, 0.911 mins elapsed.
2025-01-15 18:56:11.159332 : Organizing Assays (1 of 1), 0.911 mins elapsed.
2025-01-15 18:56:12.569337 : Constructing SummarizedExperiment, 0.935 mins elapsed.
2025-01-15 18:56:18.166288 : Finished Matrix Creation, 1.028 mins elapsed.

2025-01-15 19:38:10.159373 : getMatrixFromProject Input-Parameters, Class = list

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


2025-01-15 19:38:15.784735 : Reading GeneIntegrationMatrix : SRR13927737(1 of 4), 0.094 mins elapsed.

2025-01-15 19:38:15.843208 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-15 19:38:19.129453 : featureDF SRR13927737, Class = DFrame
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

2025-01-15 19:38:19.16081 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927737.arrow, 0.055 mins elapsed.

2025-01-15 19:38:59.430947 : mat SRR13927737, Class = dgCMatrix
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


2025-01-15 19:38:59.890836 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927737.arrow, 0.734 mins elapsed.

2025-01-15 19:39:00.363205 : se SRR13927737, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 4439 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(4439): SRR13927737#GTCACGGGTCAACAGG-1 SRR13927737#TGCCTGTTCTTACTCA-1 ...
  SRR13927737#TTGCTATAGATGTTCC-1 SRR13927737#GCTCACTAGTTAGCGG-1
colData names(21): BlacklistRatio nDiFrags ... ReadsInPeaks FRIP

2025-01-15 19:39:00.368087 : Completed GeneIntegrationMatrix : SRR13927737(1 of 4), 0.837 mins elapsed.
2025-01-15 19:39:00.368683 : Reading GeneIntegrationMatrix : SRR13927738(2 of 4), 0.837 mins elapsed.

2025-01-15 19:39:00.430362 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-15 19:39:03.508071 : featureDF SRR13927738, Class = DFrame
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

2025-01-15 19:39:03.521142 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927738.arrow, 0.052 mins elapsed.

2025-01-15 19:39:36.194774 : mat SRR13927738, Class = dgCMatrix
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


2025-01-15 19:39:36.435571 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927738.arrow, 0.6 mins elapsed.

2025-01-15 19:39:36.740073 : se SRR13927738, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 3736 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(3736): SRR13927738#TAGCATGAGGCTAAAT-1 SRR13927738#TGTTAGGTCGGTCAGC-1 ...
  SRR13927738#GCTCCTATCCGTGCGA-1 SRR13927738#GGCACGTAGTCGATAA-1
colData names(21): BlacklistRatio nDiFrags ... ReadsInPeaks FRIP

2025-01-15 19:39:36.744202 : Completed GeneIntegrationMatrix : SRR13927738(2 of 4), 1.443 mins elapsed.
2025-01-15 19:39:36.744787 : Reading GeneIntegrationMatrix : SRR13927735(3 of 4), 1.443 mins elapsed.

2025-01-15 19:39:36.863299 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-15 19:39:39.674414 : featureDF SRR13927735, Class = DFrame
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

2025-01-15 19:39:39.686795 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927735.arrow, 0.047 mins elapsed.

2025-01-15 19:40:06.345423 : mat SRR13927735, Class = dgCMatrix
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


2025-01-15 19:40:06.734702 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927735.arrow, 0.498 mins elapsed.

2025-01-15 19:40:07.070079 : se SRR13927735, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 3813 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(3813): SRR13927735#TTATGTCTCCAGGTAT-1 SRR13927735#TATTGCTCATCAGAAA-1 ...
  SRR13927735#GATCATGGTGTTGTTG-1 SRR13927735#GTCCATCCATTTCACT-1
colData names(21): BlacklistRatio nDiFrags ... ReadsInPeaks FRIP

2025-01-15 19:40:07.073942 : Completed GeneIntegrationMatrix : SRR13927735(3 of 4), 1.949 mins elapsed.
2025-01-15 19:40:07.074484 : Reading GeneIntegrationMatrix : SRR13927736(4 of 4), 1.949 mins elapsed.

2025-01-15 19:40:07.209244 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-15 19:40:10.102265 : featureDF SRR13927736, Class = DFrame
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

2025-01-15 19:40:10.114392 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927736.arrow, 0.048 mins elapsed.

2025-01-15 19:40:19.294117 : mat SRR13927736, Class = dgCMatrix
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


2025-01-15 19:40:19.69999 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927736.arrow, 0.208 mins elapsed.

2025-01-15 19:40:20.014909 : se SRR13927736, Class = SummarizedExperiment
class: SummarizedExperiment 
dim: 23525 3534 
metadata(0):
assays(1): GeneIntegrationMatrix
rownames: NULL
rowData names(6): seqnames start ... name idx
colnames(3534): SRR13927736#GCTGTTCTCCTTGACC-1 SRR13927736#AGTTACGAGTAGGTCG-1 ...
  SRR13927736#TAAACCGCAAGTCTGT-1 SRR13927736#TGCATTTGTTCTACCC-1
colData names(21): BlacklistRatio nDiFrags ... ReadsInPeaks FRIP

2025-01-15 19:40:20.019362 : Completed GeneIntegrationMatrix : SRR13927736(4 of 4), 2.164 mins elapsed.
2025-01-15 19:40:20.020888 : Organizing colData, 2.164 mins elapsed.
2025-01-15 19:40:20.220388 : Organizing rowData, 2.168 mins elapsed.
2025-01-15 19:40:20.223362 : Organizing rowRanges, 2.168 mins elapsed.
2025-01-15 19:40:20.228166 : Organizing Assays (1 of 1), 2.168 mins elapsed.
2025-01-15 19:40:24.264148 : Constructing SummarizedExperiment, 2.235 mins elapsed.
2025-01-15 19:40:25.056994 : Finished Matrix Creation, 2.247 mins elapsed.

2025-01-15 19:40:25.069407 : getMatrixFromProject Input-Parameters, Class = list

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


2025-01-15 19:40:32.69433 : Reading PeakMatrix : SRR13927737(1 of 4), 0.127 mins elapsed.

2025-01-15 19:40:32.73519 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-15 19:40:35.111775 : featureDF SRR13927737, Class = DFrame
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

2025-01-15 19:40:35.140963 : Getting PeakMatrix from ArrowFile : SRR13927737.arrow, 0.04 mins elapsed.

2025-01-15 19:40:57.957421 : mat SRR13927737, Class = dgCMatrix
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


2025-01-15 19:40:58.156281 : Organizing SE PeakMatrix from ArrowFile : SRR13927737.arrow, 0.424 mins elapsed.

2025-01-15 19:41:01.820286 : se SRR13927737, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 4439 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(4439): SRR13927737#GTCACGGGTCAACAGG-1 SRR13927737#TGCCTGTTCTTACTCA-1 ...
  SRR13927737#TTGCTATAGATGTTCC-1 SRR13927737#GCTCACTAGTTAGCGG-1
colData names(21): BlacklistRatio nDiFrags ... ReadsInPeaks FRIP

2025-01-15 19:41:01.827174 : Completed PeakMatrix : SRR13927737(1 of 4), 0.613 mins elapsed.
2025-01-15 19:41:01.827762 : Reading PeakMatrix : SRR13927738(2 of 4), 0.613 mins elapsed.

2025-01-15 19:41:02.011222 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-15 19:41:05.211815 : featureDF SRR13927738, Class = DFrame
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

2025-01-15 19:41:05.243803 : Getting PeakMatrix from ArrowFile : SRR13927738.arrow, 0.054 mins elapsed.

2025-01-15 19:41:30.930406 : mat SRR13927738, Class = dgCMatrix
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


2025-01-15 19:41:31.242285 : Organizing SE PeakMatrix from ArrowFile : SRR13927738.arrow, 0.487 mins elapsed.

2025-01-15 19:41:31.559924 : se SRR13927738, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 3736 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(3736): SRR13927738#TAGCATGAGGCTAAAT-1 SRR13927738#TGTTAGGTCGGTCAGC-1 ...
  SRR13927738#GCTCCTATCCGTGCGA-1 SRR13927738#GGCACGTAGTCGATAA-1
colData names(21): BlacklistRatio nDiFrags ... ReadsInPeaks FRIP

2025-01-15 19:41:31.565202 : Completed PeakMatrix : SRR13927738(2 of 4), 1.108 mins elapsed.
2025-01-15 19:41:31.565813 : Reading PeakMatrix : SRR13927735(3 of 4), 1.108 mins elapsed.

2025-01-15 19:41:31.653501 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-15 19:41:40.837307 : featureDF SRR13927735, Class = DFrame
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

2025-01-15 19:41:40.86786 : Getting PeakMatrix from ArrowFile : SRR13927735.arrow, 0.154 mins elapsed.

2025-01-15 19:41:58.137003 : mat SRR13927735, Class = dgCMatrix
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


2025-01-15 19:41:58.315795 : Organizing SE PeakMatrix from ArrowFile : SRR13927735.arrow, 0.444 mins elapsed.

2025-01-15 19:41:58.836026 : se SRR13927735, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 3813 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(3813): SRR13927735#TTATGTCTCCAGGTAT-1 SRR13927735#TATTGCTCATCAGAAA-1 ...
  SRR13927735#GATCATGGTGTTGTTG-1 SRR13927735#GTCCATCCATTTCACT-1
colData names(21): BlacklistRatio nDiFrags ... ReadsInPeaks FRIP

2025-01-15 19:41:58.84118 : Completed PeakMatrix : SRR13927735(3 of 4), 1.563 mins elapsed.
2025-01-15 19:41:58.841782 : Reading PeakMatrix : SRR13927736(4 of 4), 1.563 mins elapsed.

2025-01-15 19:41:59.031921 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-15 19:42:11.751915 : featureDF SRR13927736, Class = DFrame
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

2025-01-15 19:42:11.781042 : Getting PeakMatrix from ArrowFile : SRR13927736.arrow, 0.212 mins elapsed.

2025-01-15 19:42:26.822701 : mat SRR13927736, Class = dgCMatrix
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


2025-01-15 19:42:26.993977 : Organizing SE PeakMatrix from ArrowFile : SRR13927736.arrow, 0.466 mins elapsed.

2025-01-15 19:42:27.270858 : se SRR13927736, Class = RangedSummarizedExperiment
class: RangedSummarizedExperiment 
dim: 210721 3534 
metadata(0):
assays(1): PeakMatrix
rownames: NULL
rowData names(1): idx
colnames(3534): SRR13927736#GCTGTTCTCCTTGACC-1 SRR13927736#AGTTACGAGTAGGTCG-1 ...
  SRR13927736#TAAACCGCAAGTCTGT-1 SRR13927736#TGCATTTGTTCTACCC-1
colData names(21): BlacklistRatio nDiFrags ... ReadsInPeaks FRIP

2025-01-15 19:42:27.27586 : Completed PeakMatrix : SRR13927736(4 of 4), 2.037 mins elapsed.
2025-01-15 19:42:27.277137 : Organizing colData, 2.037 mins elapsed.
2025-01-15 19:42:27.470744 : Organizing rowData, 2.04 mins elapsed.
2025-01-15 19:42:27.476422 : Organizing rowRanges, 2.04 mins elapsed.
2025-01-15 19:42:27.484324 : Organizing Assays (1 of 1), 2.04 mins elapsed.
2025-01-15 19:42:28.359997 : Constructing SummarizedExperiment, 2.055 mins elapsed.
2025-01-15 19:42:38.516148 : Finished Matrix Creation, 2.17 mins elapsed.

2025-01-15 21:36:33.71355 : getMatrixFromProject Input-Parameters, Class = list

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


2025-01-15 21:36:46.631726 : Reading GeneIntegrationMatrix : SRR13927737(1 of 4), 0.215 mins elapsed.

2025-01-15 21:36:46.700676 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-15 21:36:49.73617 : featureDF SRR13927737, Class = DFrame
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

2025-01-15 21:36:49.756434 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927737.arrow, 0.051 mins elapsed.

2025-01-15 21:37:08.229643 : mat SRR13927737, Class = dgCMatrix
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


2025-01-15 21:37:08.583917 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927737.arrow, 0.365 mins elapsed.

2025-01-15 21:37:08.975924 : se SRR13927737, Class = SummarizedExperiment
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

2025-01-15 21:37:08.998612 : Completed GeneIntegrationMatrix : SRR13927737(1 of 4), 0.588 mins elapsed.
2025-01-15 21:37:08.99913 : Reading GeneIntegrationMatrix : SRR13927738(2 of 4), 0.588 mins elapsed.

2025-01-15 21:37:09.066797 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-15 21:37:12.443145 : featureDF SRR13927738, Class = DFrame
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

2025-01-15 21:37:12.454365 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927738.arrow, 0.056 mins elapsed.

2025-01-15 21:37:22.686836 : mat SRR13927738, Class = dgCMatrix
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


2025-01-15 21:37:22.966164 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927738.arrow, 0.232 mins elapsed.

2025-01-15 21:37:23.274007 : se SRR13927738, Class = SummarizedExperiment
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

2025-01-15 21:37:23.2795 : Completed GeneIntegrationMatrix : SRR13927738(2 of 4), 0.826 mins elapsed.
2025-01-15 21:37:23.280019 : Reading GeneIntegrationMatrix : SRR13927735(3 of 4), 0.826 mins elapsed.

2025-01-15 21:37:23.406806 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-15 21:37:26.440152 : featureDF SRR13927735, Class = DFrame
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

2025-01-15 21:37:26.451963 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927735.arrow, 0.051 mins elapsed.

2025-01-15 21:37:37.910875 : mat SRR13927735, Class = dgCMatrix
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


2025-01-15 21:37:38.225269 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927735.arrow, 0.247 mins elapsed.

2025-01-15 21:37:38.620096 : se SRR13927735, Class = SummarizedExperiment
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

2025-01-15 21:37:38.624619 : Completed GeneIntegrationMatrix : SRR13927735(3 of 4), 1.082 mins elapsed.
2025-01-15 21:37:38.625111 : Reading GeneIntegrationMatrix : SRR13927736(4 of 4), 1.082 mins elapsed.

2025-01-15 21:37:38.686842 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-15 21:37:41.256316 : featureDF SRR13927736, Class = DFrame
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

2025-01-15 21:37:41.266172 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927736.arrow, 0.043 mins elapsed.

2025-01-15 21:37:50.625235 : mat SRR13927736, Class = dgCMatrix
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


2025-01-15 21:37:50.896865 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927736.arrow, 0.204 mins elapsed.

2025-01-15 21:37:51.058316 : se SRR13927736, Class = SummarizedExperiment
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

2025-01-15 21:37:51.062064 : Completed GeneIntegrationMatrix : SRR13927736(4 of 4), 1.289 mins elapsed.
2025-01-15 21:37:51.063288 : Organizing colData, 1.289 mins elapsed.
2025-01-15 21:37:51.157805 : Organizing rowData, 1.291 mins elapsed.
2025-01-15 21:37:51.160542 : Organizing rowRanges, 1.291 mins elapsed.
2025-01-15 21:37:51.164808 : Organizing Assays (1 of 1), 1.291 mins elapsed.
2025-01-15 21:37:53.188852 : Constructing SummarizedExperiment, 1.325 mins elapsed.
2025-01-15 21:37:53.590729 : Finished Matrix Creation, 1.331 mins elapsed.

2025-01-15 21:37:53.592141 : getMatrixFromProject Input-Parameters, Class = list

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


2025-01-15 21:38:05.252897 : Reading PeakMatrix : SRR13927737(1 of 4), 0.194 mins elapsed.

2025-01-15 21:38:05.346115 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-15 21:38:09.097579 : featureDF SRR13927737, Class = DFrame
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

2025-01-15 21:38:09.125518 : Getting PeakMatrix from ArrowFile : SRR13927737.arrow, 0.063 mins elapsed.

2025-01-15 21:38:16.217868 : mat SRR13927737, Class = dgCMatrix
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


2025-01-15 21:38:16.41933 : Organizing SE PeakMatrix from ArrowFile : SRR13927737.arrow, 0.185 mins elapsed.

2025-01-15 21:38:16.612773 : se SRR13927737, Class = RangedSummarizedExperiment
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

2025-01-15 21:38:16.619235 : Completed PeakMatrix : SRR13927737(1 of 4), 0.384 mins elapsed.
2025-01-15 21:38:16.619766 : Reading PeakMatrix : SRR13927738(2 of 4), 0.384 mins elapsed.

2025-01-15 21:38:16.736688 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-15 21:38:19.565295 : featureDF SRR13927738, Class = DFrame
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

2025-01-15 21:38:19.593984 : Getting PeakMatrix from ArrowFile : SRR13927738.arrow, 0.048 mins elapsed.

2025-01-15 21:38:26.048344 : mat SRR13927738, Class = dgCMatrix
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


2025-01-15 21:38:26.268525 : Organizing SE PeakMatrix from ArrowFile : SRR13927738.arrow, 0.159 mins elapsed.

2025-01-15 21:38:26.483618 : se SRR13927738, Class = RangedSummarizedExperiment
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

2025-01-15 21:38:26.488129 : Completed PeakMatrix : SRR13927738(2 of 4), 0.548 mins elapsed.
2025-01-15 21:38:26.488621 : Reading PeakMatrix : SRR13927735(3 of 4), 0.548 mins elapsed.

2025-01-15 21:38:26.574287 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-15 21:38:30.892493 : featureDF SRR13927735, Class = DFrame
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

2025-01-15 21:38:30.925351 : Getting PeakMatrix from ArrowFile : SRR13927735.arrow, 0.073 mins elapsed.

2025-01-15 21:38:36.979413 : mat SRR13927735, Class = dgCMatrix
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


2025-01-15 21:38:37.157728 : Organizing SE PeakMatrix from ArrowFile : SRR13927735.arrow, 0.176 mins elapsed.

2025-01-15 21:38:37.50311 : se SRR13927735, Class = RangedSummarizedExperiment
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

2025-01-15 21:38:37.507963 : Completed PeakMatrix : SRR13927735(3 of 4), 0.732 mins elapsed.
2025-01-15 21:38:37.508493 : Reading PeakMatrix : SRR13927736(4 of 4), 0.732 mins elapsed.

2025-01-15 21:38:37.591969 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-15 21:38:40.443665 : featureDF SRR13927736, Class = DFrame
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

2025-01-15 21:38:40.4722 : Getting PeakMatrix from ArrowFile : SRR13927736.arrow, 0.048 mins elapsed.

2025-01-15 21:38:46.263871 : mat SRR13927736, Class = dgCMatrix
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


2025-01-15 21:38:46.436159 : Organizing SE PeakMatrix from ArrowFile : SRR13927736.arrow, 0.147 mins elapsed.

2025-01-15 21:38:46.648208 : se SRR13927736, Class = RangedSummarizedExperiment
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

2025-01-15 21:38:46.652832 : Completed PeakMatrix : SRR13927736(4 of 4), 0.884 mins elapsed.
2025-01-15 21:38:46.654066 : Organizing colData, 0.884 mins elapsed.
2025-01-15 21:38:46.734579 : Organizing rowData, 0.886 mins elapsed.
2025-01-15 21:38:46.740593 : Organizing rowRanges, 0.886 mins elapsed.
2025-01-15 21:38:46.747854 : Organizing Assays (1 of 1), 0.886 mins elapsed.
2025-01-15 21:38:48.031828 : Constructing SummarizedExperiment, 0.907 mins elapsed.
2025-01-15 21:38:53.872608 : Finished Matrix Creation, 1.005 mins elapsed.

2025-01-15 22:33:51.699515 : getMatrixFromProject Input-Parameters, Class = list

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


2025-01-15 22:34:16.724296 : Reading GeneIntegrationMatrix : SRR13927737(1 of 4), 0.417 mins elapsed.

2025-01-15 22:34:16.79496 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-15 22:34:19.96159 : featureDF SRR13927737, Class = DFrame
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

2025-01-15 22:34:19.981536 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927737.arrow, 0.053 mins elapsed.

2025-01-15 22:34:38.188608 : mat SRR13927737, Class = dgCMatrix
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


2025-01-15 22:34:38.551891 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927737.arrow, 0.363 mins elapsed.

2025-01-15 22:34:38.909972 : se SRR13927737, Class = SummarizedExperiment
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

2025-01-15 22:34:38.914389 : Completed GeneIntegrationMatrix : SRR13927737(1 of 4), 0.787 mins elapsed.
2025-01-15 22:34:38.914853 : Reading GeneIntegrationMatrix : SRR13927738(2 of 4), 0.787 mins elapsed.

2025-01-15 22:34:38.974627 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-15 22:34:42.41198 : featureDF SRR13927738, Class = DFrame
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

2025-01-15 22:34:42.423738 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927738.arrow, 0.057 mins elapsed.

2025-01-15 22:34:52.15456 : mat SRR13927738, Class = dgCMatrix
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


2025-01-15 22:34:52.46032 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927738.arrow, 0.225 mins elapsed.

2025-01-15 22:34:52.66779 : se SRR13927738, Class = SummarizedExperiment
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

2025-01-15 22:34:52.671033 : Completed GeneIntegrationMatrix : SRR13927738(2 of 4), 1.016 mins elapsed.
2025-01-15 22:34:52.67147 : Reading GeneIntegrationMatrix : SRR13927735(3 of 4), 1.016 mins elapsed.

2025-01-15 22:34:52.956146 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-15 22:34:55.837122 : featureDF SRR13927735, Class = DFrame
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

2025-01-15 22:34:55.847504 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927735.arrow, 0.048 mins elapsed.

2025-01-15 22:35:06.284218 : mat SRR13927735, Class = dgCMatrix
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


2025-01-15 22:35:06.595172 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927735.arrow, 0.227 mins elapsed.

2025-01-15 22:35:06.903121 : se SRR13927735, Class = SummarizedExperiment
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

2025-01-15 22:35:06.906448 : Completed GeneIntegrationMatrix : SRR13927735(3 of 4), 1.253 mins elapsed.
2025-01-15 22:35:06.906904 : Reading GeneIntegrationMatrix : SRR13927736(4 of 4), 1.253 mins elapsed.

2025-01-15 22:35:06.973091 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-15 22:35:09.439218 : featureDF SRR13927736, Class = DFrame
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

2025-01-15 22:35:09.449307 : Getting GeneIntegrationMatrix from ArrowFile : SRR13927736.arrow, 0.041 mins elapsed.

2025-01-15 22:35:18.421908 : mat SRR13927736, Class = dgCMatrix
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


2025-01-15 22:35:18.689032 : Organizing SE GeneIntegrationMatrix from ArrowFile : SRR13927736.arrow, 0.195 mins elapsed.

2025-01-15 22:35:18.851962 : se SRR13927736, Class = SummarizedExperiment
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

2025-01-15 22:35:18.855426 : Completed GeneIntegrationMatrix : SRR13927736(4 of 4), 1.453 mins elapsed.
2025-01-15 22:35:18.856524 : Organizing colData, 1.453 mins elapsed.
2025-01-15 22:35:18.940121 : Organizing rowData, 1.454 mins elapsed.
2025-01-15 22:35:18.942316 : Organizing rowRanges, 1.454 mins elapsed.
2025-01-15 22:35:18.945889 : Organizing Assays (1 of 1), 1.454 mins elapsed.
2025-01-15 22:35:20.8752 : Constructing SummarizedExperiment, 1.486 mins elapsed.
2025-01-15 22:35:21.278247 : Finished Matrix Creation, 1.493 mins elapsed.

2025-01-15 22:35:21.279493 : getMatrixFromProject Input-Parameters, Class = list

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


2025-01-15 22:35:33.699196 : Reading PeakMatrix : SRR13927737(1 of 4), 0.207 mins elapsed.

2025-01-15 22:35:33.788708 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-15 22:35:36.72677 : featureDF SRR13927737, Class = DFrame
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

2025-01-15 22:35:36.754665 : Getting PeakMatrix from ArrowFile : SRR13927737.arrow, 0.049 mins elapsed.

2025-01-15 22:35:43.779405 : mat SRR13927737, Class = dgCMatrix
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


2025-01-15 22:35:43.97718 : Organizing SE PeakMatrix from ArrowFile : SRR13927737.arrow, 0.17 mins elapsed.

2025-01-15 22:35:44.173182 : se SRR13927737, Class = RangedSummarizedExperiment
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

2025-01-15 22:35:44.17941 : Completed PeakMatrix : SRR13927737(1 of 4), 0.382 mins elapsed.
2025-01-15 22:35:44.179872 : Reading PeakMatrix : SRR13927738(2 of 4), 0.382 mins elapsed.

2025-01-15 22:35:44.267517 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-15 22:35:47.087235 : featureDF SRR13927738, Class = DFrame
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

2025-01-15 22:35:47.11633 : Getting PeakMatrix from ArrowFile : SRR13927738.arrow, 0.047 mins elapsed.

2025-01-15 22:35:53.540618 : mat SRR13927738, Class = dgCMatrix
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


2025-01-15 22:35:53.749349 : Organizing SE PeakMatrix from ArrowFile : SRR13927738.arrow, 0.158 mins elapsed.

2025-01-15 22:35:54.142228 : se SRR13927738, Class = RangedSummarizedExperiment
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

2025-01-15 22:35:54.146556 : Completed PeakMatrix : SRR13927738(2 of 4), 0.548 mins elapsed.
2025-01-15 22:35:54.147002 : Reading PeakMatrix : SRR13927735(3 of 4), 0.548 mins elapsed.

2025-01-15 22:35:54.232421 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-15 22:35:57.241785 : featureDF SRR13927735, Class = DFrame
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

2025-01-15 22:35:57.270054 : Getting PeakMatrix from ArrowFile : SRR13927735.arrow, 0.051 mins elapsed.

2025-01-15 22:36:02.968998 : mat SRR13927735, Class = dgCMatrix
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


2025-01-15 22:36:03.140349 : Organizing SE PeakMatrix from ArrowFile : SRR13927735.arrow, 0.148 mins elapsed.

2025-01-15 22:36:03.366094 : se SRR13927735, Class = RangedSummarizedExperiment
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

2025-01-15 22:36:03.370315 : Completed PeakMatrix : SRR13927735(3 of 4), 0.702 mins elapsed.
2025-01-15 22:36:03.370772 : Reading PeakMatrix : SRR13927736(4 of 4), 0.702 mins elapsed.

2025-01-15 22:36:03.446211 : getMatrixFromArrow Input-Parameters, Class = list

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



2025-01-15 22:36:06.129455 : featureDF SRR13927736, Class = DFrame
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

2025-01-15 22:36:06.158121 : Getting PeakMatrix from ArrowFile : SRR13927736.arrow, 0.045 mins elapsed.

2025-01-15 22:36:11.437416 : mat SRR13927736, Class = dgCMatrix
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


2025-01-15 22:36:11.607237 : Organizing SE PeakMatrix from ArrowFile : SRR13927736.arrow, 0.136 mins elapsed.

2025-01-15 22:36:11.797392 : se SRR13927736, Class = RangedSummarizedExperiment
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

2025-01-15 22:36:11.801747 : Completed PeakMatrix : SRR13927736(4 of 4), 0.842 mins elapsed.
2025-01-15 22:36:11.802854 : Organizing colData, 0.842 mins elapsed.
2025-01-15 22:36:11.880943 : Organizing rowData, 0.843 mins elapsed.
2025-01-15 22:36:11.885433 : Organizing rowRanges, 0.843 mins elapsed.
2025-01-15 22:36:11.891576 : Organizing Assays (1 of 1), 0.844 mins elapsed.
2025-01-15 22:36:13.138072 : Constructing SummarizedExperiment, 0.864 mins elapsed.
2025-01-15 22:36:18.680176 : Finished Matrix Creation, 0.957 mins elapsed.
