## Aprena scRNAseq Analysis

## libraries 
library(Seurat)
library(harmony)
library(dplyr)
library(ggplot2)

sessionInfo()
# R version 4.1.0 (2021-05-18)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS Big Sur 10.16
# 
# Matrix products: default
# LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] ggplot2_3.3.5      dplyr_1.0.7        harmony_0.1.0      Rcpp_1.0.7        
# [5] SeuratObject_4.0.2 Seurat_4.0.5      
# 
# loaded via a namespace (and not attached):
#   [1] nlme_3.1-152          matrixStats_0.61.0    spatstat.sparse_2.0-0 RcppAnnoy_0.0.19     
# [5] RColorBrewer_1.1-2    httr_1.4.2            sctransform_0.3.2     tools_4.1.0          
# [9] utf8_1.2.2            R6_2.5.1              irlba_2.3.3           rpart_4.1-15         
# [13] KernSmooth_2.23-20    uwot_0.1.10           mgcv_1.8-35           lazyeval_0.2.2       
# [17] colorspace_2.0-2      withr_2.4.2           tidyselect_1.1.1      gridExtra_2.3        
# [21] compiler_4.1.0        plotly_4.10.0         scales_1.1.1          lmtest_0.9-38        
# [25] spatstat.data_2.1-0   ggridges_0.5.3        pbapply_1.5-0         goftest_1.2-3        
# [29] stringr_1.4.0         digest_0.6.28         spatstat.utils_2.2-0  pkgconfig_2.0.3      
# [33] htmltools_0.5.2       parallelly_1.28.1     fastmap_1.1.0         htmlwidgets_1.5.4    
# [37] rlang_0.4.11          shiny_1.7.1           generics_0.1.0        zoo_1.8-9            
# [41] jsonlite_1.7.2        ica_1.0-2             magrittr_2.0.1        patchwork_1.1.1      
# [45] Matrix_1.3-3          munsell_0.5.0         fansi_0.5.0           abind_1.4-5          
# [49] reticulate_1.22       lifecycle_1.0.1       stringi_1.7.5         MASS_7.3-54          
# [53] Rtsne_0.15            plyr_1.8.6            grid_4.1.0            parallel_4.1.0       
# [57] listenv_0.8.0         promises_1.2.0.1      ggrepel_0.9.1         crayon_1.4.1         
# [61] miniUI_0.1.1.1        deldir_1.0-5          lattice_0.20-44       cowplot_1.1.1        
# [65] splines_4.1.0         tensor_1.5            pillar_1.6.3          igraph_1.2.6         
# [69] spatstat.geom_2.3-0   future.apply_1.8.1    reshape2_1.4.4        codetools_0.2-18     
# [73] leiden_0.3.9          glue_1.4.2            data.table_1.14.2     png_0.1-7            
# [77] vctrs_0.3.8           httpuv_1.6.3          gtable_0.3.0          RANN_2.6.1           
# [81] purrr_0.3.4           spatstat.core_2.3-0   polyclip_1.10-0       tidyr_1.1.4          
# [85] scattermore_0.7       future_1.22.1         mime_0.12             xtable_1.8-4         
# [89] later_1.3.0           survival_3.2-11       viridisLite_0.4.0     tibble_3.1.5         
# [93] cluster_2.1.2         globals_0.14.0        fitdistrplus_1.1-6    ellipsis_0.3.2       
# [97] ROCR_1.0-11    

## define directories
datadir <- "/Volumes/blue/ferrallm/ferrallm/Moffitt-CICPT-3181-Sallman-Amy-10x"
outdir <- paste(datadir,"/analysis",sep="")

## IMPORT INDIVIDUAL SAMPLES
S1.data <- Read10X(data.dir=paste(datadir,"/S1_ACR_001_005_SCR_v2/outs/filtered_feature_bc_matrix/",sep=""))
S1 <- CreateSeuratObject(counts=S1.data)
Idents(object=S1) <- "S1"
S1$Tx <- 'control'
S1$Rsp <- 'responder'
S1$timept <- 0 

S2.data <- Read10X(data.dir=paste(datadir,"/S2_A-SD_001_040_SCR/outs/filtered_feature_bc_matrix/",sep=""))
S2 <- CreateSeuratObject(counts=S2.data)
Idents(object=S2) <- "S2"
S2$Tx <- 'control'
S2$Rsp <- 'nonresponder'
S2$timept <- 0

S3.data <- Read10X(data.dir=paste(datadir,"/S3_ASD_001_058_SCR/outs/filtered_feature_bc_matrix/",sep=""))
S3 <- CreateSeuratObject(counts=S3.data)
Idents(object=S3) <- "S3"
S3$Tx <- 'control'
S3$Rsp <- 'nonresponder'
S3$timept <- 0

S4.data <- Read10X(data.dir=paste(datadir,"/S4_CCR_001_008_SCR/outs/filtered_feature_bc_matrix/",sep=""))
S4 <- CreateSeuratObject(counts=S4.data)
Idents(object=S4) <- "S4"
S4$Tx <- 'combo'
S4$Rsp <- 'responder'
S4$timept <- 0

S5.data <- Read10X(data.dir=paste(datadir,"/S5_CCR_001_008_mCRwHI/outs/filtered_feature_bc_matrix/",sep=""))
S5 <- CreateSeuratObject(counts=S5.data)
Idents(object=S5) <- "S5"
S5$Tx <- 'combo'
S5$Rsp <- 'responder'
S5$timept <- 2

S6.data <- Read10X(data.dir=paste(datadir,"/S6_CCR_002_019_SCR/outs/filtered_feature_bc_matrix/",sep=""))
S6 <- CreateSeuratObject(counts=S6.data)
Idents(object=S6) <- "S6"
S6$Tx <- 'combo'
S6$Rsp <- 'responder'
S6$timept <- 0

S7.data <- Read10X(data.dir=paste(datadir,"/S7_CCR_002_019_CR/outs/filtered_feature_bc_matrix/",sep=""))
S7 <- CreateSeuratObject(counts=S7.data)
Idents(object=S7) <- "S7"
S7$Tx <- 'combo'
S7$Rsp <- 'responder'
S7$timept <- 1

S8.data <- Read10X(data.dir=paste(datadir,"/S8_CCR_333_136_SCR/outs/filtered_feature_bc_matrix/",sep=""))
S8 <- CreateSeuratObject(counts=S8.data)
Idents(object=S8) <- "S8"
S8$Tx <- 'combo'
S8$Rsp <- 'responder'
S8$timept <- 0

S9.data <- Read10X(data.dir=paste(datadir,"/S9_CCR_333_136_CR/outs/filtered_feature_bc_matrix/",sep=""))
S9 <- CreateSeuratObject(counts=S9.data)
Idents(object=S9) <- "S9"
S9$Tx <- 'combo'
S9$Rsp <- 'responder'
S9$timept <- 1

S10.data <- Read10X(data.dir=paste(datadir,"/S10_CSD_008_127_SCR/outs/filtered_feature_bc_matrix/",sep=""))
S10 <- CreateSeuratObject(counts=S10.data)
Idents(object=S10) <- "S10"
S10$Tx <- 'combo'
S10$Rsp <- 'nonresponder'
S10$timept <- 0 

S11.data <- Read10X(data.dir=paste(datadir,"/S11_CSD_008_137_SCR/outs/filtered_feature_bc_matrix/",sep=""))
S11 <- CreateSeuratObject(counts=S11.data)
Idents(object=S11) <- "S11"
S11$Tx <- 'combo'
S11$Rsp <- 'nonresponder'
S11$timept <- 0 

S12.data <- Read10X(data.dir=paste(datadir,"/S12_CSD_011_082_SCR/outs/filtered_feature_bc_matrix/",sep=""))
S12 <- CreateSeuratObject(counts=S12.data)
Idents(object=S12) <- "S12"
S12$Tx <- 'combo'
S12$Rsp <- 'nonresponder'
S12$timept <- 0 

### MERGING SAMPLES













