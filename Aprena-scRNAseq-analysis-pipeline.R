######===========================================######
#### Aprena scRNAseq Analysis
######===========================================######

##############################################
#### QUICK REFERENCE
##############################################
# Seurat Commands: https://satijalab.org/seurat/articles/essential_commands.html
# Harmony Integration: https://satijalab.org/signac/0.2/articles/integration.html


##############################################
#### LOAD LIBRARIES 
##############################################
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

##############################################
#### DEFINE DIRECTORIES
##############################################
datadir <- "/Volumes/blue/ferrallm/ferrallm/Moffitt-CICPT-3181-Sallman-Amy-10x"
outdir <- paste(datadir,"/analysis",sep="")

##############################################
#### IMPORT INDIVIDUAL SAMPLES
##############################################
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

##############################################
#### MERGING SAMPLES
##############################################
aprenaCohort <- merge(x=S1, y=list(S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12), add.cell.ids=c("S1","S2","S3","S4","S5","S6","S7","S8","S9","S10","S11","S12"))

saveRDS(aprenaCohort, paste(outdir,"/Aprena-scRNAseq-loaded+merged-noQC_2021-10-21.rds",sep=""))

##############################################
#### QUALITY CONTROL
##############################################
aprenaCohort[["percent.mt"]] <- PercentageFeatureSet(aprenaCohort, pattern = "^MT-")
VlnPlot(object=aprenaCohort, features="percent.mt")+geom_hline(yintercept = 25)

summary(aprenaCohort@meta.data$percent.mt)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.000   5.289   9.243  17.346  19.920  97.655

perMitCutoff <- 25

pdf(paste(outdir,'/QC_Violin_PerMito-wCutoff_seurat_pipeline_preHarmony_2021-10-21.pdf',sep=""))
  VlnPlot(object=aprenaCohort, features="percent.mt")+geom_hline(yintercept = perMitCutoff)
dev.off()

pdf(paste(outdir,'/QC_Violin_nFeat+nCount+PerMito_seurat_pipeline_preHarmony_2021-10-21.pdf',sep=""))
  VlnPlot(aprenaCohort, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

summary(aprenaCohort@meta.data$nFeature_RNA)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   26     596    1307    1716    2131   10305 

nFeatLowerN <- 200
nFeatUpperN <- mean(aprenaCohort@meta.data$nFeature_RNA, na.rm=TRUE) + 2*sd(aprenaCohort@meta.data$nFeature_RNA, na.rm=TRUE) 
# 4858.044

pdf(paste(outdir,'/QC_Violin_nFeat-wCutoff_seurat_pipeline_preHarmony_2021-10-21.pdf',sep=""))
  VlnPlot(object=aprenaCohort, features="nFeature_RNA")+geom_hline(yintercept = nFeatLowerN)+geom_hline(yintercept = nFeatUpperN)
dev.off()

dim(aprenaCohort)
# [1]  36601 148986

aprenaCohort <- subset(aprenaCohort, subset = nFeature_RNA > nFeatLowerN & nFeature_RNA < nFeatUpperN & percent.mt < perMitCutoff)

dim(aprenaCohort)
# [1]  36601 106751 (lost ~28% of cells with this filtering)

saveRDS(aprenaCohort, paste(outdir,"/Aprena-scRNAseq-loaded+merged-postQCsubset_2021-10-21.rds",sep=""))

##############################################
#### LOG NORMALIZE SAMPLES
##############################################
aprenaCohort <- NormalizeData(aprenaCohort, normalization.method = "LogNormalize", scale.factor = 10000)

##############################################
#### FIND VARIABLE FEATURES
##############################################
aprenaCohort <- FindVariableFeatures(aprenaCohort, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes and plot
top10 <- head(VariableFeatures(aprenaCohort), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(aprenaCohort)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
pdf(paste(outdir,'/HVG_Top10-labeled_seurat_pipeline_preHarmony_2021-10-21.pdf',sep=""))
  plot2
dev.off()

##############################################
#### SCALE DATA
##############################################
# (only use HVG), regressing out effects of nCountRNA and percent.mito
aprenaCohort <- ScaleData(aprenaCohort, features = VariableFeatures(aprenaCohort), vars.to.regress = c("nCount_RNA","percent.mt"))

##############################################
#### DIMENSION REDUCTION: PCA
##############################################
aprenaCohort <- RunPCA(aprenaCohort, features = VariableFeatures(object = aprenaCohort))

# PC_ 1 
# Positive:  IFI30, CD68, CLEC7A, TYROBP, PLAUR, VIM, S100A11, BRI3, FCN1, PSAP 
# CTSD, MAFB, FCER1G, COTL1, TIMP1, CXCL8, LGALS1, CEBPD, CTSS, SERPINA1 
# NAMPT, LYZ, CST3, SAT1, C5AR1, RAB31, FTL, ATP2B1-AS1, CSTA, ANXA5 
# Negative:  AHSP, SLC4A1, CA1, HBA1, HBD, HEMGN, ANK1, HBM, PRDX2, HMBS 
# SPTA1, GYPA, GYPB, ALAS2, RHAG, FECH, SNCA, CA2, SOX6, OSBP2 
# SLC2A1, UROD, SLC25A37, TLCD4, TFRC, GLRX5, NFIA, CENPF, ABCB10, SELENBP1 
# PC_ 2 
# Positive:  LTB, IL7R, IL32, GAS5, FKBP11, CCL5, TRBC1, BIRC3, ABLIM1, DUSP2 
# CD27, SELENOM, CTSW, NKG7, KLRB1, SEC11C, JUN, CST7, CD79A, CD8A 
# KLRD1, GZMA, MZB1, IL2RB, SSR4, ZBP1, GZMH, GNLY, SNHG7, CD79B 
# Negative:  CD36, CA2, MKI67, GYPA, SPTA1, CENPF, RHAG, TFRC, HEMGN, HMBS 
# ANK1, BLVRB, SLC4A1, NUSAP1, NFIA, TLCD4, SOX6, ABCB10, SLC25A37, UROD 
# GYPB, ASPM, AHSP, TFDP1, H1F0, CPOX, TOP2A, PRC1, AQP1, CA1 
# PC_ 3 
# Positive:  HSP90AA1, SNHG29, STMN1, TUBA1B, DUT, H2AFZ, ATP5IF1, HSPD1, RANBP1, NUCB2 
# IGFBP7, GSTP1, TIMM13, HSPE1, ATP5MC1, GAPDH, GIHCG, TMEM14C, SYNGR1, HIST1H4C 
# TUBB, PCLAF, PPP1R14B, CDK6, TXN, MYB, NME1, CCT6A, PRSS57, IMPDH2 
# Negative:  ALAS2, BPGM, TRIM58, IFIT1B, SLC4A1, TMCC2, SLC25A37, DCAF12, SNCA, CPEB4 
# SELENBP1, XPO7, HBA1, HBM, TRAK2, LGALS3, SLC14A1, TNS1, OSBP2, SLC2A1 
# ACSL6, HBD, AC100835.2, FECH, SOX6, TBCEL, TSPAN5, EPB42, SMOX, ARHGEF12 
# PC_ 4 
# Positive:  MNDA, MS4A6A, AC020656.1, LST1, AIF1, FGL2, LYZ, RNASE2, S100A4, TNFSF13B 
# LGALS2, S100A8, AP1S2, CYBB, ADA2, CTSS, SAMHD1, RNASE6, MARCH1, HLA-DRA 
# GSTP1, JAML, CLEC12A, KCTD12, GRN, HLA-DRB5, CSTA, S100A9, NCF2, FCN1 
# Negative:  CXCL12, CCDC80, FSTL1, DCN, EPAS1, COL6A1, FN1, PLPP3, VCAM1, CDH11 
# COL6A2, LUM, CP, IL1R1, CHL1, C1S, ANGPTL4, MGP, LEPR, FBN1 
# CALD1, FRMD6, CFH, C11orf96, COL6A3, COL1A1, IGFBP4, C1R, TF, COL1A2 
# PC_ 5 
# Positive:  THBS1, MET, SDC2, PHLDA1, FTH1, UPP1, SERPINB2, C15orf48, ITGB8, SEMA6B 
# CCL7, CXCL3, AQP9, CXCL2, CD109, CXCL16, PLIN2, NINJ1, G0S2, CXCL8 
# IFITM10, IRAK1, FNIP2, NEU4, FABP5, JARID2, ITGAX, CCL5, GAPDH, NR4A3 
# Negative:  MNDA, VCAN, MS4A6A, FGL2, MARCH1, TNFSF13B, LGALS2, AC020656.1, RNASE2, KCTD12 
# CD302, ADA2, CYBB, RNASE6, CLEC12A, HLA-DRA, LST1, AIF1, HLA-DRB5, JAML 
# HLA-DMB, AP1S2, TRIM58, DPYSL2, ALDH2, CTSS, S100A8, HLA-DRB1, TMEM176B, GRN 

pdf(paste(outdir,'/PCA_seurat_pipeline_preHarmony_2021-10-21.pdf',sep=""))
  DimPlot(aprenaCohort, reduction = "pca", pt.size = 0.0001)
dev.off()

# Determine dimensionality of dataset
pdf(paste(outdir,'/PCA_Elbow-Plot_seurat_pipeline_preHarmony_2021-10-21.pdf',sep=""))
  ElbowPlot(aprenaCohort, ndims = 50)
dev.off()

saveRDS(aprenaCohort, paste(outdir,"/Aprena-scRNAseq-loaded+merged-postPCA_2021-10-21.rds",sep=""))

############ TO DO
### follow back up with JackStraw analysis

##############################################
#### Harmony to remove batch effects
##############################################
aprenaCohort <- RunHarmony(aprenaCohort, group.by.vars="Tx")

# Harmony 1/10
# 0%   10   20   30   40   50   60   70   80   90   100%
#   [----|----|----|----|----|----|----|----|----|----|
#      **************************************************|
#      Harmony 2/10
#    0%   10   20   30   40   50   60   70   80   90   100%
#      [----|----|----|----|----|----|----|----|----|----|
#         **************************************************|
#         Harmony 3/10
#       0%   10   20   30   40   50   60   70   80   90   100%
#         [----|----|----|----|----|----|----|----|----|----|
#            **************************************************|
#            Harmony 4/10
#          0%   10   20   30   40   50   60   70   80   90   100%
#            [----|----|----|----|----|----|----|----|----|----|
#               **************************************************|
#               Harmony 5/10
#             0%   10   20   30   40   50   60   70   80   90   100%
#               [----|----|----|----|----|----|----|----|----|----|
#                  **************************************************|
#                  Harmony 6/10
#                0%   10   20   30   40   50   60   70   80   90   100%
#                  [----|----|----|----|----|----|----|----|----|----|
#                     **************************************************|
#                     Harmony 7/10
#                   0%   10   20   30   40   50   60   70   80   90   100%
#                     [----|----|----|----|----|----|----|----|----|----|
#                        **************************************************|
#                        Harmony 8/10
#                      0%   10   20   30   40   50   60   70   80   90   100%
#                        [----|----|----|----|----|----|----|----|----|----|
#                           **************************************************|
#                           Harmony converged after 8 iterations
#                         Warning: Invalid name supplied, making object name syntactically valid. New object name is Seurat..ProjectDim.RNA.harmony; see ?make.names for more details on syntax validity

saveRDS(aprenaCohort, paste(outdir,"/Aprena-scRNAseq-loaded+merged-postHarmony_2021-10-21.rds",sep=""))

##############################################
#### Find Neighbors, Cluster and Visualize with UMAP & Seurat Defaults
##############################################
aprenaCohort <- FindNeighbors(aprenaCohort, dims = 1:30)
aprenaCohort <- FindClusters(aprenaCohort, resolution = 0.8, verbose = FALSE)
aprenaCohort <- RunUMAP(aprenaCohort, dims = 1:30)

# Warning: The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
# To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
# This message will be shown once per session
# 17:01:32 UMAP embedding parameters a = 0.9922 b = 1.112
# 17:01:32 Read 106751 rows and found 30 numeric columns
# 17:01:32 Using Annoy for neighbor search, n_neighbors = 30
# 17:01:32 Building Annoy index with metric = cosine, n_trees = 50
# 0%   10   20   30   40   50   60   70   80   90   100%
#   [----|----|----|----|----|----|----|----|----|----|
#      **************************************************|
#      17:01:41 Writing NN index file to temp file /var/folders/hb/_pjzcr2n6ds4dpjtd4fxckp41d_t7x/T//RtmpUdXr8O/file22511a45a968
#    17:01:41 Searching Annoy index using 1 thread, search_k = 3000
#    17:02:09 Annoy recall = 100%
#    17:02:10 Commencing smooth kNN distance calibration using 1 thread
#    17:02:14 Initializing from normalized Laplacian + noise
#    17:02:25 Commencing optimization for 200 epochs, with 4951980 positive edges
#    0%   10   20   30   40   50   60   70   80   90   100%
#      [----|----|----|----|----|----|----|----|----|----|
#         **************************************************|
#         17:03:20 Optimization finished
#       Warning message:
#         In UseMethod("depth") :
#         no applicable method for 'depth' applied to an object of class "NULL"

saveRDS(aprenaCohort, paste(outdir,"/Aprena-scRNAseq-loaded+merged-postHarmony+Clustering+UMAP_2021-10-21.rds",sep=""))
 
##############################################
#### VISUALIZE DATA
##############################################

pdf(paste(outdir,'/UMAP_Clusters_seurat_pipeline_postHarmony_2021-10-21.pdf',sep=""))
  DimPlot(aprenaCohort, label = TRUE)
dev.off()

pdf(paste(outdir,'/UMAP_Group-by-Tx_seurat_pipeline_postHarmony_2021-10-21.pdf',sep=""))
  DimPlot(aprenaCohort, label = TRUE, group.by="Tx")
dev.off()

pdf(paste(outdir,'/UMAP_Group-by-ResponderStatus_seurat_pipeline_postHarmony_2021-10-21.pdf',sep=""))
  DimPlot(aprenaCohort, label = TRUE, group.by="Rsp")
dev.off()

pdf(paste(outdir,'/UMAP_Group-by-TimePoint_seurat_pipeline_postHarmony_2021-10-21.pdf',sep=""))
  DimPlot(aprenaCohort, label = TRUE, group.by="timept")
dev.off()


##############################################
#### IF HARMONY BY TIMEPOINT RATHER THAN TREATMENT
##############################################
aprenaCohortAlt <- readRDS(paste(outdir,"/Aprena-scRNAseq-loaded+merged-postPCA_2021-10-21.rds",sep=""))
aprenaCohortAlt <- RunHarmony(aprenaCohortAlt, group.by.vars="timept")

# Harmony 1/10
# 0%   10   20   30   40   50   60   70   80   90   100%
#   [----|----|----|----|----|----|----|----|----|----|
#      **************************************************|
#      Harmony 2/10
#    0%   10   20   30   40   50   60   70   80   90   100%
#      [----|----|----|----|----|----|----|----|----|----|
#         **************************************************|
#         Harmony 3/10
#       0%   10   20   30   40   50   60   70   80   90   100%
#         [----|----|----|----|----|----|----|----|----|----|
#            **************************************************|
#            Harmony converged after 3 iterations
#          Warning: Invalid name supplied, making object name syntactically valid. New object name is Seurat..ProjectDim.RNA.harmony; see ?make.names for more details on syntax validity

saveRDS(aprenaCohortAlt, paste(outdir,"/Aprena-scRNAseq-loaded+merged-postHarmony_ALT-byTimept_2021-10-21.rds",sep=""))

aprenaCohortAlt <- FindNeighbors(aprenaCohortAlt, dims = 1:30)
aprenaCohortAlt <- FindClusters(aprenaCohortAlt, resolution = 0.8, verbose = FALSE)
aprenaCohortAlt <- RunUMAP(aprenaCohortAlt, dims = 1:30)

# 17:24:03 UMAP embedding parameters a = 0.9922 b = 1.112
# 17:24:03 Read 106751 rows and found 30 numeric columns
# 17:24:03 Using Annoy for neighbor search, n_neighbors = 30
# 17:24:03 Building Annoy index with metric = cosine, n_trees = 50
# 0%   10   20   30   40   50   60   70   80   90   100%
#   [----|----|----|----|----|----|----|----|----|----|
#      **************************************************|
#      17:24:12 Writing NN index file to temp file /var/folders/hb/_pjzcr2n6ds4dpjtd4fxckp41d_t7x/T//RtmpUdXr8O/file2251ff08db0
#    17:24:12 Searching Annoy index using 1 thread, search_k = 3000
#    17:24:41 Annoy recall = 100%
#    17:24:42 Commencing smooth kNN distance calibration using 1 thread
#    17:24:46 Initializing from normalized Laplacian + noise
#    17:24:58 Commencing optimization for 200 epochs, with 4951980 positive edges
#    0%   10   20   30   40   50   60   70   80   90   100%
#      [----|----|----|----|----|----|----|----|----|----|
#         **************************************************|
#         17:25:54 Optimization finished
#       Warning message:
#         In grid.Call.graphics(C_upviewport, as.integer(n)) :
#         cannot pop the top-level viewport ('grid' and 'graphics' output mixed?)

saveRDS(aprenaCohortAlt, paste(outdir,"/Aprena-scRNAseq-loaded+merged-postHarmony+Clustering+UMAP_ALT-byTimept_2021-10-21.rds",sep=""))

pdf(paste(outdir,'/UMAP_Clusters_seurat_pipeline_postHarmony_ALT-byTimept_2021-10-21.pdf',sep=""))
  DimPlot(aprenaCohortAlt, label = TRUE)
dev.off()

pdf(paste(outdir,'/UMAP_Group-by-Tx_seurat_pipeline_postHarmony_ALT-byTimept_2021-10-21.pdf',sep=""))
  DimPlot(aprenaCohortAlt, label = TRUE, group.by="Tx")
dev.off()

pdf(paste(outdir,'/UMAP_Group-by-ResponderStatus_seurat_pipeline_postHarmony_ALT-byTimept_2021-10-21.pdf',sep=""))
  DimPlot(aprenaCohortAlt, label = TRUE, group.by="Rsp")
dev.off()

pdf(paste(outdir,'/UMAP_Group-by-TimePoint_seurat_pipeline_postHarmony_ALT-byTimept_2021-10-21.pdf',sep=""))
  DimPlot(aprenaCohortAlt, label = TRUE, group.by="timept")
dev.off()


##############################################
#### EXPORT DATA FOR DIVERSITY ANALYSIS
##############################################

datadir <- "/Volumes/blue/ferrallm/ferrallm/Moffitt-CICPT-3181-Sallman-Amy-10x"
outdir <- paste(datadir,"/analysis",sep="")
aprenaCohortAlt <- readRDS(paste(outdir,"/Aprena-scRNAseq-loaded+merged-postHarmony+Clustering+UMAP_ALT-byTimept_2021-10-21.rds",sep=""))

library(data.table)

date <- "_2021-11-04"
res <- "0.8"

file1 <- paste(outdir,"tmp1_names",res,date,".csv", sep="")
file2 <- paste(outdir,"tmp2_data",res,date,".csv", sep="")
divout <- paste(outdir,"/CellBreakdown_PerClusterPerType_res",res,date,".csv", sep="")

cat(aprenaCohortAlt@meta.data$RNA_snn_res.0.8, file=file2, sep=",\n")
aprenaCohortAlt$res0p8 <- aprenaCohortAlt@meta.data$RNA_snn_res.0.8

listNames <- aprenaCohortAlt@assays$RNA@data@Dimnames[2]
write.table(data.frame(listNames),
            row.names=FALSE,
            col.names = FALSE, 
            file = file1,
            sep=",")

mydat1 <- read.csv(file2)
mydat2 <- read.csv(file1)
fulldat <- cbind(mydat2[1],mydat1[1])
fulltab <- as.data.table(fulldat)
# name table columns
names(fulltab)[1] <- paste("UMI")
names(fulltab)[2] <- paste("cluster")
# group data based on clusters
group_by(fulltab, cluster)
# create a table counting unqiue UMIs/cells per cluster
tabPerClus <- fulltab %>% group_by(cluster) %>% count()
type <- sub("\\_.*","",fulltab$UMI)
fulltab <- cbind(fulltab, type)

A <- fulltab %>% group_by(cluster) %>% count(type)
write.csv(A, file=divout)


##############################################
#### DIFFERNETIAL GENE EXPRESSION - BY CONDITIONS OF INTEREST
##############################################

## rename samples with Patient IDs
aprenaCohortAlt@meta.data$PatientID <- unlist(tstrsplit(rownames(aprenaCohortAlt@meta.data),"_")[1])
unique(aprenaCohortAlt@meta.data$PatientID)

## set PatientID as active.ident for grouping
aprenaCohortAlt <- SetIdent(aprenaCohortAlt, value = aprenaCohortAlt@meta.data$PatientID)

## group cells
ctrls <- WhichCells(aprenaCohortAlt, idents = c("S1", "S2", "S3"))
combo <- WhichCells(aprenaCohortAlt, idents = c("S4", "S5", "S6","S7","S8","S9","S10","S11","S12"))

pdf(paste(outdir,"/UMAP_res0.8_Highlighted_Combo-Tx",date,".pdf",sep=""))
dPlot <- DimPlot(aprenaCohortAlt, reduction="umap", cells.highlight=combo) + theme(legend.position = "none")
print(dPlot)
dev.off()

ctrlvcombo.de.markers <- FindMarkers(aprenaCohortAlt, ident.1 = ctrls, ident.2 = combo)
#   |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=16m 18s
write.csv(ctrlvcombo.de.markers, paste(outdir,"/DiffExpGeneList_res0.8_1-Ctrl-vs-2-Combo_viaWilcoxonRankTest",date,".csv",sep=""))

rspd <- WhichCells(aprenaCohortAlt, idents = c("S1", "S4", "S5", "S6","S7","S8"))
nonrspd <- WhichCells(aprenaCohortAlt, idents = c("S2", "S3", "S10", "S11","S12"))

pdf(paste(outdir,"/UMAP_res0.8_Highlighted_Responders",date,".pdf",sep=""))
dPlot <- DimPlot(aprenaCohortAlt, reduction="umap", cells.highlight=rspd) + theme(legend.position = "none")
print(dPlot)
dev.off()

RespVNonResp.de.markers <- FindMarkers(aprenaCohortAlt, ident.1 = nonrspd, ident.2 = rspd)
#   |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=05m 47s
write.csv(RespVNonResp.de.markers, paste(outdir,"/DiffExpGeneList_res0.8_1-Nonresp-vs-2-Responder_viaWilcoxonRankTest",date,".csv",sep=""))

time0 <- WhichCells(aprenaCohortAlt, idents = c("S1", "S2", "S3", "S4", "S6","S8", "S10", "S11", "S12"))
time1p <- WhichCells(aprenaCohortAlt, idents = c("S5", "S7", "S9"))

pdf(paste(outdir,"/UMAP_res0.8_Highlighted_LaterTimePoint",date,".pdf",sep=""))
dPlot <- DimPlot(aprenaCohortAlt, reduction="umap", cells.highlight=time1p) + theme(legend.position = "none")
print(dPlot)
dev.off()

overtime.de.markers <- FindMarkers(aprenaCohortAlt, ident.1 = time0, ident.2 = time1p)
# |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=08m 10s
write.csv(overtime.de.markers, paste(outdir,"/DiffExpGeneList_res0.8_1-time-0-vs-2-later-time_viaWilcoxonRankTest",date,".csv",sep=""))

## additional combinations to compare - after completion of the cluster differential gene expression 
## on combo, responder vs non-responder
## on control, responder vs non-responder
## early versus late time point for a given individual

aprenaCohortAlt <- SetIdent(aprenaCohortAlt, value = aprenaCohortAlt@meta.data$PatientID)

## control, resp vs non-responder
ctrlRspd <- WhichCells(aprenaCohortAlt, idents = c("S1"))
ctrlNonrspd <- WhichCells(aprenaCohortAlt, idents = c("S2", "S3"))
ctrl.RvNR.de.markers <- FindMarkers(aprenaCohortAlt, ident.1 = ctrlNonrspd, ident.2 = ctrlRspd)
# For a more efficient implementation of the Wilcoxon Rank Sum Test,
# (default method for FindMarkers) please install the limma package
# --------------------------------------------
#   install.packages('BiocManager')
# BiocManager::install('limma')
# --------------------------------------------
#   After installation of limma, Seurat will automatically use the more 
# efficient implementation (no further action necessary).
# This message will be shown once per session
# |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=02m 44s
write.csv(ctrl.RvNR.de.markers, paste(outdir,"/DiffExpGeneList_res0.8_Ctrl-1-Nonresp-vs-2-Resp_viaWilcoxonRankTest",date,".csv",sep=""))

## combo, resp vs non-responder
comboRspd <- WhichCells(aprenaCohortAlt, idents = c("S4","S5","S6","S7","S8","S9"))
comboNonrspd <- WhichCells(aprenaCohortAlt, idents = c("S10", "S11", "S12"))
combo.RvNR.de.markers <- FindMarkers(aprenaCohortAlt, ident.1 = comboNonrspd, ident.2 = comboRspd)
# |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=08m 25s
write.csv(combo.RvNR.de.markers, paste(outdir,"/DiffExpGeneList_res0.8_Combo-1-Nonresp-vs-2-Resp_viaWilcoxonRankTest",date,".csv",sep=""))

## early versus late time point for a given individual
CCR002019early <- WhichCells(aprenaCohortAlt, idents = c("S6"))
CCR002019late <- WhichCells(aprenaCohortAlt, idents = c("S7"))
CCR002019.de.markers <- FindMarkers(aprenaCohortAlt, ident.1 = CCR002019early, ident.2 = CCR002019late)
# |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=01m 31s
write.csv(CCR002019.de.markers, paste(outdir,"/DiffExpGeneList_res0.8_C-CR_002-019-1-time0-vs-2-later-time_viaWilcoxonRankTest",date,".csv",sep=""))

CCR333136early <- WhichCells(aprenaCohortAlt, idents = c("S8"))
CCR333136late <- WhichCells(aprenaCohortAlt, idents = c("S9"))
CCR333136.de.markers <- FindMarkers(aprenaCohortAlt, ident.1 = CCR333136early, ident.2 = CCR333136late)
# |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=02m 02s
write.csv(CCR333136.de.markers, paste(outdir,"/DiffExpGeneList_res0.8_C-CR_333-136-1-time0-vs-2-later-time_viaWilcoxonRankTest",date,".csv",sep=""))

CCR001008early <- WhichCells(aprenaCohortAlt, idents = c("S4"))
CCR001008late <- WhichCells(aprenaCohortAlt, idents = c("S5"))
CCR001008.de.markers <- FindMarkers(aprenaCohortAlt, ident.1 = CCR001008early, ident.2 = CCR001008late)
# |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=07m 22s
write.csv(CCR001008.de.markers, paste(outdir,"/DiffExpGeneList_res0.8_C-CR_001-008-1-time0-vs-2-later-time_viaWilcoxonRankTest",date,".csv",sep=""))


##############################################
#### DIFFERNETIAL GENE EXPRESSION - BY CLUSTER
##############################################

## setting cluster assignment as active identity
aprenaCohortAlt <- SetIdent(aprenaCohortAlt, value = aprenaCohortAlt$res0p8)

levels(aprenaCohortAlt)
# clusters 0 through 37

clus0 <- WhichCells(aprenaCohortAlt, idents = c("0"))
clus1 <- WhichCells(aprenaCohortAlt, idents = c("1"))
clus2 <- WhichCells(aprenaCohortAlt, idents = c("2"))
clus3 <- WhichCells(aprenaCohortAlt, idents = c("3"))
clus4 <- WhichCells(aprenaCohortAlt, idents = c("4"))
clus5 <- WhichCells(aprenaCohortAlt, idents = c("5"))
clus6 <- WhichCells(aprenaCohortAlt, idents = c("6"))
clus7 <- WhichCells(aprenaCohortAlt, idents = c("7"))
clus8 <- WhichCells(aprenaCohortAlt, idents = c("8"))
clus9 <- WhichCells(aprenaCohortAlt, idents = c("9"))
clus10 <- WhichCells(aprenaCohortAlt, idents = c("10"))
clus11 <- WhichCells(aprenaCohortAlt, idents = c("11"))
clus12 <- WhichCells(aprenaCohortAlt, idents = c("12"))
clus13 <- WhichCells(aprenaCohortAlt, idents = c("13"))
clus14 <- WhichCells(aprenaCohortAlt, idents = c("14"))
clus15 <- WhichCells(aprenaCohortAlt, idents = c("15"))
clus16 <- WhichCells(aprenaCohortAlt, idents = c("16"))
clus17 <- WhichCells(aprenaCohortAlt, idents = c("17"))
clus18 <- WhichCells(aprenaCohortAlt, idents = c("18"))
clus19 <- WhichCells(aprenaCohortAlt, idents = c("19"))
clus20 <- WhichCells(aprenaCohortAlt, idents = c("20"))
clus21 <- WhichCells(aprenaCohortAlt, idents = c("21"))
clus22 <- WhichCells(aprenaCohortAlt, idents = c("22"))
clus23 <- WhichCells(aprenaCohortAlt, idents = c("23"))
clus24 <- WhichCells(aprenaCohortAlt, idents = c("24"))
clus25 <- WhichCells(aprenaCohortAlt, idents = c("25"))
clus26 <- WhichCells(aprenaCohortAlt, idents = c("26"))
clus27 <- WhichCells(aprenaCohortAlt, idents = c("27"))
clus28 <- WhichCells(aprenaCohortAlt, idents = c("28"))
clus29 <- WhichCells(aprenaCohortAlt, idents = c("29"))
clus30 <- WhichCells(aprenaCohortAlt, idents = c("30"))
clus31 <- WhichCells(aprenaCohortAlt, idents = c("31"))
clus32 <- WhichCells(aprenaCohortAlt, idents = c("32"))
clus33 <- WhichCells(aprenaCohortAlt, idents = c("33"))
clus34 <- WhichCells(aprenaCohortAlt, idents = c("34"))
clus35 <- WhichCells(aprenaCohortAlt, idents = c("35"))
clus36 <- WhichCells(aprenaCohortAlt, idents = c("36"))
clus37 <- WhichCells(aprenaCohortAlt, idents = c("37"))

clus0.de.markers <- FindMarkers(aprenaCohortAlt, ident.1 = clus0)
#   |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=24m 38s
write.csv(clus0.de.markers, paste(outdir,"/DiffExpGeneList_res0.8_1-Clus0-vs-2-All-Others_viaWilcoxonRankTest",date,".csv",sep=""))

clus1.de.markers <- FindMarkers(aprenaCohortAlt, ident.1 = clus1)
#   |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=21m 19s
write.csv(clus1.de.markers, paste(outdir,"/DiffExpGeneList_res0.8_1-Clus1-vs-2-All-Others_viaWilcoxonRankTest",date,".csv",sep=""))

clus2.de.markers <- FindMarkers(aprenaCohortAlt, ident.1 = clus2)
# |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=01h 29m 22s
write.csv(clus2.de.markers, paste(outdir,"/DiffExpGeneList_res0.8_1-Clus2-vs-2-All-Others_viaWilcoxonRankTest",date,".csv",sep=""))

clus3.de.markers <- FindMarkers(aprenaCohortAlt, ident.1 = clus3)
# |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=01h 04m 05s
write.csv(clus3.de.markers, paste(outdir,"/DiffExpGeneList_res0.8_1-Clus3-vs-2-All-Others_viaWilcoxonRankTest",date,".csv",sep=""))

clus4.de.markers <- FindMarkers(aprenaCohortAlt, ident.1 = clus4)
# |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=29m 33s
write.csv(clus4.de.markers, paste(outdir,"/DiffExpGeneList_res0.8_1-Clus4-vs-2-All-Others_viaWilcoxonRankTest",date,".csv",sep=""))

clus5.de.markers <- FindMarkers(aprenaCohortAlt, ident.1 = clus5)
# |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=29m 45s
write.csv(clus5.de.markers, paste(outdir,"/DiffExpGeneList_res0.8_1-Clus5-vs-2-All-Others_viaWilcoxonRankTest",date,".csv",sep=""))

clus6.de.markers <- FindMarkers(aprenaCohortAlt, ident.1 = clus6)
# |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=41m 12s
write.csv(clus6.de.markers, paste(outdir,"/DiffExpGeneList_res0.8_1-Clus6-vs-2-All-Others_viaWilcoxonRankTest",date,".csv",sep=""))

clus7.de.markers <- FindMarkers(aprenaCohortAlt, ident.1 = clus7)
# |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=53m 07s
write.csv(clus7.de.markers, paste(outdir,"/DiffExpGeneList_res0.8_1-Clus7-vs-2-All-Others_viaWilcoxonRankTest",date,".csv",sep=""))

clus8.de.markers <- FindMarkers(aprenaCohortAlt, ident.1 = clus8)
# |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=44m 54s
write.csv(clus8.de.markers, paste(outdir,"/DiffExpGeneList_res0.8_1-Clus8-vs-2-All-Others_viaWilcoxonRankTest",date,".csv",sep=""))

clus9.de.markers <- FindMarkers(aprenaCohortAlt, ident.1 = clus9)
# |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=40m 34s
write.csv(clus9.de.markers, paste(outdir,"/DiffExpGeneList_res0.8_1-Clus9-vs-2-All-Others_viaWilcoxonRankTest",date,".csv",sep=""))

clus10.de.markers <- FindMarkers(aprenaCohortAlt, ident.1 = clus10)
# |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=39m 07s
write.csv(clus10.de.markers, paste(outdir,"/DiffExpGeneList_res0.8_1-Clus10-vs-2-All-Others_viaWilcoxonRankTest",date,".csv",sep=""))

clus11.de.markers <- FindMarkers(aprenaCohortAlt, ident.1 = clus11)
# |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=23m 56s
write.csv(clus11.de.markers, paste(outdir,"/DiffExpGeneList_res0.8_1-Clus11-vs-2-All-Others_viaWilcoxonRankTest",date,".csv",sep=""))

clus12.de.markers <- FindMarkers(aprenaCohortAlt, ident.1 = clus12)
# |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=24m 05s
write.csv(clus12.de.markers, paste(outdir,"/DiffExpGeneList_res0.8_1-Clus12-vs-2-All-Others_viaWilcoxonRankTest",date,".csv",sep=""))

clus13.de.markers <- FindMarkers(aprenaCohortAlt, ident.1 = clus13)
#|++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=27m 35s
write.csv(clus13.de.markers, paste(outdir,"/DiffExpGeneList_res0.8_1-Clus13-vs-2-All-Others_viaWilcoxonRankTest",date,".csv",sep=""))

clus14.de.markers <- FindMarkers(aprenaCohortAlt, ident.1 = clus14)
# |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=21m 17s
write.csv(clus14.de.markers, paste(outdir,"/DiffExpGeneList_res0.8_1-Clus14-vs-2-All-Others_viaWilcoxonRankTest",date,".csv",sep=""))

clus15.de.markers <- FindMarkers(aprenaCohortAlt, ident.1 = clus15)
# |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=01h 20m 50s
write.csv(clus15.de.markers, paste(outdir,"/DiffExpGeneList_res0.8_1-Clus15-vs-2-All-Others_viaWilcoxonRankTest",date,".csv",sep=""))

clus16.de.markers <- FindMarkers(aprenaCohortAlt, ident.1 = clus16)
# |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=34m 53s
write.csv(clus16.de.markers, paste(outdir,"/DiffExpGeneList_res0.8_1-Clus16-vs-2-All-Others_viaWilcoxonRankTest",date,".csv",sep=""))

clus17.de.markers <- FindMarkers(aprenaCohortAlt, ident.1 = clus17)
# |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=32m 53s
write.csv(clus17.de.markers, paste(outdir,"/DiffExpGeneList_res0.8_1-Clus17-vs-2-All-Others_viaWilcoxonRankTest",date,".csv",sep=""))

clus18.de.markers <- FindMarkers(aprenaCohortAlt, ident.1 = clus18)
# |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=26m 57s
write.csv(clus18.de.markers, paste(outdir,"/DiffExpGeneList_res0.8_1-Clus18-vs-2-All-Others_viaWilcoxonRankTest",date,".csv",sep=""))

clus19.de.markers <- FindMarkers(aprenaCohortAlt, ident.1 = clus19)
# |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=33m 26s
write.csv(clus19.de.markers, paste(outdir,"/DiffExpGeneList_res0.8_1-Clus19-vs-2-All-Others_viaWilcoxonRankTest",date,".csv",sep=""))

clus20.de.markers <- FindMarkers(aprenaCohortAlt, ident.1 = clus20)
# |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=44m 36s
write.csv(clus20.de.markers, paste(outdir,"/DiffExpGeneList_res0.8_1-Clus20-vs-2-All-Others_viaWilcoxonRankTest",date,".csv",sep=""))

clus21.de.markers <- FindMarkers(aprenaCohortAlt, ident.1 = clus21)
# |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=22m 59s
write.csv(clus21.de.markers, paste(outdir,"/DiffExpGeneList_res0.8_1-Clus21-vs-2-All-Others_viaWilcoxonRankTest",date,".csv",sep=""))

clus22.de.markers <- FindMarkers(aprenaCohortAlt, ident.1 = clus22)
# |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=29m 17s
write.csv(clus22.de.markers, paste(outdir,"/DiffExpGeneList_res0.8_1-Clus22-vs-2-All-Others_viaWilcoxonRankTest",date,".csv",sep=""))

clus23.de.markers <- FindMarkers(aprenaCohortAlt, ident.1 = clus23)
# |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=30m 12s
write.csv(clus23.de.markers, paste(outdir,"/DiffExpGeneList_res0.8_1-Clus23-vs-2-All-Others_viaWilcoxonRankTest",date,".csv",sep=""))

clus24.de.markers <- FindMarkers(aprenaCohortAlt, ident.1 = clus24)
# |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=26m 04s
write.csv(clus24.de.markers, paste(outdir,"/DiffExpGeneList_res0.8_1-Clus24-vs-2-All-Others_viaWilcoxonRankTest",date,".csv",sep=""))

clus25.de.markers <- FindMarkers(aprenaCohortAlt, ident.1 = clus25)
# |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=26m 06s
write.csv(clus25.de.markers, paste(outdir,"/DiffExpGeneList_res0.8_1-Clus25-vs-2-All-Others_viaWilcoxonRankTest",date,".csv",sep=""))

date <- "_2021-11-04"

clus26.de.markers <- FindMarkers(aprenaCohortAlt, ident.1 = clus26)
# |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=46m 41s
write.csv(clus26.de.markers, paste(outdir,"/DiffExpGeneList_res0.8_1-Clus26-vs-2-All-Others_viaWilcoxonRankTest",date,".csv",sep=""))

clus27.de.markers <- FindMarkers(aprenaCohortAlt, ident.1 = clus27)
# |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=27m 12s
write.csv(clus27.de.markers, paste(outdir,"/DiffExpGeneList_res0.8_1-Clus27-vs-2-All-Others_viaWilcoxonRankTest",date,".csv",sep=""))

clus28.de.markers <- FindMarkers(aprenaCohortAlt, ident.1 = clus28)
# |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=16m 32s
write.csv(clus28.de.markers, paste(outdir,"/DiffExpGeneList_res0.8_1-Clus28-vs-2-All-Others_viaWilcoxonRankTest",date,".csv",sep=""))

clus29.de.markers <- FindMarkers(aprenaCohortAlt, ident.1 = clus29)
# |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=34m 11s
write.csv(clus29.de.markers, paste(outdir,"/DiffExpGeneList_res0.8_1-Clus29-vs-2-All-Others_viaWilcoxonRankTest",date,".csv",sep=""))

clus30.de.markers <- FindMarkers(aprenaCohortAlt, ident.1 = clus30)
# |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=32m 39s
write.csv(clus30.de.markers, paste(outdir,"/DiffExpGeneList_res0.8_1-Clus30-vs-2-All-Others_viaWilcoxonRankTest",date,".csv",sep=""))

clus31.de.markers <- FindMarkers(aprenaCohortAlt, ident.1 = clus31)
# |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=29m 10s
write.csv(clus31.de.markers, paste(outdir,"/DiffExpGeneList_res0.8_1-Clus31-vs-2-All-Others_viaWilcoxonRankTest",date,".csv",sep=""))

clus32.de.markers <- FindMarkers(aprenaCohortAlt, ident.1 = clus32)
# |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=33m 24s
write.csv(clus32.de.markers, paste(outdir,"/DiffExpGeneList_res0.8_1-Clus32-vs-2-All-Others_viaWilcoxonRankTest",date,".csv",sep=""))

clus33.de.markers <- FindMarkers(aprenaCohortAlt, ident.1 = clus33)
# |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=40m 45s
write.csv(clus33.de.markers, paste(outdir,"/DiffExpGeneList_res0.8_1-Clus33-vs-2-All-Others_viaWilcoxonRankTest",date,".csv",sep=""))

clus34.de.markers <- FindMarkers(aprenaCohortAlt, ident.1 = clus34)
# |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=42m 40s
write.csv(clus34.de.markers, paste(outdir,"/DiffExpGeneList_res0.8_1-Clus34-vs-2-All-Others_viaWilcoxonRankTest",date,".csv",sep=""))

clus35.de.markers <- FindMarkers(aprenaCohortAlt, ident.1 = clus35)
# |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=44m 09s
write.csv(clus35.de.markers, paste(outdir,"/DiffExpGeneList_res0.8_1-Clus35-vs-2-All-Others_viaWilcoxonRankTest",date,".csv",sep=""))

clus36.de.markers <- FindMarkers(aprenaCohortAlt, ident.1 = clus36)
# |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=16m 52s
write.csv(clus36.de.markers, paste(outdir,"/DiffExpGeneList_res0.8_1-Clus36-vs-2-All-Others_viaWilcoxonRankTest",date,".csv",sep=""))

clus37.de.markers <- FindMarkers(aprenaCohortAlt, ident.1 = clus37)
# |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=39m 53s
write.csv(clus37.de.markers, paste(outdir,"/DiffExpGeneList_res0.8_1-Clus37-vs-2-All-Others_viaWilcoxonRankTest",date,".csv",sep=""))

#<------------------------------------------- HERE






