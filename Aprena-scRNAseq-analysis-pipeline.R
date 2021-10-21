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

##############################################
#### DEFINE DIRECTORIES
##############################################
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
# [1]  36601 106751 (losst ~28% of cells with this filtering)

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
#<------------------------------------------- HERE
##############################################
#### Harmony to remove batch effects
##############################################
aprenaCohort <- RunHarmony(aprenaCohort, group.by.vars="Tx")

saveRDS(aprenaCohort, paste(outdir,"/Aprena-scRNAseq-loaded+merged-postHarmony_2021-10-21.rds",sep=""))

##############################################
#### Find Neighbors, Cluster and Visualize with UMAP & Seurat Defaults
##############################################
aprenaCohort <- FindNeighbors(aprenaCohort, dims = 1:30)
aprenaCohort <- FindClusters(aprenaCohort, resolution = 0.8, verbose = FALSE)
aprenaCohort <- RunUMAP(aprenaCohort, dims = 1:30)

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
#### DIFFERNETIAL GENE EXPRESSION - BY CONDITIONS OF INTEREST
##############################################


##############################################
#### EXPORT DATA FOR DIVERSITY ANALYSIS
##############################################














