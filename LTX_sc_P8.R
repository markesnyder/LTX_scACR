#Code for figure 4 - does not include analysis of circulating T cells which was performed 
#using ImmunoSeq software 

#Loading libraries, some of these are superfluous for this specific figure

library(cowplot)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(sctransform)
library(ggpubr)
library(scRepertoire)
library(ggrepel)

#creating seruat objects and formating for merge together and with TCR data
acr.data <- Read10X(data.dir = "pathway-to-P8-acr/filtered_feature_bc_matrix/")
aP8 <- CreateSeuratObject(counts = acr.data, project = "acr_P8", min.cells = 3, min.features = 200)
aP8 <- PercentageFeatureSet(aP8, pattern = "^MT-", col.name = "percent.mt")
aP8 <- SCTransform(aP8, vars.to.regress = "percent.mt", verbose = FALSE)

r_aP8 <- RenameCells(object = aP8, add.cell.id = "2_ACR_P8")

late.data <- Read10X(data.dir = "pathway-to-P8-treated/filtered_feature_bc_matrix/")
tP8 <- CreateSeuratObject(counts = late.data, project = "tx_P8", min.cells = 3, min.features = 200)
tP8 <- PercentageFeatureSet(tP8, pattern = "^MT-", col.name = "percent.mt")
tP8 <- SCTransform(tP8, vars.to.regress = "percent.mt", verbose = FALSE)

r_tP8 <- RenameCells(object = tP8, add.cell.id = "3_Treated_P8")


p8.list <- list(r_aP8, r_tP8)

#options(future.globals.maxSize= 891289600)
p8.features <- SelectIntegrationFeatures(object.list = p8.list, nfeatures = 3000)
p8.list <- PrepSCTIntegration(object.list = p8.list, anchor.features = p8.features, 
                                verbose = FALSE)
p8.anchors <- FindIntegrationAnchors(object.list = p8.list, normalization.method = "SCT", 
                                       anchor.features = p8.features, verbose = FALSE)
p8.integrated <- IntegrateData(anchorset = p8.anchors, normalization.method = "SCT", 
                                 verbose = FALSE)

#Subset to isolate just the T cells
p8 <- subset(p8.integrated, subset = CD68 < 0.2)

#Generating UMAP
p8 <- RunPCA(p8, verbose = FALSE)
p8 <- RunUMAP(p8, dims = 1:10, verbose = FALSE)

p8 <- FindNeighbors(p8, dims = 1:10, verbose = FALSE)
p8 <- FindClusters(p8, verbose = FALSE)
DimPlot(p8, label = TRUE) + NoLegend()
DimPlot(object = p8, group.by = "orig.ident")

#Saving work that has thusfar been accomplished.  
saveRDS(p8, file = "p8ACRvTnew.gex.rds")

#Loading TCR data for just ACR and Treatment
atcr <- read.csv("path-to-P8/ACR/TCR/outs/filtered_contig_annotations.csv")
ttcr <- read.csv("path-to-P8/Treated/TCR/outs/filtered_contig_annotations.csv")

tcr_list <- list(atcr, ttcr)

tail(tcr_list[[1]])

combined <- combineTCR(tcr_list, 
                       samples = c("2_ACR", "3_Treated"), 
                       ID = c("P8", "P8"), 
                       cells ="T-AB")

#view TCR data
quantContig(combined, cloneCall="gene+nt", scale = T)
quantContig_output <- quantContig(combined, cloneCall="gene+nt", 
                                  scale = T, exportTable = T)
quantContig_output
abundanceContig(combined, cloneCall = "gene", scale = F)
clonalHomeostasis(combined, cloneCall = "gene")
clonalProportion(combined, cloneCall = "gene") 
clonalOverlap(combined, cloneCall = "gene+nt", method = "morisita")
clonesizeDistribution(combined, cloneCall = "gene+nt", method="ward.D2")
clonalDiversity(combined, cloneCall = "gene", group = "samples")
scatterClonotype(combined, cloneCall ="gene", 
                 x.axis = "1_Early", y.axis = "2_ACR",
                 graph = "proportion")
lengthContig(combined, cloneCall="aa", chains = "combined") 
lengthContig(combined, cloneCall="nt", chains = "single") 

#Combining clonotype data with seurat object using scRepertoire
seurat <- readRDS(file = "p8ACRvTnew.gex.rds")
DimPlot(seurat, label = T) + NoLegend()
table(Idents(seurat))

#seurat is our seurat object, combined is our clonotypes
seurat <- combineExpression(combined, seurat, 
                            cloneCall="gene", groupBy = "sample", proportion = FALSE, 
                            cloneTypes=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))

saveRDS(seurat, file = "p8ACRvTnew.gexTCR.rds")

#Start analysis 
seurat <- readRDS(file = "p8ACRvTnew.gexTCR.rds")

colorblind_vector <- colorRampPalette(c("#FF4B20", "#FFB433", 
                                        "#C6FDEC", "#7AC5FF", "#0348A6"))
DimPlot(seurat, group.by = "orig.ident") + NoLegend() +
  scale_color_manual(values=colorblind_vector(3))


DimPlot(seurat, group.by = "cloneType") +
  scale_color_manual(values = colorblind_vector(5), na.value="grey")

clonalDiversity(combined, cloneCall = "gene", group = "sample")

#Bar graph of clonal expansion by cluster
occupiedscRepertoire(seurat, x.axis = "cluster")

#Showing clonal overalap by cluster and timepoint
alluvialClonotypes(seurat, cloneCall = "gene", 
                   y.axes = c("orig.ident", "cluster"), 
                   color = "TRAV12-2.TRAJ42.TRAC_TRBV20-1.TRBJ2-3.TRBD2.TRBC2") + 
  scale_fill_manual(values = c("grey", colorblind_vector(1)))

#Figure 4a
Fig_4a <- DimPlot(seurat, label = F) 

#Figure 4b
Fig_4b <- FeaturePlot(seurat, features = c("CD8A"), max.cutoff = 1, min.cutoff = -1, pt.size = 0.4,#Defines
                      cols = c("grey", "navyblue"))
#Figure 4c
#Isolating the top expanded clones at the time of rejection
seurat <- highlightClonotypes(seurat, cloneCall= "aa", sequence = c("CAVGASKAAGNKLTF_CASRRTGRNQPQHF", "NA_CASRRTGRNQPQHF", "CAVPEREGSNYKLTF_CAWSVELGGSNQPQHF", "NA_CAWSVELGGSNQPQHF", "CAVGESKAAGNKLTF_CASRPTGRNQPQHF", "NA_CASRPTGRNQPQHF", "CAYRSLIIIQGAQKLVF_CASTNMRGRGYTF", "NA_CASTNMRGRGYTF"))
Fig_c <- DimPlot(seurat, group.by = "highlight", cols = c('red', 'red', 'navyblue', 'navyblue', 'green3', 'green3', 'violetred4', 'violetred4'), pt.size = 0.6, split.by = 'orig.ident')


############################################################
# Figure 4G - Generating figures for BAL circulating clonotypes
############################################################

#For this figure we identfied top TCRB expanded circulating clones at the time of ACR
#We then identify these clones within our seurat object of the BAL 
#We include both TRA_TRB pairs in the BAL as well as orphaned TRB for each clone. 


#Isolating the top expanded clones at the time of rejection
seurat <- highlightClonotypes(seurat, cloneCall= "aa", sequence = c("CAVGASKAAGNKLTF_CASRRTGRNQPQHF", "NA_CASRRTGRNQPQHF", "CAVPEREGSNYKLTF_CAWSVELGGSNQPQHF", "NA_CAWSVELGGSNQPQHF", "CAVGESKAAGNKLTF_CASRPTGRNQPQHF", "NA_CASRPTGRNQPQHF", "CAYRSLIIIQGAQKLVF_CASTNMRGRGYTF", "NA_CASTNMRGRGYTF"))
Fig_c <- DimPlot(seurat, group.by = "highlight", cols = c('red', 'red', 'navyblue', 'navyblue', 'green3', 'green3', 'violetred4', 'violetred4'), pt.size = 0.6, split.by = 'orig.ident')

#Isolating the top expanded clones in the periphery 
seurat <- highlightClonotypes(seurat, cloneCall= "aa", sequence = c("CALSETALGDMRF_CASGLGQPGQGDTGELFF", "NA_CASSLWTGASNTGELFF", "CALDMSEIKAAGNKLTF_CASGGRYGYTF", "NA_CASGGRYGYTF", "CAVRARMDSSYKLIF_CASSPPSGAGATNEKLFF", "NA_CASSPPSGAGATNEKLFF", "CALSEPQGGSEKLVF_CASSPWTPGDEQYF", "NA_CASSPWTPGDEQYF", "NA_CASSSHEYRGLLTEAFF", "NA_CASSTRDRGSYEQYF"))
Fig_d <- DimPlot(seurat, group.by = "highlight", cols = c('red', 'red', 'violetred4', 'violetred4', 'navyblue', 'navyblue','green3', 'green3', 'orange', 'orange','orange'), pt.size = 0.6, split.by = 'orig.ident')


#Most abundant clones across all samples
seurat <- highlightClonotypes(seurat, cloneCall= "aa", sequence = c("CALSETALGDMRF_CASGLGQPGQGDTGELFF", "CAVRARMDSSYKLIF_CASSPPSGAGATNEKLFF", "NA_CASSTRDRGSYEQYF", "NA_CASSSHEYRGLLTEAFF" ))
Fig_d <- DimPlot(seurat, group.by = "highlight", cols = c('red', 'red', 'red', 'red'), pt.size = 0.6, split.by = 'orig.ident')


#isolating what is CD8 and what is CD4
new_Fig_b <- FeaturePlot(seurat, features = c("CD8A", "CD4"), max.cutoff = 1, min.cutoff = -1, pt.size = 0.4,#Defines
                         cols = c("grey", "navyblue"))

#Isolating the top expanded clones at the time of rejection
seurat <- highlightClonotypes(seurat, cloneCall= "aa", sequence = c("CAVGASKAAGNKLTF_CASRRTGRNQPQHF", "NA_CASRRTGRNQPQHF", "CAVPEREGSNYKLTF_CAWSVELGGSNQPQHF", "NA_CAWSVELGGSNQPQHF", "CAVGESKAAGNKLTF_CASRPTGRNQPQHF", "NA_CASRPTGRNQPQHF", "CAYRSLIIIQGAQKLVF_CASTNMRGRGYTF", "NA_CASTNMRGRGYTF"))
new_Fig_c <- DimPlot(seurat, group.by = "highlight", cols = c('red', 'red', 'navyblue', 'navyblue', 'green3', 'green3', 'violetred4', 'violetred4'), pt.size = 0.6, split.by = 'orig.ident')


figure4ac <- ggarrange(new_Fig_b,new_Fig_c,
                       labels = c("B","C"), 
                       ncol = 2, nrow = 1)


