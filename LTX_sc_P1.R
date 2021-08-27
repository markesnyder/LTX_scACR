#Code used to generate Figure 3

library(cowplot)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(sctransform)
library(ggpubr)
library(scRepertoire)
library(ggrepel)

acr.data <- Read10X(data.dir = "/Users/MES374/Documents/Seurat/TXR00030/acr/filtered_feature_bc_matrix/")
aP1 <- CreateSeuratObject(counts = acr.data, project = "acrP1", min.cells = 3, min.features = 200)
aP1 <- PercentageFeatureSet(aP1, pattern = "^MT-", col.name = "percent.mt")
aP1 <- SCTransform(aP1, vars.to.regress = "percent.mt", verbose = FALSE)
head(aP1@meta.data, 5)

r_aP1 <- RenameCells(object = aP1, add.cell.id = "2_ACR_P1")

late.data <- Read10X(data.dir = "/Users/MES374/Documents/Seurat/TXR00030/late/outs/filtered_feature_bc_matrix/")
tP1 <- CreateSeuratObject(counts = late.data, project = "treatedP1", min.cells = 3, min.features = 200)
tP1 <- PercentageFeatureSet(tP1, pattern = "^MT-", col.name = "percent.mt")
tP1 <- SCTransform(tP1, vars.to.regress = "percent.mt", verbose = FALSE)

r_tP1 <- RenameCells(object = tP1, add.cell.id = "3_Treated_P1")

p1.list <- list(r_aP1, r_tP1)

#options(future.globals.maxSize= 891289600)
p1.features <- SelectIntegrationFeatures(object.list = p1.list, nfeatures = 3000)
p1.list <- PrepSCTIntegration(object.list = p1.list, anchor.features = p1.features, 
                               verbose = FALSE)
p1.anchors <- FindIntegrationAnchors(object.list = p1.list, normalization.method = "SCT", 
                                      anchor.features = p1.features, verbose = FALSE)
p1.integrated <- IntegrateData(anchorset = p1.anchors, normalization.method = "SCT", 
                                verbose = FALSE)
#Subset to isolate just the T cells
p1 <- subset(p1.integrated, subset = CD68 < 0.2)

p1 <- RunPCA(p1, verbose = FALSE)
p1 <- RunUMAP(p1, dims = 1:10, verbose = FALSE)

p1 <- FindNeighbors(p1, dims = 1:10, verbose = FALSE)
p1 <- FindClusters(p1, verbose = FALSE)
DimPlot(p1, label = TRUE) + NoLegend()
DimPlot(object = p1, group.by = "orig.ident")

seurat2.markers <- FindAllMarkers(seurat2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
seurat2.markers %>%
  group_by(cluster) %>%
  top_n(n = 2, wt = avg_log2FC)


top10 <- seurat2.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)
DoHeatmap(seurat2, features = top10$gene) + NoLegend()

#Saving work that has thusfar been accomplished.  
saveRDS(p1, file = "p1ACRvTnew.gex.rds")

#Loading the TCR data which we will later integrate into our seurat object using scRepertoire
atcr <- read.csv("/Users/MES374/Documents/scLTX/TXR00030/LateACR/TCR/outs/filtered_contig_annotations.csv")
ttcr <- read.csv("/Users/MES374/Documents/scLTX/TXR00030/Treated/TCR/outs/filtered_contig_annotations.csv")

tcr_list <- list(atcr, ttcr)

combined2 <- combineTCR(tcr_list, 
                        samples = c("2_ACR", "3_Treated"), 
                        ID = c("P1", "P1"), 
                        cells ="T-AB")


#Visualizing our clonotype data
quantContig(combined2, cloneCall="gene+nt", scale = T)

quantContig_output <- quantContig(combined2, cloneCall="gene+nt", 
                                  scale = T, exportTable = T)
quantContig_output

abundanceContig(combined2, cloneCall = "gene", scale = F)

compareClonotypes(combined2, numbers = 20, samples = c("2_ACR_P1", "3_Treated_P1"), 
                  cloneCall="aa", graph = "alluvial")

compareClonotypes(combined2, numbers = 20, samples = c("2_ACR_P1", "3_Treated_P1"), 
                  cloneCall="gene", graph = "alluvial")

clonalHomeostasis(combined2, cloneCall = "gene")

clonalProportion(combined2, cloneCall = "gene") 

clonalOverlap(combined2, cloneCall = "gene+nt", method = "morisita")

clonesizeDistribution(combined2, cloneCall = "gene+nt", method="ward.D2")

clonalDiversity(combined2, cloneCall = "gene", group = "samples")

scatterClonotype(combined2, cloneCall ="gene", 
                 x.axis = "1_Early", y.axis = "2_ACR",
                 graph = "proportion")

clonalOverlap(combined2, cloneCall="aa", method="overlap")

lengthContig(combined2, cloneCall="aa", chains = "combined") 
lengthContig(combined2, cloneCall="nt", chains = "single") 

#Now we need to integrate the VDJ data (combined2) with our seurat object using scRepertoire
#load the seurat object which has already been combined and clustered wiht SCTransform
seurat2 <- readRDS(file = "p1ACRvTnew.gex.rds")
DimPlot(seurat2, label = T) + NoLegend()
table(Idents(seurat2))
head(seurat2@meta.data, 5)
tail(seurat2@meta.data, 5)

head(combined2[[1]])

#Contig data does not match: seurat objec has an additional _1, _2, _3, per list
seurat2 <- combineExpression(combined2, seurat2, 
                             cloneCall="gene", groupBy = "sample", proportion = FALSE, 
                             cloneTypes=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))

saveRDS(seurat2, file = "cp1.Seurat2.GexTCR.rds")

#Now we perform our downstream analysis
#First visualizing our clonal expansion in relation to our previous clusters
DimPlot(seurat2, group.by = "orig.ident") + NoLegend() +
  scale_color_manual(values=colorblind_vector(3))


DimPlot(seurat2, group.by = "cloneType") +
  scale_color_manual(values = colorblind_vector(5), na.value="grey")

clonalDiversity(seurat2, cloneCall = "gene", group = "sample")

occupiedscRepertoire(seurat2, x.axis = "cluster")


alluvialClonotypes(seurat2, cloneCall = "gene", 
                   y.axes = c("orig.ident", "cluster"), 
                   color = "TRAV12-2.TRAJ42.TRAC_TRBV20-1.TRBJ2-3.TRBD2.TRBC2") + 
  scale_fill_manual(values = c("grey", colorblind_vector(1)))

#Generating figures
Fig3a <- DimPlot(seurat2, reduction = "umap")

Fig3b <- clonalOverlay(seurat2, reduction = "umap", 
                       freq.cutpoint = 20, bins = 25, facet = "orig.ident") + 
  guides(color = FALSE)

#Looking at clonal expansion by cluster
occupiedscRepertoire(seurat2, x.axis = "cluster")

########################
#Figure 3 C and D, identifying the top expanded clones at the time of ACR then repeat for top expanded clones 
#at the time of successful treatment

seurat2 <- highlightClonotypes(seurat2, cloneCall= "aa", sequence = c("CAVNFGNEKLTF_CASSRGGAQAFF", "NA_CASSRGGAQAFF", "CVVSDLAGKLIF_CASSLSTDLSNQPQHF", "NA_CASSLSTDLSNQPQHF", "CAVREDTGFQKLVF_CASSDHRETGANVLTF", "NA_CASSDHRETGANVLTF", "CAESPNYGGSQGNLIF_CASSFVSGINYEQYF", "NA_CASSFVSGINYEQYF"))
Fig2c <- DimPlot(seurat2, group.by = "highlight", cols = c('red', 'red', 'navyblue', 'navyblue', 'green3', 'green3', 'violetred4', 'violetred4'), pt.size = 0.6, split.by = 'orig.ident')

seurat2 <- highlightClonotypes(seurat2, cloneCall= "aa", sequence = c("CAMTPSGSARQLTF_CSARSAGRLSYNEQFF", "NA_CSARSAGRLSYNEQFF", "CAETVFGNEKLTF_CASSLGPNNQPQHF", "NA_CASSLGPNNQPQHF", "CAVNFGNEKLTF_CASSRGGAQAFF", "NA_CASSRGGAQAFF", "CAESPNYGGSQGNLIF_CASSFVSGINYEQYF", "NA_CASSFVSGINYEQYF"))
Fig2d <- DimPlot(seurat2, group.by = "highlight", cols = c('red', 'red', 'navyblue', 'navyblue', 'green3', 'green3', 'violetred4','violetred4'), pt.size = 0.6, split.by = 'orig.ident')

head(seurat2@meta.data, 5)

DimPlot(seurat2, reduction = "umap")

figure3cd <- ggarrange(Fig2c, Fig2d,
                       labels = c("B","C"), 
                       ncol = 1, nrow = 2)


#generating feature plots for figure 3E

Fig2d <- FeaturePlot(seurat2, features = c("CD8A"), max.cutoff = 1, min.cutoff = -1, pt.size = 0.4,#Defines
                     cols = c("grey", "navyblue"))

Fig2e <- FeaturePlot(seurat2, features = c("TNFRSF4"), max.cutoff = 1, min.cutoff = -1, pt.size = 0.4,#Defines
                     cols = c("grey", "navyblue"))

Fig2f <- FeaturePlot(seurat2, features = c("CCR7"), max.cutoff = 1, min.cutoff = -1, pt.size = 0.4,#Defines
                     cols = c("grey", "navyblue"))

Fig2f1 <- FeaturePlot(seurat2, features = c("S1PR1"), max.cutoff = 1, min.cutoff = -1, pt.size = 0.4,#Defines
                      cols = c("grey", "navyblue"))

Fig2g <- FeaturePlot(seurat2, features = c("TCF7"), max.cutoff = 1, min.cutoff = -1, pt.size = 0.4,#Defines
                     cols = c("grey", "navyblue"))

Fig2h <- FeaturePlot(seurat2, features = c("ITGAE"), max.cutoff = 1, min.cutoff = -1, pt.size = 0.4,#Defines
                     cols = c("grey", "navyblue"))

Fig2h1 <- FeaturePlot(seurat2, features = c("ITGA1"), max.cutoff = 1, min.cutoff = -1, pt.size = 0.4,#Defines
                      cols = c("grey", "navyblue"))

Fig2i <- FeaturePlot(seurat2, features = c("PRDM1"), max.cutoff = 1, min.cutoff = -1, pt.size = 0.4,#Defines
                     cols = c("grey", "navyblue"))

Fig2j <- FeaturePlot(seurat2, features = c("LAG3"), max.cutoff = 1, min.cutoff = -1, pt.size = 0.4,#Defines
                     cols = c("grey", "navyblue"))

Fig2k <- FeaturePlot(seurat2, features = c("CXCR6"), max.cutoff = 1, min.cutoff = -1, pt.size = 0.4,#Defines
                     cols = c("grey", "navyblue"))

Fig2l <- FeaturePlot(seurat2, features = c("KLRC1"), max.cutoff = 1, min.cutoff = -1, pt.size = 0.4,#Defines
                     cols = c("grey", "navyblue"))

Fig2l1 <- FeaturePlot(seurat2, features = c("CRIP1"), max.cutoff = 1, min.cutoff = -1, pt.size = 0.4,#Defines
                      cols = c("grey", "navyblue"))

Fig2m <- FeaturePlot(seurat2, features = c("GZMB"), max.cutoff = 1, min.cutoff = -1, pt.size = 0.4,#Defines
                     cols = c("grey", "navyblue"))

Fig2n <- FeaturePlot(seurat2, features = c("GZMK"), max.cutoff = 1, min.cutoff = -1, pt.size = 0.4,#Defines
                     cols = c("grey", "navyblue"))

Fig2o <- FeaturePlot(seurat2, features = c("PRF1"), max.cutoff = 1, min.cutoff = -1, pt.size = 0.4,#Defines
                     cols = c("grey", "navyblue"))

Fig2o1 <- FeaturePlot(seurat2, features = c("IFNG"), max.cutoff = 1, min.cutoff = -1, pt.size = 0.4,#Defines
                      cols = c("grey", "navyblue"))

fig_3e <- ggarrange(Fig2d, Fig2e, Fig2f, Fig2f1, Fig2g,
                           Fig2h, Fig2h1, Fig2i, Fig2j, Fig2k, 
                           Fig2l, Fig2l1, Fig2m, Fig2n, Fig2o, Fig2o1,
                           #labels = c("D","E","F","G", "H", "I", "J", "J", "K"), 
                           ncol = 4, nrow = 4)

#Figure 3F
#Volcano plot of top clones at time of ACR
levels(seurat2)
head(seurat2@meta.data, 5)
#Set top clones at the time of rejection
top2_acr_clones <- c("CAVNFGNEKLTF_CASSRGGAQAFF", "NA_CASSRGGAQAFF", "CVVSDLAGKLIF_CASSLSTDLSNQPQHF", "NA_CASSLSTDLSNQPQHF", "CAVREDTGFQKLVF_CASSDHRETGANVLTF", "NA_CASSDHRETGANVLTF", "CAESPNYGGSQGNLIF_CASSFVSGINYEQYF", "NA_CASSFVSGINYEQYF")

#subset seurat object to include only those clones
seurat2_top_clones <- subset(seurat2, cells = WhichCells(object = seurat2, expression = CTaa == "CAVNFGNEKLTF_CASSRGGAQAFF" | CTaa == "NA_CASSRGGAQAFF" |CTaa == "CVVSDLAGKLIF_CASSLSTDLSNQPQHF" | CTaa == "NA_CASSLSTDLSNQPQHF" 
                                                         | CTaa =="CAVREDTGFQKLVF_CASSDHRETGANVLTF" | CTaa == "NA_CASSDHRETGANVLTF" | CTaa =="CAESPNYGGSQGNLIF_CASSFVSGINYEQYF" | CTaa == "NA_CASSFVSGINYEQYF"))

head(seurat2_top_clones@meta.data, 5)
tail(seurat2_top_clones@meta.data, 5)
Idents(seurat2_top_clones) <- 'orig.ident'
acr2.response <- FindMarkers(seurat2_top_clones, assay = "SCT", ident.1 = "acr342", ident.2 = "late342", verbose = FALSE)


acr2.response['genes'] <- row.names(acr2.response) 
head(acr2.response)
write.csv(acr2.response, "TopCloneDEG2.csv")

lfc = 0.8
alpha = 0.05

tab2 = data.frame(genes = acr2.response$genes, logFCDP = acr2.response$avg_log2FC, DPpval = -log10(acr2.response$p_val))
signGenes = (abs(tab2$logFCDP) > lfc & tab2$DPpval > -log10(alpha))
signGenesUp = (tab2$logFCDP) > lfc & tab2$DPpval > -log10(alpha)
signGenesDown = (tab2$logFCDP) < -lfc & tab2$DPpval > -log10(alpha)

#Setup genes to list on the graph
tab2$genelabels <- factor(tab2$genes, levels = c("CRIP1", "CCL4", "VIM", "PTPRCAP","GZMB","PRF1",
                                                 "KLRD1", "VCAM1", "CXCL13","LAG3", "MTRNR2L12",
                                                 "TAP1","TUBA1B", "EEF1G", "RPL41", "KLRC1"))

head(tab2)

fig2p <- ggplot(data = tab2, mapping = aes(x=logFCDP, y=DPpval, label = genelabels))+
  geom_point() +
  xlim(-3,3) +
  # ylim(0,20) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("ACR vs Treated") +
  ylab("-log10 pvalue") +  
  geom_vline(xintercept = 0.8, linetype = "dotted") +
  geom_vline(xintercept = -0.8, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "red") +
  geom_hline(yintercept = -log10(alpha), linetype = "dotted") +
  geom_point(color = "blue", data = tab2[signGenes,]) +
  geom_label_repel(col = "black", size = 6, na.rm = TRUE, max.overlaps = 20, box.padding = unit(0.45, "lines"), hjust = 1)


