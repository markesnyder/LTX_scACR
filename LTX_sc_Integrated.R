#Code that went into creating Figure 6
#first we integrated just the three ACR files and compared Top expanded clones versus all others
#next we integrated all 6 files (3 ACR and 3 Treatment), subsetted only the top expanded clones
#and compared the gene expression between clones at the time of ACR and after treatment
#For this aspect of the analysis, we just embedded the TCR data into the metadata manually

library(cowplot)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(sctransform)
library(ggpubr)
library(scRepertoire)
library(ggrepel)

###############################################################################
# Creating the seurat object and adding VDJ data to seurat object metadata manually
###############################################################################

#Step 1 & 2: Creat Seurat objects to integrate and run SCTransform individually

#P8 data
acr.data <- Read10X(data.dir = "/Users/MES374/Documents/Seurat/TXR00342/acr/filtered_feature_bc_matrix/")
aP8 <- CreateSeuratObject(counts = acr.data, project = "acrP8", min.cells = 3, min.features = 200)
aP8 <- PercentageFeatureSet(aP8, pattern = "^MT-", col.name = "percent.mt")
aP8 <- SCTransform(aP8, vars.to.regress = "percent.mt", verbose = FALSE)
head(aP8@meta.data, 5)

#adding clonotype data
atcr <- read.csv("/Users/MES374/Documents/scLTX/TXR00342/LateACR/TCR/outs/filtered_contig_annotations.csv")
atcr <- atcr[!duplicated(atcr$barcode), ]
atcr <- atcr[,c("barcode", "raw_clonotype_id")]
names(atcr)[names(atcr) == "raw_clonotype_id"] <- "clonotype_id"
clono <- read.csv("/Users/MES374/Documents/scLTX/TXR00342/LateACR/TCR/outs/clonotypes.csv")
atcr <- merge(atcr, clono[, c("clonotype_id", "cdr3s_aa")])
names(atcr)[names(atcr) == "cdr3s_aa"] <- "cdr3s_aa"
atcr <- atcr[, c(2,1,3)]
rownames(atcr) <- atcr[,1]
atcr[,1] <- NULL
head(atcr)
aP8 <- AddMetaData(object=aP8, metadata=atcr)

#include only those with a clonotype
aP8 <- subset(aP8, cells = colnames(aP8)[(!is.na(aP8$clonotype_id))])
head(aP8@meta.data, 5) 


late.data <- Read10X(data.dir = "/Users/MES374/Documents/Seurat/TXR00342/late/filtered_feature_bc_matrix/")
tP8 <- CreateSeuratObject(counts = late.data, project = "treatedP8", min.cells = 3, min.features = 200)
tP8 <- PercentageFeatureSet(tP8, pattern = "^MT-", col.name = "percent.mt")
tP8 <- SCTransform(tP8, vars.to.regress = "percent.mt", verbose = FALSE)


#adding clonotype data
ltcr <- read.csv("/Users/MES374/Documents/scLTX/TXR00342/Treated/TCR/outs/filtered_contig_annotations.csv")
ltcr <- ltcr[!duplicated(ltcr$barcode), ]
ltcr <- ltcr[,c("barcode", "raw_clonotype_id")]
names(ltcr)[names(ltcr) == "raw_clonotype_id"] <- "clonotype_id"
clono <- read.csv("/Users/MES374/Documents/scLTX/TXR00342/Treated/TCR/outs/clonotypes2.csv")
ltcr <- merge(ltcr, clono[, c("clonotype_id", "cdr3s_aa")])
names(ltcr)[names(ltcr) == "cdr3s_aa"] <- "cdr3s_aa"
ltcr <- ltcr[, c(2,1,3)]
rownames(ltcr) <- ltcr[,1]
ltcr[,1] <- NULL
head(ltcr)
tP8 <- AddMetaData(object=tP8, metadata=ltcr)

tP8 <- subset(tP8, cells = colnames(tP8)[(!is.na(tP8$clonotype_id))])
head(l342@meta.data, 5)

#P1 data
acr.data <- Read10X(data.dir = "/Users/MES374/Documents/Seurat/TXR00030/acr/filtered_feature_bc_matrix/")
aP1 <- CreateSeuratObject(counts = acr.data, project = "acrP1", min.cells = 3, min.features = 200)
aP1 <- PercentageFeatureSet(aP1, pattern = "^MT-", col.name = "percent.mt")
aP1 <- SCTransform(aP1, vars.to.regress = "percent.mt", verbose = FALSE)
head(aP1@meta.data, 5)

#adding clonotype data
atcr <- read.csv("/Users/MES374/Documents/scLTX/TXR00030/LateACR/TCR/outs/filtered_contig_annotations.csv")
atcr <- atcr[!duplicated(atcr$barcode), ]
atcr <- atcr[,c("barcode", "raw_clonotype_id")]
names(atcr)[names(atcr) == "raw_clonotype_id"] <- "clonotype_id"
clono <- read.csv("/Users/MES374/Documents/scLTX/TXR00030/LateACR/TCR/outs/clonotypes.csv")
atcr <- merge(atcr, clono[, c("clonotype_id", "cdr3s_aa")])
names(atcr)[names(atcr) == "cdr3s_aa"] <- "cdr3s_aa"
atcr <- atcr[, c(2,1,3)]
rownames(atcr) <- atcr[,1]
atcr[,1] <- NULL
head(atcr)
aP1 <- AddMetaData(object=aP1, metadata=atcr)

#include only those with a clonotype
aP1 <- subset(aP1, cells = colnames(aP1)[(!is.na(aP1$clonotype_id))])
head(aP1@meta.data, 5) 

late.data <- Read10X(data.dir = "/Users/MES374/Documents/Seurat/TXR00030/late/outs/filtered_feature_bc_matrix/")
tP1 <- CreateSeuratObject(counts = late.data, project = "treatedP1", min.cells = 3, min.features = 200)
tP1 <- PercentageFeatureSet(tP1, pattern = "^MT-", col.name = "percent.mt")
tP1 <- SCTransform(tP1, vars.to.regress = "percent.mt", verbose = FALSE)


#adding clonotype data
ltcr <- read.csv("/Users/MES374/Documents/scLTX/TXR00030/Treated/TCR/outs/filtered_contig_annotations.csv")
ltcr <- ltcr[!duplicated(ltcr$barcode), ]
ltcr <- ltcr[,c("barcode", "raw_clonotype_id")]
names(ltcr)[names(ltcr) == "raw_clonotype_id"] <- "clonotype_id"
clono <- read.csv("/Users/MES374/Documents/scLTX/TXR00030/Treated/TCR/outs/clonotypes2.csv")
ltcr <- merge(ltcr, clono[, c("clonotype_id", "cdr3s_aa")])
names(ltcr)[names(ltcr) == "cdr3s_aa"] <- "cdr3s_aa"
ltcr <- ltcr[, c(2,1,3)]
rownames(ltcr) <- ltcr[,1]
ltcr[,1] <- NULL
head(ltcr)
tP1 <- AddMetaData(object=tP1, metadata=ltcr)

tP1 <- subset(tP1, cells = colnames(tP1)[(!is.na(tP1$clonotype_id))])
head(tP1@meta.data, 5)

#P3 data 
acr.data <- Read10X(data.dir = "/Users/MES374/Documents/scLTX/TXR00235/LateACR/filtered_feature_bc_matrix/")
aP3 <- CreateSeuratObject(counts = acr.data, project = "acrP3", min.cells = 3, min.features = 200)
aP3 <- PercentageFeatureSet(aP3, pattern = "^MT-", col.name = "percent.mt")
aP3 <- SCTransform(aP3, vars.to.regress = "percent.mt", verbose = FALSE)
head(aP3@meta.data, 5)

#adding clonotype data
atcr <- read.csv("/Users/MES374/Documents/scLTX/TXR00235/LateACR/TCR/outs/filtered_contig_annotations.csv")
atcr <- atcr[!duplicated(atcr$barcode), ]
atcr <- atcr[,c("barcode", "raw_clonotype_id")]
names(atcr)[names(atcr) == "raw_clonotype_id"] <- "clonotype_id"
clono <- read.csv("/Users/MES374/Documents/scLTX/TXR00235/LateACR/TCR/outs/clonotypes.csv")
atcr <- merge(atcr, clono[, c("clonotype_id", "cdr3s_aa")])
names(atcr)[names(atcr) == "cdr3s_aa"] <- "cdr3s_aa"
atcr <- atcr[, c(2,1,3)]
rownames(atcr) <- atcr[,1]
atcr[,1] <- NULL
head(atcr)
aP3 <- AddMetaData(object=aP3, metadata=atcr)

#include only those with a clonotype
aP3 <- subset(aP3, cells = colnames(aP3)[(!is.na(aP3$clonotype_id))])
head(aP3@meta.data, 5) 

late.data <- Read10X(data.dir = "/Users/MES374/Documents/scLTX/TXR00235/Treated/filtered_feature_bc_matrix/")
tP3 <- CreateSeuratObject(counts = late.data, project = "treatedP3", min.cells = 3, min.features = 200)
tP3 <- PercentageFeatureSet(tP3, pattern = "^MT-", col.name = "percent.mt")
tP3 <- SCTransform(tP3, vars.to.regress = "percent.mt", verbose = FALSE)


#adding clonotype data
ltcr <- read.csv("/Users/MES374/Documents/scLTX/TXR00235/Treated/TCR/outs/filtered_contig_annotations.csv")
ltcr <- ltcr[!duplicated(ltcr$barcode), ]
ltcr <- ltcr[,c("barcode", "raw_clonotype_id")]
names(ltcr)[names(ltcr) == "raw_clonotype_id"] <- "clonotype_id"
clono <- read.csv("/Users/MES374/Documents/scLTX/TXR00235/Treated/TCR/outs/clonotypes2.csv")
ltcr <- merge(ltcr, clono[, c("clonotype_id", "cdr3s_aa")])
names(ltcr)[names(ltcr) == "cdr3s_aa"] <- "cdr3s_aa"
ltcr <- ltcr[, c(2,1,3)]
rownames(ltcr) <- ltcr[,1]
ltcr[,1] <- NULL
head(ltcr)
tP3 <- AddMetaData(object=tP3, metadata=ltcr)

tP3 <- subset(tP3, cells = colnames(tP3)[(!is.na(tP3$clonotype_id))])
head(tP3@meta.data, 5)


###############################################################################
# Integrate the data
###############################################################################
options(future.globals.maxSize = 1000 * 1024^2)

#Step 3: creating list to integrate and running PrepSCTIntegration
concat.list <- list(aP8, tP8, aP1, tP1, aP3, tP3)

concat.features <- SelectIntegrationFeatures(object.list = concat.list, nfeatures = 3000)
concat.list <- PrepSCTIntegration(object.list = concat.list, anchor.features = concat.features, 
                                  verbose = FALSE)

#Step 4: Integrate the datasets
concat.anchors <- FindIntegrationAnchors(object.list = concat.list, normalization.method = "SCT", 
                                         anchor.features = concat.features, verbose = FALSE)
concat.integrated <- IntegrateData(anchorset = concat.anchors, normalization.method = "SCT", 
                                   verbose = FALSE)

#Subset to isolate just the T cells
concat.integrated <- subset(concat.integrated, subset = CD68 < 0.2)

#Next we begin the downstream analysis, but first we need to normalize by RNA
concat.integrated <- RunPCA(concat.integrated, verbose = FALSE)
concat.integrated <- RunUMAP(concat.integrated, dims = 1:30)
head(concat.integrated@meta.data, 5)

DefaultAssay(concat.integrated) <- "RNA"
# Normalize RNA data for visualization purposes
concat.integrated <- NormalizeData(concat.integrated, verbose = FALSE)

#Making sure our data looks ok
DimPlot(concat.integrated, group.by = "orig.ident", combine = FALSE)
FeaturePlot(concat.integrated, features = c("GZMB"), split.by = "orig.ident", max.cutoff = 3, 
            cols = c("grey", "red"))
FeaturePlot(concat.integrated, features = c("CD8A"), max.cutoff = 3, 
            cols = c("grey", "red"))

FeaturePlot(concat.integrated, features = c("GZMB"), max.cutoff = 3, 
            cols = c("grey", "red"))

FeaturePlot(concat.integrated, features = c("ITGAE"), max.cutoff = 3, 
            cols = c("grey", "red"))


#Just reducing to has clonotypeID
table(!is.na(concat.integrated$clonotype_id))

RidgePlot(concat.integrated, features="nCount_RNA")
DimPlot(concat.integrated, reduction = "umap", label = TRUE)

#Saving work that has thusfar been accomplished.  
saveRDS(concat.integrated, file = "concat.integrated.rds")


###############################################################################
# Can start analysis from here, downloading work from prior analysis
###############################################################################
concatTCR <- readRDS(file = "concat.integrated.rds")
head(concatTCR@meta.data, 5)

#Need to label rejection versus treatment 
concatTCR$condition <- plyr::mapvalues(
  x = concatTCR$orig.ident, 
  from = c("acrP8", "treatedP8", "acrP1", "treatedP1", "acrP3", "treatedP3"), 
  to = c("Rejection", "Treated", "Rejection", "Treated", "Rejection", "Treated")
)

head(concatTCR@meta.data, 5)

#Generate UMAP showing top clones at ACR and Treatment
DimPlot(object = concatTCR, cells.highlight = WhichCells(object = concatTCR, expression = cdr3s_aa == "TRA:CAVNFGNEKLTF;TRB:CASSRGGAQAFF" | cdr3s_aa =="TRB:CASSRGGAQAFF"| cdr3s_aa =="TRA:CVVSDLAGKLIF;TRB:CASSLSTDLSNQPQHF"| cdr3s_aa =="TRB:CASSLSTDLSNQPQHF"
                                                    | cdr3s_aa =="TRA:CAVREDTGFQKLVF;TRB:CASSDHRETGANVLTF" | cdr3s_aa =="TRB:CASSDHRETGANVLTF" | cdr3s_aa =="TRA:CAESPNYGGSQGNLIF;TRB:CASSFVSGINYEQYF" | cdr3s_aa =="TRB:CASSFVSGINYEQYF"
                                                    | cdr3s_aa =="TRA:CAVGASKAAGNKLTF;TRB:CASRRTGRNQPQHF" |cdr3s_aa =="TRB:CASRRTGRNQPQHF" |cdr3s_aa =="TRA:CAVPEREGSNYKLTF;TRB:CAWSVELGGSNQPQHF" | cdr3s_aa =="TRB:CAWSVELGGSNQPQHF"
                                                    | cdr3s_aa =="TRA:CAYRSLIIIQGAQKLVF;TRB:CASTNMRGRGYTF"| cdr3s_aa =="TRB:CASTNMRGRGYTF" | cdr3s_aa =="TRA:CAVGESKAAGNKLTF;TRB:CASRPTGRNQPQHF"| cdr3s_aa == "TRB:CASRPTGRNQPQHF"
                                                    | cdr3s_aa =="TRA:CSPMNSGGYQKVTF;TRB:CASKQTERSSYNEQFF"| cdr3s_aa =="TRB:CASKQTERSSYNEQFF"| cdr3s_aa =="TRA:CAASEGNARLMF;TRB:CATSRLIGLGQGIPGELFF"| cdr3s_aa =="TRB:CATSRLIGLGQGIPGELFF"
                                                    | cdr3s_aa =="TRA:CAYRSFNAGKSTF;TRB:CASSGNSGGTPYEQYF" | cdr3s_aa =="TRB:CASSGNSGGTPYEQYF"| cdr3s_aa =="TRA:CAASTGTASKLTF;TRB:CASSLRRSEQFF" | cdr3s_aa =="TRB:CASSLRRSEQFF"), split.by = "condition")


###############################################################################
# Generate DEG for top 4 clones from ACR versus Treatment 
###############################################################################

#subset seurat object to include only those clones
concat_top_clones <- subset(concatTCR, cells = WhichCells(object = concatTCR, expression = cdr3s_aa == "TRA:CAVNFGNEKLTF;TRB:CASSRGGAQAFF" | cdr3s_aa == "TRB:CASSRGGAQAFF" |cdr3s_aa == "TRA:CVVSDLAGKLIF;TRB:CASSLSTDLSNQPQHF" | cdr3s_aa == "TRB:CASSLSTDLSNQPQHF" #TXR00030
                                                     | cdr3s_aa =="TRA:CAVREDTGFQKLVF;TRB:CASSDHRETGANVLTF" | cdr3s_aa == "TRB:CASSDHRETGANVLTF" | cdr3s_aa =="TRA:CAESPNYGGSQGNLIF;TRB:CASSFVSGINYEQYF" | cdr3s_aa == "TRB:CASSFVSGINYEQYF" #TXR00030
                                                     | cdr3s_aa =="TRA:CAVGASKAAGNKLTF;TRB:CASRRTGRNQPQHF" | cdr3s_aa == "TRB:CASRRTGRNQPQHF" | cdr3s_aa =="TRA:CAVPEREGSNYKLTF;TRB:CAWSVELGGSNQPQHF" | cdr3s_aa == "TRB:CAWSVELGGSNQPQHF" #TXR00342
                                                     | cdr3s_aa =="TRA:CAYRSLIIIQGAQKLVF;TRB:CASTNMRGRGYTF" | cdr3s_aa == "TRB:CASTNMRGRGYTF" | cdr3s_aa =="TRA:CAVGESKAAGNKLTF;TRB:CASRPTGRNQPQHF" | cdr3s_aa == "TRB:CASRPTGRNQPQHF"  #TXR00342
                                                     | cdr3s_aa =="TRA:CSPMNSGGYQKVTF;TRB:CASKQTERSSYNEQFF" | cdr3s_aa == "TRB:CASKQTERSSYNEQFF" | cdr3s_aa =="TRA:CAASEGNARLMF;TRB:CATSRLIGLGQGIPGELFF" | cdr3s_aa == "TRB:CATSRLIGLGQGIPGELFF" #TXR00235
                                                     | cdr3s_aa =="TRA:CAYRSFNAGKSTF;TRB:CASSGNSGGTPYEQYF" | cdr3s_aa == "TRB:CASSGNSGGTPYEQYF" | cdr3s_aa =="TRA:CAASTGTASKLTF;TRB:CASSLRRSEQFF" | cdr3s_aa == "TRB:CASSLRRSEQFF")) #TXR00235


head(concat_top_clones@meta.data, 5)
tail(concat_top_clones@meta.data, 5)
Idents(concat_top_clones) <- 'condition'
ConTop4.ACRvTx.response <- FindMarkers(concat_top_clones, assay = "SCT", ident.1 = "Rejection", ident.2 = "Treated", verbose = FALSE)


ConTop4.ACRvTx.response['genes'] <- row.names(ConTop4.ACRvTx.response) 
head(ConTop4.ACRvTx.response)
write.csv(ConTop4.ACRvTx.response, "ConTop4.ACRvTx.DEG.csv")

head(concat_top_clones@meta.data, 5)
Idents(concat_top_clones) <- 'orig.ident'
cluster.averages <- AverageExpression(concat_top_clones, return.seurat = TRUE)
head(cluster.averages@meta.data, 5)

Idents(concat_top_clones) <- 'orig.ident'

DoHeatmap(cluster.averages, features = unlist(TopFeatures(concat_top_clones[["pca"]], balanced = TRUE)), size = 3, 
          draw.lines = FALSE)

levels(cluster.averages)
levels(cluster.averages) <- c("acrP8", "treatedP8", "acrP1", "treatedP1", "acrP3", "treatedP3")

gene.list <- c("GZMB","GZMK","TBX21","EOMES","ZFP36","PRF1","IFNG","LAMP1","CCL5","LAG3", "PRDM1","TNF",  "JUN", "TNFAIP3","MCL1","ZNF683",   "ITGAE", "CD69", "PDCD1", "IL32", "DUSP4","KLF6","PTPRC","CREM","CXCR6", "NFKB2", "RUNX3","CRIP1", "NME2")

DoHeatmap(cluster.averages, features = gene.list,  size = 3, disp.min = -1, disp.max = 1,
          draw.lines = FALSE)




###############################################################################
# Concat only the ACR samples and compare top clones versus not top clones in DEG
###############################################################################

#Generate concat file from just the ACR files
options(future.globals.maxSize = 1000 * 1024^2)

#Step 3: creating list to integrate and running PrepSCTIntegration
concat.list <- list(aP8, aP1, aP3)

concat.features <- SelectIntegrationFeatures(object.list = concat.list, nfeatures = 3000)
concat.list <- PrepSCTIntegration(object.list = concat.list, anchor.features = concat.features, 
                                  verbose = FALSE)

#Step 4: Integrate the datasets
concat.anchors <- FindIntegrationAnchors(object.list = concat.list, normalization.method = "SCT", 
                                         anchor.features = concat.features, verbose = FALSE)
concat.integrated <- IntegrateData(anchorset = concat.anchors, normalization.method = "SCT", 
                                   verbose = FALSE)

#Subset to isolate just the T cells
concat.integrated <- subset(concat.integrated, subset = CD68 < 0.2)

#Next we begin the downstream analysis, but first we need to normalize by RNA
concat.integrated <- RunPCA(concat.integrated, verbose = FALSE)
concat.integrated <- RunUMAP(concat.integrated, dims = 1:30)
head(concat.integrated@meta.data, 5)

DefaultAssay(concat.integrated) <- "RNA"
# Normalize RNA data for visualization purposes
concat.integrated <- NormalizeData(concat.integrated, verbose = FALSE)

#Making sure our data looks ok
DimPlot(concat.integrated, group.by = "orig.ident", combine = FALSE)

#Saving work that has thusfar been accomplished.  
saveRDS(concat.integrated, file = "concat.ACRonly.rds")

###############################################################################
# Start here for DEG analysis of Top Clones versus not top clones for ACR only
###############################################################################
#load concat ACR file
cACR <- readRDS(file = "concat.ACRonly.rds")
head(cACR@meta.data, 5)

#See what file looks like
DimPlot(cACR, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

FeaturePlot(cACR, features = c("S1PR1"), max.cutoff = 1, min.cutoff = -1,
            cols = c("grey", "red"))

#Plot Concat file highlighting top clones we wish to compare: 
DimPlot(object = cACR, cells.highlight = WhichCells(object = cACR, expression = cdr3s_aa == "TRA:CAVNFGNEKLTF;TRB:CASSRGGAQAFF" | cdr3s_aa =="TRB:CASSRGGAQAFF"| cdr3s_aa =="TRA:CVVSDLAGKLIF;TRB:CASSLSTDLSNQPQHF"| cdr3s_aa =="TRB:CASSLSTDLSNQPQHF"
                                                    | cdr3s_aa =="TRA:CAVREDTGFQKLVF;TRB:CASSDHRETGANVLTF" | cdr3s_aa =="TRB:CASSDHRETGANVLTF" | cdr3s_aa =="TRA:CAESPNYGGSQGNLIF;TRB:CASSFVSGINYEQYF" | cdr3s_aa =="TRB:CASSFVSGINYEQYF"
                                                    | cdr3s_aa =="TRA:CAVGASKAAGNKLTF;TRB:CASRRTGRNQPQHF" |cdr3s_aa =="TRB:CASRRTGRNQPQHF" |cdr3s_aa =="TRA:CAVPEREGSNYKLTF;TRB:CAWSVELGGSNQPQHF" | cdr3s_aa =="TRB:CAWSVELGGSNQPQHF"
                                                    | cdr3s_aa =="TRA:CAYRSLIIIQGAQKLVF;TRB:CASTNMRGRGYTF"| cdr3s_aa =="TRB:CASTNMRGRGYTF" | cdr3s_aa =="TRA:CAVGESKAAGNKLTF;TRB:CASRPTGRNQPQHF"| cdr3s_aa == "TRB:CASRPTGRNQPQHF"
                                                    | cdr3s_aa =="TRA:CSPMNSGGYQKVTF;TRB:CASKQTERSSYNEQFF"| cdr3s_aa =="TRB:CASKQTERSSYNEQFF"| cdr3s_aa =="TRA:CAASEGNARLMF;TRB:CATSRLIGLGQGIPGELFF"| cdr3s_aa =="TRB:CATSRLIGLGQGIPGELFF"
                                                    | cdr3s_aa =="TRA:CAYRSFNAGKSTF;TRB:CASSGNSGGTPYEQYF" | cdr3s_aa =="TRB:CASSGNSGGTPYEQYF"| cdr3s_aa =="TRA:CAASTGTASKLTF;TRB:CASSLRRSEQFF" | cdr3s_aa =="TRB:CASSLRRSEQFF"))


#Add two new columns in metadata - we will attempt to make this Yes, if within the top clones
cACR$top4 <- "No"
cACR$top4m <- cACR$cdr3s_aa
head(cACR@meta.data, 5)

#Replace the No with Yes for all samples where top clone is present

# Generating list of top 4 clones from each sample at the time of ACR
top4_acr_clones <- c("TRA:CAVNFGNEKLTF;TRB:CASSRGGAQAFF", "TRB:CASSRGGAQAFF", "TRA:CVVSDLAGKLIF;TRB:CASSLSTDLSNQPQHF", "TRB:CASSLSTDLSNQPQHF", #TXR00030
                     "TRA:CAVREDTGFQKLVF;TRB:CASSDHRETGANVLTF", "TRB:CASSDHRETGANVLTF", "TRA:CAESPNYGGSQGNLIF;TRB:CASSFVSGINYEQYF", "TRB:CASSFVSGINYEQYF", #TXR00030
                     "TRA:CAVGASKAAGNKLTF;TRB:CASRRTGRNQPQHF", "TRB:CASRRTGRNQPQHF", "TRA:CAVPEREGSNYKLTF;TRB:CAWSVELGGSNQPQHF", "TRB:CAWSVELGGSNQPQHF", #TXR00342
                     "TRA:CAYRSLIIIQGAQKLVF;TRB:CASTNMRGRGYTF", "TRB:CASTNMRGRGYTF", "TRA:CAVGESKAAGNKLTF;TRB:CASRPTGRNQPQHF", "TRB:CASRPTGRNQPQHF",  #TXR00342
                     "TRA:CSPMNSGGYQKVTF;TRB:CASKQTERSSYNEQFF", "TRB:CASKQTERSSYNEQFF", "TRA:CAASEGNARLMF;TRB:CATSRLIGLGQGIPGELFF", "TRB:CATSRLIGLGQGIPGELFF", #TXR00235
                     "TRA:CAYRSFNAGKSTF;TRB:CASSGNSGGTPYEQYF" , "TRB:CASSGNSGGTPYEQYF", "TRA:CAASTGTASKLTF;TRB:CASSLRRSEQFF" , "TRB:CASSLRRSEQFF") #TXR00235

cACR$top4 <- "Yes"

#Subset those cells that are top clones
cACR_top <- subset(cACR, cells = WhichCells(object = cACR, expression = cdr3s_aa == "TRA:CAVNFGNEKLTF;TRB:CASSRGGAQAFF" | cdr3s_aa == "TRB:CASSRGGAQAFF" |cdr3s_aa == "TRA:CVVSDLAGKLIF;TRB:CASSLSTDLSNQPQHF" | cdr3s_aa == "TRB:CASSLSTDLSNQPQHF" #TXR00030
                                            | cdr3s_aa =="TRA:CAVREDTGFQKLVF;TRB:CASSDHRETGANVLTF" | cdr3s_aa == "TRB:CASSDHRETGANVLTF" | cdr3s_aa =="TRA:CAESPNYGGSQGNLIF;TRB:CASSFVSGINYEQYF" | cdr3s_aa == "TRB:CASSFVSGINYEQYF" #TXR00030
                                            | cdr3s_aa =="TRA:CAVGASKAAGNKLTF;TRB:CASRRTGRNQPQHF" | cdr3s_aa == "TRB:CASRRTGRNQPQHF" | cdr3s_aa =="TRA:CAVPEREGSNYKLTF;TRB:CAWSVELGGSNQPQHF" | cdr3s_aa == "TRB:CAWSVELGGSNQPQHF" #TXR00342
                                            | cdr3s_aa =="TRA:CAYRSLIIIQGAQKLVF;TRB:CASTNMRGRGYTF" | cdr3s_aa == "TRB:CASTNMRGRGYTF" | cdr3s_aa =="TRA:CAVGESKAAGNKLTF;TRB:CASRPTGRNQPQHF" | cdr3s_aa == "TRB:CASRPTGRNQPQHF"  #TXR00342
                                            | cdr3s_aa =="TRA:CSPMNSGGYQKVTF;TRB:CASKQTERSSYNEQFF" | cdr3s_aa == "TRB:CASKQTERSSYNEQFF" | cdr3s_aa =="TRA:CAASEGNARLMF;TRB:CATSRLIGLGQGIPGELFF" | cdr3s_aa == "TRB:CATSRLIGLGQGIPGELFF" #TXR00235
                                            | cdr3s_aa =="TRA:CAYRSFNAGKSTF;TRB:CASSGNSGGTPYEQYF" | cdr3s_aa == "TRB:CASSGNSGGTPYEQYF" | cdr3s_aa =="TRA:CAASTGTASKLTF;TRB:CASSLRRSEQFF" | cdr3s_aa == "TRB:CASSLRRSEQFF")) #TXR00235

cACR_top$top4 <- "Yes"
head(cACR_top@meta.data, 5)

#Subset those that are teh bottom clones
cACR_bottom <- subset(cACR, cells = WhichCells(object = cACR, expression = cdr3s_aa != "TRA:CAVNFGNEKLTF;TRB:CASSRGGAQAFF" & cdr3s_aa != "TRB:CASSRGGAQAFF" & cdr3s_aa != "TRA:CVVSDLAGKLIF;TRB:CASSLSTDLSNQPQHF" & cdr3s_aa != "TRB:CASSLSTDLSNQPQHF" #TXR00030
                                               & cdr3s_aa !="TRA:CAVREDTGFQKLVF;TRB:CASSDHRETGANVLTF" & cdr3s_aa != "TRB:CASSDHRETGANVLTF" & cdr3s_aa !="TRA:CAESPNYGGSQGNLIF;TRB:CASSFVSGINYEQYF" & cdr3s_aa != "TRB:CASSFVSGINYEQYF" #TXR00030
                                               & cdr3s_aa !="TRA:CAVGASKAAGNKLTF;TRB:CASRRTGRNQPQHF" & cdr3s_aa != "TRB:CASRRTGRNQPQHF" & cdr3s_aa !="TRA:CAVPEREGSNYKLTF;TRB:CAWSVELGGSNQPQHF" & cdr3s_aa != "TRB:CAWSVELGGSNQPQHF" #TXR00342
                                               & cdr3s_aa !="TRA:CAYRSLIIIQGAQKLVF;TRB:CASTNMRGRGYTF" & cdr3s_aa != "TRB:CASTNMRGRGYTF" & cdr3s_aa !="TRA:CAVGESKAAGNKLTF;TRB:CASRPTGRNQPQHF" & cdr3s_aa != "TRB:CASRPTGRNQPQHF"  #TXR00342
                                               & cdr3s_aa !="TRA:CSPMNSGGYQKVTF;TRB:CASKQTERSSYNEQFF" & cdr3s_aa != "TRB:CASKQTERSSYNEQFF" & cdr3s_aa != "TRA:CAASEGNARLMF;TRB:CATSRLIGLGQGIPGELFF" & cdr3s_aa != "TRB:CATSRLIGLGQGIPGELFF" #TXR00235
                                               & cdr3s_aa != "TRA:CAYRSFNAGKSTF;TRB:CASSGNSGGTPYEQYF" & cdr3s_aa != "TRB:CASSGNSGGTPYEQYF" & cdr3s_aa != "TRA:CAASTGTASKLTF;TRB:CASSLRRSEQFF" & cdr3s_aa != "TRB:CASSLRRSEQFF")) #TXR00235

cACR_bottom$top4 <- "No"
head(cACR_bottom@meta.data, 5)

#Merge the two files
combo <- merge(x = cACR_top, y = cACR_bottom)
head(combo@meta.data, 5)
tail(combo@meta.data, 5)

# Now set up DEG

Idents(combo) <- 'top4'
ACR.Top4vAll.response <- FindMarkers(combo, assay = "SCT", ident.1 = "Yes", ident.2 = "No", verbose = FALSE)

ACR.Top4vAll.response['genes'] <- row.names(ACR.Top4vAll.response) 
head(ACR.Top4vAll.response)
write.csv(ACR.Top4vAll.response, "ACR.Top4vAll.response.csv")

##########################################################################
#Volcano Plot of Top Clones versus All other clones
##########################################################################

lfc = 0.5
alpha = 0.05

tab2 = data.frame(genes = ACR.Top4vAll.response$genes, logFCDP = ACR.Top4vAll.response$avg_log2FC, DPpval = -log10(ACR.Top4vAll.response$p_val))
signGenes = (abs(tab2$logFCDP) > lfc & tab2$DPpval > -log10(alpha))
signGenesUp = (tab2$logFCDP) > lfc & tab2$DPpval > -log10(alpha)
signGenesDown = (tab2$logFCDP) < -lfc & tab2$DPpval > -log10(alpha)

#Setup genes to list on the graph
tab2$genelabels <- factor(tab2$genes, levels = c("GZMB", "GZMH", "CCL5", "HLA-DRB1","NKG7","ZNF683",
                                                 "XCL2", "XCL1", "PRF1","ITGAE", "KLRB1",
                                                 "LAG3","KLRC1", "FTH1", "S1PR1", "CCR7",
                                                 "FOS","IFITM1","CD8A","TNFRSF4","HOPX"))

head(tab2)

ggplot(data = tab2, mapping = aes(x=logFCDP, y=DPpval, label = genelabels))+
  geom_point() +
  xlim(-1,2) +
  ylim(0,180) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("Top 4 clones from each sample versus all others") +
  ylab("-log10 pvalue") +  
  geom_vline(xintercept = 0.5, linetype = "dotted") +
  geom_vline(xintercept = -0.5, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "red") +
  geom_hline(yintercept = -log10(alpha), linetype = "dotted") +
  geom_point(color = "blue", data = tab2[signGenes,]) +
  geom_label_repel(col = "black", size = 6, na.rm = TRUE, max.overlaps = 20, box.padding = unit(0.45, "lines"), hjust = 1)











