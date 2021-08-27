#Figure 2 - clonotype overlap over time for P1, P3, and P8 

#Load the TCR data for P8
etcr <- read.csv("path-to-P1/Early/TCR/outs/filtered_contig_annotations.csv")
atcr <- read.csv("path-to-P8/ACR/TCR/outs/filtered_contig_annotations.csv")
ttcr <- read.csv("path-to-P8/Treated/TCR/outs/filtered_contig_annotations.csv")

tcr_list <- list(etcr, atcr, ttcr)

tail(tcr_list[[1]])

combined <- combineTCR(tcr_list, 
                       samples = c("1_Early","2_ACR", "3_Treated"), 
                       ID = c("P8","P8", "P8"), 
                       cells ="T-AB")

head(combined[[1]])

quantContig(combined, cloneCall="gene+nt", scale = T)

quantContig_output <- quantContig(combined, cloneCall="gene+nt", 
                                  scale = T, exportTable = T)
quantContig_output

abundanceContig(combined, cloneCall = "gene", scale = F)

compareClonotypes(combined, numbers = NULL, samples = c("1_Early_P8", "2_ACR_P8", "3_Treated_P8"), 
                  cloneCall="aa", graph = "alluvial")

compareClonotypes(combined, numbers = NULL, samples = c("2_ACR_P8", "3_Treated_P8"), 
                  cloneCall="aa", graph = "alluvial")

Fig2a <- compareClonotypes(combined, numbers = NULL, samples = c("1_Early_T342", "2_ACR_T342", "3_Treated_T342"), 
                           cloneCall="aa", graph = "alluvial")

ggsave("Fig2a.pdf", width = 200, height = 50, units = "cm", scale = 8,limitsize = FALSE)

#Repeating for P1
atcrP1 <- read.csv("path-to-P1/ACR/TCR/outs/filtered_contig_annotations.csv")
ttcrP1 <- read.csv("path-to-P1/Treated/TCR/outs/filtered_contig_annotations.csv")


tcrP1_list <- list(atcrP1, ttcrP1)

tail(tcrP1_list[[1]])

combinedP1 <- combineTCR(tcrP1_list, 
                         samples = c("1_ACR", "2_Treated"), 
                         ID = c("P1", "P1"), 
                         cells ="T-AB")

head(combinedP1[[1]])

Figure2b <- compareClonotypes(combinedP1, numbers = NULL, samples = c("1_ACR_P1", "2_Treated_P1"), 
                              cloneCall="aa", graph = "alluvial")

ggsave("Fig2b.pdf", width = 150, height = 50, units = "cm", scale = 6,limitsize = FALSE)

#And again for P3
atcrP3 <- read.csv("path-to-P3/ACR/TCR/outs/filtered_contig_annotations.csv")
ttcrP3 <- read.csv("path-to-P3/Treated/TCR/outs/filtered_contig_annotations.csv")


tcrP3_list <- list(atcrP3, ttcrP3)

tail(tcrP3_list[[1]])

combinedP3 <- combineTCR(tcrP3_list, 
                          samples = c("1_ACR", "2_Treated"), 
                          ID = c("P3", "P3"), 
                          cells ="T-AB")

head(combinedP3[[1]])

Figure2c <- compareClonotypes(combinedP3, numbers = NULL, samples = c("1_ACR_P3", "2_Treated_P3"), 
                              cloneCall="aa", graph = "alluvial")

ggsave("Fig2c.pdf", width = 210, height = 50, units = "cm", scale = 7,limitsize = FALSE)

#######################
# Supplemental figure 1 - clonal diversity estimation
########################################################
# Clonal Tracking of P8, P1, and P3
########################################################

tcr_list <- list(atcr, atcrP1, atcrP3, ttcr, ttcrP1, ttcrP3)

tail(tcr_list[[1]])



combined <- combineTCR(tcr_list, 
                       samples = c("ACR_P8","ACR_P1", "ACR_P3","Tx_P8","Tx_P1","Tx_P3"), 
                       ID = c("ACR","ACR", "ACR", "Tx","Tx", "Tx"), 
                       cells ="T-AB")

SupplFig1 <- clonalDiversity(combined, cloneCall = "gene", 
                             n.boots = 100)

