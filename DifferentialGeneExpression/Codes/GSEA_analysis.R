## Gene Set Enrichment Analysis.
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(ggplot2)
library(DOSE)
GapmervsUT <- read.csv("DESeq2_GapmervsUT.csv", row.names = 1)
Gapmer_controlvsUT <- read.csv("DESeq2_Gapmer_controlvsUT.csv", row.names = 1)

MOEvsUT <- read.csv("DESeq2_MOEvsUT.csv", row.names = 1)
MOE_GTCvsUT <- read.csv("DESeq2_MOE_GTCvsUT.csv", row.names = 1)


TMOvsUT <- read.csv("DESeq2_TMOvsUT.csv", row.names = 1)
TMO_GTCvsUT <- read.csv("DESeq2_TMO_GTCvsUT.csv", row.names = 1)


PMOvsUT <- read.csv("DESeq2_PMOvsUT.csv", row.names = 1)
PMO_GTCvsUT <- read.csv("DESeq2_PMO_GTCvsUT.csv", row.names = 1)



#Extract log2FC, and sort genes based on the log2FC values

#Getting log2 fold change
Gapmer_gene_list <- GapmervsUT$log2FoldChange
GapmerControl_gene_list <- Gapmer_controlvsUT$log2FoldChange

MOE_gene_list <- MOEvsUT$log2FoldChange
MOE_GTC_gene_list <- MOE_GTCvsUT$log2FoldChange


TMO_gene_list <- TMOvsUT$log2FoldChange
TMO_GTC_gene_list <- TMO_GTCvsUT$log2FoldChange


PMO_gene_list <- PMOvsUT$log2FoldChange
PMO_GTC_gene_list <- PMO_GTCvsUT$log2FoldChange



# name the vector
names(Gapmer_gene_list) <- rownames(GapmervsUT)
names(GapmerControl_gene_list) <- rownames(Gapmer_controlvsUT)

names(MOE_gene_list) <- rownames(MOEvsUT)
names(MOE_GTC_gene_list) <- rownames(MOE_GTCvsUT)


names(TMO_gene_list) <- rownames(TMOvsUT)
names(TMO_GTC_gene_list) <- rownames(TMO_GTCvsUT)
names(TMO_SCR_gene_list) <- rownames(TMO_SCRvsUT)

names(PMO_gene_list) <- rownames(PMOvsUT)
names(PMO_GTC_gene_list) <- rownames(PMO_GTCvsUT)


# omit any NA values 
Gapmer_gene_list<-na.omit(Gapmer_gene_list)
GapmerControl_gene_list<-na.omit(GapmerControl_gene_list)

MOE_gene_list<-na.omit(MOE_gene_list)
MOE_GTC_gene_list<-na.omit(MOE_GTC_gene_list)


TMO_gene_list<-na.omit(TMO_gene_list)
TMO_GTC_gene_list<-na.omit(TMO_GTC_gene_list)



PMO_gene_list<-na.omit(PMO_gene_list)
PMO_GTC_gene_list<-na.omit(PMO_GTC_gene_list)




# sort the list in decreasing order
Gapmer_gene_list = sort(Gapmer_gene_list, decreasing = TRUE)
GapmerControl_gene_list = sort(GapmerControl_gene_list, decreasing = TRUE)

MOE_gene_list = sort(MOE_gene_list, decreasing = TRUE)
MOE_GTC_gene_list = sort(MOE_GTC_gene_list, decreasing = TRUE)


TMO_gene_list = sort(TMO_gene_list, decreasing = TRUE)
TMO_GTC_gene_list = sort(TMO_GTC_gene_list, decreasing = TRUE)


PMO_gene_list = sort(PMO_gene_list, decreasing = TRUE)
PMO_GTC_gene_list = sort(PMO_GTC_gene_list, decreasing = TRUE)


##############################################################
set.seed(123)

gse_Gapmer <- gseGO(geneList=Gapmer_gene_list, 
                    ont ="ALL", 
                    keyType = "SYMBOL", 
                    minGSSize = 3, 
                    maxGSSize = 800, 
                    pvalueCutoff = 0.05, 
                    verbose = TRUE, 
                    OrgDb = "org.Hs.eg.db")  

gse_GapmerControl <- gseGO(geneList=GapmerControl_gene_list, #no term enriched under specific pvalueCutoff
                    ont ="ALL", 
                    keyType = "SYMBOL", 
                    minGSSize = 3, 
                    maxGSSize = 800, 
                    pvalueCutoff = 0.05, 
                    verbose = TRUE, 
                    OrgDb = "org.Hs.eg.db") 
########################################################

gse_MOE <- gseGO(geneList=MOE_gene_list, ##no term enriched under specific pvalueCutoff
                    ont ="ALL", 
                    keyType = "SYMBOL", 
                    minGSSize = 3, 
                    maxGSSize = 800, 
                    pvalueCutoff = 0.05, 
                    verbose = TRUE, 
                    OrgDb = "org.Hs.eg.db") 


gse_MOE_GTC <- gseGO(geneList=MOE_GTC_gene_list, 
                 ont ="ALL", 
                 keyType = "SYMBOL", 
                 minGSSize = 3, 
                 maxGSSize = 800, 
                 pvalueCutoff = 0.05, 
                 verbose = TRUE, 
                 OrgDb = "org.Hs.eg.db") 


######################################################

gse_TMO <- gseGO(geneList=TMO_gene_list, 
                 ont ="ALL", 
                 keyType = "SYMBOL", 
                 minGSSize = 3, 
                 maxGSSize = 800, 
                 pvalueCutoff = 0.05, 
                 verbose = TRUE, 
                 OrgDb = "org.Hs.eg.db")  


gse_TMO_GTC <- gseGO(geneList=TMO_GTC_gene_list, ##no term enriched under specific pvalueCutoff
                     ont ="ALL", 
                     keyType = "SYMBOL", 
                     minGSSize = 3, 
                     maxGSSize = 800, 
                     pvalueCutoff = 0.05, 
                     verbose = TRUE, 
                     OrgDb = "org.Hs.eg.db") 


################################
gse_PMO <- gseGO(geneList=PMO_gene_list, #no term enriched under specific pvalueCutoff
                 ont ="ALL", 
                 keyType = "SYMBOL", 
                 minGSSize = 3, 
                 maxGSSize = 800, 
                 pvalueCutoff = 0.05, 
                 verbose = TRUE, 
                 OrgDb = "org.Hs.eg.db") 

gse_PMO_GTC <- gseGO(geneList=PMO_GTC_gene_list, 
                     ont ="ALL", 
                     keyType = "SYMBOL", 
                     minGSSize = 3, 
                     maxGSSize = 800, 
                     pvalueCutoff = 0.05, 
                     verbose = TRUE, 
                     OrgDb = "org.Hs.eg.db") 


######################
write.csv(gse_Gapmer@result, file = "GSE_Gapmer.csv")
write.csv(gse_MOE_GTC@result, file = "GSE_MOE_GTC.csv")
write.csv(gse_TMO@result, file = "GSE_TMO.csv")
write.csv(gse_PMO_GTC@result, file = "GSE_PMO_GTC.csv")


dotplot(gse_Gapmer, font.size = 10, label_format = 40, title = "" , split=".sign") + 
  facet_grid(.~.sign)+ scale_x_continuous(name = "GeneRatio", limits = c(0, 1), breaks = seq(0.25, 1, by = 0.25)) +
  theme(plot.margin = margin(5, 5, 5, 5, "mm"))

dotplot(gse_MOE_GTC, font.size = 10, label_format = 40, title = "" , split=".sign") + 
  facet_grid(.~.sign)+ scale_x_continuous(name = "GeneRatio", limits = c(0, 1), breaks = seq(0.25, 1, by = 0.25)) +
  theme(plot.margin = margin(5, 5, 5, 5, "mm"))


dotplot(gse_TMO, font.size = 10, label_format = 40, title = "" , split=".sign") + 
  facet_grid(.~.sign)+ scale_x_continuous(name = "GeneRatio", limits = c(0, 1), breaks = seq(0.25, 1, by = 0.25)) +
  theme(plot.margin = margin(5, 5, 5, 5, "mm"))

dotplot(gse_PMO_GTC, font.size = 10, label_format = 40, title = "" , split=".sign") + 
  facet_grid(.~.sign)+ scale_x_continuous(name = "GeneRatio", limits = c(0, 1), breaks = seq(0.25, 1, by = 0.25)) +
  theme(plot.margin = margin(5, 5, 5, 5, "mm"))

###########################################################################
