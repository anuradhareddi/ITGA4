## Gene Set Enrichment Analysis.
library(clusterProfiler)

GapmervsUT <- read.csv("DESeq2_GapmervsUT.csv", row.names = 1)
Gapmer_controlvsUT <- read.csv("DESeq2_Gapmer_controlvsUT.csv", row.names = 1)

MOEvsUT <- read.csv("DESeq2_MOEvsUT.csv", row.names = 1)
MOE_controlvsUT <- read.csv("DESeq2_MOE_controlvsUT.csv", row.names = 1)


nASO2vsUT <- read.csv("DESeq2_nASO2vsUT.csv", row.names = 1)
nASO2_controlvsUT <- read.csv("DESeq2_nASO2_controlvsUT.csv", row.names = 1)


PMOvsUT <- read.csv("DESeq2_PMOvsUT.csv", row.names = 1)
PMO_controlvsUT <- read.csv("DESeq2_PMO_controlvsUT.csv", row.names = 1)



#Extract log2FC, and sort genes based on the log2FC values

#Getting log2 fold change
Gapmer_gene_list <- GapmervsUT$log2FoldChange
Gapmer_control_gene_list <- Gapmer_controlvsUT$log2FoldChange

MOE_gene_list <- MOEvsUT$log2FoldChange
MOE_control_gene_list <- MOE_controlvsUT$log2FoldChange


nASO2_gene_list <- nASO2vsUT$log2FoldChange
nASO2_control_gene_list <- nASO2_controlvsUT$log2FoldChange


PMO_gene_list <- PMOvsUT$log2FoldChange
PMO_control_gene_list <- PMO_controlvsUT$log2FoldChange



# name the vector
names(Gapmer_gene_list) <- rownames(GapmervsUT)
names(Gapmer_control_gene_list) <- rownames(Gapmer_controlvsUT)

names(MOE_gene_list) <- rownames(MOEvsUT)
names(MOE_control_gene_list) <- rownames(MOE_controlvsUT)


names(nASO2_gene_list) <- rownames(nASO2vsUT)
names(nASO2_control_gene_list) <- rownames(nASO2_controlvsUT)


names(PMO_gene_list) <- rownames(PMOvsUT)
names(PMO_control_gene_list) <- rownames(PMO_controlvsUT)


# omit any NA values 
Gapmer_gene_list<-na.omit(Gapmer_gene_list)
Gapmer_control_gene_list<-na.omit(Gapmer_control_gene_list)

MOE_gene_list<-na.omit(MOE_gene_list)
MOE_control_gene_list<-na.omit(MOE_control_gene_list)


nASO2_gene_list<-na.omit(nASO2_gene_list)
nASO2_control_gene_list<-na.omit(nASO2_control_gene_list)



PMO_gene_list<-na.omit(PMO_gene_list)
PMO_control_gene_list<-na.omit(PMO_control_gene_list)




# sort the list in decreasing order
Gapmer_gene_list = sort(Gapmer_gene_list, decreasing = TRUE)
Gapmer_control_gene_list = sort(Gapmer_control_gene_list, decreasing = TRUE)

MOE_gene_list = sort(MOE_gene_list, decreasing = TRUE)
MOE_control_gene_list = sort(MOE_control_gene_list, decreasing = TRUE)


nASO2_gene_list = sort(nASO2_gene_list, decreasing = TRUE)
nASO2_control_gene_list = sort(nASO2_control_gene_list, decreasing = TRUE)


PMO_gene_list = sort(PMO_gene_list, decreasing = TRUE)
PMO_control_gene_list = sort(PMO_control_gene_list, decreasing = TRUE)


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

gse_Gapmer_control <- gseGO(geneList=Gapmer_control_gene_list, #no term enriched under specific pvalueCutoff
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


gse_MOE_control <- gseGO(geneList=MOE_control_gene_list, 
                 ont ="ALL", 
                 keyType = "SYMBOL", 
                 minGSSize = 3, 
                 maxGSSize = 800, 
                 pvalueCutoff = 0.05, 
                 verbose = TRUE, 
                 OrgDb = "org.Hs.eg.db") 


######################################################

gse_nASO2 <- gseGO(geneList=nASO2_gene_list, 
                 ont ="ALL", 
                 keyType = "SYMBOL", 
                 minGSSize = 3, 
                 maxGSSize = 800, 
                 pvalueCutoff = 0.05, 
                 verbose = TRUE, 
                 OrgDb = "org.Hs.eg.db")  


gse_nASO2_control <- gseGO(geneList=nASO2_control_gene_list, ##no term enriched under specific pvalueCutoff
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

gse_PMO_control <- gseGO(geneList=PMO_control_gene_list, 
                     ont ="ALL", 
                     keyType = "SYMBOL", 
                     minGSSize = 3, 
                     maxGSSize = 800, 
                     pvalueCutoff = 0.05, 
                     verbose = TRUE, 
                     OrgDb = "org.Hs.eg.db") 


######################
write.csv(gse_Gapmer@result, file = "GSE_Gapmer.csv")
write.csv(gse_MOE_control@result, file = "GSE_MOE_control.csv")
write.csv(gse_nASO2@result, file = "GSE_nASO2.csv")
write.csv(gse_PMO_control@result, file = "GSE_PMO_control.csv")


dotplot(gse_Gapmer, font.size = 10, label_format = 40, title = "" , split=".sign") + 
  facet_grid(.~.sign)+ scale_x_continuous(name = "GeneRatio", limits = c(0, 1), breaks = seq(0.25, 1, by = 0.25)) +
  theme(plot.margin = margin(5, 5, 5, 5, "mm"))

dotplot(gse_MOE_control, font.size = 10, label_format = 40, title = "" , split=".sign") + 
  facet_grid(.~.sign)+ scale_x_continuous(name = "GeneRatio", limits = c(0, 1), breaks = seq(0.25, 1, by = 0.25)) +
  theme(plot.margin = margin(5, 5, 5, 5, "mm"))


dotplot(gse_nASO2, font.size = 10, label_format = 40, title = "" , split=".sign") + 
  facet_grid(.~.sign)+ scale_x_continuous(name = "GeneRatio", limits = c(0, 1), breaks = seq(0.25, 1, by = 0.25)) +
  theme(plot.margin = margin(5, 5, 5, 5, "mm"))

dotplot(gse_PMO_control, font.size = 10, label_format = 40, title = "" , split=".sign") + 
  facet_grid(.~.sign)+ scale_x_continuous(name = "GeneRatio", limits = c(0, 1), breaks = seq(0.25, 1, by = 0.25)) +
  theme(plot.margin = margin(5, 5, 5, 5, "mm"))

###########################################################################
