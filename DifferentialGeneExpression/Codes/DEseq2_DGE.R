setwd("D:/ITGA4_STAR_aligment/DifferentialGeneExpression")

# read counts data
Counts <- read.csv("D:/ITGA4_STAR_aligment/DifferentialGeneExpression/RawCounts.csv", row.names = 1)
sample_info <- read.csv("D:/ITGA4_STAR_aligment/DifferentialGeneExpression/sample_info.csv", row.names = 1)

#################################
# making sure the row names in colData matches to column names in counts_data
all(colnames(Counts) %in% rownames(sample_info))

# are they in the same order?
all(colnames(Counts) == rownames(sample_info))
##########################################
library(DESeq2)
library(ggplot2)
library(EnhancedVolcano)

############################
############################


dds <- DESeqDataSetFromMatrix(countData = Counts,
                              colData = sample_info,
                              design =  ~ Condition)

#################################################################

keep <- rowSums(counts(dds) >= 10) >=3
dds <- dds[keep,]


####################################################

deseq <- DESeq(dds)

##################
## set the factor level(here we are comparing with untreated)
relevel(deseq$Condition, ref = "Untreated")

######################################DIFFERENTAL EXPRESSION RESULTS

# Extract the results for each contrast
results_GapmervsUT <- results(deseq, contrast = c("Condition", "Gapmer", "Untreated"))
results_Gapmer_controlvsUT <- results(deseq, contrast = c("Condition", "GapmerControl", "Untreated"))

###
results_MOEvsUT <- results(deseq, contrast = c("Condition", "MOE", "Untreated"))
results_MOE_GTCvsUT <- results(deseq, contrast = c("Condition", "MOE_GTC", "Untreated"))



##
results_TMOvsUT <- results(deseq, contrast = c("Condition", "TMO", "Untreated"))
results_TMO_GTCvsUT <- results(deseq, contrast = c("Condition", "TMO_GTC", "Untreated"))

###
results_PMOvsUT <- results(deseq, contrast = c("Condition", "PMO", "Untreated"))
results_PMO_GTCvsUT <- results(deseq, contrast = c("Condition", "PMO_GTC", "Untreated"))



##

summary(results_GapmervsUT)
mcols(results_GapmervsUT)$description

########################################################################
############################
#Sort by p-value

results_GapmervsUT <- results_GapmervsUT[order(results_GapmervsUT$padj),]
results_Gapmer_controlvsUT <- results_Gapmer_controlvsUT[order(results_Gapmer_controlvsUT$padj),]

results_MOEvsUT <- results_MOEvsUT[order(results_MOEvsUT$padj),]
results_MOE_GTCvsUT <- results_MOE_GTCvsUT[order(results_MOE_GTCvsUT$padj),]



results_TMOvsUT <- results_TMOvsUT[order(results_TMOvsUT$padj),]
results_TMO_GTCvsUT <- results_TMO_GTCvsUT[order(results_TMO_GTCvsUT$padj),]

results_PMOvsUT <- results_PMOvsUT[order(results_PMOvsUT$padj),]
results_PMO_GTCvsUT <- results_PMO_GTCvsUT[order(results_PMO_GTCvsUT$padj),]

##############################################      VOLCANO 

results_GapmervsUT.df <- as.data.frame(results_GapmervsUT)
results_Gapmer_controlvsUT.df <- as.data.frame(results_Gapmer_controlvsUT)

results_MOEvsUT.df <- as.data.frame(results_MOEvsUT)
results_MOE_GTCvsUT.df <- as.data.frame(results_MOE_GTCvsUT)

results_TMOvsUT.df <- as.data.frame(results_TMOvsUT)
results_TMO_GTCvsUT.df <- as.data.frame(results_TMO_GTCvsUT)

results_PMOvsUT.df <- as.data.frame(results_PMOvsUT)
results_PMO_GTCvsUT.df <- as.data.frame(results_PMO_GTCvsUT)

###########
write.csv(results_GapmervsUT.df, file = "DESeq2_GapmervsUT.csv")
write.csv(results_Gapmer_controlvsUT.df, file = "DESeq2_Gapmer_controlvsUT.csv")

write.csv(results_MOEvsUT.df, file = "DESeq2_MOEvsUT.csv")
write.csv(results_MOE_GTCvsUT.df, file = "DESeq2_MOE_GTCvsUT.csv")


write.csv(results_TMOvsUT.df, file = "DESeq2_TMOvsUT.csv")
write.csv(results_TMO_GTCvsUT.df, file = "DESeq2_TMO_GTCvsUT.csv")


write.csv(results_PMOvsUT.df, file = "DESeq2_PMOvsUT.csv")
write.csv(results_PMO_GTCvsUT.df, file = "DESeq2_PMO_GTCvsUT.csv")


########


library(EnhancedVolcano)

GapmervsUT_volc <- EnhancedVolcano(
  results_GapmervsUT.df,
  lab = rownames(results_GapmervsUT.df),
  x = 'log2FoldChange',
  y = 'padj',
  xlim = c(-10, 10),
  ylim = c(0, 10),
  pCutoff = 0.05,
  FCcutoff = 0.2,
  col = c("black", "black", "black", "red"),
  colAlpha = 1,
  labSize = 3,
  legendPosition = "",
  legendLabSize = 10,
  caption = NULL,
  gridlines.major = TRUE,
  gridlines.minor = TRUE,
  title = "Gapmer",
  subtitle = "",
  border = 'full',
  borderWidth = 0.8,
  borderColour = 'black')


GapmervsUT_volc <- GapmervsUT_volc + 
  xlab("log2FoldChange") +
  ylab("-log10(padj)") +
  theme(plot.title = element_text(margin = margin(b = -18)))

print(GapmervsUT_volc)

####
Gapmer_controlvsUT_volc <- EnhancedVolcano(
  results_Gapmer_controlvsUT.df,
  lab = rownames(results_Gapmer_controlvsUT.df),
  x = 'log2FoldChange',
  y = 'padj',
  xlim = c(-10, 10),
  ylim = c(0, 20),
  pCutoff = 0.05,
  FCcutoff = 0.2,
  col = c("black", "black", "black", "red"),
  colAlpha = 1,
  labSize = 3,
  legendPosition = "",
  legendLabSize = 10,
  caption = NULL,
  gridlines.major = TRUE,
  gridlines.minor = TRUE,
  title = "Gapmer Control",
  subtitle = "",
  border = 'full',
  borderWidth = 0.8,
  borderColour = 'black')

Gapmer_controlvsUT_volc <- Gapmer_controlvsUT_volc + 
  xlab("log2FoldChange") +
  ylab("-log10(padj)") +
  theme(plot.title = element_text(margin = margin(b = -18)))

print(Gapmer_controlvsUT_volc)
###
MOEvsUT_volc <- EnhancedVolcano(
  results_MOEvsUT.df,
  lab = rownames(results_MOEvsUT.df),
  x = 'log2FoldChange',
  y = 'padj',
  xlim = c(-10, 10),
  ylim = c(0, 10),
  pCutoff = 0.05,
  FCcutoff = 0.2,
  col = c("black", "black", "black", "red"),
  colAlpha = 1,
  labSize = 3,
  legendPosition = "",
  legendLabSize = 10,
  caption = NULL,
  gridlines.major = TRUE,
  gridlines.minor = TRUE,
  title = "MOE",
  subtitle = "",
  border = 'full',
  borderWidth = 0.8,
  borderColour = 'black')


MOEvsUT_volc <- MOEvsUT_volc + 
  xlab("log2FoldChange") +
  ylab("-log10(padj)") +
  theme(plot.title = element_text(margin = margin(b = -18)))

print(MOEvsUT_volc)
##

MOE_GTCvsUT_volc <- EnhancedVolcano(
  results_MOE_GTCvsUT.df,
  lab = rownames(results_MOE_GTCvsUT.df),
  x = 'log2FoldChange',
  y = 'padj',
  xlim = c(-10, 10),
  ylim = c(0, 20),
  pCutoff = 0.05,
  FCcutoff = 0.2,
  col = c("black", "black", "black", "red"),
  colAlpha = 1,
  labSize = 3,
  legendPosition = "",
  legendLabSize = 10,
  caption = NULL,
  gridlines.major = TRUE,
  gridlines.minor = TRUE,
  title = "MOE GTC",
  subtitle = "",
  border = 'full',
  borderWidth = 0.8,
  borderColour = 'black')


MOE_GTCvsUT_volc <- MOE_GTCvsUT_volc + 
  xlab("log2FoldChange") +
  ylab("-log10(padj)") +
  theme(plot.title = element_text(margin = margin(b = -18)))

print(MOE_GTCvsUT_volc)
##


##########################################

TMOvsUT_volc <- EnhancedVolcano(
  results_TMOvsUT.df,
  lab = rownames(results_TMOvsUT.df),
  x = 'log2FoldChange',
  y = 'padj',
  xlim = c(-10, 10),
  ylim = c(0, 10),
  pCutoff = 0.05,
  FCcutoff = 0.2,
  col = c("black", "black", "black", "red"),
  colAlpha = 1,
  labSize = 3,
  legendPosition = "",
  legendLabSize = 10,
  caption = NULL,
  gridlines.major = TRUE,
  gridlines.minor = TRUE,
  title = "TMO",
  subtitle = "",
  border = 'full',
  borderWidth = 0.8,
  borderColour = 'black')


TMOvsUT_volc <- TMOvsUT_volc + 
  xlab("log2FoldChange") +
  ylab("-log10(padj)") +
  theme(plot.title = element_text(margin = margin(b = -18)))

print(TMOvsUT_volc)
##

TMO_GTCvsUT_volc <- EnhancedVolcano(
  results_TMO_GTCvsUT.df,
  lab = rownames(results_TMO_GTCvsUT.df),
  x = 'log2FoldChange',
  y = 'padj',
  xlim = c(-10, 10),
  ylim = c(0, 20),
  pCutoff = 0.05,
  FCcutoff = 0.2,
  col = c("black", "black", "black", "red"),
  colAlpha = 1,
  labSize = 3,
  legendPosition = "",
  legendLabSize = 10,
  caption = NULL,
  gridlines.major = TRUE,
  gridlines.minor = TRUE,
  title = "TMO GTC",
  subtitle = "",
  border = 'full',
  borderWidth = 0.8,
  borderColour = 'black')


TMO_GTCvsUT_volc <- TMO_GTCvsUT_volc + 
  xlab("log2FoldChange") +
  ylab("-log10(padj)") +
  theme(plot.title = element_text(margin = margin(b = -18)))

print(TMO_GTCvsUT_volc)
##

######################################################


PMOvsUT_volc <- EnhancedVolcano(
  results_PMOvsUT.df,
  lab = rownames(results_PMOvsUT.df),
  x = 'log2FoldChange',
  y = 'padj',
  xlim = c(-10, 10),
  ylim = c(0, 10),
  pCutoff = 0.05,
  FCcutoff = 0.2,
  col = c("black", "black", "black", "red"),
  colAlpha = 1,
  labSize = 3,
  legendPosition = "",
  legendLabSize = 10,
  caption = NULL,
  gridlines.major = TRUE,
  gridlines.minor = TRUE,
  title = "PMO",
  subtitle = "",
  border = 'full',
  borderWidth = 0.8,
  borderColour = 'black')


PMOvsUT_volc <- PMOvsUT_volc + 
  xlab("log2FoldChange") +
  ylab("-log10(padj)") +
  theme(plot.title = element_text(margin = margin(b = -18)))

print(PMOvsUT_volc)
##

PMO_GTCvsUT_volc <- EnhancedVolcano(
  results_PMO_GTCvsUT.df,
  lab = rownames(results_PMO_GTCvsUT.df),
  x = 'log2FoldChange',
  y = 'padj',
  xlim = c(-10, 10),
  ylim = c(0, 20),
  pCutoff = 0.05,
  FCcutoff = 0.2,
  col = c("black", "black", "black", "red"),
  colAlpha = 1,
  labSize = 3,
  legendPosition = "",
  legendLabSize = 10,
  caption = NULL,
  gridlines.major = TRUE,
  gridlines.minor = TRUE,
  title = "PMO GTC",
  subtitle = "",
  border = 'full',
  borderWidth = 0.8,
  borderColour = 'black')


PMO_GTCvsUT_volc <- PMO_GTCvsUT_volc + 
  xlab("log2FoldChange") +
  ylab("-log10(padj)") +
  theme(plot.title = element_text(margin = margin(b = -18)))

print(PMO_GTCvsUT_volc)
####

####

library(ggpubr)

ggarrange(GapmervsUT_volc, PMOvsUT_volc, TMOvsUT_volc, MOEvsUT_volc, ncol = 4, nrow = 1)
ggarrange(Gapmer_controlvsUT_volc, PMO_GTCvsUT_volc, TMO_GTCvsUT_volc, MOE_GTCvsUT_volc, ncol = 4, nrow = 1)
ggarrange(PMO_GTCvsUT_volc, TMO_GTCvsUT_volc, MOE_GTCvsUT_volc, ncol = 3, nrow = 1)

###########################VENN DIAGRAMS
#https://mkempenaar.github.io/gene_expression_analysis/chapter-5.html


## Data preparation
pval_threshold <- 0.05
logFC_threshold <- 0.2

Gapmer_DEGS <- rownames(results_GapmervsUT.df[which(results_GapmervsUT.df$padj < pval_threshold &
                                                abs(results_GapmervsUT.df$log2FoldChange) >= logFC_threshold), ])
Gapmer_DEGS

MOE_DEGS <- rownames(results_MOEvsUT.df[which(results_MOEvsUT.df$padj < pval_threshold &
                                                abs(results_MOEvsUT.df$log2FoldChange) >= logFC_threshold), ])
MOE_DEGS

TMO_DEGS <- rownames(results_TMOvsUT.df[which(results_TMOvsUT.df$padj < pval_threshold &
                                                abs(results_TMOvsUT.df$log2FoldChange) >= logFC_threshold), ])
TMO_DEGS

PMO_DEGS <- rownames(results_PMOvsUT.df[which(results_PMOvsUT.df$padj < pval_threshold &
                                                abs(results_PMOvsUT.df$log2FoldChange) >= logFC_threshold), ])
PMO_DEGS
########################################
Gapmer_control_DEGS <- rownames(results_Gapmer_controlvsUT.df[which(results_Gapmer_controlvsUT.df$padj < pval_threshold &
                                                      abs(results_Gapmer_controlvsUT.df$log2FoldChange) >= logFC_threshold), ])
Gapmer_control_DEGS

MOE_GTC_DEGS <- rownames(results_MOE_GTCvsUT.df[which(results_MOE_GTCvsUT.df$padj < pval_threshold &
                                                abs(results_MOE_GTCvsUT.df$log2FoldChange) >= logFC_threshold), ])
MOE_GTC_DEGS

TMO_GTC_DEGS <- rownames(results_TMO_GTCvsUT.df[which(results_TMO_GTCvsUT.df$padj < pval_threshold &
                                                abs(results_TMO_GTCvsUT.df$log2FoldChange) >= logFC_threshold), ])
TMO_GTC_DEGS

PMO_GTC_DEGS <- rownames(results_PMO_GTCvsUT.df[which(results_PMO_GTCvsUT.df$padj < pval_threshold &
                                                        abs(results_PMO_GTCvsUT.df$log2FoldChange) >= logFC_threshold), ])
PMO_GTC_DEGS


#################################

##Venn diagrams
library(VennDiagram)
# Define custom fill colors using predefined color names
custom_fill <- c("salmon1", "palegreen", "lightskyblue", "plum1")

# Create the Venn diagram with custom fill colors
venn.plot1 <- venn.diagram(
  x = list(Gapmer = Gapmer_DEGS, PMO = PMO_DEGS, TMO = TMO_DEGS, MOE = MOE_DEGS),
  category.names = c("Gapmer", "PMO", "TMO", "MOE"),
  filename = NULL,
  fill = custom_fill

)

# Plot the Venn diagram
grid.draw(venn.plot1)

#####################
venn.plot2 <- venn.diagram(
  x = list(GapmerControl = Gapmer_control_DEGS, PMO_GTC = PMO_GTC_DEGS, TMO_GTC = TMO_GTC_DEGS, MOE_GTC = MOE_GTC_DEGS),
  category.names = c("GapmerControl", "PMO_GTC", "TMO_GTC", "MOE_GTC"),
  filename = NULL,
  fill = custom_fill
  
)

# Plot the Venn diagram
grid.draw(venn.plot2)

#######################

#############################################################################
############################## UPSET PLOT
library(ComplexHeatmap)
library(UpSetR)

# Define your list of gene sets
x <- list(
  Gapmer_ITGA4 = Gapmer_DEGS,
  PMO_ITGA4 = PMO_DEGS,
  TMO_ITGA4 = TMO_DEGS,
  MOE_ITGA4 = MOE_DEGS,
  Gapmer_Control = Gapmer_control_DEGS,
  PMO_GTC = PMO_GTC_DEGS,
  TMO_GTC = TMO_GTC_DEGS,
  MOE_GTC = MOE_GTC_DEGS
)

# Create the combination matrix
m <- make_comb_mat(x)

ss = set_size(m)
ss
cs = comb_size(m)
ht = UpSet(m, 
           set_order = order(ss),
           comb_order = order(comb_degree(m), -cs),
           top_annotation = HeatmapAnnotation(
             "Intersection Size" = anno_barplot(cs, 
                                                ylim = c(0, max(cs)*1.1),
                                                border = FALSE, 
                                                gp = gpar(fill = "black"), 
                                                height = unit(4, "cm")
             ), 
             annotation_name_side = "left", 
             annotation_name_rot = 90),
           left_annotation = rowAnnotation(
             "Set Size" = anno_barplot(-ss, 
                                       baseline = 0,
                                       axis_param = list(
                                         at = c(0, -1000, -2000, -3000, -4000),
                                         labels = c(0, 1000, 2000, 3000, 4000),
                                         labels_rot = 0),
                                       border = FALSE, 
                                       gp = gpar(fill = "black"), 
                                       width = unit(4, "cm")
             ),
             set_name = anno_text(set_name(m), 
                                  location = 0.5, 
                                  just = "center",
                                  width = max_text_width(set_name(m)) + unit(4, "mm"))
           ), 
           right_annotation = NULL,
           show_row_names = FALSE)
ht = draw(ht)
od = column_order(ht)
decorate_annotation("Intersection Size", {
  grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"), 
            default.units = "native", just = c("center", "bottom"), 
            gp = gpar(fontsize = 8, col = "black"), rot = 0)
})


