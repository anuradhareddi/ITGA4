
# read counts data
Counts <- read.csv("DifferentialGeneExpression/RawCounts.csv", row.names = 1)
sample_info <- read.csv("DifferentialGeneExpression/sample_info.csv", row.names = 1)

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
results_Gapmer_controlvsUT <- results(deseq, contrast = c("Condition", "Gapmer_control", "Untreated"))

###
results_MOEvsUT <- results(deseq, contrast = c("Condition", "MOE", "Untreated"))
results_MOE_controlvsUT <- results(deseq, contrast = c("Condition", "MOE_control", "Untreated"))



##
results_nASO2vsUT <- results(deseq, contrast = c("Condition", "nASO2", "Untreated"))
results_nASO2_controlvsUT <- results(deseq, contrast = c("Condition", "nASO2_control", "Untreated"))

###
results_PMOvsUT <- results(deseq, contrast = c("Condition", "PMO", "Untreated"))
results_PMO_controlvsUT <- results(deseq, contrast = c("Condition", "PMO_control", "Untreated"))



##

summary(results_GapmervsUT)
mcols(results_GapmervsUT)$description

########################################################################
############################
#Sort by p-value

results_GapmervsUT <- results_GapmervsUT[order(results_GapmervsUT$padj),]
results_Gapmer_controlvsUT <- results_Gapmer_controlvsUT[order(results_Gapmer_controlvsUT$padj),]

results_MOEvsUT <- results_MOEvsUT[order(results_MOEvsUT$padj),]
results_MOE_controlvsUT <- results_MOE_controlvsUT[order(results_MOE_controlvsUT$padj),]



results_nASO2vsUT <- results_nASO2vsUT[order(results_nASO2vsUT$padj),]
results_nASO2_controlvsUT <- results_nASO2_controlvsUT[order(results_nASO2_controlvsUT$padj),]

results_PMOvsUT <- results_PMOvsUT[order(results_PMOvsUT$padj),]
results_PMO_controlvsUT <- results_PMO_controlvsUT[order(results_PMO_controlvsUT$padj),]

##############################################      VOLCANO 

results_GapmervsUT.df <- as.data.frame(results_GapmervsUT)
results_Gapmer_controlvsUT.df <- as.data.frame(results_Gapmer_controlvsUT)

results_MOEvsUT.df <- as.data.frame(results_MOEvsUT)
results_MOE_controlvsUT.df <- as.data.frame(results_MOE_controlvsUT)

results_nASO2vsUT.df <- as.data.frame(results_nASO2vsUT)
results_nASO2_controlvsUT.df <- as.data.frame(results_nASO2_controlvsUT)

results_PMOvsUT.df <- as.data.frame(results_PMOvsUT)
results_PMO_controlvsUT.df <- as.data.frame(results_PMO_controlvsUT)

###########
write.csv(results_GapmervsUT.df, file = "DESeq2_GapmervsUT.csv")
write.csv(results_Gapmer_controlvsUT.df, file = "DESeq2_Gapmer_controlvsUT.csv")

write.csv(results_MOEvsUT.df, file = "DESeq2_MOEvsUT.csv")
write.csv(results_MOE_controlvsUT.df, file = "DESeq2_MOE_controlvsUT.csv")


write.csv(results_nASO2vsUT.df, file = "DESeq2_nASO2vsUT.csv")
write.csv(results_nASO2_controlvsUT.df, file = "DESeq2_nASO2_controlvsUT.csv")


write.csv(results_PMOvsUT.df, file = "DESeq2_PMOvsUT.csv")
write.csv(results_PMO_controlvsUT.df, file = "DESeq2_PMO_controlvsUT.csv")


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



##########################################

nASO2vsUT_volc <- EnhancedVolcano(
  results_nASO2vsUT.df,
  lab = rownames(results_nASO2vsUT.df),
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
  title = "nASO2",
  subtitle = "",
  border = 'full',
  borderWidth = 0.8,
  borderColour = 'black')


nASO2vsUT_volc <- nASO2vsUT_volc + 
  xlab("log2FoldChange") +
  ylab("-log10(padj)") +
  theme(plot.title = element_text(margin = margin(b = -18)))

print(nASO2vsUT_volc)
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


###########################VENN DIAGRAMS

pval_threshold <- 0.05
logFC_threshold <- 0.2

Gapmer_DEGS <- rownames(results_GapmervsUT.df[which(results_GapmervsUT.df$padj < pval_threshold &
                                                abs(results_GapmervsUT.df$log2FoldChange) >= logFC_threshold), ])
Gapmer_DEGS

MOE_DEGS <- rownames(results_MOEvsUT.df[which(results_MOEvsUT.df$padj < pval_threshold &
                                                abs(results_MOEvsUT.df$log2FoldChange) >= logFC_threshold), ])
MOE_DEGS

nASO2_DEGS <- rownames(results_nASO2vsUT.df[which(results_nASO2vsUT.df$padj < pval_threshold &
                                                abs(results_nASO2vsUT.df$log2FoldChange) >= logFC_threshold), ])
nASO2_DEGS

PMO_DEGS <- rownames(results_PMOvsUT.df[which(results_PMOvsUT.df$padj < pval_threshold &
                                                abs(results_PMOvsUT.df$log2FoldChange) >= logFC_threshold), ])
PMO_DEGS
########################################
Gapmer_control_DEGS <- rownames(results_Gapmer_controlvsUT.df[which(results_Gapmer_controlvsUT.df$padj < pval_threshold &
                                                      abs(results_Gapmer_controlvsUT.df$log2FoldChange) >= logFC_threshold), ])
Gapmer_control_DEGS

MOE_control_DEGS <- rownames(results_MOE_controlvsUT.df[which(results_MOE_controlvsUT.df$padj < pval_threshold &
                                                abs(results_MOE_controlvsUT.df$log2FoldChange) >= logFC_threshold), ])
MOE_control_DEGS

nASO2_control_DEGS <- rownames(results_nASO2_controlvsUT.df[which(results_nASO2_controlvsUT.df$padj < pval_threshold &
                                                abs(results_nASO2_controlvsUT.df$log2FoldChange) >= logFC_threshold), ])
nASO2_control_DEGS

PMO_control_DEGS <- rownames(results_PMO_controlvsUT.df[which(results_PMO_controlvsUT.df$padj < pval_threshold &
                                                        abs(results_PMO_controlvsUT.df$log2FoldChange) >= logFC_threshold), ])
PMO_control_DEGS


#################################

##Venn diagrams
library(VennDiagram)
# Define custom fill colors using predefined color names
custom_fill <- c("salmon1", "palegreen", "lightskyblue", "plum1")

# Create the Venn diagram with custom fill colors
venn.plot1 <- venn.diagram(
  x = list(Gapmer = Gapmer_DEGS, PMO = PMO_DEGS, nASO2 = nASO2_DEGS, MOE = MOE_DEGS),
  category.names = c("Gapmer", "PMO", "nASO2", "MOE"),
  filename = NULL,
  fill = custom_fill

)

# Plot the Venn diagram
grid.draw(venn.plot1)



#######################

#############################################################################
############################## UPSET PLOT
library(UpSetR)

# Define your list of gene sets
x <- list(
  Gapmer_ITGA4 = Gapmer_DEGS,
  PMO_ITGA4 = PMO_DEGS,
  nASO2_ITGA4 = nASO2_DEGS,
  MOE_ITGA4 = MOE_DEGS,
  Gapmer_Control = Gapmer_control_DEGS,
  PMO_control = PMO_control_DEGS,
  nASO2_control = nASO2_control_DEGS,
  MOE_control = MOE_control_DEGS
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


