setwd("D:/ITGA4_STAR_aligment/splicewiz")
library(SpliceWiz)
pb_path <- file.path("D:/ITGA4_STAR_aligment/splicewiz/pbOutput")

expr <- findSpliceWizOutput(pb_path)
expr

MOE_expr <- expr[c(7:9,34:36),]
PMO_expr <- expr[c(16:18,34:36),]
TMO_expr <- expr[c(25:27,34:36),]


NxtSE_MOE <- file.path("D:/ITGA4_STAR_aligment/splicewiz/NxtSE/NxtSE_MOE")
NxtSE_PMO <- file.path("D:/ITGA4_STAR_aligment/splicewiz/NxtSE/NxtSE_PMO")
NxtSE_TMO <- file.path("D:/ITGA4_STAR_aligment/splicewiz/NxtSE/NxtSE_TMO")


se_MOE <- makeSE(NxtSE_MOE)
se_PMO <- makeSE(NxtSE_PMO)
se_TMO <- makeSE(NxtSE_TMO)


colData(se_MOE)$condition <- c("MOE", "MOE", "MOE", "UT", "UT", "UT")
colData(se_PMO)$condition <- c("PMO", "PMO", "PMO", "UT", "UT", "UT")
colData(se_TMO)$condition <- c("TMO", "TMO", "TMO", "UT", "UT", "UT")


se_MOE.filtered <- se_MOE[applyFilters(se_MOE),]
se_PMO.filtered <- se_PMO[applyFilters(se_PMO),]
se_TMO.filtered <- se_TMO[applyFilters(se_TMO),]

library(DESeq2)

deseq_MOE <- ASE_DESeq(
  se = se_MOE.filtered,
  test_factor = "condition",
  test_nom = "MOE",
  test_denom = "UT",
  IRmode = "annotated"
)

deseq_PMO <- ASE_DESeq(
  se = se_PMO.filtered,
  test_factor = "condition",
  test_nom = "PMO",
  test_denom = "UT",
  IRmode = "annotated"
)

deseq_TMO <- ASE_DESeq(
  se = se_TMO.filtered,
  test_factor = "condition",
  test_nom = "TMO",
  test_denom = "UT",
  IRmode = "annotated"
)
####################################

MOE <- read.csv("D:/ITGA4_STAR_aligment/splicewiz/DESeq2/deseq_MOE_sigs.csv")
PMO <- read.csv("D:/ITGA4_STAR_aligment/splicewiz/DESeq2/deseq_PMO_sigs.csv")
TMO <- read.csv("D:/ITGA4_STAR_aligment/splicewiz/DESeq2/deseq_TMO_sigs.csv")
####################



go_MOE <- goASE(
  enrichedEventNames = MOE$EventName,
  universeEventNames = NULL,
  se = se_MOE.filtered
)
plotGO(go_MOE, filter_n_terms = 10)

###############################

go_PMO <- goASE(
  enrichedEventNames = PMO$EventName,
  universeEventNames = NULL,
  se = se_PMO.filtered
)
plotGO(go_PMO, filter_n_terms = 10)
#######################

go_TMO <- goASE(
  enrichedEventNames = TMO$EventName,
  universeEventNames = NULL,
  se = se_TMO.filtered
)

plotGO(go_TMO, filter_n_terms = 10)
############################################

library(data.table)

# Convert the list column to a character column by concatenating list elements
go_MOE$overlapGenes <- sapply(go_MOE$overlapGenes, function(x) paste(x, collapse = ";"))
go_PMO$overlapGenes <- sapply(go_PMO$overlapGenes, function(x) paste(x, collapse = ";"))
go_TMO$overlapGenes <- sapply(go_TMO$overlapGenes, function(x) paste(x, collapse = ";"))

# Save the modified data table to a CSV file
fwrite(go_MOE, file = "D:/ITGA4_STAR_aligment/splicewiz/DESeq2/go_MOE.csv")
fwrite(go_PMO, file = "D:/ITGA4_STAR_aligment/splicewiz/DESeq2/go_PMO.csv")
fwrite(go_TMO, file = "D:/ITGA4_STAR_aligment/splicewiz/DESeq2/go_TMO.csv")

