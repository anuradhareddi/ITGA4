library(SpliceWiz)

MOE <- read.csv("splicewiz/DESeq2/deseq_MOE_sigs.csv")
PMO <- read.csv("splicewiz/DESeq2/deseq_PMO_sigs.csv")

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




