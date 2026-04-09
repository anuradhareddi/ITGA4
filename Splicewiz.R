library(SpliceWiz)
#spliceWiz()
##############################################################################

ref_path <- file.path("Reference_gencode_v44")

######################################################################################

pb_path <- file.path("splicewiz/pbOutput")

##############################################################################


#find samples

expr <- findSpliceWizOutput(pb_path)
expr

Gapmer_expr <- expr[c(1:3,34:36),]

MOE_expr <- expr[c(7:9,34:36),]

PMO_expr <- expr[c(16:18,34:36),]

nASO2_expr <- expr[c(25:27,34:36),]


###################################################################
#Collate the experiment
NxtSE_Gapmer <- file.path("splicewiz/NxtSE/NxtSE_Gapmer")

NxtSE_MOE <- file.path("splicewiz/NxtSE/NxtSE_MOE")

NxtSE_PMO <- file.path("splicewiz/NxtSE/NxtSE_PMO")

NxtSE_nASO2 <- file.path("splicewiz/NxtSE/NxtSE_nASO2")



####
collateData(
  Experiment = Gapmer_expr,
  reference_path = ref_path,
  output_path = NxtSE_Gapmer,
  n_threads = 8
 )


#######################################


collateData(
  Experiment = MOE_expr,
  reference_path = ref_path,
  output_path = NxtSE_MOE,
  n_threads = 8
)


#######################################################

collateData(
  Experiment = PMO_expr,
  reference_path = ref_path,
  output_path = NxtSE_PMO,
  n_threads = 8
)



#######################################################
collateData(
  Experiment = nASO2_expr,
  reference_path = ref_path,
  output_path = NxtSE_nASO2,
  n_threads = 8
)


#######################################################


se_Gapmer <- makeSE(NxtSE_Gapmer)

se_MOE <- makeSE(NxtSE_MOE)

se_PMO <- makeSE(NxtSE_PMO)

se_nASO2 <- makeSE(NxtSE_nASO2)




# Assigning annotations to samples

colData(se_Gapmer)$condition <- c("Gapmer", "Gapmer", "Gapmer", "UT", "UT", "UT")

colData(se_MOE)$condition <- c("MOE", "MOE", "MOE", "UT", "UT", "UT")

colData(se_PMO)$condition <- c("PMO", "PMO", "PMO", "UT", "UT", "UT")

colData(se_nASO2)$condition <- c("nASO2", "nASO2", "nASO2", "UT", "UT", "UT")


# Filtering with default filters
se_Gapmer.filtered <- se_Gapmer[applyFilters(se_Gapmer),]
se_MOE.filtered <- se_MOE[applyFilters(se_MOE),]
se_PMO.filtered <- se_PMO[applyFilters(se_PMO),]
se_nASO2.filtered <- se_nASO2[applyFilters(se_nASO2),]


#Differential splicing
library(DESeq2)
deseq_Gapmer <- ASE_DESeq(
  se = se_Gapmer.filtered,
  test_factor = "condition",
  test_nom = "Gapmer",
  test_denom = "UT",
  IRmode = "annotated"
)


deseq_MOE <- ASE_DESeq(
  se = se_MOE.filtered,
  test_factor = "condition",
  test_nom = "MOE",
  test_denom = "UT",
  IRmode = "annotated"
)


#############
deseq_PMO <- ASE_DESeq(
  se = se_PMO.filtered,
  test_factor = "condition",
  test_nom = "PMO",
  test_denom = "UT",
  IRmode = "annotated"
)


##################

deseq_nASO2 <- ASE_DESeq(
  se = se_nASO2.filtered,
  test_factor = "condition",
  test_nom = "nASO2",
  test_denom = "UT",
  IRmode = "annotated"
)


write.csv(deseq_Gapmer, file = "deseq_Gapmer.csv")
write.csv(deseq_MOE, file = "deseq_MOE.csv")
write.csv(deseq_PMO, file = "deseq_PMO.csv")
write.csv(deseq_nASO2, file = "deseq_nASO2.csv")


