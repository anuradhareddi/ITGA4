setwd("D:/ITGA4_STAR_aligment/splicewiz")
library(SpliceWiz)
#spliceWiz()
##############################################################################

ref_path <- file.path("C:/Users/20150013/Downloads/RNAseq_Reference/Gencode/Reference_gencode_v44")

######################################################################################

pb_path <- file.path("D:/ITGA4_STAR_aligment/splicewiz/pbOutput")

##############################################################################


#find samples

expr <- findSpliceWizOutput(pb_path)
expr

Gapmer_expr <- expr[c(1:3,34:36),]
GapmerControl_expr <- expr[c(4:6,34:36),]

MOE_expr <- expr[c(7:9,34:36),]
MOE_GTC_expr <- expr[c(10:12,34:36),]

PMO_expr <- expr[c(16:18,34:36),]
PMO_GTC_expr <- expr[c(19:21,34:36),]

TMO_expr <- expr[c(25:27,34:36),]
TMO_GTC_expr <- expr[c(28:30,34:36),]


###################################################################
#Collate the experiment
NxtSE_Gapmer <- file.path("D:/ITGA4_STAR_aligment/splicewiz/NxtSE/NxtSE_Gapmer")
NxtSE_GapmerControl <- file.path("D:/ITGA4_STAR_aligment/splicewiz/NxtSE/NxtSE_GapmerControl")

NxtSE_MOE <- file.path("D:/ITGA4_STAR_aligment/splicewiz/NxtSE/NxtSE_MOE")
NxtSE_MOE_GTC <- file.path("D:/ITGA4_STAR_aligment/splicewiz/NxtSE/NxtSE_MOE_GTC")

NxtSE_PMO <- file.path("D:/ITGA4_STAR_aligment/splicewiz/NxtSE/NxtSE_PMO")
NxtSE_PMO_GTC <- file.path("D:/ITGA4_STAR_aligment/splicewiz/NxtSE/NxtSE_PMO_GTC")

NxtSE_TMO <- file.path("D:/ITGA4_STAR_aligment/splicewiz/NxtSE/NxtSE_TMO")
NxtSE_TMO_GTC <- file.path("D:/ITGA4_STAR_aligment/splicewiz/NxtSE/NxtSE_TMO_GTC")



####
collateData(
  Experiment = Gapmer_expr,
  reference_path = ref_path,
  output_path = NxtSE_Gapmer,
  n_threads = 8
 )

collateData(
  Experiment = GapmerControl_expr,
  reference_path = ref_path,
  output_path = NxtSE_GapmerControl,
  n_threads = 8
  
)
#######################################


collateData(
  Experiment = MOE_expr,
  reference_path = ref_path,
  output_path = NxtSE_MOE,
  n_threads = 8
)

collateData(
  Experiment = MOE_GTC_expr,
  reference_path = ref_path,
  output_path = NxtSE_MOE_GTC,
  n_threads = 8
)


#######################################################

collateData(
  Experiment = PMO_expr,
  reference_path = ref_path,
  output_path = NxtSE_PMO,
  n_threads = 8
)

collateData(
  Experiment = PMO_GTC_expr,
  reference_path = ref_path,
  output_path = NxtSE_PMO_GTC,
  n_threads = 8
)



#######################################################
collateData(
  Experiment = TMO_expr,
  reference_path = ref_path,
  output_path = NxtSE_TMO,
  n_threads = 8
)

collateData(
  Experiment = TMO_GTC_expr,
  reference_path = ref_path,
  output_path = NxtSE_TMO_GTC,
  n_threads = 8
)


#######################################################


se_Gapmer <- makeSE(NxtSE_Gapmer)
se_GapmerControl <- makeSE(NxtSE_GapmerControl)

se_MOE <- makeSE(NxtSE_MOE)
se_MOE_GTC <- makeSE(NxtSE_MOE_GTC)


se_PMO <- makeSE(NxtSE_PMO)
se_PMO_GTC <- makeSE(NxtSE_PMO_GTC)


se_TMO <- makeSE(NxtSE_TMO)
se_TMO_GTC <- makeSE(NxtSE_TMO_GTC)



#Step1: Assigning annotations to samples

colData(se_Gapmer)$condition <- c("Gapmer", "Gapmer", "Gapmer", "UT", "UT", "UT")
colData(se_GapmerControl)$condition <- c("GapmerControl", "GapmerControl", "GapmerControl", "UT", "UT", "UT")

colData(se_MOE)$condition <- c("MOE", "MOE", "MOE", "UT", "UT", "UT")
colData(se_MOE_GTC)$condition <- c("MOE_GTC", "MOE_GTC", "MOE_GTC", "UT", "UT", "UT")


colData(se_PMO)$condition <- c("PMO", "PMO", "PMO", "UT", "UT", "UT")
colData(se_PMO_GTC)$condition <- c("PMO_GTC", "PMO_GTC", "PMO_GTC", "UT", "UT", "UT")


colData(se_TMO)$condition <- c("TMO", "TMO", "TMO", "UT", "UT", "UT")
colData(se_TMO_GTC)$condition <- c("TMO_GTC", "TMO_GTC", "TMO_GTC", "UT", "UT", "UT")



#Step2: Filtering with default filters
se_Gapmer.filtered <- se_Gapmer[applyFilters(se_Gapmer),]
se_GapmerControl.filtered <- se_GapmerControl[applyFilters(se_GapmerControl),]

se_MOE.filtered <- se_MOE[applyFilters(se_MOE),]
se_MOE_GTC.filtered <- se_MOE_GTC[applyFilters(se_MOE_GTC),]

se_PMO.filtered <- se_PMO[applyFilters(se_PMO),]
se_PMO_GTC.filtered <- se_PMO_GTC[applyFilters(se_PMO_GTC),]


se_TMO.filtered <- se_TMO[applyFilters(se_TMO),]
se_TMO_GTC.filtered <- se_TMO_GTC[applyFilters(se_TMO_GTC),]


#Differential splicing
library(DESeq2)
deseq_Gapmer <- ASE_DESeq(
  se = se_Gapmer.filtered,
  test_factor = "condition",
  test_nom = "Gapmer",
  test_denom = "UT",
  IRmode = "annotated"
)

deseq_GapmerControl <- ASE_DESeq(
  se = se_GapmerControl.filtered,
  test_factor = "condition",
  test_nom = "GapmerControl",
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

deseq_MOE_GTC <- ASE_DESeq(
  se = se_MOE_GTC.filtered,
  test_factor = "condition",
  test_nom = "MOE_GTC",
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

deseq_PMO_GTC <- ASE_DESeq(
  se = se_PMO_GTC.filtered,
  test_factor = "condition",
  test_nom = "PMO_GTC",
  test_denom = "UT",
  IRmode = "annotated"
)


##################

deseq_TMO <- ASE_DESeq(
  se = se_TMO.filtered,
  test_factor = "condition",
  test_nom = "TMO",
  test_denom = "UT",
  IRmode = "annotated"
)

deseq_TMO_GTC <- ASE_DESeq(
  se = se_TMO_GTC.filtered,
  test_factor = "condition",
  test_nom = "TMO_GTC",
  test_denom = "UT",
  IRmode = "annotated"
)


write.csv(deseq_Gapmer, file = "D:/ITGA4_STAR_aligment/splicewiz/DESeq2/deseq_Gapmer.csv")
write.csv(deseq_GapmerControl, file = "D:/ITGA4_STAR_aligment/splicewiz/DESeq2/deseq_GapmerControl.csv")

write.csv(deseq_MOE, file = "D:/ITGA4_STAR_aligment/splicewiz/DESeq2/deseq_MOE.csv")
write.csv(deseq_MOE_GTC, file = "D:/ITGA4_STAR_aligment/splicewiz/DESeq2/deseq_MOE_GTC.csv")

write.csv(deseq_PMO, file = "D:/ITGA4_STAR_aligment/splicewiz/DESeq2/deseq_PMO.csv")
write.csv(deseq_PMO_GTC, file = "D:/ITGA4_STAR_aligment/splicewiz/DESeq2/deseq_PMO_GTC.csv")

write.csv(deseq_TMO, file = "D:/ITGA4_STAR_aligment/splicewiz/DESeq2/deseq_TMO.csv")
write.csv(deseq_TMO_GTC, file = "D:/ITGA4_STAR_aligment/splicewiz/DESeq2/deseq_TMO_GTC.csv")

