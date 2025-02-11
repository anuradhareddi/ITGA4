setwd("D:/ITGA4_STAR_aligment/splicewiz")

pb_path <- file.path("D:/ITGA4_STAR_aligment/splicewiz/pbOutput")

expr <- findSpliceWizOutput(pb_path)
expr

Gapmer_expr <- expr[c(1:3,34:36),]
MOE_expr <- expr[c(7:9,34:36),]
PMO_expr <- expr[c(16:18,34:36),]
TMO_expr <- expr[c(25:27,34:36),]

NxtSE_Gapmer <- file.path("D:/ITGA4_STAR_aligment/splicewiz/NxtSE/NxtSE_Gapmer")
NxtSE_MOE <- file.path("D:/ITGA4_STAR_aligment/splicewiz/NxtSE/NxtSE_MOE")
NxtSE_PMO <- file.path("D:/ITGA4_STAR_aligment/splicewiz/NxtSE/NxtSE_PMO")
NxtSE_TMO <- file.path("D:/ITGA4_STAR_aligment/splicewiz/NxtSE/NxtSE_TMO")

se_Gapmer <- makeSE(NxtSE_Gapmer)
se_MOE <- makeSE(NxtSE_MOE)
se_PMO <- makeSE(NxtSE_PMO)
se_TMO <- makeSE(NxtSE_TMO)

colData(se_Gapmer)$condition <- c("Gapmer", "Gapmer", "Gapmer", "UT", "UT", "UT")
colData(se_MOE)$condition <- c("MOE", "MOE", "MOE", "UT", "UT", "UT")
colData(se_PMO)$condition <- c("PMO", "PMO", "PMO", "UT", "UT", "UT")
colData(se_TMO)$condition <- c("TMO", "TMO", "TMO", "UT", "UT", "UT")

se_Gapmer.filtered <- se_Gapmer[applyFilters(se_Gapmer),]
se_MOE.filtered <- se_MOE[applyFilters(se_MOE),]
se_PMO.filtered <- se_PMO[applyFilters(se_PMO),]
se_TMO.filtered <- se_TMO[applyFilters(se_TMO),]

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

######################## Sashimi plot

plotCoverage(
  se_Gapmer.filtered,
  Event = "SE:ITGA4-201-exon4;ITGA4-209-int2",
  Gene = "ITGA4",
  seqname = "chr2",
  start = 181475067,
  end = 181478756,
  strand = "*",
  zoom_factor = 0.05,
  bases_flanking = 100,
  tracks = c("Gapmer", "UT"),
  condition = "condition",
  ribbon_mode = "none",
  selected_transcripts = c("ITGA4-201", "ITGA4-209"),
  reverseGenomeCoords = FALSE,
  plotJunctions = TRUE,
  junctionThreshold = 0.2,
  plot_key_isoforms = TRUE,
  condense_tracks = FALSE,
  stack_tracks = FALSE,
  t_test = FALSE,
  norm_event = "SE:ITGA4-201-exon4;ITGA4-209-int2",
  usePlotly = FALSE
)

plotCoverage(
  se_MOE.filtered,
  Event = "SE:ITGA4-201-exon4;ITGA4-209-int2",
  Gene = "ITGA4",
  seqname = "chr2",
  start = 181475067,
  end = 181478756,
  strand = "*",
  zoom_factor = 0.05,
  bases_flanking = 100,
  tracks = c("MOE", "UT"),
  condition = "condition",
  ribbon_mode = "none",
  selected_transcripts = c("ITGA4-201", "ITGA4-209"),
  reverseGenomeCoords = FALSE,
  plotJunctions = TRUE,
  junctionThreshold = 0.1,
  plot_key_isoforms = TRUE,
  condense_tracks = FALSE,
  stack_tracks = FALSE,
  t_test = FALSE,
  norm_event = "SE:ITGA4-201-exon4;ITGA4-209-int2",
  usePlotly = FALSE
)

plotCoverage(
  se_TMO.filtered,
  Event = "SE:ITGA4-201-exon4;ITGA4-209-int2",
  Gene = "ITGA4",
  seqname = "chr2",
  start = 181475067,
  end = 181478756,
  strand = "*",
  zoom_factor = 0.05,
  bases_flanking = 100,
  tracks = c("TMO", "UT"),
  condition = "condition",
  ribbon_mode = "none",
  selected_transcripts = c("ITGA4-201", "ITGA4-209"),
  reverseGenomeCoords = FALSE,
  plotJunctions = TRUE,
  junctionThreshold = 0.2,
  plot_key_isoforms = TRUE,
  condense_tracks = FALSE,
  stack_tracks = FALSE,
  t_test = FALSE,
  norm_event = "SE:ITGA4-201-exon4;ITGA4-209-int2",
  usePlotly = FALSE
)

plotCoverage(
  se_PMO.filtered,
  Event = "SE:ITGA4-201-exon4;ITGA4-209-int2",
  Gene = "ITGA4",
  seqname = "chr2",
  start = 181475067,
  end = 181478756,
  strand = "*",
  zoom_factor = 0.05,
  bases_flanking = 100,
  tracks = c("PMO", "UT"),
  condition = "condition",
  ribbon_mode = "none",
  selected_transcripts = c("ITGA4-201", "ITGA4-209"),
  reverseGenomeCoords = FALSE,
  plotJunctions = TRUE,
  junctionThreshold = 0.2,
  plot_key_isoforms = TRUE,
  condense_tracks = FALSE,
  stack_tracks = FALSE,
  t_test = FALSE,
  norm_event = "SE:ITGA4-201-exon4;ITGA4-209-int2",
  usePlotly = FALSE
)
