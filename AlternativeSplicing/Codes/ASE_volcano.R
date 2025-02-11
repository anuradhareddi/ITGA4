Gapmer <- read.csv("D:/ITGA4_STAR_aligment/splicewiz/DESeq2/deseq_Gapmer.csv")
GapmerControl <- read.csv("D:/ITGA4_STAR_aligment/splicewiz/DESeq2/deseq_GapmerControl.csv")

MOE <- read.csv("D:/ITGA4_STAR_aligment/splicewiz/DESeq2/deseq_MOE.csv")
MOE_GTC <- read.csv("D:/ITGA4_STAR_aligment/splicewiz/DESeq2/deseq_MOE_GTC.csv")

PMO <- read.csv("D:/ITGA4_STAR_aligment/splicewiz/DESeq2/deseq_PMO.csv")
PMO_GTC <- read.csv("D:/ITGA4_STAR_aligment/splicewiz/DESeq2/deseq_PMO_GTC.csv")

TMO <- read.csv("D:/ITGA4_STAR_aligment/splicewiz/DESeq2/deseq_TMO.csv")
TMO_GTC <- read.csv("D:/ITGA4_STAR_aligment/splicewiz/DESeq2/deseq_TMO_GTC.csv")


########################################## VISUALIZATION

### VOLCANO PLOTS

############################ VOLCANO
library(ggrepel)

extract_gene_symbol <- function(event_name) {
  if (grepl(":", event_name)) {
    gene_symbol <- strsplit(event_name, ":")[[1]][2]
    gene_symbol <- strsplit(gene_symbol, "-")[[1]][1]
  } else if (grepl("/", event_name)) {
    gene_symbol <- strsplit(event_name, "/")[[1]][1]
  } else {
    gene_symbol <- event_name
  }
  return(gene_symbol)
}
######
Gapmer$Gene_Symbol <- sapply(Gapmer$EventName,extract_gene_symbol)
GapmerControl$Gene_Symbol <- sapply(GapmerControl$EventName,extract_gene_symbol)

MOE$Gene_Symbol <- sapply(MOE$EventName,extract_gene_symbol)
MOE_GTC$Gene_Symbol <- sapply(MOE_GTC$EventName,extract_gene_symbol)

PMO$Gene_Symbol <- sapply(PMO$EventName,extract_gene_symbol)
PMO_GTC$Gene_Symbol <- sapply(PMO_GTC$EventName,extract_gene_symbol)

TMO$Gene_Symbol <- sapply(TMO$EventName,extract_gene_symbol)
TMO_GTC$Gene_Symbol <- sapply(TMO_GTC$EventName,extract_gene_symbol)

########################################################################
library(ggplot2)
library(dplyr)


filtered_Gapmer <- Gapmer %>%
  filter(EventType != "SE")


ggplot(filtered_Gapmer, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(data = subset(filtered_Gapmer, -log10(padj) > 1.3 & abs(log2FoldChange) > 0.2),
             aes(color = "red"), size = 3) +
  geom_point(data = subset(filtered_Gapmer, !(-log10(padj) > 1.3 & abs(log2FoldChange) > 0.2)),
             size = 3) +
  facet_wrap(vars(EventType), ncol = 6, strip.position = "bottom") +
  labs(title = "Gapmer", x = "Log2-fold change", y = "-log10(padj)") +
  geom_text_repel(data = subset(filtered_Gapmer, -log10(padj) > 1.3 & abs(log2FoldChange) > 0.2),
                  aes(label = Gene_Symbol), size = 2, show.legend = FALSE, max.overlaps = 30) +
  scale_color_manual(values = c("red")) +  # Define color for significant points
  scale_y_continuous(limits = c(0, 5)) +  # Setting y-axis limits
  scale_x_continuous(limits = c(-20, 20)) +  # Setting y-axis limits
  theme(axis.text = element_text(size = 8),
        axis.text.x = element_text(size = 8),
        axis.title = element_text(size = 8, face = "bold"),
        legend.position = "none",
        strip.text = element_text(size = 8))
##########
filtered_MOE <- MOE %>%
  filter(EventType != "SE")


ggplot(filtered_MOE, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(data = subset(filtered_MOE, -log10(padj) > 1.3 & abs(log2FoldChange) > 0.2),
             aes(color = "red"), size = 3) +
  geom_point(data = subset(filtered_MOE, !(-log10(padj) > 1.3 & abs(log2FoldChange) > 0.2)),
             size = 3) +
  facet_wrap(vars(EventType), ncol = 6, strip.position = "bottom") +
  labs(title = "MOE", x = "Log2-fold change", y = "-log10(padj)") +
  geom_text_repel(data = subset(filtered_MOE, -log10(padj) > 1.3 & abs(log2FoldChange) > 0.2),
                  aes(label = Gene_Symbol), size = 2, show.legend = FALSE, max.overlaps = 30) +
  scale_color_manual(values = c("red")) +  # Define color for significant points
  scale_y_continuous(limits = c(0, 5)) +  # Setting y-axis limits
  scale_x_continuous(limits = c(-20, 20)) +  # Setting y-axis limits
  theme(axis.text = element_text(size = 8),
        axis.text.x = element_text(size = 8),
        axis.title = element_text(size = 8, face = "bold"),
        legend.position = "none",
        strip.text = element_text(size = 8))
##########
filtered_PMO <- PMO %>%
  filter(EventType != "SE")


ggplot(filtered_PMO, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(data = subset(filtered_PMO, -log10(padj) > 1.3 & abs(log2FoldChange) > 0.2),
             aes(color = "red"), size = 3) +
  geom_point(data = subset(filtered_PMO, !(-log10(padj) > 1.3 & abs(log2FoldChange) > 0.2)),
             size = 3) +
  facet_wrap(vars(EventType), ncol = 6, strip.position = "bottom") +
  labs(title = "PMO", x = "Log2-fold change", y = "-log10(padj)") +
  geom_text_repel(data = subset(filtered_PMO, -log10(padj) > 1.3 & abs(log2FoldChange) > 0.2),
                  aes(label = Gene_Symbol), size = 2, show.legend = FALSE, max.overlaps = 30) +
  scale_color_manual(values = c("red")) +  # Define color for significant points
  scale_y_continuous(limits = c(0, 5)) +  # Setting y-axis limits
  scale_x_continuous(limits = c(-20, 20)) +  # Setting y-axis limits
  theme(axis.text = element_text(size = 8),
        axis.text.x = element_text(size = 8),
        axis.title = element_text(size = 8, face = "bold"),
        legend.position = "none",
        strip.text = element_text(size = 8))
##########
filtered_TMO <- TMO %>%
  filter(EventType != "SE")


ggplot(filtered_TMO, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(data = subset(filtered_TMO, -log10(padj) > 1.3 & abs(log2FoldChange) > 0.2),
             aes(color = "red"), size = 3) +
  geom_point(data = subset(filtered_TMO, !(-log10(padj) > 1.3 & abs(log2FoldChange) > 0.2)),
             size = 3) +
  facet_wrap(vars(EventType), ncol = 6, strip.position = "bottom") +
  labs(title = "TMO", x = "Log2-fold change", y = "-log10(padj)") +
  geom_text_repel(data = subset(filtered_TMO, -log10(padj) > 1.3 & abs(log2FoldChange) > 0.2),
                  aes(label = Gene_Symbol), size = 2, show.legend = FALSE, max.overlaps = 30) +
  scale_color_manual(values = c("red")) +  # Define color for significant points
  scale_y_continuous(limits = c(0, 5)) +  # Setting y-axis limits
  scale_x_continuous(limits = c(-20, 20)) +  # Setting y-axis limits
  theme(axis.text = element_text(size = 8),
        axis.text.x = element_text(size = 8),
        axis.title = element_text(size = 8, face = "bold"),
        legend.position = "none",
        strip.text = element_text(size = 8))
##########

ggplot(GapmerControl, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(data = subset(GapmerControl, -log10(padj) > 1.3 & abs(log2FoldChange) > 0.2),
             aes(color = "red"), size = 3) +
  geom_point(data = subset(GapmerControl, !(-log10(padj) > 1.3 & abs(log2FoldChange) > 0.2)),
             size = 3) +
  facet_wrap(vars(EventType), ncol = 7, strip.position = "bottom") +
  labs(title = "GapmerControl", x = "Log2-fold change", y = "-log10(padj)") +
  geom_text_repel(data = subset(GapmerControl, -log10(padj) > 1.3 & abs(log2FoldChange) > 0.2),
                  aes(label = Gene_Symbol), size = 2, show.legend = FALSE, max.overlaps = 30) +
  scale_color_manual(values = c("red")) +  # Define color for significant points
  scale_y_continuous(limits = c(0, 5)) +  # Setting y-axis limits
  scale_x_continuous(limits = c(-20, 20)) +  # Setting y-axis limits
  theme(axis.text = element_text(size = 8),
        axis.text.x = element_text(size = 8),
        axis.title = element_text(size = 8, face = "bold"),
        legend.position = "none",
        strip.text = element_text(size = 8))
###

ggplot(MOE_GTC, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(data = subset(MOE_GTC, -log10(padj) > 1.3 & abs(log2FoldChange) > 0.2),
             aes(color = "red"), size = 3) +
  geom_point(data = subset(MOE_GTC, !(-log10(padj) > 1.3 & abs(log2FoldChange) > 0.2)),
             size = 3) +
  facet_wrap(vars(EventType), ncol = 7, strip.position = "bottom") +
  labs(title = "MOE_GTC", x = "Log2-fold change", y = "-log10(padj)") +
  geom_text_repel(data = subset(MOE_GTC, -log10(padj) > 1.3 & abs(log2FoldChange) > 0.2),
                  aes(label = Gene_Symbol), size = 2, show.legend = FALSE, max.overlaps = 30) +
  scale_color_manual(values = c("red")) +  # Define color for significant points
  scale_y_continuous(limits = c(0, 5)) +  # Setting y-axis limits
  scale_x_continuous(limits = c(-20, 20)) +  # Setting y-axis limits
  theme(axis.text = element_text(size = 8),
        axis.text.x = element_text(size = 8),
        axis.title = element_text(size = 8, face = "bold"),
        legend.position = "none",
        strip.text = element_text(size = 8))

##

ggplot(PMO_GTC, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(data = subset(PMO_GTC, -log10(padj) > 1.3 & abs(log2FoldChange) > 0.2),
             aes(color = "red"), size = 3) +
  geom_point(data = subset(PMO_GTC, !(-log10(padj) > 1.3 & abs(log2FoldChange) > 0.2)),
             size = 3) +
  facet_wrap(vars(EventType), ncol = 7, strip.position = "bottom") +
  labs(title = "PMO_GTC", x = "Log2-fold change", y = "-log10(padj)") +
  geom_text_repel(data = subset(PMO_GTC, -log10(padj) > 1.3 & abs(log2FoldChange) > 0.2),
                  aes(label = Gene_Symbol), size = 2, show.legend = FALSE, max.overlaps = 30) +
  scale_color_manual(values = c("red")) +  # Define color for significant points
  scale_y_continuous(limits = c(0, 5)) +  # Setting y-axis limits
  scale_x_continuous(limits = c(-20, 20)) +  # Setting y-axis limits
  theme(axis.text = element_text(size = 8),
        axis.text.x = element_text(size = 8),
        axis.title = element_text(size = 8, face = "bold"),
        legend.position = "none",
        strip.text = element_text(size = 8))
###

##
ggplot(TMO_GTC, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(data = subset(TMO_GTC, -log10(padj) > 1.3 & abs(log2FoldChange) > 0.2),
             aes(color = "red"), size = 3) +
  geom_point(data = subset(TMO_GTC, !(-log10(padj) > 1.3 & abs(log2FoldChange) > 0.2)),
             size = 3) +
  facet_wrap(vars(EventType), ncol = 7, strip.position = "bottom") +
  labs(title = "TMO_GTC", x = "Log2-fold change", y = "-log10(padj)") +
  geom_text_repel(data = subset(TMO_GTC, -log10(padj) > 1.3 & abs(log2FoldChange) > 0.2),
                  aes(label = Gene_Symbol), size = 2, show.legend = FALSE, max.overlaps = 30) +
  scale_color_manual(values = c("red")) +  # Define color for significant points
  scale_y_continuous(limits = c(0, 5)) +  # Setting y-axis limits
  scale_x_continuous(limits = c(-20, 20)) +  # Setting y-axis limits
  theme(axis.text = element_text(size = 8),
        axis.text.x = element_text(size = 8),
        axis.title = element_text(size = 8, face = "bold"),
        legend.position = "none",
        strip.text = element_text(size = 8))
###



