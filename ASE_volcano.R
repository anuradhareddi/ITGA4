Gapmer <- read.csv("splicewiz/DESeq2/deseq_Gapmer.csv")

MOE <- read.csv("splicewiz/DESeq2/deseq_MOE.csv")

PMO <- read.csv("splicewiz/DESeq2/deseq_PMO.csv")

nASO2 <- read.csv("splicewiz/DESeq2/deseq_nASO2.csv")



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

MOE$Gene_Symbol <- sapply(MOE$EventName,extract_gene_symbol)

PMO$Gene_Symbol <- sapply(PMO$EventName,extract_gene_symbol)

nASO2$Gene_Symbol <- sapply(nASO2$EventName,extract_gene_symbol)

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
filtered_nASO2 <- nASO2 %>%
  filter(EventType != "SE")


ggplot(filtered_nASO2, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(data = subset(filtered_nASO2, -log10(padj) > 1.3 & abs(log2FoldChange) > 0.2),
             aes(color = "red"), size = 3) +
  geom_point(data = subset(filtered_nASO2, !(-log10(padj) > 1.3 & abs(log2FoldChange) > 0.2)),
             size = 3) +
  facet_wrap(vars(EventType), ncol = 6, strip.position = "bottom") +
  labs(title = "nASO2", x = "Log2-fold change", y = "-log10(padj)") +
  geom_text_repel(data = subset(filtered_nASO2, -log10(padj) > 1.3 & abs(log2FoldChange) > 0.2),
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

