# Load libraries
library(recount)
library(DESeq2)
library(tidyverse)
library(ggpubr)
library(EnhancedVolcano)
library(EnsDb.Hsapiens.v86)
library(msigdbr)
library(clusterProfiler)
library(ggplot2)
library(pheatmap)

###############################
### Download and Load Data ###
###############################
exp <-  'SRP043043'

url <- download_study(exp)

# Load data
load(file.path(exp, 'rse_gene.Rdata'))

# Examine the read count matrix
cts <- assay(rse_gene)

# Examine the sample metadata to see if it's ready for DESeq2 analysis
cd <- as.data.frame(colData(rse_gene))

# Add condition status, modify 'title' to indicate Replicates
rse_gene$cell = gsub('^[[:alpha:]]*[0-9]*\\_', '', rse_gene$title)
# rse_gene$condition <- factor(ifelse(grepl("siScram", colData(rse_gene)$title), "Untreated", "Treated"))
# 
# # change ensembl id to remove vers
rownames(rse_gene) <- gsub('\\.[0-9]*$', '', rownames(rse_gene))
rse_gene$condition <- c(rep("siSCR", 3), rep("siZNF217", 3))

#####################################
### Principal Component Analysis ###
####################################
# Create DEseq dataset: design formula variable of interest is treatment status
dds <- DESeqDataSet(rse_gene, ~condition)
dds$condition <- relevel(dds$condition, ref = "siSCR")

# Create PCA to inspect batch effect
# Perform regularized logarithm transformation (rlog) on the data; VST (>30 samples)
vsd<- vst(dds)

# Create PCA plot
pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color = condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() + labs(color='Treatment')

####################
### DEG Analysis ###
####################
# Analyse
dds <- DESeq(dds)

# Pre-filter low counts
#keep <- rowSums(counts(dds)) >= 10
#dds <- dds[keep,]

# Obtain results, contrast to indicate which is ref.
res <- results(dds, contrast = c("condition", "siZNF217", "siSCR"))


# LFC shrink
normRes <- lfcShrink(dds, res = res, coef= "condition_siZNF217_vs_siSCR", type= "apeglm")

# MA Plot
plotMA(normRes)


############################
### Volcano Plot of DEG ###
###########################
# Make a df out of the res
resdf <- normRes %>%
  data.frame() %>%
  rownames_to_column(var="ensembl_id") 

# Creating a tibble of ENSEMBL (GENEID), ENTREZID, SYMBOL, GENENAME using EnsDb.Hsapiens.v86 db
anno <- AnnotationDbi::select(EnsDb.Hsapiens.v86,resdf$ensembl_id, 
                                        columns=c("GENEID", "SYMBOL"), 
                                        keytype="GENEID")

# left_join the res and the gene information
finRes <- resdf %>%
  left_join(anno, by = c("ensembl_id" ="GENEID"))

# Plot Volcano
EnhancedVolcano( finRes, lab = finRes$SYMBOL, 
                 x = 'log2FoldChange', y = 'padj',
                 xlim = c(-8, 8), title = 'Treated vs Untreated',
                 pCutoff = 1e-100, FCcutoff = 1.25, 
                 pointSize = 2.0, 
                 labSize = 3.0,
                 border = "full", borderWidth = 1.5, borderColour = "black", 
                 gridlines.major = FALSE, gridlines.minor = FALSE)


#########################
### Extract sig DEGs ###
########################
# Significant genes
sigRes <- finRes[which(finRes$padj < 0.05 & abs(finRes$log2FoldChange) >= 1.5 & finRes$baseMean >= 20), ]

# Heatmap of all DEGs
mat <- assay(vsd)  
idx <- sigRes$ensembl_id 
DEgenes <- mat[idx,]

annotation <- as.data.frame(colData(vsd)[, c("cell", "condition")])

pheatmap(DEgenes, scale = "row", show_rownames = FALSE, clustering_distance_rows = "correlation", annotation_col = annotation, main="Differentially Expressed Genes")

# Extract Top 20 Over and Under-expressed Genes
# Up
top20_up <- sigRes %>%
  dplyr::filter(log2FoldChange > 1.5 & padj < 0.05) %>%
  arrange(desc(log2FoldChange)) %>%
  head(20)

# Down
top20_down <- sigRes %>%
  dplyr::filter(log2FoldChange < 1.5 & padj < 0.05) %>%
  arrange(log2FoldChange) %>%
 head(20)

# Plot Heatmaps
# Up HM
top20_up_hm <- mat[top20_up$ensembl_id,]
rownames(top20_up_hm) <- top20_up$SYMBOL

pheatmap(top20_up_hm, scale = "row", clustering_distance_rows = "correlation", annotation_col = annotation, main="Top Over Expressed Genes")

# Under HM
top20_down_hm <- mat[top20_down$ensembl_id,]
rownames(top20_down_hm) <- top20_down$SYMBOL

pheatmap(top20_down_hm, scale = "row", clustering_distance_rows = "correlation", annotation_col = annotation, main="Top 20 Under Expressed Genes")

###############
### Enrichr ###
###############

# Extract over expressed genes 
over_expressed_genes <- finRes %>%
  dplyr::filter(padj < .05 & log2FoldChange > 1.5) %>%
  pull(SYMBOL)

OEgenes <- finRes %>%
  dplyr::filter(padj < .05 & log2FoldChange > 1.5) %>%
  write_csv("OEgenes.csv")

# Create gene set for ref
gene_sets <- msigdbr(species = "Homo sapiens", category = "C5")
gene_sets <- gene_sets %>%
  dplyr::select(gs_name, gene_symbol)

# Use enrichr to look at OE genes (compared to GO)
egmt <- enricher(gene = over_expressed_genes,
                 TERM2GENE = gene_sets)
edf <- as.data.frame(egmt)


#############
### GSEA ###
############
# Adding a score for GSEA
gsea_df <- finRes%>%
  arrange(padj) %>%
  mutate(gsea_metric = -log10(padj) * sign(log2FoldChange))

# Deal with inf (swap up the inf with the highest pvalue)
gsea_df <- gsea_df %>%
  mutate(padj = case_when(padj == 0 ~ .Machine$double.xmin,
                          TRUE ~ padj)) %>%
  mutate(gsea_metric = -log10(padj) * sign(log2FoldChange)) 

# Remove NAs and order by GSEA
gsea_df <- gsea_df  %>%
  dplyr::filter(!is.na(gsea_metric)) %>%
  arrange(desc(gsea_metric)) # needed to run GSEA later

# GSEA value histogram
hist(gsea_df$gsea_metric, breaks = 100)

# Get the ranked GSEA vector
ranks <- gsea_df  %>%
  dplyr::select(SYMBOL, gsea_metric) %>%
  distinct(SYMBOL, .keep_all = TRUE) %>%
  deframe()

# Run GSEA
gseares <- GSEA(geneList = ranks, 
                TERM2GENE = gene_sets,
                pvalueCutoff = 1)

gsearesdf <- as.data.frame(gseares)

View(gsearesdf)

# Plot GSEA results
# Top 5 Over
# Make GSEA plot for top and bottom results
top_pathways <- gsearesdf %>%
  top_n(n = 5, wt = NES) %>%
  pull(ID)

#  Make gseaplot for each and return as list
top_pathway_plots <- lapply(top_pathways, function(pathway) {
  gseaplot(gseares, geneSetID = pathway, title = pathway)
})

# Arrange with labels as a multi-panel plot
top_pathway_plot <- ggarrange(plotlist = top_pathway_plots,
                              ncol = 2, nrow = 3, labels = "AUTO")


#  Save it
ggsave(top_pathway_plot, filename = "top_GSEA_up.png",
       height = 11, width = 18)

# Top 5 Under
bottom_pathways <- gsearesdf %>%
  top_n(n = 5, wt = -NES) %>%
  pull(ID)

bottom_pathway_plots <- lapply(bottom_pathways, function(pathway) {
  gseaplot(gseares, geneSetID = pathway, title = pathway)
})

bottom_pathway_plot <- ggarrange(plotlist = bottom_pathway_plots,
                                 ncol = 2, nrow = 3, labels = "AUTO")

ggsave(bottom_pathway_plot, filename = "top_GSEA_down.png",
       height = 11, width = 18)
