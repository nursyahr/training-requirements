# Load libraries
suppressMessages( library(recount))
suppressMessages( library(DESeq2))
suppressMessages( library(dplyr))
suppressMessages( library(tibble) )
suppressMessages( library(ggplot2))
suppressMessages( library(plotly))
suppressMessages( library(org.Hs.eg.db) )
suppressMessages( library( AnnotationDbi ) )
suppressMessages( library( apeglm) )
suppressMessages(library(genefilter))
suppressMessages(library(pheatmap))
suppressMessages(library(PoiClaClu))
suppressMessages(library(RColorBrewer))
suppressMessages(library(EnhancedVolcano))
suppressMessages(library(EnsDb.Hsapiens.v86))
################################
#### DOWNLOAD AND LOAD DATA ####
################################

# AML3 cells were treated with Azacytidine and compared against untreated cells O
# Overall design: We used RNA-Seq to detail the global programme of gene expression 
# in human untreated and Azacytidine treated AML3 cells

exp <- 'SRP038101'

# Download data
url <- download_study(exp)

# Load data
load(file.path(exp, 'rse_gene.Rdata'))

# Looking at rse
rse_gene

# Sample phenotype provided by recount
colData(rse_gene)

## At the gene level, the row data includes the gene Gencode ids, the gene
## symbols and the sum of the disjoint exons widths, which can be used for
## taking into account the gene length.
rowData(rse_gene)

## At the exon level, you can get the gene Gencode ids from the names of:
# rowRanges(rse_exon)


# Thresholds
pThr <- 0.01
logFCThr <- 1
baseMeanThr <-20
cpmThr <- 1

## PHENOTYPE INFO
# Browse study project at SRA
browse_study(exp)

# View GEO ids
colData(rse_gene)$geo_accession

# Extract sample characteristics
geochar <- lapply(split(colData(rse_gene), seq_len(nrow(colData(rse_gene)))), geo_characteristics)

# Sample information
sample_info <- data.frame(
  run = colData(rse_gene)$run,
  group = ifelse(grepl("Untreated", colData(rse_gene)$title), "untrt", "trt"),
  cell.line = "AML3"
)

# Scale counts to read counts (needed for recount data)
rse <- scale_counts(rse_gene)

#######
## Add sample information for DE analysis
colData(rse)$group <- sample_info$group

# Create DEseq dataset: design formula variable of interest is treatment status
dds <- DESeqDataSet(rse, ~ group)

# Pre-filter to remove low read counts: Make file size smaller and computationally faster to run downstream
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Specify untreated as the control (ref) group
colData(dds)$group <- factor(colData(dds)$group, levels = c("untrt","trt"))

# Normalize counts (used for downstream visualization) and save file as .txt
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
write.table(normalized_counts, file="normalized_counts.txt", sep="\t", quote=F, col.names=NA)

###############################
#### EXPLORING DE RESULTS I ####
##############################

# Perform DE analysis
dds <- DESeq(dds)

# Extract results
res <- results(dds, alpha = 0.05)

# Shrink LFC
## Create unshrunk one to compare
unshrunk_res <- res

# Apply LFC shrink
res <- lfcShrink(dds, coef = "group_trt_vs_untrt", type = "apeglm")

# Summarize results to see no of genes tested, genes up/down de 
summary(res, alpha = 0.05)

# Create a tibble of results
res_tb <- res %>%
  data.frame() %>%
  rownames_to_column(var="ensembl_id") %>% 
  as_tibble() %>%
  # Edit ENSMEBL IDs to remove version (digits after decimal)
  mutate(ensembl_id= gsub('\\.[0-9]*$', '', ensembl_id))

##################################
#### OBTAINING SIG GENES LIST ####
##################################

# Creating a tibble of ENSEMBL (GENEID), ENTREZID, SYMBOL, GENENAME using EnsDb.Hsapiens.v86 db
anno <- as_tibble(AnnotationDbi::select(EnsDb.Hsapiens.v86,res_tb$ensembl_id, 
                                        columns=c("GENEID", "ENTREZID", "SYMBOL", "GENENAME"), 
                                        keytype="GENEID"))

# Filter acc to thresholds set earlier to obtain list of sig. genes. 
sigRes <- res_tb %>%
  dplyr::filter(res$padj <= pThr & 
           abs( res$log2FoldChange ) >= logFCThr & 
           res$baseMean >= baseMeanThr ) %>%
  # Edit ENSMEBL IDs to remove version (digits after decimal)
  mutate(ensembl_id= gsub('\\.[0-9]*$', '', ensembl_id)) %>%
  # Use <<left_join>> to merge `anno` and `sigRes` 
  left_join(anno, by = c("ensembl_id" ="GENEID"))


# Use <<left_join>> to merge `anno` and `res` 
wholeRes <- res_tb %>%
  left_join(anno, by = c("ensembl_id" ="GENEID"))

#############################
#### DATA VISUALIZATION ####
#############################
## PRE-DESEQ2
# Heatmap of poisson distances between samples. Input: non-normalized raw count data
# Euclidean distance for normalized
# Purpose: To check similarity between samples. Ensure that samples are differentiated by treatment cond
poisd <- PoissonDistance(t(counts(dds)))
samplePoisDistMatrix <- as.matrix( poisd$dd)
rownames(samplePoisDistMatrix) <- dds$title
colnames(samplePoisDistMatrix) <- NULL
colors = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors)

## MA Plot to see avg expression VS log2fc
#Sh runk
plotMA(res, ylim=c(-2,2))

# Unshrunk
plotMA(unshrunk_res, ylim=c(-2,2))

####

# Histogram of adj. p-values of sigRes (prev alrdy filtered out genes with low count with alpha = 0.05)
hist(res$pvalue, breaks= 20, col="grey50", border="white" )

## Volcano Plot with whole gene
#Default cutoffs are log2FC > |2| and adjusted P-value < 0.05
# Good to show this to illustrate your definition of up/downreg genes bc others can see 
# roughly the thresholds you use? compared to just stating how many is up/down
EnhancedVolcano( wholeRes, lab = wholeRes$SYMBOL, 
                 x = 'log2FoldChange', y = 'padj',
                 xlim = c(-8, 8), title = 'Treated vs Untreated',
                 pCutoff = 0.05, FCcutoff = 0.58, pointSize = 2.0, 
                 labSize = 3.0,
                 border = "full", borderWidth = 1.5, borderColour = "black", 
                 gridlines.major = FALSE, gridlines.minor = FALSE)

# PCA Plot: Normalized data
# Perform regularized logarithm transformation (rlog) on the data; VST (>30 samples)
rld <- rlog(dds)
rownames(rld) <- gsub('\\.[0-9]*$', '', rownames(rld))

# Plot PCA: How much of the variation in the data is attributed to treatment condn? Look at PC1.
pcaData <- plotPCA(rld, intgroup = "group", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color = group)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() + labs(color='Treatment')

######
# SigRes




##Later transfer on top###
dds$cell = as.factor(dds$title)
dds$cell = gsub('^[[:alpha:]]*\\.', 'Replicate ', dds$cell)
dds$cell = as.factor(dds$cell)
#######
# gene with largest positive log2 Fold Change
large_GeneCounts <- plotCounts(dds, gene=which.max(sigRes$log2FoldChange), intgroup=c("group", "cell"), returnData = TRUE)
large_Gene <- sigRes$GENENAME[which.max(sigRes$log2FoldChange)]
ggplot(large_GeneCounts, aes(x = group, y = count, color = cell, group = cell)) + geom_point(size = 4) + geom_line(size = 1,linetype = 2)  + theme_classic()

######
# Heatmap of Top 20 Up/Down reg
mat <- assay(rld)   # Use the rlog normalized counts
idx <- sigRes$ensembl_id 
DEgenes <- mat[idx,]

mat$cell = gsub('^[[:alpha:]]*\\.', 'Replicate ', dds$cell)
mat$cell = factor(dds$cell, levels = c("Replicate 1", "Replicate 2", "Replicate 3"))

annotation <- as.data.frame(colData(rld)[, c("title", "group")])
annotation$title = as.factor(gsub('^[[:alpha:]]*\\.', 'Replicate ', annotation$title))

pheatmap(DEgenes, scale = "row", show_rownames = FALSE, clustering_distance_rows = "correlation", annotation_col = annotation, main="Differentially Expressed Genes", returnData = TRUE)

# View top 20 sigres by log2fc
top_20_DE <- sigRes %>%
  mutate(absLFC = abs(log2FoldChange)) %>%
  arrange(desc(absLFC)) %>%
  top_n(20)

mat <- assay(rld) 
top_20_DE_hm <- mat[top_20_DE$ensembl_id,]
rownames(top_20_DE_hm) <- top_20_DE$SYMBOL

pheatmap(top_20_DE_hm, scale = "row", clustering_distance_rows = "correlation", annotation_col = annotation, main="Top 20 Differentially Expressed genes")


