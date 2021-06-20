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


### DOWNLOAD AND LOAD DATA
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

# Scale counts to read countså
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

# MA Plot
## Shrunk
plotMA(res, ylim=c(-2,2))

## Unshrunk
plotMA(unshrunk_res, ylim=c(-2,2))

# Summarize results to see no of genes tested, genes up/down de 
summary(res, alpha = 0.05)

# Create a tibble of results
res_tb <- res %>%
  data.frame() %>%
  rownames_to_column(var="ensembl_id") %>% 
  as_tibble()

##################################
#### OBTAINING SIG GENES LIST ####
##################################

# Filter acc to thresholds set earlier to obtain list of sig. genes. 
sigRes <- res_tb %>%
  filter(res$padj <= pThr & 
           abs( res$log2FoldChange ) >= logFCThr & 
           res$baseMean >= baseMeanThr ) %>%
  # Edit ENSMEBL IDs to remove version (digits after decimal)
  mutate(ensembl_id= gsub('\\.[0-9]*$', '', ensembl_id))


# Creating a tibble of ENSEMBL, ENTREZID, SYMBOL, GENENAME 
anno <- as_tibble(AnnotationDbi::select(org.Hs.eg.db, sigRes$gene, 
                              columns=c("ENSEMBL", "ENTREZID", "SYMBOL", "GENENAME"), 
                              keytype="ENSEMBL"))

# Use <<left_join>> to merge `anno` and `sigRes` 
sigRes <- sigRes %>%
  left_join(anno, by = c("ensembl_id" ="ENSEMBL"))








