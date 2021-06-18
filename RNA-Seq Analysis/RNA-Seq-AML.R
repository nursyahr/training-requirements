Load libraries
library(recount)
library(DESeq2)
library(dplyr)


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

# Perform DE analysis
dds <- DESeq(dds)

# Extract results
res <- results(dds)

# Select genes that are significant according to thresholds set previously
# Get the index position of genes that meet these threshold
