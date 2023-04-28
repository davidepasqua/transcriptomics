#======================================================================

# Loading Packages
#======================================================================

Load necessary libraries for analysis
library("DESeq2") # DESeq2: a package for statistical analysis of RNASeq data
library("factoextra") # factoextra: a package for multivariate analysis (PCA)
library("EnhancedVolcano") # EnhancedVolcano: a package for creating volcano plots
library("gprofiler2") # gprofiler2: a package for functional enrichment analysis

#======================================================================

# Loading Data
#======================================================================

# setwd - set working directory this is the directory where input files are located and where output # data will be saved

setwd("/Users/davide/Desktop/transcriptomics/")

# path of input files (count matrix and metadata)
input_file_data <- "input_genes_DESeq2.txt"
input_file_meta <- "PRJNA867309/files/metadata.tsv"

# load dataset and metadata using read.table function
# sep="\t", indicates that the separator in the file is a tab
# row.names=1, indicates that the first column contains row names
# header=TRUE, indicates that the first row contains column names the data is stored in a data frame
dataset <- read.table(input_file_data, sep="\t", row.names=1, header=TRUE)
metadata <- read.table(input_file_meta, sep="\t", row.names=1, header=TRUE)

#======================================================================

# Multivariate Analysis of Data
#======================================================================

# for PCA the data frame must be transposed using the t() function
my.data <- t(dataset)

# remove genes that are never identified (the sum of the columns is 0!)
my.data = my.data[,colSums(my.data) > 0]

prcomp performs a transformation of the data so that it has mean=0 and standard deviation=1
pca <- stats::prcomp(my.data, scale = TRUE, center = TRUE)

# PCA plot
O <- fviz_pca_ind(pca, label="all", habillage=as.factor(metadata$Phenotype), addEllipses=TRUE, ellipse.level=0.95, ellipse.type="convex")
O

# plot percentage of variance explained by each axis
A <- fviz_screeplot(pca, addlabels = TRUE, ylim = c(0, 60))

# plot contribution of first 10 eigenvectors on first two principal components
P <- fviz_pca_biplot(pca, label="var", habillage=as.factor(metadata$Phenotype),
addEllipses=TRUE, geom.var="text", ellipse.type="convex", ellipse.level=0.95, select.var=list(contrib=10))

# contribution for each axis
P1 <- fviz_contrib(pca, choice="var", axes = 1 , top=10)
P2 <- fviz_contrib(pca, choice="var", axes = 2 , top=10)
P3 <- fviz_contrib(pca, choice="var", axes = 3 , top=10)
P

#======================================================================

Analysis with DESeq2
#======================================================================

# create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = dataset, colData = metadata, design = ~ genotype cell type)

# perform analysis
des <- DESeq( dds )

# get results
res <- results( des )

# save results in a file
write.table(res,"results_deseq2.txt",quote=FALSE,dec=",",sep="\t")

#======================================================================
Volcano plot
#======================================================================

EnhancedVolcano(res, lab = NA, x = 'log2FoldChange', y = 'padj', xlim = c(-5, 5), title = 'Volcano plot', pCutoff = 0.05)
#======================================================================

# Functional enrichment
#======================================================================

# Save results with padj FDR<=0.05 and fold-change >= |2| (absolute fold-change)
genes <- subset(res, padj<=0.05 & abs(log2FoldChange)>=2.0)

# Save gene names
gene_names <- rownames(genes)

# Perform enrichment
gostres <- gost(query = gene_names,
organism = "hsapiens", ordered_query = FALSE,
multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
measure_underrepresentation = FALSE, evcodes = TRUE,
user_threshold = 0.05, correction_method = "g_SCS",
domain_scope = "annotated", custom_bg = NULL,
numeric_ns = "", sources = NULL, as_short_link = FALSE)

# Show interactive plot
gostplot(gostres, capped = TRUE, interactive = TRUE)

# Save results in a matrix
gostres_res <- gostres$result
gostres_res <- as.matrix(gostres_res)
write.table(gostres_res, "results_gos_NF_no_fdr.txt", quote = FALSE, dec = ",", sep = "\t")