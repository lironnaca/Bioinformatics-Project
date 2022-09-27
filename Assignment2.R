# Required Libraries
library(DESeq2)
library(ggplot2)
library(magrittr)
library(M3C)
library(reshape2)
library(devtools)
library(ComplexHeatmap)
library(clusterProfiler)
library(enrichplot)
library(ggridges)
library(ggnewscale)
library(gprofiler2)
library(GenomicSuperSignature)
library(bcellViper)
library("org.Hs.eg.db", character.only = TRUE)

# Required installation
install_github("jokergoo/ComplexHeatmap")
if (!require("BiocManager"))
  install.packages("BiocManager")
BiocManager::install("clusterProfiler", version = "3.15")
BiocManager::install("pathview")
BiocManager::install("enrichplot")
BiocManager::install("org.Hs.eg.db", character.only = TRUE)
BiocManager::install("GenomicSuperSignature")
BiocManager::install("bcellViper")
# Move to the working directory of the project

setwd("G:/.shortcut-targets-by-id/1tOZ9FLIG74AOVRYD5wR_OLB7btxxeZ_X/student exchange/courses/bioinformatics/Project/Bio2")

#### Define folders for later use####
if (!dir.exists("data")) {
  dir.create("data")
}

plots_dir <- "plots"

if (!dir.exists(plots_dir)) {
  dir.create(plots_dir)
}

results_dir <- "results"

if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}

# Define paths 
data_file <- file.path("data", "data.csv")

metadata_file <- file.path("data", "metadata.csv")

# Setting seed to be able to keep track of results
set.seed(12345)

# Retrieving the Data and the Metadata files
metadata <- readr::read_csv(metadata_file)

expression_df <- readr::read_csv(data_file) %>%
  tibble::column_to_rownames("Gene")

metadata <- metadata %>%
  dplyr::mutate(type = factor(type, levels = c("critical", "non-critical")))

gene_matrix <- round(expression_df)

# Starting the differential expression analysis
ddset <- DESeqDataSetFromMatrix(countData = gene_matrix,
                                colData = metadata,
                                design = ~type)

deseq_object <- DESeq(ddset)

deseq_results <- results(deseq_object)

deseq_results <- lfcShrink(deseq_object, coef = 2, res = deseq_results)

deseq_df <- deseq_results %>% as.data.frame() %>%
            tibble::rownames_to_column("Gene") %>%
            dplyr::mutate(threshold = padj < 0.05) %>%
            dplyr::arrange(dplyr::desc(log2FoldChange))

# plotCounts(ddset, gene = "DPM1", intgroup = "type")

### Task 1 ###

# Making the density plots
mini <- apply(expression_df, 1, min)
mini <- data.frame(mini)
colnames(mini)[1] <- "count"
mini_plot <- ggplot(mini, aes(x=count))+
  geom_density(color="darkblue", fill="lightblue")
mini_plot

max <- apply(expression_df, 1, max)
max <- data.frame(max)
colnames(max)[1] <- "count"
max_plot <- ggplot(max, aes(x=count))+
  geom_density(color="darkblue", fill="lightblue")
max_plot

median <- apply(expression_df, 1, median)
median <- data.frame(median)
colnames(median)[1] <- "count"
median_plot <- ggplot(median, aes(x=count))+
  geom_density(color="darkblue", fill="lightblue")
median_plot

### Task 2 ###

# PCA
rld <- rlog(ddset)
plotPCA(rld, intgroup=c("type"))

# T-SNE
tsne(gene_matrix, labels=as.factor(metadata$type))

### Task 3 ###

volcano_plot <- EnhancedVolcano::EnhancedVolcano(
  deseq_df,
  lab = deseq_df$Gene,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.05 
)
volcano_plot

### Task 4 ###
filterd_data <- deseq_df[deseq_df[5] < 0.01 , ]
filterd_data <- data.frame(filterd_data)

filterd_data <- filterd_data[filterd_data[3] < -1 | filterd_data[3] > 1, ]
rownames(filterd_data) <- filterd_data[,1]

filterd_data <- filterd_data[,-1]
mat <- counts(ddset)[rownames(filterd_data), ]
mat <- t(apply(mat, 1, scale))
coldata <- metadata %>% tibble::column_to_rownames("id")
colnames(mat) <- rownames(coldata)

map <- Heatmap(mat, cluster_rows = T, cluster_columns = T,
               column_labels = colnames(mat), name="try")

map
### Task 5 ###

# clustProfiler

original_gene_list <- deseq_df$log2FoldChange

names(original_gene_list) <- deseq_df$Gene

gene_list<-na.omit(original_gene_list)

gene_list = sort(gene_list, decreasing = TRUE)


gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "SYMBOL", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = "org.Hs.eg.db", 
             pAdjustMethod = "none")
require(DOSE)
clustProfiler_plot_1 <- dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)
clustProfiler_plot_1

# gProfiler2
gostres2 <- gost(query = deseq_df$Gene, 
                 organism = "hsapiens")
gProfiler2_plot <- gostplot(gostres2, capped = FALSE, interactive = FALSE)
gProfiler2_plot

# GenomicSuperSignature
RAVmodel <- getModel("PLIERpriors", load=TRUE)
val_all <- validate(ddset, RAVmodel)
GenomicSuperSignature_plot <- plotValidate(val_all, interactive = FALSE)
GenomicSuperSignature_plot
