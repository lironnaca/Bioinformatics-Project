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
library(stringr)
library(data.table)
library(topGO)
library(genefilter)
library("writexl")
library(GO.db)
library(biomaRt)
library(Rgraphviz)
library(factoextra)
library(cluster)
library(ggalluvial)

# # Required installation
# install_github("jokergoo/ComplexHeatmap")
# if (!require("BiocManager"))
#   install.packages("BiocManager")
# BiocManager::install("clusterProfiler", version = "3.15")
# BiocManager::install("pathview")
# BiocManager::install("enrichplot")
# BiocManager::install("org.Hs.eg.db", character.only = TRUE)
# BiocManager::install("GenomicSuperSignature")
# BiocManager::install("bcellViper")
# BiocManager::install("topGO")
# BiocManager::install("GO.db")
# BiocManager::install("biomaRt")
# BiocManager::install("Rgraphviz")
# BiocManager::install("genefilter")
# install.packages("factoextra")

# Move to the working directory of the project

# # Liron
setwd("G:/.shortcut-targets-by-id/1tOZ9FLIG74AOVRYD5wR_OLB7btxxeZ_X/student exchange/courses/bioinformatics/Project/Assignment3_Rproject")

# Seggev
# setwd("G:/My Drive/Folder With Liron/student exchange/courses/bioinformatics/Project/Assignment3_Rproject")

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

metadata <- metadata %>%
  dplyr::mutate(type = factor(type, levels = c("critical", "non-critical")))

expression_df <- readr::read_csv(data_file)

# preprocsessing the data
gene_names <- str_split_fixed(expression_df$Gene, "_", 2)

expression_df$Gene <- gene_names[,-1]

no_dup_data <- as.data.table(expression_df)

keys <- c('Gene')

no_dup_data <- no_dup_data[,lapply(.SD, mean),keys]

expression_df <- as.data.frame(no_dup_data) %>% tibble::column_to_rownames("Gene")

Min <- min(expression_df)

func <- function(x){
  return(x-Min)
}
expression_df <- apply(expression_df, 2, func)

gene_matrix <- round(expression_df)




### 2.a ###

vari <- apply(expression_df, 1, var)
vari <- data.frame(vari)
colnames(vari)[1] <- "varience"
vari <- cbind(expression_df, vari)
vari <- vari[order(-vari$varience),]
most_varied <- vari[1:5000,-70]
most_varied <- t(most_varied)

df <- scale(most_varied)

### task 2 - K-means ###

set.seed(1)

km <- kmeans(df, centers = 5, nstart = 25)

fviz_cluster(km, data = df)


### task 2 - Hierarchical clustering (hclust) ###

distance_mat <- dist(most_varied, method = 'euclidean')

set.seed(1)

Hierar_cl <- hclust(distance_mat, method = "average")

plot(Hierar_cl)

fit <- cutree(Hierar_cl, k = 5 )

rect.hclust(Hierar_cl, k = 5, border = "red")

fviz_cluster(list(data = most_varied, cluster = fit))


### task 2 - PAM clustering ###

set.seed(1)

fviz_nbclust(df, pam, method ="silhouette")+theme_minimal()
pamResult <- pam(df, k = 5)

fviz_cluster(pamResult, data=df)


### task 2 - Alluvial plot ###

final <- cbind(metadata$type, as.data.frame(km$cluster), as.data.frame(fit), as.data.frame(pamResult$clustering))

colnames(final) <- c('type', 'km', 'hc', 'pam')
freq_table <- as.data.frame(table(final))

ggplot(as.data.frame(freq_table),
       aes(y = Freq,
           axis1 = km, axis2 = hc, axis3 = pam)) +
  geom_alluvium(aes(fill = type),
                width = 0, knot.pos = 0, reverse = FALSE) +
  guides(fill = FALSE) +
  geom_stratum(width = 1/8, reverse = FALSE) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),
            reverse = FALSE) +
  scale_x_continuous(breaks = 1:3, labels = c("km", "hc", "pam")) +
  coord_flip() +
  ggtitle("hara")



##### task 3 #####
library(circlize)
ha <- HeatmapAnnotation(km = km$cluster, hc = fit, pam = pamResult$clustering,
                        col=list(km=colorRamp2(c(1, 3, 5), c("brown", "white", "orange")),
                                 hc=colorRamp2(c(1, 3, 5), c("brown", "white", "orange")),
                                 pam=colorRamp2(c(1, 3, 5), c("brown", "white", "orange"))))
map <- Heatmap(t(most_varied), cluster_rows = T, cluster_columns = F,
               column_labels = colnames(t(most_varied)), name = "Heat Bar", show_row_name = FALSE, bottom_annotation = ha)

map











# # Starting the differential expression analysis
# ddset <- DESeqDataSetFromMatrix(countData = gene_matrix,
#                                 colData = metadata,
#                                 design = ~type)
# 
# deseq_object <- DESeq(ddset)
# 
# deseq_results <- results(deseq_object)
# 
# deseq_results <- lfcShrink(deseq_object, coef = 2, res = deseq_results)
# 
# deseq_df <- deseq_results %>% as.data.frame() %>%
#   tibble::rownames_to_column("Gene") %>%
#   dplyr::mutate(threshold = padj < 0.05) %>%
#   dplyr::arrange(dplyr::desc(log2FoldChange))
# 
# write_xlsx(deseq_df,"results/deseq.xlsx")

# plotCounts(ddset, gene = "DPM1", intgroup = "type")


# ### Task 4 ###
# filterd_data <- deseq_df[deseq_df[5] < 0.01 , ]
# filterd_data <- data.frame(filterd_data)
# 
# filterd_data <- filterd_data[filterd_data[3] < -1 | filterd_data[3] > 1, ]
# rownames(filterd_data) <- filterd_data[,1]
# 
# filterd_data <- filterd_data[,-1]
# mat <- counts(ddset)[rownames(filterd_data), ]
# mat <- t(apply(mat, 1, scale))
# coldata <- metadata %>% tibble::column_to_rownames("id")
# colnames(mat) <- rownames(coldata)
# 
# map <- Heatmap(t(most_varied), cluster_rows = T, cluster_columns = F,
#                column_labels = colnames(t(most_varied)), name = "Heat Bar")
# 
# map


