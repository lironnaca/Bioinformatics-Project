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

# Move to the working directory of the project

# # Liron
# setwd("G:/.shortcut-targets-by-id/1tOZ9FLIG74AOVRYD5wR_OLB7btxxeZ_X/student exchange/courses/bioinformatics/Project/Bio2")

# Seggev
setwd("G:/My Drive/Folder With Liron/student exchange/courses/bioinformatics/Project/Bio2")

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

write_xlsx(deseq_df,"results/deseq.xlsx")

# plotCounts(ddset, gene = "DPM1", intgroup = "type")

### Task 1 ###

# Making the density plots
mini <- apply(expression_df, 1, min)
mini <- data.frame(mini)
colnames(mini)[1] <- "count"
mini_plot <- ggplot(mini, aes(x=count))+
  geom_density(color="darkblue", fill="lightblue") + ggtitle("Minimum density plot")
mini_plot

max <- apply(expression_df, 1, max)
max <- data.frame(max)
colnames(max)[1] <- "count"
max_plot <- ggplot(max, aes(x=count))+
  geom_density(color="darkblue", fill="lightblue") + ggtitle("Maximum density plot") 
max_plot

median <- apply(expression_df, 1, median)
median <- data.frame(median)
colnames(median)[1] <- "count"
median_plot <- ggplot(median, aes(x=count))+
  geom_density(color="darkblue", fill="lightblue") + ggtitle("Median density plot")
median_plot

### Task 2 ###

# PCA
rld <- rlog(ddset)
plotPCA(rld, intgroup=c("type")) + ggtitle("PCA")

# T-SNE
tsne(gene_matrix, labels=as.factor(metadata$type))

### Task 3 ###

volcano_plot <- EnhancedVolcano::EnhancedVolcano(
  deseq_df,
  lab = deseq_df$Gene,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.05,
  title = "Volcano Plot"
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

map <- Heatmap(mat, cluster_rows = T, cluster_columns = F,
               column_labels = colnames(mat), name = "Heat Bar")

map
### Task 5 ###

# topGo

exp_data= read.table('data/data.csv', header=TRUE, sep=',')
bg_genes=as.character(exp_data[,1])

candidate_list =read.table('data/data.csv', header=TRUE, sep=',')
candidate_list= as.character(candidate_list[,1])

length(bg_genes)
head(bg_genes)

length(candidate_list)
head(candidate_list)

db= useMart('ENSEMBL_MART_ENSEMBL',dataset='hsapiens_gene_ensembl', host="www.ensembl.org")
go_ids= getBM(attributes=c('go_id', 'external_gene_name', 'namespace_1003'), filters='external_gene_name', values=bg_genes, mart=db)

gene_2_GO=unstack(go_ids[,c(1,2)])

keep = candidate_list %in% go_ids[,2]
keep =which(keep==TRUE)
candidate_list=candidate_list[keep]

geneList=factor(as.integer(bg_genes %in% candidate_list))
names(geneList)= bg_genes

GOdata=new('topGOdata', ontology='BP', allGenes = geneList, annot = annFUN.gene2GO, gene2GO = gene_2_GO)

classic_fisher_result=runTest(GOdata, algorithm='classic', statistic='fisher')

weight_fisher_result=runTest(GOdata, algorithm='weight01', statistic='fisher')

allGO=usedGO(GOdata)
all_res=GenTable(GOdata, weightFisher=weight_fisher_result, orderBy='weightFisher', topNodes=length(allGO))

p.adj=round(p.adjust(all_res$weightFisher,method="BH"),digits = 0)

all_res_final=cbind(all_res,p.adj)
all_res_final=all_res_final[order(all_res_final$p.adj),]

results.table.p= all_res_final[which(all_res_final$weightFisher<=0.001),]

results.table.bh=all_res_final[which(all_res_final$p.adj<=0.05),]

df_topGO <- as.data.frame(all_res_final)
write_xlsx(df_topGO,"results/topGO.xlsx")

pdf(file='plots/topGO_subgraph.pdf', height=12, width=12, paper='special', pointsize=18)
showSigOfNodes(GOdata, score(weight_fisher_result), useInfo = "none", sigForAll=FALSE, firstSigNodes=2,.NO.CHAR=50)
dev.off()

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

df_gse <- as.data.frame(gse)
write_xlsx(df_gse,"results/clustprofiler.xlsx")

# gProfiler2
gostres2 <- gost(query = deseq_df$Gene, 
                 organism = "hsapiens")
gProfiler2_plot <- gostplot(gostres2, capped = FALSE, interactive = FALSE)
gProfiler2_plot

df_gostres2 <- as.data.frame(gostres2$result)
write_xlsx(df_gostres2,"results/gprofiler2.xlsx")

# GenomicSuperSignature
RAVmodel <- getModel("PLIERpriors", load=TRUE)
val_all <- validate(ddset, RAVmodel)
GenomicSuperSignature_plot <- plotValidate(val_all, interactive = FALSE)
GenomicSuperSignature_plot

val_all$id <- rownames(val_all)
write_xlsx(val_all,"results/genomicsupersignature.xlsx")

