# Required Libraries
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
library(Biobase)
library(mclust)
library(circlize)
library(GGally)

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
# BiocManager::install("ConsensusClusterPlus")

# Move to the working directory of the project

# # # Liron
# setwd("G:/.shortcut-targets-by-id/1tOZ9FLIG74AOVRYD5wR_OLB7btxxeZ_X/student exchange/courses/bioinformatics/Project/Assignment3_Rproject")
# 
# # Seggev
# # setwd("G:/My Drive/Folder With Liron/student exchange/courses/bioinformatics/Project/Assignment3_Rproject")

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
vari <- vari[,-70]

all_the_clusters <- as.data.frame(metadata$type)

colnames(all_the_clusters)[1] <- 'metadata'

arr <- c('km_10_k_2', 'hc_10_k_2', 'pam_10_k_2', 'mc_10_k_2',
         'km_100_k_2', 'hc_100_k_2', 'pam_100_k_2', 'mc_100_k_2',
         'km_5000_k_2', 'hc_5000_k_2', 'pam_5000_k_2', 'mc_5000_k_2',
         'km_10000_k_2', 'hc_10000_k_2', 'pam_10000_k_2', 'mc_10000_k_2',
         'km_5000_k_4', 'hc_5000_k_4', 'pam_5000_k_4', 'mc_5000_k_4',
         'km_5000_k_6', 'hc_5000_k_6', 'pam_5000_k_6', 'mc_5000_k_6',
         'km_5000_k_10', 'hc_5000_k_10', 'pam_5000_k_10', 'mc_5000_k_10')

km_cluster <- function(df, gene_num, k){
  
  set.seed(1)
  
  km <- kmeans(df, centers = k, nstart = 25)
  
  all_the_clusters <<- cbind(all_the_clusters, km$cluster)
  
  colnames(all_the_clusters)[i] <<- paste("km_", gene_num, "_k_", k, sep="")
  
  i <<- i+1

  name_1 <- paste("Cluster Method:KM, K is", k,
                "Number of Genes:",
                gene_num, sep =" ") 
  name_2 <- paste(".\\plots\\KM_k_", k,"_Gene_num_", gene_num, ".png", sep="")
  
  png(filename = name_2, width = 900, height = 600)
  
  plot(fviz_cluster(km, data = df, main= name_1 ))
  
  dev.off()
}

hc_cluter <- function(most_varied, df, gene_num, k){
  
  distance_mat <- dist(most_varied, method = 'euclidean')
  
  set.seed(1)
  
  Hierar_cl <- hclust(distance_mat, method = "average")
  
  # plot(Hierar_cl)
  
  fit <- cutree(Hierar_cl, k = k )
  
  all_the_clusters <<- cbind(all_the_clusters, as.data.frame(fit))
  
  colnames(all_the_clusters)[i] <<- paste("hc_", gene_num, "_k_", k, sep="")
  
  i <<- i+1
  
  # rect.hclust(Hierar_cl, k = k, border = "red")
  
  name_1 <- paste("Cluster Method:HC, K is", k,
                  "Number of Genes:",
                  gene_num, sep =" ") 
  
  name_2 <- paste(".\\plots\\HC_k_", k,"_Gene_num_", gene_num, ".png", sep="")
  
  png(filename = name_2, width = 900, height = 600)
  
  plot(fviz_cluster(list(data = most_varied, cluster = fit), data = df, main= name_1 ))
  
  dev.off()
}

pam_cluster <- function(df, gene_num, k){
  set.seed(1)
  
  pamResult <- pam(df, k = k)
  
  all_the_clusters <<- cbind(all_the_clusters, pamResult$clustering)
  
  colnames(all_the_clusters)[i] <<- paste("pam_", gene_num, "_k_", k, sep="")
  
  i <<- i+1
  
  name_1 <- paste("Cluster Method:PAM, K is", k,
                  "Number of Genes:",
                  gene_num, sep =" ") 
  name_2 <- paste(".\\plots\\PAM_k_", k,"_Gene_num_", gene_num,".png", sep="")
  
  png(filename = name_2, width = 900, height = 600)
  
  plot(fviz_cluster(pamResult, data = df, main= name_1 ))
  
  dev.off()
  
}

mc_cluster <- function(most_varied,df, gene_num, k){
  
  set.seed(1)
  
  mc <- Mclust(most_varied,G=k)
  
  all_the_clusters <<- cbind(all_the_clusters, mc$classification)
  
  colnames(all_the_clusters)[i] <<- paste("mc_", gene_num, "_k_", k, sep="")
  
  i <<- i+1
  
  name_1 <- paste("Cluster Method:MC, K is", k,
                  "Number of Genes:",
                  gene_num, sep =" ") 
  name_2 <- paste(".\\plots\\MC_k_", k,"_Gene_num_", gene_num, ".png", sep="")
  
  png(filename = name_2, width = 900, height = 600)
  
  plot(fviz_cluster(mc, data = df, main= name_1 ))
  
  dev.off()
  
}

create_four_clusters <- function(gene_num, k){
  most_varied <- vari[1:gene_num,]
  most_varied <- t(most_varied)
  
  df <- scale(most_varied)
  

  print(km_cluster(df, gene_num, k))
  
  print(hc_cluter(most_varied, df, gene_num, k))

  print(pam_cluster(df, gene_num, k))
  
  print(mc_cluster(most_varied, df, gene_num, k))
  
}

i <- 2

for (gene_num in c(10,100,1000,5000,10000)){
  create_four_clusters(gene_num, 2)
}

for (k in c(4, 6, 10)){
  create_four_clusters(5000, k)
}


### task 3 - Alluvial plot ###

final <- cbind(metadata$type, as.data.frame(all_the_clusters$km_10_k_2),
               as.data.frame(all_the_clusters$km_100_k_2),
               as.data.frame(all_the_clusters$km_1000_k_2),
               as.data.frame(all_the_clusters$km_10000_k_2))

freq_table <- as.data.frame(table(final))
colnames(freq_table) <- c('type', 'KM_10', 'KM_100', 'KM_1000', 'KM_10000', 'Freq')


ggplot(as.data.frame(freq_table),
       aes(y = Freq,
           axis1 = KM_10, axis2 = KM_100, axis3 = KM_1000, axis4= KM_10000)) +
  geom_alluvium(aes(fill = type),
                width = 0, knot.pos = 0, reverse = FALSE) +
  guides(fill = FALSE) +
  geom_stratum(width = 1/8, reverse = FALSE) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),
            reverse = FALSE) +
  scale_x_continuous(breaks = 1:4, labels = c('KM_10', 'KM_100', 'KM_1000', 'KM_10000')) +
  coord_flip() +
  ggtitle("Alluvium Plot- KM")


final <- cbind(metadata$type, as.data.frame(all_the_clusters$hc_10_k_2),
               as.data.frame(all_the_clusters$hc_100_k_2),
               as.data.frame(all_the_clusters$hc_1000_k_2),
               as.data.frame(all_the_clusters$hc_10000_k_2))

freq_table <- as.data.frame(table(final))
colnames(freq_table) <- c('type', 'HC_10', 'HC_100', 'HC_1000', 'HC_10000', 'Freq')


ggplot(as.data.frame(freq_table),
       aes(y = Freq,
           axis1 = HC_10, axis2 = HC_100, axis3 = HC_1000, axis4= HC_10000)) +
  geom_alluvium(aes(fill = type),
                width = 0, knot.pos = 0, reverse = FALSE) +
  guides(fill = FALSE) +
  geom_stratum(width = 1/8, reverse = FALSE) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),
            reverse = FALSE) +
  scale_x_continuous(breaks = 1:4, labels = c('HC_10', 'HC_100', 'HC_1000', 'HC_10000')) +
  coord_flip() +
  ggtitle("Alluvium Plot- HC")


final <- cbind(metadata$type, as.data.frame(all_the_clusters$pam_10_k_2),
               as.data.frame(all_the_clusters$pam_100_k_2),
               as.data.frame(all_the_clusters$pam_1000_k_2),
               as.data.frame(all_the_clusters$pam_10000_k_2))

freq_table <- as.data.frame(table(final))
colnames(freq_table) <- c('type', 'PAM_10', 'PAM_100', 'PAM_1000', 'PAM_10000', 'Freq')


ggplot(as.data.frame(freq_table),
       aes(y = Freq,
           axis1 = PAM_10, axis2 = PAM_100, axis3 = PAM_1000, axis4= PAM_10000)) +
  geom_alluvium(aes(fill = type),
                width = 0, knot.pos = 0, reverse = FALSE) +
  guides(fill = FALSE) +
  geom_stratum(width = 1/8, reverse = FALSE) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),
            reverse = FALSE) +
  scale_x_continuous(breaks = 1:4, labels = c('PAM_10', 'PAM_100', 'PAM_1000', 'PAM_10000')) +
  coord_flip() +
  ggtitle("Alluvium Plot- PAM")


final <- cbind(metadata$type, as.data.frame(all_the_clusters$mc_10_k_2),
               as.data.frame(all_the_clusters$mc_100_k_2),
               as.data.frame(all_the_clusters$mc_1000_k_2),
               as.data.frame(all_the_clusters$mc_10000_k_2))

freq_table <- as.data.frame(table(final))
colnames(freq_table) <- c('type', 'MC_10', 'MC_100', 'MC_1000', 'MC_10000', 'Freq')


ggplot(as.data.frame(freq_table),
       aes(y = Freq,
           axis1 = MC_10, axis2 = MC_100, axis3 = MC_1000, axis4= MC_10000)) +
  geom_alluvium(aes(fill = type),
                width = 0, knot.pos = 0, reverse = FALSE) +
  guides(fill = FALSE) +
  geom_stratum(width = 1/8, reverse = FALSE) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),
            reverse = FALSE) +
  scale_x_continuous(breaks = 1:4, labels = c('MC_10', 'MC_100', 'MC_1000', 'MC_10000')) +
  coord_flip() +
  ggtitle("Alluvium Plot- MC")

##### task 3 #####

most_varied <- vari[1:5000,]
most_varied <- t(most_varied)

ha <- HeatmapAnnotation(km = all_the_clusters$km_5000_k_4,
                        hc = all_the_clusters$hc_5000_k_4,
                        pam = all_the_clusters$pam_5000_k_4,
                        mc = all_the_clusters$mc_5000_k_4,
                        col=list(km=colorRamp2(c(1, 2, 4), c("brown", "white", "orange")),
                                 hc=colorRamp2(c(1, 2, 4), c("brown", "white", "orange")),
                                 pam=colorRamp2(c(1, 2, 4), c("brown", "white", "orange")),
                                 mc=colorRamp2(c(1, 2, 4), c("brown", "white", "orange"))))
map <- Heatmap(t(most_varied), cluster_rows = T, cluster_columns = F,
               column_labels = colnames(t(most_varied)), name = "Heat Bar", show_row_name = FALSE, bottom_annotation = ha)

map


#### task 4.a ####

options(warn=0)

real_vs_km_table <- table(metadata$type, all_the_clusters$km_5000_k_2)
real_vs_km <- chisq.test(metadata$type, all_the_clusters$km_5000_k_2, correct = FALSE)

real_vs_hc_table <- table(metadata$type, all_the_clusters$hc_5000_k_2)
real_vs_hc <- chisq.test(metadata$type, all_the_clusters$hc_5000_k_2, correct = FALSE)

real_vs_pam_table <- table(metadata$type, all_the_clusters$pam_5000_k_2)
real_vs_pam <- chisq.test(metadata$type, all_the_clusters$pam_5000_k_2, correct = FALSE)

real_vs_mc_table <- table(metadata$type, all_the_clusters$mc_5000_k_2)
real_vs_mc <- chisq.test(metadata$type, all_the_clusters$mc_5000_k_2, correct = FALSE)


#### task 4.b ####

couples <- data.frame()
arr <- c('metadata','km_10_k_2', 'hc_10_k_2', 'pam_10_k_2', 'mc_10_k_2',
         'km_100_k_2', 'hc_100_k_2', 'pam_100_k_2', 'mc_100_k_2',
         'km_1000_k_2', 'hc_1000_k_2', 'pam_1000_k_2', 'mc_1000_k_2',
         'km_5000_k_2', 'hc_5000_k_2', 'pam_5000_k_2', 'mc_5000_k_2',
         'km_10000_k_2', 'hc_10000_k_2', 'pam_10000_k_2', 'mc_10000_k_2')

for (i in 1:20){
  for (j in (i+1):21){
    result <- chisq.test(table(cbind(all_the_clusters[i], all_the_clusters[j])), correct = FALSE)
    couples <- rbind(couples,data.frame(Method1 = arr[i],Method2 = arr[j], p_value=result$p.value, x_square=result$statistic))
  }
}

couples <- cbind(couples, p.adjust(couples$p_value))
colnames(couples)[5] <- "p.adjust"

little_couples <- couples[13:16,]
little_couples <- little_couples[-4]



#### task 4.e ####

ggpairs(couples, c(1,2,3,5), cardinality_threshold=20)
