# figures for unsupervised clustering of tcga samples 
library(data.table)
library(ggplot2)
library(ggpubr)
library(ggrepel)
ggplot2::theme_set(ggpubr::theme_pubclean())

source('/fast/AG_Akalin/buyar/flexynesis_manuscript_work/flexynesis_manuscript/src/utils.R')
folder <- '/fast/AG_Akalin/buyar/flexynesis_manuscript_work/analyses/unsupervised_cancertype/'

get_plot <- function(df_tsne, factors, label_size = 5) {
  centroids <- data.frame(get_basis_matrix(df_tsne, factors))
  colnames(centroids) <- c('x', 'y')
  df_tsne$factor <- factors
  centroids$group <- as.factor(rownames(centroids))
  p <- ggplot(df_tsne, aes(x = tSNE1, y = tSNE2)) +
    geom_point(aes(color = factor), alpha = 0.6) +
    geom_text_repel(data = centroids, aes(x = x, y = y, 
                                    label = group), 
              size = label_size) + 
  theme(legend.position = 'none', text = element_text(size = 18)) 
  return(p)
}

E <- read.csv(file.path(folder, 'cancertype.embeddings_train.csv'), row.names = 1)
df_tsne <- plot_tsne(E, returnData = T)

clusters <- read.csv(file.path(folder, 'clusters.csv'), header = T, row.names = 1)
rownames(clusters) <- clusters$sample
clusters$cohort <- gsub("TCGA-", "", clusters$label)

p1 <- get_plot(df_tsne, as.factor(clusters[rownames(df_tsne), 'cluster']), label_size = 5)
p1 
p2 <- get_plot(df_tsne, as.factor(clusters[rownames(df_tsne), 'cohort']), label_size = 5)
p2
p3 <- plot_cluster_comparison(clusters[,'cluster',drop=F], clusters[,'cohort',drop=F])
p3

p <- cowplot::plot_grid(cowplot::plot_grid(p1, p2, ncol = 1, labels = c('A', 'B')), 
                   p3, nrow = 1, rel_widths = c(1.5, 1), labels = c('A', 'C'))

ggsave(filename = 'tcga_cancertype_clustering.pdf', 
       plot = p, width = 10, height = 8)





