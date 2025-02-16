# figures for unsupervised clustering of tcga samples 
library(data.table)
library(ggplot2)
library(ggpubr)
library(ggrepel)
ggplot2::theme_set(ggpubr::theme_pubclean())

args <- commandArgs(trailingOnly = T)

utils_script <- args[1]
folder <- args[2]

source(utils_script)

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

# plot correspondence between cluster labels and cancer types 
m <- table(clusters$cohort, clusters$cluster)
m <- t(apply(m, 1, function(x) round(x/sum(x)*100,1)))
colors <- colorRampPalette(c("white",'darkblue', "red"))(100)
ami <- aricode::AMI(clusters$cohort, clusters$cluster)
p3 <- pheatmap::pheatmap(m, silent = F, cluster_rows = T, cluster_cols = T, 
                   display_numbers = F, color = colors, cellwidth = 20, cellheight = 10,
                   treeheight_row = 0, treeheight_col = 0, 
                   main = paste0("Correspondence between cluster labels and cancer types\nAdjusted Mutual Information: ",round(ami,2)))
p <- cowplot::plot_grid(cowplot::plot_grid(p1, p2, ncol = 2, labels = c('A', 'B')), 
                   p3$gtable, nrow = 2, labels = c('A', 'C'))

ggsave(filename = 'tcga_cancertype_clustering.pdf', 
       plot = p, width = 10, height = 8)







