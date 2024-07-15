# analyse cross-modality prediction results from depmap

args <- commandArgs(trailingOnly = T)

input_path <- args[1] # path to depmap data 
output_path <- args[2] # path to flexynesis output 

library(data.table)
library(ggpubr)
library(ggplot2)

ggplot2::theme_set(ggpubr::theme_pubclean())
data.table::setDTthreads(4)
# compute average correlation between known and predicted dependency scores
compute_cors <- function(obs, pred, ori = 1) {
  cols <- intersect(colnames(obs), colnames(pred))
  rows <- intersect(rownames(obs), rownames(pred))
  obs <- obs[rows, cols]
  pred <- pred[rows, cols]
  if(ori == 1) {
    cors <- pbapply::pbsapply(rows, function(x) {
      cor(obs[x,], pred[x,])
    })
  } else if (ori == 2) {
    cors <- pbapply::pbsapply(cols, function(x) {
      cor(obs[,x], pred[,x])
    })
  }
  
  print(mean(cors, na.rm = T))
  return(cors)
}
getm <- function(dt) {
  as.matrix(data.frame(dt[,-1], row.names = dt[[1]], check.names = F))
}

f1_score <- function(var1, var2) {
  cont_table <- table(var1, var2)
  # Assuming 'TRUE' and 'FALSE' are the levels, we extract counts
  TP <- cont_table["TRUE", "TRUE"]
  TN <- cont_table["FALSE", "FALSE"]
  FP <- cont_table["FALSE", "TRUE"]
  FN <- cont_table["TRUE", "FALSE"]
  # Calculate Precision and Recall
  Precision <- TP / (TP + FP)
  Recall <- TP / (TP + FN)
  # Calculate F1 Score
  F1_Score <- 2 * (Precision * Recall) / (Precision + Recall)
  # Print the F1 Score
  return(F1_Score)
}

compute_f1s <- function(obs, pred, thr = 0.5) {
  cols <- intersect(colnames(obs), colnames(pred))
  rows <- intersect(rownames(obs), rownames(pred))
  obs <- obs[rows, cols]
  pred <- pred[rows, cols]
  
  scores <- pbapply::pbsapply(rows, function(x) {
    f1_score(obs[x, cols] > thr, pred[x, cols] > thr)
  })
  return(scores)
}

main <- function(prefix, layer) {
  train <- getm(data.table::fread(file.path(input_path, 'train', 
                                            paste0(layer, '.csv'))))
  test <- getm(data.table::fread(file.path(input_path, 'test', 
                                           paste0(layer, '.csv'))))
  
  pred_train <- getm(data.table::fread(file.path(output_path, 
                                                 paste0(prefix,'.train_decoded.', layer, '.csv'))))
  pred_test <- getm(data.table::fread(file.path(output_path, 
                                                paste0(prefix,'.test_decoded.', layer, '.csv'))))
  
  cors_train <- compute_cors(train, pred_train, ori = 1)
  cors_test <- compute_cors(test, pred_test, ori = 1)
  
  res <- data.table(rbind(data.frame('sample' = names(cors_train), 'cor' = as.numeric(cors_train),  'split' = 'train', 'N' = ncol(pred_train)),
        data.frame('sample' = names(cors_test), 'cor' = as.numeric(cors_test), 'split' = 'test', 'N' = ncol(pred_test))))
  return(res)
}

results <- do.call(rbind, lapply(c('gex', 'gex_pt', 'gex_pt_dp'), function(x) {
  message(date(), " => Analysing ",x,"\n")
  dt <- main(prefix = x, layer = 'crispr')
  dt$analysis <- x
  return(dt)
}))

labels <- list('gex' = 'gex_only', 
               'gex_pt' = 'gex\n+protein embeddings', 
               'gex_pt_dp' = 'gex\n+protein embeddings\n+describeProt')

results$group <- as.factor(ifelse(results$split == 'train', 
                        paste0("Training samples: N=",results[split == 'train']$N[1]), 
                        paste0("Test samples: N=",results[split == 'test']$N[1])))
results$group <- factor(results$group, levels = rev(levels(results$group)))

results$features <- as.character(labels[results$analysis])

p1 <- ggboxplot(results, x = 'features', y = 'cor', fill = 'features', 
          facet.by = 'group') + 
  theme(axis.text.x = element_text(angle = 15, hjust = 1, vjust = 1)) +
  #ggtitle(label = "Predicted vs Known Gene Dependency Scores",
   #       subtitle = paste("Pearson correlation distribution for",length(unique(results$sample)),"cell lines")) + 
  labs(y = 'Pearson Correlation') +
  theme(text = element_text(size = 12), legend.position = 'none')


stats <- do.call(rbind, lapply(c('gex', 'gex_pt', 'gex_pt_dp'), function(x) {
  message(date(), " => Analysing ",x,"\n")
  dt <- data.table::fread(file.path(output_path, 
                                    paste0(x, '.stats.csv')))
  dt$analysis <- x
  return(dt)
}))

stats$features <- as.character(labels[stats$analysis])
p2 <- ggbarplot(stats[metric == 'pearson_corr'], x = 'features', y = 'value', fill = 'features') +
  theme(text = element_text(size = 12),
        axis.text.x = element_text(angle = 15, hjust = 1), 
        legend.position = 'none') + 
  labs(y = "Pearson Correlation\n(predicted vs known hubness)")

p <- cowplot::plot_grid(p1, p2, labels = 'AUTO', rel_widths = c(2, 1.2))

ggsave(filename = 'depmap_gene_dependency_prediction.pdf', 
       plot = p, width = 210, height = 140, units = 'mm')

message(date()," => Finished making the figures")



