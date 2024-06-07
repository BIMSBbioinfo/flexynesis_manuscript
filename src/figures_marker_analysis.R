# drug response markers 

args <- commandArgs(trailingOnly = T)

utils_script <- args[1]
workdir <- args[2]

source(utils_script)

library(ggplot2)
library(ggpubr)
library(data.table)
ggplot2::theme_set(ggpubr::theme_pubclean())
data.table::setDTthreads(12)

files <- dir(workdir, "stats.csv$", recursive = T)
stats <- do.call(rbind, pbapply::pblapply(files[!grepl('baseline', files)], function(x) {
  dt <- data.table::fread(file.path(workdir, x))
  dt$prefix <- gsub(".stats.csv", '', basename(x))
  return(dt)
}))

# combine with analysis table
analysis_table <- data.table::fread(file.path(workdir, 'analysis_table.csv'))[,-1]

# merge stats with analysis table
stats <- merge.data.table(analysis_table, stats, by = 'prefix')

#best <- stats[metric == 'pearson_corr',.SD[which.max(value)],by = c('target', 'tool')]
best <- stats[metric == 'pearson_corr',.SD[which.max(value)],by = var]




imp <- do.call(rbind, sapply(simplify = F, best$prefix, function(x) {
  f <- file.path(workdir, 'results', paste0(x, '.feature_importance.csv'))
  cat(f,"\n")
  if(file.exists(f)) {
    dt <- data.table::fread(f)
    return(dt)
  }
}))

ens2hgnc <- readRDS('/data/local/buyar/datasets/ens2hgnc.RDS')

imp$genename <- imp$name
imp[grep('ENSG', name)]$genename <- ens2hgnc[match(imp[grep('ENSG', name)]$name, ref_gene_id)]$hgnc_symbol

civic <- data.table::fread('/fast/AG_Akalin/buyar/flexynesis_manuscript_work/analyses/marker_analysis/01-Jan-2023-ClinicalEvidenceSummaries.tsv')

p1 <- ggplot(best, aes(x = reorder(var, value), y = value)) + 
  geom_bar(stat = 'identity', position = 'dodge', fill = 'red', alpha = 0.5) +
  geom_text(aes(label = paste(data_types, tool)), y = 0.01, hjust = 0, size = 4) + 
  theme(axis.title.y = element_blank()) + 
  labs(y = 'Best Pearson Correlation') + coord_flip()

plots <- list()
plots[[1]] <- p1

plots <- c(plots, lapply(unique(imp$target_variable), function(x) {
  dt <- imp[target_variable == x][order(importance, decreasing = T)]
  dt$importance <- dt$importance/max(dt$importance)
  dt <- dt[1:10]
  # get genes from civic table
  dt$in_civic_db <- ifelse(dt$gene %in% unique(civic[drugs == x][evidence_direction == 'Supports']$gene),
                           'civic', "")
  dt$label <- paste(dt$genename, dt$layer)

  ggplot(dt, aes(x = reorder(label, importance), y = importance)) + 
    geom_bar(stat = 'identity', aes(fill = layer), position = 'dodge', show.legend = F) +
    geom_text(aes(label = in_civic_db), y = 0.01, hjust = 0, size = 4) + 
    scale_fill_manual(values = list('rna' = 'lightgreen', 'mutation' = 'lightblue')) + 
    labs(title = x) + coord_flip() + theme(axis.title.y = element_blank())
}))

p <- cowplot::plot_grid(plotlist = plots, nrow = 3, labels = 'AUTO')

ggsave(filename = 'marker_analysis.pdf', 
       plot = p, width = 11, height = 8)


# # import drug response values and see if we can see different distributions in marker combinations 
# clin <- data.table::fread(file.path(workdir, 'data', 'ccle_vs_gdsc', 'train', 'clin.csv'))
# mut <- data.table::fread(file.path(workdir, 'data', 'ccle_vs_gdsc', 'train', 'mutation.csv'))
# 
# mut <- melt.data.table(mut, id.vars = 'V1')
# 
# 
# d <- 'Selumetinib'
# dt <- data.table('sample' = clin$V1, 'aac' = clin[[d]])
# im <- imp[target_variable == d][order(importance, decreasing = T)][layer == 'mutation'][1:2]
# dt <- cbind(dt, sapply(im$name, function(x) {mut[V1 == x][match(dt$sample, variable)]$value}))
# dt <- na.omit(dt)
# dt$comb <- apply(dt[,-c(1:2)], 1, function(x) {
#   cm <- paste(names(x), x, collapse = ' ')
#   cm <- gsub("1", '(mut)', cm)
#   cm <- gsub("0", "(wt)", cm)
#   return(cm)
# })
# ggboxplot(dt, x = 'comb', y = 'aac', add = 'jitter', color = 'comb') + 
#   labs(y = 'AAC: high AAC = high response', 
#        title = paste(d, " Top 2 Mutation Markers")) + coord_flip()  
# 
# dt$order <- tapply(dt$aac, dt$comb, median)
# ggplot(dt, aes(x = reorder(comb,), y = aac)) + 
#   geom_boxplot()
# 
# 





