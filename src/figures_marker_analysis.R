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
analysis_table <- data.table::fread(file.path(workdir, 'output', 'analysis_table.csv'))[,-1]

# merge stats with analysis table
stats <- merge.data.table(analysis_table, stats, by = 'prefix')

#best <- stats[metric == 'pearson_corr',.SD[which.max(value)],by = c('target', 'tool')]
best <- stats[metric == 'pearson_corr',.SD[which.max(value)],by = var]




imp <- do.call(rbind, sapply(simplify = F, best$prefix, function(x) {
  f <- file.path(workdir, 'output', 'results', paste0(x, '.feature_importance.IntegratedGradients.csv'))
  cat(f,"\n")
  if(file.exists(f)) {
    dt <- data.table::fread(f)
    return(dt)
  }
}))

ens2hgnc <- readRDS(file.path(workdir, 'ens2hgnc.RDS'))

imp$genename <- imp$name
imp[grep('ENSG', name)]$genename <- ens2hgnc[match(imp[grep('ENSG', name)]$name, ref_gene_id)]$hgnc_symbol

civic <- data.table::fread(file.path(workdir, '01-Jan-2023-ClinicalEvidenceSummaries.tsv'))

best$label <- paste0(best$data_types, "\n", best$tool)

p1 <- ggplot(best, aes(x = reorder(var, value), y = value)) + 
  geom_bar(stat = 'identity', position = 'dodge', aes(fill = value), alpha = 0.5, width = 0.5) +
  geom_text(aes(label = label, y = value + 0.05), size = 4) + 
  theme(axis.title.x = element_blank(), 
        text = element_text(size = 16), legend.position = 'none') + 
  labs(y = 'Best Pearson Correlation') + 
  scale_fill_gradient(high = 'darkred', low = 'blue')

p1
plots <- list()
plots[[1]] <- p1
plots <- c(plots, lapply(unique(imp$target_variable), function(x) {
  dt <- imp[target_variable == x][order(importance, decreasing = T)]
  dt$importance <- dt$importance/max(dt$importance)
  dt <- dt[1:10]
  # get genes from civic table
  dt$in_civic_db <- ifelse(dt$gene %in% unique(civic[drugs == x][evidence_direction == 'Supports']$gene),
                           'Known Target (CIVIC)', "")
  dt$label <- paste(dt$genename, dt$layer)

  ggplot(dt, aes(x = reorder(label, importance), y = importance)) + 
    geom_bar(stat = 'identity', aes(fill = layer), position = 'dodge', alpha = 0.7, show.legend = F) +
    geom_text(aes(label = in_civic_db), y = 0.01, hjust = 0, size = 3.5) + 
    scale_fill_brewer(type = 'qual', palette = 3) + 
    #scale_fill_manual(values = list('rna' = 'lightgreen', 'mutation' = 'lightblue')) + 
    labs(title = x) + coord_flip() + theme(axis.title.y = element_blank())
}))

p <- cowplot::plot_grid(cowplot::plot_grid(cowplot::ggdraw(), plots[[1]], cowplot::ggdraw(), 
                                      nrow = 1, rel_widths  = c(1, 10, 1)), 
                   cowplot::plot_grid(plotlist = plots[2:9], nrow = 2),
                   ncol = 1, rel_heights = c(1, 2), labels = c('A', 'B'))

ggsave(filename = 'marker_analysis.pdf', 
       plot = p, width = 12, height = 12)





