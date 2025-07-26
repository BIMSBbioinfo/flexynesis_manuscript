# drug response markers 

args <- commandArgs(trailingOnly = T)

utils_script <- args[1]
workdir <- args[2]

source(utils_script)

library(ggplot2)
library(ggpubr)
library(data.table)
library(openxlsx)
ggplot2::theme_set(ggpubr::theme_pubclean())
data.table::setDTthreads(12)

files <- dir(file.path(workdir, 'output', 'results'), "stats.csv$", full.names = T)
files <- files[!grepl('baseline', files)]
stats <- do.call(rbind, pbapply::pblapply(files, function(x) {
  dt <- data.table::fread(x)
  dt$prefix <- gsub(".stats.csv", '', basename(x))
  return(dt)
}))

# combine with analysis table
analysis_table <- data.table::fread(file.path(workdir, 'output', 'analysis_table.csv'))[,-1]

# merge stats with analysis table
stats <- merge.data.table(analysis_table, stats, by = 'prefix')

#best <- stats[metric == 'pearson_corr',.SD[which.max(value)],by = c('target', 'tool')]
best <- stats[metric == 'pearson_corr',.SD[which.max(value)],by = target]


imp <- do.call(rbind, sapply(simplify = F, best$prefix, function(x) {
    dt1 <- data.table::fread(file.path(workdir, 'output', 'results', paste0(x, '.feature_importance.IntegratedGradients.csv')))
    dt2 <- data.table::fread(file.path(workdir, 'output', 'results', paste0(x, '.feature_importance.GradientShap.csv')))
    dt <- merge.data.table(dt1, dt2, by = c('target_variable', 'target_class', 'target_class_label', 'layer', 'name'))
    return(dt)
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
        text = element_text(size = 16)) + 
  labs(y = 'Best Pearson Correlation') + 
  scale_fill_gradient(high = 'darkred', low = 'blue')

p1
plots <- list()
plots[[1]] <- p1
plots <- c(plots, lapply(unique(imp$target_variable), function(x) {
  dt <- imp[target_variable == x][order(importance.x, decreasing = T)]
  dt$importance <- dt$importance.x/max(dt$importance.x)
  dt <- dt[1:10]
  # get genes from civic table
  dt$in_civic_db <- ifelse(dt$gene %in% unique(civic[drugs == x][evidence_direction == 'Supports']$gene),
                           'Known Target (CIVIC)', "")
  dt$label <- paste(dt$genename, dt$layer)

  ggplot(dt, aes(x = reorder(label, importance), y = importance)) + 
    geom_bar(stat = 'identity', aes(fill = layer), position = 'dodge', alpha = 0.7, show.legend = F) +
    geom_text(aes(label = in_civic_db), y = 0.01, hjust = 0, size = 3.5) + 
    scale_fill_brewer(type = 'qual', palette = 3) + 
    labs(title = x) + coord_flip() + theme(axis.title.y = element_blank())
}))

p <- cowplot::plot_grid(cowplot::plot_grid(cowplot::ggdraw(), plots[[1]], cowplot::ggdraw(), 
                                      nrow = 1, rel_widths  = c(1, 10, 1)), 
                   cowplot::plot_grid(plotlist = plots[2:9], nrow = 2),
                   ncol = 1, rel_heights = c(1, 2), labels = c('A', 'B'))

fig8b_source_data <- do.call(rbind, lapply(unique(imp$target_variable), function(x) {
  dt <- imp[target_variable == x][order(importance.x, decreasing = T)]
  dt$importance <- dt$importance.x/max(dt$importance.x)
  dt <- dt[1:10]
  # get genes from civic table
  dt$in_civic_db <- ifelse(dt$gene %in% unique(civic[drugs == x][evidence_direction == 'Supports']$gene),
                           'Known Target (CIVIC)', "")
  dt$label <- paste(dt$genename, dt$layer)
  return(dt)
}))

fig8_source_data <- list('Figure8a' = plots[[1]]$data, 
                         'Figure8b' = fig8b_source_data)
# print figure source data 
lapply(names(fig8_source_data), function(x) {
  write.csv(fig8_source_data[[x]], file = paste0(x, ".source_data.tsv"))
})

ggsave(filename = 'Figure8.pdf', 
       plot = p, width = 12, height = 12)

# make a plot about the correlation between integrated gradients and gradientshap
p <- ggscatter(imp, x = 'importance.x', y = 'importance.y', add = 'reg.line',  facet.by = 'target_variable', scales = 'free', 
          color = 'layer', alpha = 0.5, size = 2, cor.coef = T, cor.method = 'spearman', 
          add.params = list(color = "red", alpha = 0.2, size = 0.5)) + 
  theme(text = element_text(size = 14)) + 
  labs(x = 'IntegratedGradients', 
       y=  'GradientSHAP', title = 'Feature Attribution Score Correlation:\nIntegratedGradients vs GradientSHAP')

ggsave(filename = 'SupplementaryFigure9.pdf', 
       plot = p, width = 12, height = 10)


# save supplementary table 
table1 <- best
table2 <- imp[order(-importance.x), .SD[1:10], by = target_variable]

wb <- createWorkbook()
addWorksheet(wb, "SupplementaryTable_Fig8A")
writeData(wb, "SupplementaryTable_Fig8A", table1)
addWorksheet(wb, "SupplementaryTable_Fig8B")
writeData(wb, "SupplementaryTable_Fig8B", table2)
saveWorkbook(wb, file = "SupplementaryTable5.xlsx", overwrite = TRUE)







