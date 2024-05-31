# Make figures to showcase single-task training examples 

args <- commandArgs(trailingOnly = T)

utils_script <- args[1]
workdir <- args[2]

source(utils_script)

library(ggplot2)
library(ggpubr)
library(data.table)
ggplot2::theme_set(ggpubr::theme_pubclean())

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

analyses <- stats[task == 'metabric'][tool == 'DirectPred'][data_types == 'gex,mut'][metric == 'kappa']

dat <- sapply(simplify = F, c('CLAUDIN_SUBTYPE', 'CHEMOTHERAPY', 'CLAUDIN_SUBTYPE,CHEMOTHERAPY'), function(x) {
  pf <- analyses[target == x]$prefix[1]
  get_data(pf, 'metabric', workdir)
})

# for each set of embeddings plot the embeddings colored by the two labels 

plots <- sapply(simplify = F, names(dat), function(run) {
  M <- dat[[run]]$E
  clin <- dat[[run]]$colData[rownames(M),]
  df <- plot_tsne(M, returnData = T) 
  df <- merge(df, clin, by = 'row.names')
  p1 <- ggplot(df, aes(x = tSNE1, y = tSNE2)) +
    geom_point(aes(color = CLAUDIN_SUBTYPE), size = 1, alpha = 0.4) +
    scale_color_brewer(type = 'qual', palette = 6)
  
  p2 <- ggplot(df, aes(x = tSNE1, y = tSNE2)) +
    geom_point(aes(color = CHEMOTHERAPY), size = 1, alpha = 0.4) +
    scale_color_brewer(type = 'qual', palette = 6)
  return(list('p1' = p1, 'p2' = p2))
})

combine_plots <- function(plotlist, labels) {
  legend1 <- get_legend(plotlist[[1]]$p1 + theme(legend.box.margin = margin(0, 0, 0, 12)))
  legend2 <- get_legend(plotlist[[1]]$p2 + theme(legend.box.margin = margin(0, 12, 0, 0)))
  legend = cowplot::plot_grid(legend1, legend2, ncol = 1)
  prow <- cowplot::plot_grid(
    plotlist = lapply(plotlist, function(x) {
      p1 <- x$p1 + theme(legend.position="none",
                         plot.title = element_blank(),
                         axis.title = element_text(size = 14))
      p2 <- x$p2 + theme(legend.position="none",
                         plot.title = element_blank(),
                         axis.title = element_text(size = 14))
      return(cowplot::plot_grid(p1, p2))
    }),
    labels = labels,
    nrow = 3, 
    scale = 0.9
  )
  return(cowplot::plot_grid(legend, prow, ncol = 1, rel_heights = c(1, 8)))
}

# reorder:
p <- combine_plots(plots, labels = 'AUTO')

# Combine all 
ggsave(filename = 'metabric_multitask_plot.pdf', plot = p, width = 210, height = 297, units = 'mm')
ggsave(filename = 'metabric_multitask_plot.jpg', plot = p, width = 210, height = 210, units = 'mm', dpi = 300, bg = 'white')


# 2. mixed multi task (reg+class+surv)

pf = stats[target == 'HISTOLOGICAL_DIAGNOSIS,AGE'][task == 'lgg_gbm'][metric == 'cindex'][order(value, decreasing = T)]$prefix[1]

dat <- get_data(pf, 'lgg_gbm', workdir)
dt <- dat$pred_labels[split == 'test']
samples <- unique(dt$sample_id)
E <- dat$E[samples,]
clin <- dat$colData[samples,]
#surv risk groups 
dt_surv <- dt[var == 'OS_STATUS'][match(samples, sample_id)]
dt_surv$y_hat <- as.numeric(dt_surv$y_hat)
dt_surv$risk_groups <- ifelse(dt_surv$y_hat < median(dt_surv$y_hat), 'low_risk', 'high_risk')
clin$risk_groups <- dt_surv[match(rownames(clin), sample_id)]$risk_groups
  
df <- plot_tsne(E, returnData = T)
df <- merge(df, clin, by = 'row.names')

p1 <- ggplot(df[!is.na(df$AGE),], aes(x = tSNE1, y = tSNE2)) +
  geom_point(aes(color = HISTOLOGICAL_DIAGNOSIS, size = AGE ), alpha = 0.4) +
  scale_color_brewer(type = 'qual', palette = 6) + facet_grid( ~ risk_groups) + 
  theme(legend.direction = 'vertical', text = element_text(size = 18)) 

message(date(), "=> Finished making the plots")

# top markers per outcome variable 
top_markers <- lapply(split(dat$Imp, dat$Imp$target_variable), function(x) { 
  unique(x[,.SD[which.max(importance)],by = c('layer', 'name')][order(importance, decreasing = T)])[1:10]
  })

p2 <- cowplot::plot_grid(plotlist = lapply(names(top_markers), function(x) {
  dt <- top_markers[[x]]
  ggplot(dt, aes(x = reorder(name, -importance), y = importance)) + 
    geom_bar(aes(fill = importance), stat = 'identity') + 
    scale_fill_gradient(low = 'gray', high = 'red') + 
    labs(x = 'gene', title = x) + 
    theme(text = element_text(size = 18), legend.position = 'none', 
          axis.text.x = element_text(angle = 45, hjust = 1))
}), nrow = 3)

p <- cowplot::plot_grid(p1, p2, labels = 'AUTO')

# Combine all 
ggsave(filename = 'lgg_gbm_multitask_plot.pdf', plot = p, width = 10, height = 8)
ggsave(filename = 'lgg_gbm_multitask_plot.jpg', plot = p, width = 10, height = 9, dpi = 300, bg = 'white')










