# Make figures to showcase single-task training examples 

args <- commandArgs(trailingOnly = T)
# .libPaths('/fast/AG_Akalin/buyar/flexynesis_manuscript_work/flexynesis_manuscript/manuscript/site-library/')
utils_script <- args[1]
workdir <- args[2]

source(utils_script)

library(survminer)
library(ggplot2)
library(ggpubr)
library(data.table)
library(pheatmap)
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


# 1. drug response prediction 
# get best performing model 
# get single-task models using both rna+cnv for models trained on ccle, evaluated on gdsc
analyses <- stats[task == 'ccle_vs_gdsc'][target %in% c('Lapatinib', 'Selumetinib')][metric == 'pearson_corr'][grep(',', target, invert = T)][data_types == 'rna,cnv'][,.SD[which.max(value)],by = c('target')]

preds <- sapply(simplify = F, analyses$prefix, function(x) {
  data.table::fread(file.path(workdir, 'results', paste0(x, '.predicted_labels.csv')))
})

# drug response
p1 <- cowplot::plot_grid(plotlist = lapply(preds, function(dt) {
  ggscatter(dt[split == 'test'], x = 'y', y = 'y_hat', add = 'reg.line', cor.coef = T, color = '66') + 
    labs(title = dt$var[1]) + theme(text = element_text(size = 14))
}), ncol = 2)

# cell type label prediction 
voi <- 'celltype_l1' 
pf = stats[task == 'singlecell_bonemarrow'][target == voi][metric == 'kappa'][order(value, decreasing = T)]$prefix[1]
dat <- get_data(pf, 'singlecell_bonemarrow', workdir)
dt <- dat$pred_labels[split == 'test']
cells <- dt$sample_id
p2.1 <- plot_tsne(dat$E[cells,], dat$colData[cells,][[voi]], show.labels = T) + 
  theme(legend.position = 'none')
# make a contingency table 
m <- as.matrix(table(dt$y, dt$y_hat))
m <- apply(m, 2, function(x) round(x/sum(x) * 100, 1))
p2.2 <- ggplot(melt(data.table(m, keep.rownames = T), measure.vars = rownames(m)), 
         aes(x = rn, y = variable)) + geom_tile(aes(fill = value, color = value), show.legend = F) +
  geom_text(aes(label = value)) + scale_fill_gradient(low = 'white', high = '66') +
  theme(axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), 
        text = element_text(size = 12)) 
p2 <- cowplot::plot_grid(p2.1, p2.2, ncol = 2)

# 3. survival task 
pf <- stats[task == 'lgg_gbm'][metric == 'cindex'][order(value, decreasing = T)][target == '']$prefix[1]
dat <- get_data(pf, 'lgg_gbm', workdir)
dt <- dat$pred_labels[split == 'test']
# split samples by predicted survival risk scores into 2 main groups
dt$risk_group <- ifelse(dt$y_hat > median(dt$y_hat), 'high_risk', 'low_risk')
samples <- dt$sample_id
p3.1 <- plot_tsne(dat$E[samples,], dt$risk_group, show.labels = T) + 
  theme(legend.position = 'none') + 
  scale_color_manual(values = c(66, 45))
# kaplan meier plot
df <- dat$colData[samples, c('OS_MONTHS', 'OS_STATUS')]
df <- df[!is.na(df$OS_STATUS),]
df$risk_group <- dt[match(rownames(df), sample_id)]$risk_group
surv_object <-  survival::Surv(time = df$OS_MONTHS, event = df$OS_STATUS)
fit <- surv_fit(surv_object ~ risk_group, df)
p3.2 <- ggsurvplot(fit, df, pval = TRUE, risk.table = F, surv.median.line = 'hv')[['plot']] +
  scale_color_manual(values = c(66, 45)) + theme(legend.direction = 'vertical')

p3 <- cowplot::plot_grid(p3.1, p3.2, ncol = 2)

# Combine all 
p <- cowplot::plot_grid(p1, p2, p3, labels = 'AUTO', ncol = 1, scale = 0.9)
ggsave(filename = 'single_task_plots.pdf', plot = p, width = 210, height = 297, units = 'mm')
p <- cowplot::plot_grid(p1, p2, p3, labels = 'AUTO', ncol = 1, scale = 0.9)
ggsave(filename = 'single_task_plots.jpg', plot = p, width = 210, height = 210, units = 'mm', dpi = 300, bg = 'white')
message(date(), "=> Finished making the plots")


