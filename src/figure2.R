# Make figures to showcase single-task training examples 

args <- commandArgs(trailingOnly = T)

utils_script <- args[1]
workdir <- args[2] # Single task anlaysis results folder 
pangiDir <- args[3] # directory with panGI - MSI analysis results

source(utils_script)

library(survminer)
library(ggplot2)
library(ggpubr)
library(data.table)
library(pheatmap)
library(pROC)
library(openxlsx)
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
plots <- lapply(preds, function(dt) {
  res <- cor.test(dt[split == 'test']$y, dt[split == 'test']$y_hat, method = "pearson")
  print(res$p.value)
  ggscatter(dt[split == 'test'], x = 'y', y = 'y_hat', add = 'reg.line', cor.coef = T, color = '66') + 
    labs(title = dt$var[1]) + theme(text = element_text(size = 14))
})
p1 <- cowplot::plot_grid(plotlist = plots, ncol = 2)

p1_source_data <- do.call(rbind, lapply(plots, function(x) x$data))

# 3. survival task 
pf <- stats[task == 'lgg_gbm'][metric == 'cindex'][order(value, decreasing = T)][target == '']$prefix[1]
dat <- get_data(pf, 'lgg_gbm', workdir)
dt <- dat$pred_labels[split == 'test']
# split samples by predicted survival risk scores into 2 main groups
dt$risk_group <- ifelse(dt$y_hat > median(dt$y_hat), 'high_risk', 'low_risk')
samples <- dt$sample_id
p3.1 <- plot_tsne(dat$E[samples,], dt$risk_group, show.labels = T) + 
  theme(legend.position = 'none') + 
  scale_color_brewer(type = 'qual', palette = 6)
# kaplan meier plot
df <- dat$colData[samples, c('OS_MONTHS', 'OS_STATUS')]
df$risk_group <- dt[match(rownames(df), sample_id)]$risk_group

p3_source_data <- cbind(p3.1$data, df)

df <- df[!is.na(df$OS_STATUS),]
surv_object <-  survival::Surv(time = df$OS_MONTHS, event = df$OS_STATUS)
fit <- surv_fit(surv_object ~ risk_group, df)
surv_diff <- survival::survdiff(surv_object ~ risk_group, data = df)
p_value <- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)
p3.2 <- ggsurvplot(fit, df, pval = p_value, risk.table = F, surv.median.line = 'hv')[['plot']] +
  scale_color_brewer(type = 'qual', palette = 6) + theme(legend.direction = 'vertical')
p3 <- cowplot::plot_grid(p3.1, p3.2, ncol = 2)

#2. classification task: MSI status prediction in pan-gastrointestinal cancers 
files <- dir(pangiDir, "stats.csv$", recursive = T)
stats <- do.call(rbind, pbapply::pblapply(files[!grepl('baseline', files)], function(x) {
  dt <- data.table::fread(file.path(pangiDir, x))
  dt$prefix <- gsub(".stats.csv", '', basename(x))
  return(dt)
}))
# combine with analysis table
analysis_table <- data.table::fread(file.path(pangiDir, 'analysis_table.csv'))[,-1]
# merge stats with analysis table
stats <- merge.data.table(analysis_table, stats, by = 'prefix')

voi <- 'msi' 
pf = stats[target == voi][data_types == 'gex,meth'][metric == 'average_auroc'][order(value, decreasing = T)]$prefix[1]
method_name <- stats[prefix == pf]$tool[1]

# get embeddings and prediction results 
E <- data.table::fread(file.path(pangiDir, 'results', paste0(pf, '.embeddings_test.csv')))
E <- as.matrix(data.frame(E[,-1], row.names = E[[1]], check.names = F))
pred <- data.table::fread(file.path(pangiDir, 'results', paste0(pf, '.predicted_labels.csv')))[split == 'test']
labels <- pred[match(rownames(E), sample_id)]$known_label
clin <- data.table::fread(file.path(pangiDir, 'data', 'ccle_vs_gdsc', 'test', 'clin.csv'))

p2.1 <- plot_tsne(E, clin[match(rownames(E), V1)]$msi, show.labels = T) + 
  theme(legend.position = 'top', text = element_text(size = 18)) + 
  scale_color_brewer(type = 'qual', palette = 6) 

# make a roc curve from predicted probabilities 
positive_class <- "MSI-High"
roc_dt <- pred[variable == "msi" & class_label == positive_class]
# Build binary true label (1 for MSI-High, 0 otherwise)
roc_dt[, y_true := as.integer(known_label == positive_class)]
roc_dt[, y_score := probability]
# Remove duplicates (keep one row per sample)
roc_dt <- unique(roc_dt, by = "sample_id")
roc_obj <- pROC::roc(roc_dt$y_true, roc_dt$y_score)

roc_df <- data.frame(
  fpr = 1 - roc_obj$specificities,
  tpr = roc_obj$sensitivities,
  thresholds = roc_obj$thresholds
)

p2.2 <- ggplot(roc_df, aes(x = fpr, y = tpr)) +
  geom_line(color = "blue") +
  geom_abline(slope = 1, intercept = 0, color = "grey") +
  labs(
    title = sprintf("ROC Curve for (MSI-High): AUC = %.3f", auc(roc_obj)),
    x = "False Positive Rate",
    y = "True Positive Rate"
  ) +
  theme_minimal(base_size = 18)

p2 <- cowplot::plot_grid(p2.1, p2.2, ncol = 2)

# also save the results from msi analysis to supplements 
table1 <- dcast.data.table(stats, tool + data_types + var ~ metric, value.var = 'value')[order(average_auroc, decreasing = T)]
wb <- createWorkbook()
addWorksheet(wb, "MSI_status_prediction")
writeData(wb, "MSI_status_prediction", table1)
saveWorkbook(wb, file = "SupplementaryTable2.xlsx", overwrite = TRUE)

# Combine all 
p <- cowplot::plot_grid(p1, p2, p3, labels = 'AUTO', ncol = 1, scale = 0.9)
ggsave(filename = 'Figure2.pdf', plot = p, width = 14, height = 14)
message(date(), "=> Finished making the plots")

source_data <- list('Figure2a' = p1_source_data, 
                    'Figure2b_tsne' = p2.1$data, 
                    'Figure2b_ROC' = p2.2$data, 
                    'Figure2c' = p3_source_data)
lapply(names(source_data), function(x) {
  write.csv(source_data[[x]], file = paste0(x, ".source_data.tsv"))
})








