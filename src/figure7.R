# make figures about the fine-tuning section

args <- commandArgs(trailingOnly = T)

utils_script <- args[1]
projDir <- args[2]

source(utils_script)

library(ggplot2)
library(ggpubr)
library(data.table)
library(pheatmap)
library(openxlsx)
ggplot2::theme_set(ggpubr::theme_pubclean())

get_pval_annotations <- function(dt, y_buffer = 0.05) {
  ref <- dt$Method[dt$performance_vs_reference == "reference"]
  cmp <- dt[dt$Method != ref, ]
  y_pos <- seq(max(dt$ci_upper) + y_buffer, length.out = nrow(cmp), by = y_buffer)
  stars <- cut(cmp$p_value,
               breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
               labels = c("***", "**", "*", "ns"))
  data.frame(
    group1 = ref,
    group2 = cmp$Method,
    y.position = y_pos,
    p.value = cmp$p_value,
    label = stars,
    stringsAsFactors = FALSE
  )
}

bootstrap_stats <- rbind(
  data.table::fread(file.path(projDir, 'tcga_to_ccle.bootstrap_stats.csv')),
  data.table::fread(file.path(projDir, 'ccle_vs_gdsc.bootstrap_stats.csv'))
)
bootstrap_stats$Method <- gsub("supervised", "SupervisedVAE", bootstrap_stats$Method)
bootstrap_stats$finetuning <- as.factor(ifelse(grepl('finetuned', bootstrap_stats$Method), 'With Fine Tuning', 'No Fine Tuning'))
bootstrap_stats <- split(bootstrap_stats, bootstrap_stats$target)

# Save summary stats as supplementary table 
wb <- createWorkbook()
addWorksheet(wb, "TCGA_vs_CCLE_CancerType")
writeData(wb, "TCGA_vs_CCLE_CancerType", bootstrap_stats$cancertype)
addWorksheet(wb, "CCLE_vs_GDSC_DrugResponse")
writeData(wb, "CCLE_vs_GDSC_DrugResponse", do.call(rbind, bootstrap_stats[2:9]))
saveWorkbook(wb, file = "SupplementaryTable4.xlsx", overwrite = TRUE)

# add bootstrap raw stats 
bootstrap_raw_stats <- sapply(simplify = F, dir(projDir, 'bootstrap_raw.csv$'), function(x) {
  dt <- data.table::fread(file.path(projDir, x))
  mdt <- melt.data.table(dt,id.vars = c('target', 'metric'))
  colnames(mdt)[3] <- 'Method'
  mdt$Method <- gsub("supervised", "SupervisedVAE", mdt$Method)
  mdt$finetuning <- as.factor(ifelse(grepl('finetuned', mdt$Method), 'With Fine Tuning', 'No Fine Tuning'))
  return(mdt)
})
names(bootstrap_raw_stats) <- unlist(lapply(strsplit(names(bootstrap_raw_stats), "\\."), 
                                            function(x) x[2]))

plots <- sapply(simplify = F, names(bootstrap_raw_stats), function(x) {
  dt <- bootstrap_raw_stats[[x]]
  summary_stats <- dt[, .(
    mean = mean(value),
    lower = quantile(value, 0.025),
    upper = quantile(value, 0.975)
  ), by = .(Method, finetuning)]
  p <- ggviolin(dt, x = 'Method', y = 'value', color = 'finetuning', add = NULL) +
    geom_jitter(aes(color = finetuning), alpha = 0.2, width = 0.25, size = 1) +
    
    # Overlay bootstrap CIs using precomputed summary
    geom_pointrange(data = summary_stats,
                    aes(x = Method, y = mean, ymin = lower, ymax = upper, group = finetuning),
                    position = position_dodge(width = 0.75),
                    color = "black", size = 0.4) +
    geom_hline(yintercept = summary_stats[1,]$mean, linetype = 'dashed') + 
    labs(y = paste0(dt$metric[1], ' (', dt$target[1], ')')) + 
    theme(text = element_text(size = 14), 
          legend.title = element_blank(),
          axis.text.x = element_text(angle = 30, hjust = 0.7),
          axis.title.x = element_blank()) + 
    scale_color_brewer(type = 'qual', palette = 6)
  return(p)
})
# Make a main plot 
p <- cowplot::plot_grid(plots$Selumetinib, plots$cancertype, labels = 'AUTO', 
                        ncol = 1)

ggsave(filename = 'Figure7.pdf', 
       plot = p, width = 10, height = 12)

# save figure source data 
write.csv(do.call(rbind, bootstrap_raw_stats), 
          file = "Figure7.source_data.tsv")

  



