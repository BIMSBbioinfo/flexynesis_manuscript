# make figures about the fine-tuning section

args <- commandArgs(trailingOnly = T)

utils_script <- args[1]
projDir <- args[2]

source(utils_script)

library(ggplot2)
library(ggpubr)
library(data.table)
library(pheatmap)
ggplot2::theme_set(ggpubr::theme_pubclean())

get_stats <- function(workdir) {
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
  
  stats$finetuning <- ifelse(stats$finetuning_samples > 0, 'With Fine-tuning', 'No Fine-tuning')
  return(stats)
}

# get stats for finetuning experiment for drug response prediction
stats1 <- get_stats(file.path(projDir, 'finetuning_drug_response', 'output'))
stats2 <- get_stats(file.path(projDir, "finetuning_tcga_to_ccle", "output"))

get_plot <- function(dt, legend = F) {
  ggplot(dt, aes(x = factor(finetuning_samples), y = value, fill = tool)) + 
    geom_bar(stat = "identity", position = "dodge", show.legend = legend) +
    facet_wrap( ~ finetuning, scales = "free_x") +
    scale_fill_brewer(type = 'qual', palette = 6) +
    theme_minimal(base_size = 16) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
    labs(x = "Finetuning Samples")
}

plots <- lapply(unique(stats1$target), function(x) {
  get_plot(stats1[metric == 'pearson_corr'][target == x], legend = T) + 
    labs(y = 'Pearson Correlation', title = paste("CCLE -> GDSC",x))
})

plots[[8]] <- get_plot(stats2[metric == 'f1_score'], legend = T) + labs(y = 'F1 Score', title = 'TCGA -> CCLE; Cancer Type')

p1 <- cowplot::plot_grid(plots[[1]], plots[[8]], labels = 'AUTO', ncol = 1)# main figure 

p2 <- cowplot::plot_grid(plotlist = plots[2:7], labels = 'AUTO', ncol = 2)

ggsave(filename = 'finetuning.pdf', 
       plot = p1, width = 12, height = 9)

ggsave(filename = 'finetuning.supp.pdf', 
       plot = p2, width = 12, height = 12)

  



