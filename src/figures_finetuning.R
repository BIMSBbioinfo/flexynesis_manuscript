# make figures about the fine-tuning section

args <- commandArgs(trailingOnly = T)

utils_script <- args[1]
projDir <- args[2]

source(utils_script)

library(survminer)
library(ggplot2)
library(ggpubr)
library(data.table)
library(pheatmap)
ggplot2::theme_set(ggpubr::theme_pubclean())

workdir <- file.path(projDir, "finetuning_drug_response", "output")
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

stats$finetuning <- ifelse(stats$finetuning_samples > 0, 'with_finetuning', 'no_finetuning')

stats[metric == 'pearson_corr', .SD[which.max(value)],by = c('target', 'finetuning', 'tool')]

# finetuning can improve model performance on drug response
p1 <- ggplot(stats[metric == 'pearson_corr'][tool == 'supervised_vae'], 
       aes(x = finetuning_samples, y = value)) + 
  geom_bar(stat = 'identity', aes(fill = finetuning)) + facet_grid(~ target) +
  scale_fill_brewer(type = 'qual', palette = 6) + 
  theme(text = element_text(size = 16), axis.text.x = element_text(angle = 30)) + 
  labs(y = 'Pearson Correlation')

# finetuning makes it possible to predict mycn status (tumor samples => cell lines)
# import stats without finetuning 


workdir <- file.path(projDir, "finetuning_neuroblastoma")
# no finetuning
dt1 <- data.table::fread(file.path(workdir, 'neur.stats.csv'))
dt1$method <- 'supervised_vae'
dt1 <- rbind(dt1, data.table::fread(file.path(workdir, 'neur.baseline.stats.csv')))
dt1$finetuning <- 'no_finetuning'
# with finetuning 
dt2 <- data.table::fread(file.path(workdir, 'neur_finetuned.stats.csv'))
dt2$method <- 'supervised_vae'
dt2$finetuning <- 'with_finetuning' 
dt <- rbind(dt1, dt2)

dt$method2 <- paste(dt$method, dt$finetuning, sep = '\n')

p2 <- ggplot(dt[metric == 'f1_score'], aes(x = reorder(method2, -value), y = value)) + 
  geom_bar(stat = 'identity', aes(fill = value), position = 'dodge', width = 0.5) + 
  theme(text = element_text(size = 16), legend.position = 'none') +
  labs(y = 'F1 Score', x = 'Method') 


p <- cowplot::plot_grid(p1, 
                   cowplot::plot_grid(cowplot::ggdraw(), p2, cowplot::ggdraw(), nrow = 1,
                                      rel_widths = c(1,5,1)), 
                   ncol = 1, rel_heights = c(2, 1), 
                   labels = 'AUTO')

ggsave(filename = 'finetuning.pdf', 
       plot = p, width = 11, height = 9)


  



