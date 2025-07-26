library(data.table)
library(ggplot2)
library(ggpubr)
library(pbapply)
library(DT)
library(knitr)
library(gridExtra)
library(openxlsx)
ggplot2::theme_set(ggpubr::theme_pubclean())

args <- commandArgs(trailingOnly = T)
workdir <- args[1]

files <- dir(workdir, "stats.csv$", recursive = T)
stats <- do.call(rbind, pbapply::pblapply(files[!grepl('baseline', files)], function(x) {
  dt <- data.table::fread(file.path(workdir, x))
  dt$prefix <- gsub(".stats.csv", '', basename(x))
  return(dt)
}))
# combine with analysis table
analysis_table <- data.table::fread(file.path(workdir, 'analysis_table.csv'))[,-1]
stats <- merge.data.table(analysis_table, stats, by = 'prefix')

# fusion type doesn't matter for single-modality data types 
stats[grep(',', data_types, invert = T)]$fusion <- 'early'

# for each task and target, scale the scores by maximum value in that task 
dt <- stats[metric %in% c('pearson_corr', 'cindex', 'f1_score')]
dt$group <- paste0(dt$task, "_", dt$var)
dt <- do.call(rbind, lapply(split(dt, dt$group), function(dt_sub) {
  dt_sub$score <- round((dt_sub$value / max(dt_sub$value)) * 100, 1)
  dt_sub$top <- ifelse(dt_sub$value > quantile(dt_sub$value, 1:100/100)[90], 
                       1, 0)
  dt_sub$ranking <- base::rank(-dt_sub$score, ties.method = 'first')
  return(dt_sub)
}))

dt[tool == 'RandomSurvivalForest']$tool <- 'RandomForest'
dt$omics <- ifelse(lengths(strsplit(dt$data_types, ",")) > 1, "multi", "single")
dt$learning <- ifelse(dt$tool %in% c('RandomForest', 'SVM'), 'classical', 'deep_learning')
dt$finetuning <- ifelse(dt$finetuning_samples > 0, 'with_finetuning', 'no_finetuning')

fig9_source_data <- list()

# main figure, top 10 per group
dt_sub <- dt[order(score, decreasing = T),.SD[1], by = c('group')]
p1 <- gridExtra::tableGrob(dt_sub[,c('task','var', 'tool', 'fusion', 'finetuning', 'data_types',
                                     'metric', 'value')]) 

fig9_source_data[['Figure9a']] <- dt_sub

# 1. is there a difference between deep learning vs off-the-shelf?
dt_sub <- dt[,.SD[which.max(score)],by = c('learning', 'group')]
p2 <- ggboxplot(dt_sub, 
          x = 'learning', y = 'score', add = 'jitter', color = 'learning') +
  scale_color_brewer(type = 'qual', palette = 6)+
  theme(legend.position = 'none') + stat_compare_means(label.y = 80)

fig9_source_data[['Figure9b']] <- p2$data

# 2. is there a difference between early/intermediate fusion?
dt_sub <- dt[omics == 'multi']
p3 <- ggboxplot(dt_sub[!tool %in% c('RandomForest', 'SVM', 'GNN'), .SD[which.max(score)], 
             by = c('fusion', 'group')], x = 'fusion', y = 'score', 
          add = 'jitter', color = 'fusion')  +
  scale_color_brewer(type = 'qual', palette = 6)+
  theme(legend.position = 'none') + stat_compare_means(label.y = 80)
fig9_source_data[['Figure9c']] <- p3$data


# 3. is there a difference between finetuning and no finetuning
dt_sub <- dt[learning == 'deep_learning']
p4 <- ggboxplot(dt_sub[,.SD[which.max(score)], by = c('finetuning', 'group')], 
          x = 'finetuning', y = 'score', add = 'jitter', color = 'finetuning')  +
  scale_color_brewer(type = 'qual', palette = 6)+
  theme(legend.position = 'none')  + stat_compare_means(label.y = 80)
fig9_source_data[['Figure9d']] <- p4$data

# 4. is there a difference between graph convolution methods?
p5 <- ggboxplot(dt[tool == 'GNN', .SD[which.max(score)], by = c('gnn_conv', 'group')], 
          x = 'gnn_conv', y = 'score', add = 'jitter', color = 'gnn_conv') + 
  scale_color_brewer(type = 'qual', palette = 6) +
  theme(legend.position = 'none')  + stat_compare_means(ref.group = 'GC', label.y = 50)
fig9_source_data[['Figure9e']] <- p5$data

# print figure source data 
lapply(names(fig9_source_data), function(x) {
  write.csv(fig9_source_data[[x]], file = paste0(x, ".source_data.tsv"))
})

p <- cowplot::plot_grid(p1, 
                   cowplot::plot_grid(p2, p3, p4, p5, nrow = 2, 
                                      labels = c('B', 'C', 'D', 'E')), 
                   ncol = 1, labels = c('A', ''), 
                   rel_heights = c(1, 1))

ggsave(filename = 'Figure9.pdf', 
       plot = p, width = 11, height = 14)

# save benchmark stats in table:
supptable <- dcast.data.table(stats, ... ~ metric, value.var = 'value')
supptable <- supptable[order(as.numeric(gsub("analysis", "", supptable$prefix)))]

# import bootstrap test results
bootstrap_stats <- data.table::fread(file.path(workdir, '..', 'comparisons.csv'))
wb <- createWorkbook()
addWorksheet(wb, "Benchmark_Results")
writeData(wb, "Benchmark_Results", supptable)
addWorksheet(wb, "Paired_Bootstrap_Stats")
writeData(wb, "Paired_Bootstrap_Stats", bootstrap_stats)
saveWorkbook(wb, file = "SupplementaryTable6.xlsx", overwrite = TRUE)



