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

# main figure, top 10 per group
dt_sub <- dt[order(score, decreasing = T),.SD[1], by = c('group')]
p1 <- gridExtra::tableGrob(dt_sub[,c('task','var', 'tool', 'fusion', 'finetuning', 
                                     'metric', 'value')]) 

# dt_sub <- dt[order(score, decreasing = T),.SD[1:3], by = c('group')]
# p1.1 <- ggplot(dt_sub, aes(x = ranking, y = group)) + 
#   geom_tile(aes(fill = tool), width = 0.9, height = 0.9) + 
#   geom_text(aes(label = tool)) +
#   #paste0(tool, " ", fusion, "\n", finetuning))) + 
#   labs(y = 'Task - target variable') + 
#   theme(legend.position = 'none')

# 0. How is each method's best performance
# ggboxplot(dt[,.SD[which.max(score)],by = c('tool', 'group')], 
#          x = 'tool', y = 'score', add = 'jitter', color = 'tool') 

# 1. is there a difference between deep learning vs off-the-shelf?
dt_sub <- dt[,.SD[which.max(score)],by = c('learning', 'group')]
p2 <- ggboxplot(dt[,.SD[which.max(score)],by = c('learning', 'group')], 
          x = 'learning', y = 'score', add = 'jitter', color = 'learning') +
  scale_color_brewer(type = 'qual', palette = 6)+
  theme(legend.position = 'none') 

# 2. is there a difference between early/intermediate fusion?
dt_sub <- dt[omics == 'multi']
p3 <- ggboxplot(dt_sub[!tool %in% c('RandomForest', 'SVM', 'GNN'), .SD[which.max(score)], 
             by = c('fusion', 'group')], x = 'fusion', y = 'score', 
          add = 'jitter', color = 'fusion')  +
  scale_color_brewer(type = 'qual', palette = 6)+
  theme(legend.position = 'none') 

# 3. is there a difference between finetuning and no finetuning
dt_sub <- dt[learning == 'deep_learning']
p4 <- ggboxplot(dt_sub[,.SD[which.max(score)], by = c('finetuning', 'group')], 
          x = 'finetuning', y = 'score', add = 'jitter', color = 'finetuning')  +
  scale_color_brewer(type = 'qual', palette = 6)+
  theme(legend.position = 'none') 

# 4. is there a difference between graph convolution methods?
p5 <- ggboxplot(dt[tool == 'GNN', .SD[which.max(score)], by = c('gnn_conv', 'group')], 
          x = 'gnn_conv', y = 'score', add = 'jitter', color = 'gnn_conv') + 
  scale_color_brewer(type = 'qual', palette = 6) +
  theme(legend.position = 'none') 


p <- cowplot::plot_grid(p1, 
                   cowplot::plot_grid(p2, p3, p4, p5, nrow = 1, 
                                      labels = c('B', 'C', 'D', 'E')), 
                   ncol = 1, labels = c('A', ''))

ggsave(filename = 'benchmark_summary.pdf', 
       plot = p, width = 13, height = 9)


# save benchmark stats in table:
supptable <- dcast.data.table(stats, ... ~ metric, value.var = 'value')
supptable <- supptable[order(as.numeric(gsub("analysis", "", supptable$prefix)))]
openxlsx::write.xlsx(supptable, file = 'SupplementaryTable2.xlsx')




