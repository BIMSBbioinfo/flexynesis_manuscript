# Compile run times from different tools  on cpu and gpu

library(data.table)
library(ggplot2)
library(ggpubr)

ggplot2::theme_set(ggpubr::theme_pubclean())

args <- commandArgs(trailingOnly = T)
workdir <- args[1]

f <- list.files(workdir, pattern = '(cpu|gpu).log$', full.names = T)
f
dt <- do.call(rbind, lapply(f, function(x) {
  logs = readLines(x)
  hpo_time <- as.numeric(sub("^.+?HPO: (.+) sec", "\\1", grep('Time spent in HPO', logs, value = T)))
  data_import_time <- as.numeric(sub("^.+?data import: (.+) sec", "\\1", grep('Time spent in data import', logs, value = T)))
  
  model_params <- as.numeric(sub("Total params: (.+) M", "\\1", grep('Total params', logs, value = T)))

  tool <- strsplit(basename(x), "_")[[1]][1]
  resource <- sub('^.+?(cpu|gpu).+$', "\\1", x)
  fusion <- sub('^.+?(early|intermediate).+$', "\\1", x)

  peak_system_ram <- as.numeric(sub("^.+?HPO: (.+) MB", "\\1", grep('CPU RAM after HPO', logs, value = T)))
  peak_gpu_ram <- 0
  if (resource == 'gpu') {
    peak_gpu_ram <- as.numeric(sub("^.+?GPU RAM allocated: (.+) MB", "\\1", grep('GPU RAM allocated', logs, value = T)))
  }
  
  # time in seconds; RAM in MB; params in millions 
  return(data.table('tool' = tool, 'resource' = resource, 'fusion' = fusion, 'model_params' = model_params,
                    'hpo_time' = hpo_time, 'data_import_time' = data_import_time, 
                    'peak_cpu_ram' = peak_system_ram, 'peak_gpu_ram' = peak_gpu_ram))
}))
dt
dt[tool == 'supervised']$tool <- 'supervised_vae'

p1 <- ggplot(dt, aes(x = tool, y = hpo_time)) + geom_bar(aes(fill = fusion), 
                                               stat = 'identity', position = 'dodge', width = 0.5)+
  scale_fill_brewer(type = 'qual', palette = 3) + 
  theme(text = element_text(size = 14), axis.text.x = element_text(angle = 45, hjust = 0.8)) + 
  facet_grid(~ resource) + 
  labs(x = '', y = 'Time (seconds) per HPO step')
  
p2 <- ggplot(dt[resource == 'cpu'], aes(x = tool, y = peak_cpu_ram)) + geom_bar(aes(fill = fusion), 
                                                   stat = 'identity', position = 'dodge', width = 0.5)+
  scale_fill_brewer(type = 'qual', palette = 3) + 
  theme(text = element_text(size = 14), axis.text.x = element_text(angle = 45, hjust = 0.8)) + 
  labs(x = '', y = 'Peak CPU RAM Usage (in MB) ')

p3 <- ggplot(dt[resource == 'gpu'], aes(x = tool, y = peak_gpu_ram)) + geom_bar(aes(fill = fusion), 
                                                                          stat = 'identity', position = 'dodge', width = 0.5)+
  scale_fill_brewer(type = 'qual', palette = 3) + 
  theme(text = element_text(size = 14), axis.text.x = element_text(angle = 45, hjust = 0.8)) + 
  labs(x = '', y = 'Peak GPU RAM Usage (in MB) ')  

p <- cowplot::plot_grid(p1, p2, p3, ncol = 1, labels = 'AUTO')

ggsave(filename = 'SupplementaryFigure10.pdf', 
       plot = p, width = 12, height = 10)

openxlsx::write.xlsx(dt, file = 'SupplementaryTable7.xlsx')


