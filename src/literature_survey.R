# Import Literature Survey of Multiomics Tools and make summary figures

library(openxlsx)
library(ggplot2)
library(ggpubr)
library(data.table)
library(scales)

ggplot2::theme_set(ggpubr::theme_pubclean())

tools <- openxlsx::read.xlsx('../tables/multiomics_tools.xlsx', sheet = 1)
colnames(tools)[2] <- 'code'

tools <- data.table(tools)
tools$Disease <- trimws(tools$Disease, which = 'both')

tools[code %in% c('Scripts', 'Notebooks')]$code <- 'Unpackaged Scripts or Notebooks'

p1 <- ggplot(tools[,length(Paper), by = c('code')], 
       aes(x = reorder(code, -V1), y = V1)) + 
  geom_bar(stat = 'identity', aes(fill = V1), show.legend = F) + 
  geom_text(aes(label = V1, y = V1+1)) +
  scale_fill_viridis_b() + 
  theme(axis.text.x = element_text(angle = 45)) + 
  scale_x_discrete(labels = scales::wrap_format(10)) + 
  labs(x = 'Code Availability', y = 'Number of Studies')

tools[,length(Paper),by = Disease]
tools[,length(Paper),by = Data]

ggsave(filename = '../figures/multiomics_tools_survey.pdf', plot = p1)


