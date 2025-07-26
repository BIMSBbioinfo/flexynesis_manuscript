# collate all source data files into Excel

library(openxlsx)

args <- commandArgs(trailingOnly = T)
workdir <- args[1]

files <- dir(workdir, 'Figure.*source_data.tsv')

wb <- createWorkbook()
for (x in files) {
  df <- read.csv(file.path(workdir, x))
  sheetname <- gsub(".source_data.tsv", "", x)
  addWorksheet(wb, sheetname)
  writeData(wb, sheetname, df)
}
saveWorkbook(wb, file = "Source_Data.xlsx", overwrite = TRUE)
