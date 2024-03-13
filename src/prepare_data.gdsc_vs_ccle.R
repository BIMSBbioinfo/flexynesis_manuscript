# https://bioconductor.org/packages/release/bioc/vignettes/PharmacoGx/inst/doc/PharmacoGx.pdf

args <- commandArgs(trailingOnly = T)
if(length(args) == 0) {
  stop("Provide path to data folder")
}
dataDir <- args[1] 
library(PharmacoGx)
library(data.table)
library(SummarizedExperiment)

get_omics <- function(pSet, layers) {
  dat_list <- sapply(simplify = F, layers, function(x) {
    message("importing ",x)
    obj <- PharmacoGx::summarizeMolecularProfiles(pSet,
                                                  mDataType = x,
                                                  cell.lines = cellNames(pSet),
                                                  summary.stat = ifelse(x == 'mutation', 'or', 'mean'),
                                                  verbose=FALSE, 
                                                  removeTreated = TRUE, 
                                                  fill.missing = FALSE)
    list('M' = SummarizedExperiment::assay(obj), 'colData' = data.frame(SummarizedExperiment::colData(obj), check.names = F))
  })
  # get tissueid for each and merge 
  dt <- unique(do.call(rbind, lapply(dat_list, function(x) {
    data.table(x[['colData']][,c('sampleid', 'tissueid')])
  })))
  colData <- data.frame(dt[,-1], row.names = dt[[1]], check.names = F)
  # also add drug response data
  drugs <- data.frame(t(summarizeSensitivityProfiles(
    pSet,
    sensitivity.measure='aac_recomputed',
    summary.stat="median",
    verbose=FALSE)), check.names = F)
  clin <- merge.data.frame(colData, drugs, by = 'row.names')
  clin <- data.frame(clin[,-1], row.names = clin[[1]], check.names = F)
  dat <- lapply(dat_list, function(x) x[['M']])
  dat[['clin']] <- clin
  return(dat)
}

print_dataset <- function(dat, outdir) {
  if(!dir.exists(outdir)) {
    dir.create(outdir)
  }
  lapply(names(dat), function(s) {
    p <- file.path(outdir, s)
    if(!dir.exists(p)) {
      dir.create(p)
    }
    lapply(names(dat[[s]]), function(f) {
      write.table(dat[[s]][[f]],
                  file = file.path(p, paste0(f, ".csv")), sep = ',')
    })
  })
}

message(date(), " => importing CCLE")
ccle <- readRDS(file.path(dataDir, 'CCLE.rds'))
ccle <- PharmacoGx::updateObject(ccle)
ccle_dat <- get_omics(ccle, c('rna', 'cnv'))

lapply(ccle_dat, dim)

message(date(), " => importing GDSC2")
gdsc <- readRDS(file.path(dataDir, 'GDSC2.rds'))
gdsc <- PharmacoGx::updateObject(gdsc)
gdsc_dat <- get_omics(gdsc, c('rna', 'cnv')) # ignoring mutation data because there is too few (68 genes; also with many NA values)

lapply(gdsc_dat, dim)

outdir <- file.path(dataDir, 'gdsc_vs_ccle')
if(!dir.exists(outdir)) {
  dir.create(outdir)
}

message(date(), " => printing datasets")
print_dataset(list('train' = ccle_dat, 'test' = gdsc_dat), 
              outdir = outdir)


