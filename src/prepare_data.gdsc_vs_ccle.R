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

convert_ids2names <- function(M, ens2hgnc) {
  df <- data.frame(M, check.names = F)
  dim(df)
  df$names <- ens2hgnc[match(rownames(M), ref_gene_id)]$hgnc_symbol
  dim(df)
  # drop non-unique and NA 
  df <- df[!is.na(df$names),]
  dim(df)
  df <- df[!BiocGenerics::duplicated(df$names),]
  rownames(df) <- df$names
  df$names <- NULL
  return(as.matrix(df))
}

# to convert gene ids to names 
ens2hgnc <- readRDS(file.path(dataDir, 'ens2hgnc.RDS'))


message(date(), " => importing CCLE")
ccle <- readRDS(file.path(dataDir,  'CCLE.rds'))
ccle <- PharmacoGx::updateObject(ccle)
ccle_dat <- get_omics(ccle, c('rna', 'cnv', 'mutation'))
# convert gene ids to names 
ccle_dat$rna <- convert_ids2names(ccle_dat$rna, ens2hgnc)
lapply(ccle_dat, dim)

message(date(), " => importing GDSC2")
gdsc <- readRDS(file.path(dataDir, 'GDSC2.rds'))
gdsc <- PharmacoGx::updateObject(gdsc)
gdsc_dat <- get_omics(gdsc, c('rna', 'cnv', 'mutation'))
gdsc_dat$rna <- convert_ids2names(gdsc_dat$rna, ens2hgnc)
lapply(gdsc_dat, dim)

outdir <- file.path(dataDir, '..', 'prepared', 'gdsc_vs_ccle_test')
if(!dir.exists(outdir)) {
  dir.create(outdir)
}

message(date(), " => printing datasets")
print_dataset(list('train' = ccle_dat, 'test' = gdsc_dat), 
              outdir = outdir)


