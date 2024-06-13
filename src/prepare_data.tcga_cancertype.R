# prepare data for unsupervised clustering experiment for tcga 

args = commandArgs(trailingOnly = T)

library(pbapply)

utils_script <- args[1]
dataDir <- args[2]

source(utils_script)

projects <- unlist(strsplit('TCGA-BLCA,TCGA-BRCA,TCGA-CESC,TCGA-COAD,TCGA-ESCA,TCGA-HNSC,TCGA-KIRC,TCGA-KIRP,TCGA-LGG,TCGA-LIHC,TCGA-LUAD,TCGA-LUSC,TCGA-PAAD,TCGA-PCPG,TCGA-PRAD,TCGA-SARC,TCGA-STAD,TCGA-TGCT,TCGA-THCA,TCGA-UCEC', 
                     ','))

# read gex+meth for each project

combine_matrices <- function(l) {
  common_features <- Reduce(intersect, lapply(l, rownames))
  if(length(common_features) > 100000) {
    common_features <- sample(common_features, 100000)
  }
  combined <- do.call(cbind, lapply(l, function(x) {
    x[common_features,]
  }))
  labels <- do.call(c, lapply(names(l), function(x) {
    rep(x, ncol(l[[x]]))
  }))
  return(list('dat' = combined, 'labels' = labels))
}

import_dat <- function(dataDir, projects, suffix, nodes = 3) {
  cl <- parallel::makeCluster(nodes)
  parallel::clusterExport(cl = cl, varlist = c('dataDir', 'projects'))
  dat <- pbapply::pbsapply(simplify = F, cl = cl, projects, function(x) {
    obj <- readRDS(file.path(dataDir, paste0(x, suffix)))
    M <- SummarizedExperiment::assay(obj)
    # update colnames to keep patient barcodes (not sample barcodes)
    colnames(M) <- sub("(^TCGA.{8}).+$", "\\1", colnames(M))
    return(M)
  })
  parallel::stopCluster(cl)
  return(dat)
}

gex <- import_dat(dataDir, projects, ".gex.fpkm.RDS", nodes = 10)
meth <- import_dat(dataDir, projects, ".meth.RDS", nodes = 10)

# find common samples in both and downsample 
for (x in projects) {
  s <- intersect(colnames(gex[[x]]), colnames(meth[[x]]))
  # downsample 
  s <- sample(s, 100)
  gex[[x]] <- gex[[x]][,s]
  meth[[x]] <- meth[[x]][,s]
}

gexC <- combine_matrices(gex)
methC <- combine_matrices(meth)

samples <- intersect(colnames(gexC$dat), colnames(methC$dat))

dat <- list('gex' = gexC$dat[,samples], 'meth' = methC$dat[,samples], 
            'clin' = data.frame('cohort' = gexC$labels[match(samples, colnames(gexC$dat))], 
                                row.names = samples))

dat_split <- split_dat(dat, ratio = 0.8)

print_dataset(dat_split, outdir = '../data/tcga_cancertype')



