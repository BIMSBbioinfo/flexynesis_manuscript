# prepare a dataset for testing batch correction approaches
# use tcga+ccle dataset to see if we can predict cancer type well 

library(data.table)
data.table::setDTthreads(8)

args <- commandArgs(trailingOnly = TRUE)

src_dir <- args[1]
source(file.path(src_dir, 'get_cbioportal_data.R'))
source(file.path(src_dir, 'utils.R'))

subset_dat <- function(dat, samples) {
  sapply(simplify = F, names(dat), function(x) {
    if(x == 'clin') {
      dat[[x]][samples,]
    } else {
      dat[[x]][,samples]
    }
  })
}

find_common_samples <- function(dat) {
  Reduce(intersect, lapply(names(dat), function(x) {
    if(x == 'clin') {
      return(rownames(dat[[x]]))
    } else {
      return(colnames(dat[[x]]))
    }
  }))
}

combine_matrices <- function(l) {
  common_features <- Reduce(intersect, lapply(l, rownames))
  if(length(common_features) > 100000) {
    common_features <- sample(common_features, 100000)
  }
  combined <- do.call(cbind, lapply(l, function(x) {
    x[common_features,]
  }))
  return(combined)
}

# get CCLE 
dat_ccle <- get_cbioportal_dat('ccle_broad_2019', files_list = list('gex' = 'data_mrna_seq_rpkm.txt', 
                                                               'mut' = 'data_mutations.txt', 
                                                               'cna' = 'data_cna.txt',
                                                               'clin' = 'data_clinical_sample.txt'))
# further subset ccle 
samples <- rownames(dat_ccle$clin[dat_ccle$clin$CANCER_TYPE %in%  c('Breast Cancer', 'Colorectal Cancer', 'Glioma'),])
dat_ccle <- subset_dat(dat_ccle, intersect(find_common_samples(dat_ccle), samples))

# get COADREAD
dat_coad <- get_cbioportal_dat('coadread_tcga_pan_can_atlas_2018', 
                               files_list = list('clin' = 'data_clinical_sample.txt', 
                                                 'gex' = 'data_mrna_seq_v2_rsem.txt',
                                                 'mut' = 'data_mutations.txt',
                                                 'cna' = 'data_cna.txt'
                                                 ))
# to match sample ids
colnames(dat_coad$gex) <- gsub("-..$", "", colnames(dat_coad$gex))
colnames(dat_coad$mut) <- gsub("-..$", "", colnames(dat_coad$mut))
colnames(dat_coad$cna) <- gsub("-..$", "", colnames(dat_coad$cna))
dat_coad <- subset_dat(dat_coad, sample(find_common_samples(dat_coad), 100)) #downsample to 150 samples 

# get BRCA
dat_brca <- get_cbioportal_dat('brca_tcga_pan_can_atlas_2018', 
                               files_list = list('clin' = 'data_clinical_sample.txt', 
                                                 'gex' = 'data_mrna_seq_v2_rsem.txt',
                                                 'mut' = 'data_mutations.txt',
                                                 'cna' = 'data_cna.txt'
                                                 ))
# match sample ids 
colnames(dat_brca$gex) <- gsub("-..$", "", colnames(dat_brca$gex))
colnames(dat_brca$mut) <- gsub("-..$", "", colnames(dat_brca$mut))
colnames(dat_brca$cna) <- gsub("-..$", "", colnames(dat_brca$cna))
dat_brca <- subset_dat(dat_brca, sample(find_common_samples(dat_brca), 100))


# get GBM
dat_gbm <- get_cbioportal_dat('gbm_tcga_pan_can_atlas_2018', 
                              files_list = list('clin' = 'data_clinical_sample.txt', 
                                                'gex' = 'data_mrna_seq_v2_rsem.txt',
                                                'mut' = 'data_mutations.txt',
                                                'cna' = 'data_cna.txt'
                                                ))
colnames(dat_gbm$gex) <- gsub("-..$", "", colnames(dat_gbm$gex))
colnames(dat_gbm$mut) <- gsub("-..$", "", colnames(dat_gbm$mut))
colnames(dat_gbm$cna) <- gsub("-..$", "", colnames(dat_gbm$cna))
dat_gbm <- subset_dat(dat_gbm, sample(find_common_samples(dat_gbm), 100))

# get KIRC
dat_kirc <- get_cbioportal_dat('kirc_tcga_pan_can_atlas_2018', 
                               files_list = list('clin' = 'data_clinical_sample.txt', 
                                                 'gex' = 'data_mrna_seq_v2_rsem.txt',
                                                 'mut' = 'data_mutations.txt',
                                                 'cna' = 'data_cna.txt'
                               ))
colnames(dat_kirc$gex) <- gsub("-..$", "", colnames(dat_kirc$gex))
colnames(dat_kirc$mut) <- gsub("-..$", "", colnames(dat_kirc$mut))
colnames(dat_kirc$cna) <- gsub("-..$", "", colnames(dat_kirc$cna))
dat_kirc <- subset_dat(dat_kirc, sample(find_common_samples(dat_kirc), 100))
 

# COMBINE 

dat_list <- list('CCLE' = dat_ccle,  
                 'TCGA-COAD' = dat_coad, 
                 'TCGA-GBM' = dat_gbm, 
                 'TCGA-BRCA' = dat_brca) 
#                 'TCGA-KIRC' = dat_kirc)

dat <- sapply(c('cna', 'gex'), function(x) {
  combine_matrices(lapply(dat_list, function(d) {
    d[[x]]
  }))
})
samples <- find_common_samples(dat)
length(samples)
# get sample labels (cancer types)
clin <- do.call(rbind, lapply(names(dat_list), function(x) {
  df <- dat_list[[x]][['clin']]
  s <- intersect(samples, rownames(df))
  data.frame('cancertype' = df[s,]$CANCER_TYPE, 'dataset' = x, 'source' = gsub("-.+?$", "", x), row.names = s)
}))

clin[clin$cancertype %in% c('Glioblastoma', 'Glioma'),]$cancertype <- 'Glioma'

dat$clin <- clin

# create a covariate matrix and provide that as an additional modality 
# covariates <- as.matrix(model.matrix(~ source + dataset, data = dat$clin))
# dat$cov <- t(covariates)

lapply(dat, dim)

dat <- subset_dat(dat, find_common_samples(dat))

# create train/test split and save to folder 
dat_split <- split_dat(dat, ratio = 0.5)

table(dat_split$train$clin$cancertype, dat_split$train$clin$source)
table(dat_split$test$clin$cancertype, dat_split$test$clin$source)

print_dataset(dat = dat_split, outdir = 'ccle_vs_tcga')

# Also print a dataset where the training data is TCGA and test data is CCLE 
dat_tcga <- sapply(c('cna', 'gex'), function(x) {
  combine_matrices(lapply(dat_list[c('TCGA-GBM', 'TCGA-BRCA', 'TCGA-COAD')], function(d) {
    d[[x]]
  }))
})
samples <- find_common_samples(dat_tcga)
length(samples)
# get sample labels (cancer types)
clin <- do.call(rbind, lapply(names(dat_list)[2:4], function(x) {
  df <- dat_list[[x]][['clin']]
  s <- intersect(samples, rownames(df))
  data.frame('cancertype' = df[s,]$CANCER_TYPE, 'dataset' = x, 'source' = gsub("-.+?$", "", x), row.names = s)
}))
clin[clin$cancertype %in% c('Glioblastoma', 'Glioma'),]$cancertype <- 'Glioma'
dat_tcga$clin <- clin

dat_ccle$clin$cancertype <- dat_ccle$clin$CANCER_TYPE
dat_ccle$clin$source <- 'CCLE'

print_dataset(list('train' = dat_tcga, 'test' = dat_ccle), outdir = 'tcga_to_ccle')








