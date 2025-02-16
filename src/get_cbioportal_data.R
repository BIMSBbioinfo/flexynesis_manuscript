# download and process data matrices and annotations from Cbioportal
library(httr)
library(utils)
library(readr)


# files_list: named list of files (e.g. list('gex' = 'data_mrna_seq_rpkm.txt', 'mut' = 'data_mutations.txt', 'clin' = 'data_clinical_sample.txt'))
get_cbioportal_dat <- function(study_id, files_list = NULL) {
  # Instantiate the CBioPortalData class
  cbio <- CBioPortalData$new(study_id = study_id) 
  # Download and extract a study archive
  archive_path <- cbio$download_study_archive() 
  study_dir <- cbio$extract_archive(archive_path)
  if(is.null(files_list)) {
    cbio$print_data_files()
    return()
  } 
  dat <- cbio$read_data(unlist(files_list))
  names(dat) <- names(files_list)
  
  if('clin' %in% names(dat)) {
    clin <- dat[['clin']]
    clin <- clin[!duplicated(clin[[1]])]
    clin <- data.frame(clin[,-1], row.names = clin[[1]], check.names = F)
    dat$clin <- clin
  }
  return(dat)
}

CBioPortalData <- R6::R6Class("CBioPortalData",
                              public = list(
                                base_url = "https://cbioportal-datahub.s3.amazonaws.com",
                                study_id = NULL,
                                data_files = NULL,
                                initialize = function(base_url = NULL, study_id = NULL, data_files = NULL) {
                                  if (!is.null(base_url)) {
                                    self$base_url <- base_url
                                  }
                                  if(!is.null(study_id)) {
                                    self$study_id <- study_id
                                  }
                                  self$data_files = data_files
                                },
                                download_study_archive = function() {
                                  url <- file.path(self$base_url, paste0(self$study_id, ".tar.gz"))
                                  dest_file <- paste0(self$study_id, ".tar.gz")
                                  if(!file.exists(dest_file)) {
                                    download.file(url, dest_file, mode = "wb")
                                  }
                                  return(dest_file)
                                },
                                
                                extract_archive = function(archive_path) {
                                  base = strsplit(archive_path, "\\.")[[1]][1]
                                  if(!dir.exists(base)) {
                                    untar(archive_path) 
                                  }
                                  self$data_files = dir(base, "^data_.*.txt")
                                  return(base)
                                }, 
                                
                                read_data = function(files = NULL) {
                                  if(is.null(files)) {
                                   files <- self$data_files 
                                  }
                                  cat("Will import data files",files,"\n")
                                  sapply(simplify = F, files, function(x) {
                                    cat(date(), "=> importing",file.path(self$study_id, x),"\n")
                                    dt <- data.table::fread(cmd = paste0("grep -v ^# ",file.path(self$study_id, x))) 
                                    if(grepl('mutations', x)) {
                                      cat("binarizing and converting to matrix", x, "\n")
                                      dt <- self$binarize_mutations(dt)
                                    } else if (!grepl('clinical|drug_treatment', x)) {
                                      cat("converting ",x," to matrix\n")
                                      dt <- self$process_dat(dt)
                                    } 
                                    return(dt)
                                  })
                                },
                                process_dat = function(dtc) {
                                  # exclude EntrezGeneID field 
                                  cols <- setdiff(colnames(dtc), c('Hugo_Symbol', 'Entrez_Gene_Id'))
                                  # remove non-unique rows 
                                  dtc <- dtc[!Hugo_Symbol %in% names(which(table(dtc[['Hugo_Symbol']]) > 1))]
                                  M <- as.matrix(data.frame(subset(dtc, select = cols), row.names = dtc[['Hugo_Symbol']], check.names = F))
                                  return(M)
                                },
                                binarize_mutations = function(dt) {
                                  # convert mutation data to binary matrix of genes vs samples 
                                  lapply(c("Hugo_Symbol", "Tumor_Sample_Barcode"), function(f) {
                                    if(!f %in% colnames(dt)){
                                      cat(colnames(dt),"\n")
                                      stop("Can't map mutations to sample ids.", f," not found")
                                    }
                                  })
                                  dt <- dt[,length(Variant_Classification), by = c('Hugo_Symbol', 'Tumor_Sample_Barcode')]
                                  dtc <- data.table::dcast.data.table(dt, Hugo_Symbol ~ Tumor_Sample_Barcode, value.var = 'V1')
                                  dtc <- self$process_dat(dtc) 
                                  dtc[is.na(dtc)] <- 0
                                  dtc[dtc > 0] <- 1
                                  return(dtc)
                                },
                                print_data_files = function() {
                                  print(knitr::kable(self$data_files))
                                }
                              )
)






