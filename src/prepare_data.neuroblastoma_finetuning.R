
# Download TARGET Neuroblastoma patient dataset and CCLE Neuroblastoma cell line data
# process and harmonize 

# Relevant Tasks 
# mapping cell lines to known neuroblastoma subtypes 
# use both control and confounding variables 

#args <- commandArgs(trailingOnly = TRUE)

# path to script that contains "CBioPortalData" class.
cbio_download_script = '/fast/AG_Akalin/buyar/flexynesis_manuscript_work/flexynesis_manuscript/src/get_cbioportal_data.R' #args[1]  #get_cbioportal_data.R
depmap_data_folder = '/data/local/buyar/arcas/collaborations/neuroblastoma/data/depmap' 
utils_script <- '/fast/AG_Akalin/buyar/flexynesis_manuscript_work/flexynesis_manuscript/src/utils.R'
source(utils_script)
source(cbio_download_script)

# 1. Download TARGET 2018 Neuroblastoma data 
# Instantiate the CBioPortalData class
cbio <- CBioPortalData$new(study_id = 'nbl_target_2018_pub') 
# Download and extract a study archive
archive_path <- cbio$download_study_archive() 
study_dir <- cbio$extract_archive(archive_path)
cbio$print_data_files()
files <- grep('_mrna_seq_rpkm.txt|_clinical|_cna|_mutations', cbio$data_files, value = T)
target <- cbio$read_data(files) 
names(target) <- gsub("data_|.txt", "", names(target))

gex <- target$mrna_seq_rpkm
clin <- target$clinical_sample
target_nbl <- list('gex' = log(gex+1), #convert to log scale
                   'clin' = clin[match(colnames(gex), SAMPLE_ID)])

# 2. Get depmap data
# import gex data from depmap (in log scale already)
gex <- data.table::fread(file.path(depmap_data_folder, 'OmicsExpressionProteinCodingGenesTPMLogp1.csv'))
# transpose, cleanup gene names, filter
gex <- t(as.matrix(data.frame(gex[,-1], row.names = gex[[1]], check.names = F)))
rownames(gex) <- gsub(" .+$", "", rownames(gex))
# get clin data for depmap
clin <- data.table::fread(file.path(depmap_data_folder, 'Model.csv'))
# subset depmap data only to neuroblastoma
s <- intersect(colnames(gex), clin[OncotreePrimaryDisease == 'Neuroblastoma']$ModelID)
gex <- gex[,s]
depmap <- list('gex' = gex,
               'clin' = clin[match(colnames(gex), ModelID)])

# gather some common and uncommon sample labels into a single object for all
# cell lines and TARGET samples.
# create mycn status labels
depmap$clin$mycn <- ifelse(depmap$clin[match(colnames(depmap$gex), ModelID)]$MolecularSubtype == 'MYCN_amp', 'MYCN_amp', 'other')
target_nbl$clin$mycn <- ifelse(target_nbl$clin[match(colnames(target_nbl$gex), SAMPLE_ID)]$MYCN == 'Amplified', 'MYCN_amp', 'other')

depmap$clin <- data.frame(depmap$clin[,-1], row.names = depmap$clin[[1]], check.names = F)
target_nbl$clin <- data.frame(target_nbl$clin[,-c(1:2)], row.names = target_nbl$clin$SAMPLE_ID, check.names = F)

print_dataset(list('train' = target_nbl, 
                   'test' = depmap), 'neuroblastoma')

