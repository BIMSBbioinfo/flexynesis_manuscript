# combine prot-trans embeddings withe describePROT features to make predictions about 
# protein function/localisation

utils_script <- args[1]

source(utils_script)

# read protein-level prot-trans embeddings
protTrans <- data.table::fread('/data/local/buyar/arcas/protein_embeddings/prot_trans/embeddings/uniprot/embeddings.protein_level.csv')
PT <- data.frame(protTrans[,-1], check.names = F, row.names = paste0('protTransE', protTrans$V1))


# read describe prot feature scores 
describeProt <- data.table::fread('/data/local/buyar/datasets/describePROT/9606_value.csv', skip = 21)

DP <- data.frame(t(data.frame(describeProt[,-c(1,2,3)], check.names = F, row.names = describeProt$ACC)))


# read subcellular localisation info from deeploc
deeploc <- data.table::fread('/data/local/buyar/collaborations/trendelina/data/deeploc/deeploc.tsv')


proteins <- Reduce(intersect, list('pt' = colnames(PT), 'dp' = colnames(DP), 'dl' = deeploc$uniprotAccession))

dat <- list('protTrans' = PT[,proteins], 'describeProt' = DP[,proteins], 
            'clin' = data.frame('localisation' = deeploc[match(proteins, uniprotAccession)]$localisation,
                                'localisation_simple' = deeploc[match(proteins, uniprotAccession)]$localisation.simple,
                                row.names = proteins))

# create train/test splits 
dat_split <- split_dat(dat)

# print dataset 
print_dataset(dat_split, 'protein_localisation')









