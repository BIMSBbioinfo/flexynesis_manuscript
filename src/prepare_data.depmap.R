# Goal: Predict CRISPR gene dependency probabilities from gene expression and other omics data in CCLE cell lines 
# We use omics data + protein language model embeddings + describe prot features 
# In this exercise, genes are "samples", cell lines are features. We use protein-level embeddings from prot-trans and describe-prot

args <- commandArgs(trailingOnly = T)

utils_script <- args[1]
library(data.table)

source(utils_script)

# 1. ### get DEPMAP data ### 
message(date(), " => reading DEPMAP data")
# read gene expression data from depmap 
depmapDir <- '/data/local/buyar/arcas/collaborations/neuroblastoma/data/depmap'
gex <- data.table::fread(file.path(depmapDir, 'OmicsExpressionProteinCodingGenesTPMLogp1.csv'))
gex <- as.matrix(data.frame(gex[,-1], row.names = gex[[1]], check.names = F))
# cleanup gene names 
colnames(gex) <- sub("^(.+?) \\(.+$", "\\1", colnames(gex))

# read gene dependency probabilities from depmap 
# see https://depmap.org/portal/download/all/?releasename=DepMap+Public+22Q2&filename=CRISPR_gene_dependency.csv 
# Gene Dependency Probabilities represent the likelihood that knocking out the 
# gene has a cell growth inhibition or death effect. These probabilities are derived from the scores 
# in CRISPR_gene_effect.csv as described in https://doi.org/10.1101/720243.
crispr <- data.table::fread(file.path(depmapDir, 'CRISPRGeneDependency.csv'))
crispr <- as.matrix(data.frame(crispr[,-1], row.names = crispr[[1]], check.names = F))
# cleanup gene names 
colnames(crispr) <- sub("^(.+?) \\(.+$", "\\1", colnames(crispr))

# 2. ### get prot-trans protein embeddings data ### 
message(date(), " => reading prot-trans protein embeddings")
# now we want to add more gene-level features from language models and describeprot
# read protein-level prot-trans embeddings
dt <- data.table::fread('/data/local/buyar/arcas/protein_embeddings/prot_trans/embeddings/uniprot/embeddings.protein_level.csv')
# convert prot-trans uniprot accessions to gene names 
uniprot2genenames <- readRDS('/data/local/buyar/datasets/uniprot2hgnc.RDS')
uniprot2genenames$geneName <- sub("^sp.+\\|(.+?)_HUMAN.+$", "\\1", uniprot2genenames$geneName)
# to have unique mapping we do some processing
mdt <- melt.data.table(dt, id.vars = 'V1')
colnames(mdt)[1] <- 'E' #embedding features
mdt$geneName <- uniprot2genenames[match(mdt$variable, uniprotAccession)]$geneName
# get mean values per gene name
mdt <- mdt[,mean(value),by = c('E', 'geneName')]
# dcast
protTrans <- dcast.data.table(mdt, E ~ geneName, value.var = 'V1')
protTrans <- data.frame(protTrans[,-1], check.names = F, row.names = paste0('protTransE', protTrans$E))

# 3. ### get describe prot protein features data ### 
message(date(), " => reading describe prot protein features")
# get DescribeProt feautures 
# read describe prot feature scores 
dt <- data.table::fread('/data/local/buyar/datasets/describePROT/9606_value.csv', skip = 21)
# convert uniprot accessions to gene names
mdt <- melt.data.table(dt, id.vars = 'ACC', measure.vars = colnames(dt)[-(1:3)])
mdt$geneName <- uniprot2genenames[match(mdt$ACC, uniprotAccession)]$geneName
mdt <- mdt[!is.na(geneName),mean(as.numeric(value)), by = c('geneName', 'variable')]

describeProt <- dcast.data.table(mdt, variable ~ geneName, value.var = 'V1')
describeProt <- as.matrix(data.frame(describeProt[,-1], row.names = describeProt[[1]], check.names = F))


# 4. ### get protein interactome and compute hubness scores 
# we will use some network measures for the genes as "clin" annotation for genes. 
# the hubness score of the gene could be useful in predicting the gene dependency scores 
message(date(), " => reading stringDB ")
stringdb <- data.table::fread('/data/local/buyar/datasets/interactome/string/9606.protein.links.v12.0.txt')
#stringdb <- stringdb[]
aliases <- data.table::fread('/data/local/buyar/datasets/interactome/string/9606.protein.aliases.v12.0.txt')[source == 'Ensembl_HGNC']
stringdb$gene1 <- aliases[match(stringdb$protein1, aliases$`#string_protein_id`)]$alias
stringdb$gene2 <- aliases[match(stringdb$protein2, aliases$`#string_protein_id`)]$alias
stringdb <- stringdb[!is.na(gene1)][!is.na(gene2)][gene1 != gene2][combined_score > 200] #subset by confidence score 
# create a network and compute hubness scores 
message(date(), " => computing hubness scores")
G <- igraph::graph_from_data_frame(d = stringdb[,4:5], directed = T)
scores <- igraph::hub_score(G)
scores <- data.frame('hubness' = as.numeric(scores$vector), row.names = names(scores$vector))

# 5. ## Create Train/Test splits and save datasets 

dat <- list('clin' = t(scores)[1,,drop=F], 'protTrans' = protTrans, 'describeProt' = describeProt, 
            'crispr' = crispr, 'gex' = gex)
common_genes <- Reduce(intersect, lapply(dat, colnames))

dat <- sapply(simplify = F, dat, function(x) x[,common_genes,drop=F])
dat$clin <- t(dat$clin) # transpose clin data

# create train/test splits and print to folder 

dat_split <- split_dat(dat, ratio = 0.8)

message(date(), "=> printing dataset")
print_dataset(dat_split, "depmap")

message(date(), " => finished preparing depmap dataset")







