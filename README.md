# flexynesis_manuscript
All publication material relevant for the manuscript describing the flexynesis software package  

# Project Folder

Accessible from Hulk/Beast/Max: `/fast/AG_Akalin/buyar/flexynesis_manuscript_work/` 



./raw folder contains original dataset downloaded from a source such as Cbioportal/TCGA/PharmacoGx/DepMAP
./prepared folder contains data prepared as input to flexynesis.


## Datasets used in the manuscript 

Below is a description of the datasets used in the manuscript and how to prepare them for analysis with flexynesis

### Downloaded Datasets 

Go to `/fast/AG_Akalin/buyar/flexynesis_manuscript_work/datasets`: 

./raw folder contains original dataset downloaded from a source such as Cbioportal/TCGA/PharmacoGx/DepMAP
./prepared folder contains data prepared as input to flexynesis.

./raw:

- CCLE.rds: downloaded from https://zenodo.org/record/3905462/files/CCLE.rds?download=1
- GDSC2.rds: downloaded from https://zenodo.org/record/3905481/files/GDSC2.rds?download=1
- lgggbm_tcga_pub.tar.gz: downloaded from https://www.cbioportal.org/study/summary?id=lgggbm_tcga_pub
- brca_metabric.tar.gz: downloaded from https://www.cbioportal.org/study/summary?id=brca_metabric
- depmap: downloaded from depmap portal https://depmap.org/portal/data_page/?tab=allData
- nbl_target_2018_pub.tar.gz: downloaded from cbioportal https://www.cbioportal.org/study/summary?id=nbl_target_2018_pub
- GDCData: TCGA cohort datasets for 33 cancer types downloaded using the TCGABiolinks package (see https://github.com/BIMSBbioinfo/uyar_et_al_multiomics_deeplearning)
- prot-trans: protein sequence embeddings obtained from prot-trans-xl-uniref50 model on uniprot sequences
- describeProt: protein level sequence/structure/function features from describeprot database (http://biomine.cs.vcu.edu/servers/DESCRIBEPROT/download_database_value/9606_value.csv)

### PREPARED datasets used as input to flexynesis

./prepared: (prepared data folder and which script+command was used to prepare the data)

- ccle_vs_gdsc: Drug response data from cell lines from CCLE and GDSC2 datasets
/opt/R/4.2/bin/Rscript /fast/AG_Akalin/buyar/flexynesis_manuscript_work/flexynesis_manuscript/src/prepare_data.gdsc_vs_ccle.R raw/

- lgggbm_tcga_pub_processed: Merged cohorts of LGG + GBM samples
/opt/R/4.2/bin/Rscript ../flexynesis_manuscript/src/prepare_data.LGG_GBM.R ../flexynesis_manuscript/src/get_cbioportal_data.R

- brca_metabric_processed: METABRIC dataset
/opt/R/4.2/bin/Rscript ../flexynesis_manuscript/src/prepare_data.metabric.R ../flexynesis_manuscript/src/get_cbioportal_data.R

- single_cell_bonemarrow: CITE-Seq dataset from Seurat (https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis)
/opt/R/4.2/bin/Rscript ../flexynesis_manuscript/src/prepare_data.cite_seq.R

- neuroblastoma_target_vs_depmap: neuroblastoma patient samples (TARGET study) and cell lines (depmap)
/opt/R/4.2/bin/Rscript ../flexynesis_manuscript/src/prepare_data.neuroblastoma_finetuning.R ../flexynesis_manuscript/src/get_cbioportal_data.R ./raw/depmap/ ../flexynesis_manuscript/src/utils.R

- tcga_cancertype: TCGA cancer cohort for ~21 cancer types 100 samples per each cohort
/opt/R/4.2/bin/Rscript ../flexynesis_manuscript/src/prepare_data.tcga_cancertype.R ../flexynesis_manuscript/src/utils.R ./raw/GDCdata

- depmap_gene_dependency: Dataset for gene-dependency prediction in cell lines. Consists of depmap gene expression + prottrans embeddings + describeprot features
/opt/R/4.2/bin/Rscript ../flexynesis_manuscript/src/prepare_data.depmap.R ../flexynesis_manuscript/src/utils.R ./raw/depmap/ ./raw/prot-trans/embeddings.protein_level.csv ./raw/uniprot2hgnc.RDS ./raw/describePROT/9606_value.csv

# Environment

## Install flexynesis 

```
git clone https://github.com/BIMSBbioinfo/flexynesis.git
cd flexynesis
# create conda env
conda create -n flexynesis --file spec-file.txt
conda activate flexynesis
# install flexynesis
pip install -e .
```

## Install other packages 

```
guix package --manifest=guix.scm --profile=./manuscript
```

## Activate environment

```
source ./manuscript/etc/profile
conda activate flexynesis
```


# Manuscript

- [Link to google doc containing manuscript text](https://docs.google.com/document/d/10Slme7TQLll7FEBAOtuug-y7ybUxnC1M-QuYMJ_lnUE/edit?usp=sharing)
- [Link to flexynesis package repo](https://github.com/BIMSBbioinfo/flexynesis)
- [Link to flexynesis-benchmarks package repo](https://github.com/BIMSBbioinfo/flexynesis-benchmarks)

- [Link to gdrive for graphical abstract afdesign](https://drive.google.com/file/d/1-R8KrQTxgo9ocdqsliEd8NC7Ntb5hp7I/view?usp=drive_link)

