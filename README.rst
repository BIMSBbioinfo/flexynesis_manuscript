.. image:: https://zenodo.org/badge/933743968.svg
  :target: https://doi.org/10.5281/zenodo.16443730

==============================
Flexynesis manuscript material
==============================

Publication material relevant for the manuscript describing the flexynesis software package. 

Our manuscript currently available at `BioRxiv <https://biorxiv.org/cgi/content/short/2024.07.16.603606v1>`_. 

See our github repository of `Flexynesis <https://github.com/BIMSBbioinfo/flexynesis>`_. 

Datasets used in the manuscript
===============================

Publicly available datasets 
---------------------------

* **CCLE.rds**: downloaded from `Zenodo <https://zenodo.org/record/3905462/files/CCLE.rds?download=1>`_.
* **GDSC2.rds**: downloaded from `Zenodo <https://zenodo.org/record/3905481/files/GDSC2.rds?download=1>`_.
* **lgggbm_tcga_pub.tar.gz**: downloaded from `cbioportal <https://www.cbioportal.org/study/summary?id=lgggbm_tcga_pub>`_.
* **brca_metabric.tar.gz**: downloaded from `cbioportal <https://www.cbioportal.org/study/summary?id=brca_metabric>`_.
* **depmap**: downloaded from `depmap portal <https://depmap.org/portal/data_page/?tab=allData>`_.
* **nbl_target_2018_pub.tar.gz**: downloaded from `cbioportal <https://www.cbioportal.org/study/summary?id=nbl_target_2018_pub>`_.
* **GDCData**: TCGA cohort datasets for 33 cancer types downloaded using the TCGABiolinks package (`See GitHub <https://github.com/BIMSBbioinfo/uyar_et_al_multiomics_deeplearning>`_).
* **prot-trans**: protein sequence embeddings obtained from prot-trans-xl-uniref50 model on uniprot sequences.
* **describeProt**: protein level sequence/structure/function features from describeprot database (`Download here <http://biomine.cs.vcu.edu/servers/DESCRIBEPROT/download_database_value/9606_value.csv>`_).
* **coadread_tcga_pan_can_atlas_2018.tar.gz**: downloaded from `cbioportal <https://www.cbioportal.org/study/summary?id=coadread_tcga_pan_can_atlas_2018>`_.
* **brca_tcga_pan_can_atlas_2018.tar.gz**: downloaded from `cbioportal <https://www.cbioportal.org/study/summary?id=brca_tcga_pan_can_atlas_2018>`_.
* **gbm_tcga_pan_can_atlas_2018.tar.gz**: downloaded from `cbioportal <https://www.cbioportal.org/study/summary?id=gbm_tcga_pan_can_atlas_2018>`_.


PREPARED datasets used as input to Flexynesis
---------------------------------------------

The datasets listed above were further processed to create train/test splits for training using Flexynesis. 
The prepared datasets can be downloaded from here: 
https://bimsbstatic.mdc-berlin.de/akalin/buyar/flexynesis_manuscript_material/datasets_prepared.tgz

The ``./prepared`` folder contains:

* **ccle_vs_gdsc**: Drug response data from cell lines from CCLE and GDSC2 datasets.

* **lgggbm_tcga_pub_processed**: Merged cohorts of LGG + GBM samples.

* **brca_metabric_processed**: METABRIC dataset processed.

* **single_cell_bonemarrow**: CITE-Seq dataset from Seurat.

* **tcga_vs_ccle**: TCGA tumors and CCLE cell lines from 3 different cancer types: lung cancer, glioma, and breast cancer 

* **tcga_cancertype**: TCGA cancer cohort for ~21 cancer types 100 samples per each cohort.

* **depmap_gene_dependency**: Dataset for gene-dependency prediction in cell lines. Consists of depmap gene expression + prottrans embeddings + describeprot features.

* **panGI_msi**: Gene expression and promoter methylation data from 7 different TCGA cohorts (gastrointestinal and gynocological cancers) with microsatellite instability (MSI) annotations: TCGA-COAD (Colon Adenocarcinoma), TCGA-ESCA (Esophageal Carcinoma), TCGA-PAAD (Pancreatic Adenocarcinoma), TCGA-READ (Rectum Adenocarcinoma), TCGA-STAD (Stomach Adenocarcinoma), TCGA-UCEC (Uterine Corpus Endometrial Carcinoma), TCGA-UCS (Uterine Carcinosarcoma). 

Flexynesis output for use-cases
-------------------------------

For the different use-cases described in the manuscript, Flexynesis output (along with the configurations used) 
can be downloaded from here: 
https://bimsbstatic.mdc-berlin.de/akalin/buyar/flexynesis_manuscript_material/manuscript_processed_data.tgz

Environment
===========

Clone the manuscript repo:
--------------------------

.. code-block:: bash 

    git clone https://github.com/BIMSBbioinfo/flexynesis_manuscript.git


Install flexynesis
-------------------

.. code-block:: bash

    mamba create -n flexynesisenv python==3.11 snakemake 
    mamba activate flexynesisenv
    pip install flexynesis 

Install other packages
----------------------

.. code-block:: bash

    guix package --manifest=guix.scm --profile=./manuscript

Activate environment
--------------------

.. code-block:: bash

    source ./manuscript/etc/profile
    mamba activate flexynesisenv


Figures
=======

Assuming the prepared datasets and Flexynesis output files are downloaded from the following locations:

- Datasets: https://bimsbstatic.mdc-berlin.de/akalin/buyar/flexynesis_manuscript_material/datasets_prepared.tgz
- Flexynesis output: https://bimsbstatic.mdc-berlin.de/akalin/buyar/flexynesis_manuscript_material/manuscript_processed_data.tgz

The figures in the manuscript can be reproduced using the following instructions: 

Unzip the Flexynesis datasets and output folders:

.. code-block:: bash

    tar -xzvf manuscript_processed_data.tgz
    tar -xzvf datasets_prepare.tgz 

Activate guix environment: 

.. code-block:: bash

    source ./flexynesis_manuscript/manuscript/etc/profile 

Change to folder with Flexynesis output data

.. code-block:: bash

    cd manuscript_processed_data

Figure 1: Flexynesis workflow 
-----------------------------
This figure was manually made.


Figure 2: single-task figures
-----------------------------

.. code-block:: bash

   Rscript ../flexynesis_manuscript/src/figure2.R ../flexynesis_manuscript/src/utils.R single_multi_experiments panGI_MSI_analysis/output

Figures 3 and 4: multi-task figures
-----------------------------------

.. code-block:: bash

   Rscript ../flexynesis_manuscript/src/figure3_and_4.R ../flexynesis_manuscript/src/utils.R single_multi_experiments

Figure 5: unsupervised clustering (TCGA cancer types)
-----------------------------------------------------

.. code-block:: bash 

   Rscript ../flexynesis_manuscript/src/figure5.R ../flexynesis_manuscript/src/utils.R ./unsupervised_cancertype/

Figure 6: cross-modality prediction of cell line dependency probabilities 
-------------------------------------------------------------------------

.. code-block:: bash 

   Rscript ../flexynesis_manuscript/src/figure6.R ../datasets_prepared/depmap_gene_dependency/ depmap_analysis/output/

Figure 7: demonstration of fine-tuning
--------------------------------------

.. code-block:: bash

   Rscript ../flexynesis_manuscript/src/figure7.R ../flexynesis_manuscript/src/utils.R finetuning/


Figure 8 and Supp. Figure 9: marker analysis 
--------------------------------------------

.. code-block:: bash 

   Rscript ../flexynesis_manuscript/src/figure8.R ../flexynesis_manuscript/src/utils.R marker_analysis/

Figure 9: benchmark summary
---------------------------

.. code-block:: bash

  Rscript ../flexynesis_manuscript/src/figure9.R benchmarks/output 

Supp Figure 10: Run Times & Resources 
----------------------------------

.. code-block:: bash

  Rscript ../flexynesis_manuscript/src/supp_figure_10.R runtimes/output

Supp Figure: Run Times & Resources 
----------------------------------

.. code-block:: bash

  Rscript ../flexynesis_manuscript/src/figures_runtimes.R runtimes/output

Collate all figure source data files
------------------------------------

.. code-block:: bash

  Rscript ../flexynesis_manuscript/src/collate_figure_source_data.R ../flexynesis_manuscript/data/Figure_Source_Data/


