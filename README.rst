======================
flexynesis_manuscript
======================

All publication material relevant for the manuscript describing the flexynesis software package

Project Folder
==============

Accessible from Hulk/Beast/Max: ``/fast/AG_Akalin/buyar/flexynesis_manuscript_work/``

The ``./raw`` folder contains the original dataset downloaded from a source such as Cbioportal/TCGA/PharmacoGx/DepMAP.
The ``./prepared`` folder contains data prepared as input to flexynesis.

Datasets used in the manuscript
-------------------------------

Below is a description of the datasets used in the manuscript and how to prepare them for analysis with flexynesis

Downloaded Datasets
^^^^^^^^^^^^^^^^^^^

Go to ``/fast/AG_Akalin/buyar/flexynesis_manuscript_work/datasets``:

The ``./raw`` folder contains:

* **CCLE.rds**: downloaded from `Zenodo <https://zenodo.org/record/3905462/files/CCLE.rds?download=1>`_.
* **GDSC2.rds**: downloaded from `Zenodo <https://zenodo.org/record/3905481/files/GDSC2.rds?download=1>`_.
* **lgggbm_tcga_pub.tar.gz**: downloaded from `cbioportal <https://www.cbioportal.org/study/summary?id=lgggbm_tcga_pub>`_.
* **brca_metabric.tar.gz**: downloaded from `cbioportal <https://www.cbioportal.org/study/summary?id=brca_metabric>`_.
* **depmap**: downloaded from `depmap portal <https://depmap.org/portal/data_page/?tab=allData>`_.
* **nbl_target_2018_pub.tar.gz**: downloaded from `cbioportal <https://www.cbioportal.org/study/summary?id=nbl_target_2018_pub>`_.
* **GDCData**: TCGA cohort datasets for 33 cancer types downloaded using the TCGABiolinks package (`See GitHub <https://github.com/BIMSBbioinfo/uyar_et_al_multiomics_deeplearning>`_).
* **prot-trans**: protein sequence embeddings obtained from prot-trans-xl-uniref50 model on uniprot sequences.
* **describeProt**: protein level sequence/structure/function features from describeprot database (`Download here <http://biomine.cs.vcu.edu/servers/DESCRIBEPROT/download_database_value/9606_value.csv>`_).

PREPARED datasets used as input to flexynesis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``./prepared`` folder contains:

* **ccle_vs_gdsc**: Drug response data from cell lines from CCLE and GDSC2 datasets.
  Command: 
.. code-block:: bash 

    /opt/R/4.2/bin/Rscript ./src/prepare_data.gdsc_vs_ccle.R raw/

* **lgggbm_tcga_pub_processed**: Merged cohorts of LGG + GBM samples.
  Command: 

.. code-block:: bash 

    /opt/R/4.2/bin/Rscript ./src/prepare_data.LGG_GBM.R ./src/get_cbioportal_data.R

* **brca_metabric_processed**: METABRIC dataset processed.
.. code-block:: bash 

    /opt/R/4.2/bin/Rscript ./src/prepare_data.metabric.R ./src/get_cbioportal_data.R

* **single_cell_bonemarrow**: CITE-Seq dataset from Seurat.
  Command: 
.. code-block:: bash 

    /opt/R/4.2/bin/Rscript ./src/prepare_data.cite_seq.R

* **neuroblastoma_target_vs_depmap**: neuroblastoma patient samples (TARGET study) and cell lines (depmap).
  Command: 
.. code-block:: bash 

    /opt/R/4.2/bin/Rscript ./src/prepare_data.neuroblastoma_finetuning.R ./src/get_cbioportal_data.R ./raw/depmap/ ./src/utils.R

* **tcga_cancertype**: TCGA cancer cohort for ~21 cancer types 100 samples per each cohort.
  Command: 
.. code-block:: bash 

    /opt/R/4.2/bin/Rscript ./src/prepare_data.tcga_cancertype.R ./src/utils.R ./raw/GDCdata

* **depmap_gene_dependency**: Dataset for gene-dependency prediction in cell lines. Consists of depmap gene expression + prottrans embeddings + describeprot features.
  Command: 
.. code-block:: bash 

    /opt/R/4.2/bin/Rscript ./src/prepare_data.depmap.R ./src/utils.R ./raw/depmap/ ./raw/prot-trans/embeddings.protein_level.csv ./raw/uniprot2hgnc.RDS ./raw/describePROT/9606_value.csv


Figures
==========

How to reproduce figures:

Go to ``/fast/AG_Akalin/buyar/flexynesis_manuscript_work/analyses``:

Figure 1: single-task figures
-------------------------------

.. code-block:: bash

   /opt/R/4.2/bin/Rscript ../flexynesis_manuscript/src/figures_single_task.R ../flexynesis_manuscript/src/utils.R ./output2


Figures 2 and 3: multi-task figures
-------------------------------

.. code-block:: bash

   /opt/R/4.2/bin/Rscript ../flexynesis_manuscript/src/figures_multitask.R ../flexynesis_manuscript/src/utils.R ./output2

Figure 4: unsupervised clustering (tcga cancer types)
-------------------------------

.. code-block:: bash 

   /opt/R/4.2/bin/Rscript ../flexynesis_manuscript/src/figures_tcga_unsupervised.R ../flexynesis_manuscript/src/utils.R ./unsupervised_cancertype/

Figure 5: cross-modality prediction of cell line dependency probabilities 
-------------------------------

.. code-block:: bash 

   /opt/R/4.2/bin/Rscript ../flexynesis_manuscript/src/figures_depmap.R ../datasets/prepared/depmap_gene_dependency/ depmap_analysis/output/


Figure 6: demonstration of fine-tuning
-------------------

.. code-block:: bash
   /opt/R/4.2/bin/Rscript ../flexynesis_manuscript/src/figures_finetuning.R ../flexynesis_manuscript/src/utils.R finetuning/


Figure 7: marker analysis 
-------------------

.. code-block:: bash 

   /opt/R/4.2/bin/Rscript ../flexynesis_manuscript/src/figures_marker_analysis.R ../flexynesis_manuscript/src/utils.R marker_analysis/output/


    
Environment
===========

Install flexynesis
-------------------

.. code-block:: bash

    git clone https://github.com/BIMSBbioinfo/flexynesis.git
    cd flexynesis
    conda create -n flexynesis --file spec-file.txt
    conda activate flexynesis
    pip install -e .

Install other packages
----------------------

.. code-block:: bash

    guix package --manifest=guix.scm --profile=./manuscript

Activate environment
--------------------

.. code-block:: bash

    source ./manuscript/etc/profile
    conda activate flexynesis

Manuscript
==========

* `Link to google doc containing manuscript text <https://docs.google.com/document/d/10Slme7TQLll7FEBAOtuug-y7ybUxnC1M-QuYMJ_lnUE/edit?usp=sharing>`_
* `Link to flexynesis package repo <https://github.com/BIMSBbioinfo/flexynesis>`_
* `Link to flexynesis-benchmarks package repo <https://github.com/BIMSBbioinfo/flexynesis-benchmarks>`_
* `Link to gdrive for graphical abstract afdesign <https://drive.google.com/file/d/1-R8KrQTxgo9ocdqsliEd8NC7Ntb5hp7I/view?usp=drive_link>`_
