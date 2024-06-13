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
  Command: ``/opt/R/4.2/bin/Rscript /fast/AG_Akalin/buyar/flexynesis_manuscript_work/flexynesis_manuscript/src/prepare_data.gdsc_vs_ccle.R raw/``

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

