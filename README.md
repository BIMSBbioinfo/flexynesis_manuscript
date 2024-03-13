# flexynesis_manuscript
All publication material relevant for the manuscript describing the flexynesis software package  

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

