# This project stores the machine learning model from the paper 'From driver genes to gene families: a computational analysis of oncogenic mutations and ubiquitination anomalies in hepatocellular carcinoma'
This repository serves as a repository for the supporting code and data derived from the research conducted in this study, intended for the broader scientific community to utilize. The code and data stored here are available for academic use free of charge; any other type of use is strictly prohibited.

✅
This repository primarily focuses on storing machine learning code and software from the paper "From driver genes to gene families: a computational analysis of oncogenic mutations and ubiquitination anomalies in hepatocellular carcinoma".

✅
` **premutation test**`:  Identifies protein domains with significant mutation burdens.
 
✅
` **Domain_Hotspot**`:  Identifies mutation hotspots in domains (requires first defining protein domains with significant mutation burdens).

✅` **DDGs**`:  Defines the process for identifying driver-dysregulated genes.
 
✅
 !!!!

# Introduction

!!!!!


# Contents of this repository
In this repository, you can find the following folders:
```
project
├───DDGs
│   └───Code
├───Domain_Hotspot
│   ├───Code
│   └───test data
│       └───Alignment_tsv.zip
└───premutation test
    ├───code
    └───test data
        ├───length_files
        └───pfam_annotation
```


- **DDGs**: Defines the process for identifying driver-dysregulated genes.
  - **Code**: Contains the relevant code for identifying driver-dysregulated genes.
  
- **Domain_Hotspot**: Identifies mutation hotspots in domains (requires first defining protein domains with significant mutation burdens).
  - **Code**: Contains the relevant code for Domain_Hotspot.
  - **test data**: Contains the test data for Domain_Hotspot.
    - **Alignment_tsv.zip**: Multi-sequence alignment results for domain sequences.

- **premutation test**: Identifies protein domains with significant mutation burdens.
  - **Code**: Contains the relevant code for premutation test.
  - **test data**: Contains the test data for premutation test.
    - **length_files**: Protein length files (data from the UniProt database).
    - **pfam_annotation**: Protein domain annotation files (data from PfamScan).

# Required packages


