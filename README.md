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

Hepatocellular carcinoma (HCC) is a widespread primary liver cancer with a high fatality rate. Despite several genes with oncogenic effects in HCC have been identified, many remain undiscovered. In this study, we conducted a comprehensive computational analysis to explore the involvement of genes within the same families as known driver genes in HCC.    Specifically, we expanded the concept beyond single-gene mutations to encompass gene families sharing homologous structures, integrating various omics data to comprehensively understand gene abnormalities in cancer. Our analysis identified 74 domains with an enriched mutation burden, 404 domain mutation hotspots, and 233 dysregulated driver genes. We observed that specific low-frequency somatic mutations may contribute to HCC occurrence, potentially overlooked by single-gene algorithms. Furthermore, we systematically analyzed how abnormalities in the ubiquitinated proteasome system (UPS) impact HCC, finding that abnormal genes in E3, E2, DUB families, and Degron genes often result in HCC by affecting the stability of oncogenic or tumor suppressor proteins. In conclusion, expanding the exploration of driver genes to include gene families with homologous structures emerges as a promising strategy for uncovering additional oncogenic alterations in HCC.


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
| Package    | Version |
|------------|---------|
| R          | 4.2.1   |
| data.table | 1.14.2  |
| dplyr      | 1.1.1   |
| ELMER      | 2.20.0  |
| ChAMP      | 2.26.0  |
| multiMiR   | 1.18.0  |
| ggplot2    | 4.2.3   |
| grid       | 4.2.1   |
| ggvenn     | 0.1.9   |


