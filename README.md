# Antidiabetic and Antioxidant Potential of *Oenocarpus bataua* Mart.

## About this Repository
This repository contains the workflow and data processing steps used for the statistical analysis of secondary metabolites present in *Oenocarpus bataua* Mart. (commonly known as Ungurahua).  

It includes raw and processed datasets, scripts for data treatment, and visualizations employed in the metabolomics analysis.  

The raw data are available at:  
[MetaboLights - Full Scan Data](https://www.ebi.ac.uk/metabolights/editor/guide/upload/REQ20250812212442)  

---

## Analysis Notebooks
The following notebooks describe the preprocessing and treatment of LC-MS/MS data in both ionization modes, which were subsequently used to create the input matrix for [MetaboAnalyst](https://www.metaboanalyst.ca/).  

- **Negative mode QC ([M-H]-)**: [NEG_QC_[M-H]-](https://github.com/IKIAM-NPLab/Antidiabetic-and-antioxidant-potential-of-Oenocarpus-bataua-Mart/blob/main/Treatment-Data/UNGURAGUA_NEGATIVE.md)  
- **Positive mode QC ([M+H]+)**: [POS_QC_[M+H]+](https://github.com/IKIAM-NPLab/Antidiabetic-and-antioxidant-potential-of-Oenocarpus-bataua-Mart/blob/main/Treatment-Data/UNGURAGUA_POSITIVE.md)  

---

## Processed Data for MetaboAnalyst
Final processed datasets ready for direct upload to MetaboAnalyst:

- [Positive ionization mode](https://github.com/IKIAM-NPLab/Antidiabetic-and-antioxidant-potential-of-Oenocarpus-bataua-Mart/tree/main/Results/Metaboloanalysis/Positive)  
- [Negative ionization mode](https://github.com/IKIAM-NPLab/Antidiabetic-and-antioxidant-potential-of-Oenocarpus-bataua-Mart/tree/main/Results/Metaboloanalysis/Negative)  

---

## PCA Analysis
Principal Component Analysis (PCA) plots were generated in **R** to visualize the clustering and quality control of samples.

- **Negative mode QC ([M-H]-)**  
  ![Figure_NEG_QC](https://github.com/IKIAM-NPLab/Antidiabetic-and-antioxidant-potential-of-Oenocarpus-bataua-Mart/blob/main/Results/Plots/unguragua_nrg.png)  

- **Positive mode QC ([M+H]+)**  
  ![Figure_POS_QC](https://github.com/IKIAM-NPLab/Antidiabetic-and-antioxidant-potential-of-Oenocarpus-bataua-Mart/blob/main/Results/Plots/unguragua_pos.png)  


