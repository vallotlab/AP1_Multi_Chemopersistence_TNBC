# AP-1 drives a recurrent drug persister state in triple negative breast cancer 

Scripts to reproduce the analysis in the paper **'AP-1 drives a recurrent drug persister state in triple negative breast cancer'** by _Baudre et al_.

## 0.0 Setup

In order to re-run the analysis from the paper you must first download this repository. Then, at the base of the repository. Download the processed data (e.g. count matrices) from GSEXXXXX / EGAXXXX, and place it in the "Input" folder. In the "Input" folder, each kind of data should be placed in the appropriate directory. The following hiearchy should be kept:

```
├── Annotations
│   └── TFs_network_CollecTRI.csv
├── Input
│   └── hg38
│       ├── Cell_lines
│       │   ├── Bulk_CutTag
│       │   └── scRNAseq
│       └── PDX
│           └── RNA
│               ├── BulkRNAseq
│               │   ├── HBCx14_AC_rawmat.csv
│               │   ├── HBCx14_Carboplatin_rawmat.csv
│               │   ├── HBCx14_Cisplatin_rawmat.csv
│               │   └── HBCx14_chemonaive_rawmat.csv
│               ├── Microarray
│               │   ├── HBCx10_AC_rawmat.csv
│               │   ├── HBCx10_chemonaive_rawmat.csv
│               │   ├── HBCx33_Capecitabin_rawmat.csv
│               │   ├── HBCx33_Cisplatin_rawmat.csv
│               │   └── HBCx33_chemonaive_rawmat.csv
│               └── scRNAseq
│                   ├── HBCx172_Capecitabin
│                   │   ├── barcodes.tsv.gz
│                   │   ├── features.tsv.gz
│                   │   ├── matrix.mtx.gz
│                   ├── HBCx172_chemonaive
│                   │   ├── barcodes.tsv.gz
│                   │   ├── features.tsv.gz
│                   │   ├── matrix.mtx.gz
│                   ├── HBCx218_AC
│                   │   ├── barcodes.tsv.gz
│                   │   ├── features.tsv.gz
│                   │   ├── matrix.mtx.gz
│                   ├── HBCx218_chemonaive
│                   │   ├── barcodes.tsv.gz
│                   │   ├── features.tsv.gz
│                   │   ├── matrix.mtx.gz
│                   ├── HBCx221_Capecitabin
│                   │   ├── barcodes.tsv.gz
│                   │   ├── features.tsv.gz
│                   │   ├── matrix.mtx.gz
│                   ├── HBCx221_chemonaive
│                   │   ├── barcodes.tsv.gz
│                   │   ├── features.tsv.gz
│                   │   ├── matrix.mtx.gz
│                   ├── HBCx39_Capecitabin
│                   │   ├── barcodes.tsv.gz
│                   │   ├── features.tsv.gz
│                   │   ├── matrix.mtx.gz
│                   ├── HBCx39_chemonaive
│                   │   ├── barcodes.tsv.gz
│                   │   ├── features.tsv.gz
│                   │   ├── matrix.mtx.gz
│                   ├── HBCx95_Capecitabin
│                   │   ├── barcodes.tsv.gz
│                   │   ├── features.tsv.gz
│                   │   ├── matrix.mtx.gz
│                   └── HBCx95_chemonaive
│                       ├── barcodes.tsv.gz
│                       ├── features.tsv.gz
│                       ├── matrix.mtx.gz
├── Output
│   ├── Figs
│   │   └── hg38
│   │       ├── Cell_lines
│   │       │   ├── Bulk_CutTag
│   │       │   └── scRNAseq
│   │       └── PDX
│   │           └── RNA
│   │               ├── BulkRNAseq
│   │               ├── Common_Analysis
│   │               ├── Microarray
│   │               └── scRNAseq
│   └── Objects
│       └── hg38
│           ├── Cell_lines
│           │   ├── Bulk_CutTag
│           │   └── scRNAseq
│           └── PDX
│               └── RNA
│                   ├── BulkRNAseq
│                   ├── Common_Analysis
│                   ├── Microarray
│                   └── scRNAseq
│                       ├── HBCx172_Capecitabin
│                       ├── HBCx172_chemonaive
│                       ├── HBCx218_AC
│                       ├── HBCx218_chemonaive
│                       ├── HBCx221_Capecitabin
│                       ├── HBCx221_chemonaive
│                       ├── HBCx39_Capecitabin
│                       ├── HBCx39_chemonaive
│                       ├── HBCx95_Capecitabin
│                       └── HBCx95_chemonaive
├── Scripts
│   ├── global_var.R
│   └── hg38
│       ├── Lignees
│       │   ├── Bulk_CutTag
│       │   │   ├── Analysis.Rmd
│       │   │   └── QCs_&_Objects.Rmd
│       │   └── scRNAseq
│       │       ├── Analysis.Rmd
│       │       └── QCs_&_Objects.Rmd
│       └── PDX
│           └── RNA
│               ├── Common_Analysis.Rmd
│               ├── QCs_&_Objects_BulkRNAseq.Rmd
│               ├── QCs_&_Objects_Microarray.Rmd
│               └── QCs_&_Objects_scRNAseq.Rmd

```
    
Please refer to the scripts if you have doubts where you should place your input files.

## 1.0 - PDX RNA Analysis

For this Analysis, please start by running each of the 3 scripts "QCs_&_Objects_***.Rmd", and then run Common_Analysis.Rmd. Scripts in this section mainly refers to what was used to produce Fig.1, Fig.2, Supp.Fig.1 & Supp.Fig.2 in the paper. 

## 2.0 - Cell_lines scRNAseq Analysis 
For this Analysis, please start by running the script "QCs_&_Objects.Rmd", and then run Analysis.Rmd. Scripts in this section mainly refers to what was used to produce Fig.3, Fig.5, Supp.Fig.3 & Supp.Fig.5 in the paper. 

## 3.0 - Cell_lines Bulk CUT&Tag Analysis 
For this Analysis, please start by running the script "QCs_&_Objects.Rmd", and then run Analysis.Rmd. Scripts in this section mainly refers to what was used to produce Fig.3, Fig.5, Supp.Fig.3 & Supp.Fig.5 in the paper.

