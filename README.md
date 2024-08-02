# AP-1 drives a recurrent drug persister state in triple negative breast cancer 

Scripts to reproduce the analysis in the paper 'AP-1 drives a recurrent drug persister state in triple negative breast cancer' by Baudre et al.

0.0 Setup
In order to re-run the analysis from the paper you must first download this repository. Then, at the base of the repository. Download the processed data (e.g. count matrices) from GSEXXXXX / EGAXXXX, and place it in the "Input" folder. In the "Input" folder, each kind of data should be placed in the appropriate directory. The following hiearchy should be kept:

 
├── Annotations
│   └── TFs_network_CollecTRI.csv
├── Input
│   └── hg38
│       ├── Lignees
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
│   │       ├── Lignees
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
│           ├── Lignees
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



    
Please refer to the scripts if you have doubts where you should place your input files.

1.0 scRNA-seq
PDX models
Note that:

BC976 stands for PDX_95/Patient_95
BC408 stands for PDX_39/Patient_39
BC1224 stands for PDX_172/Patient_172
MDA-MB-468 model
There are 3 sub-analyses in this directory: the main analysis ('1.Persister') and two sub analyses in response to epigenomic drug treatment ('2.UNC' and '3.UNC_5FU'). First run the 1.0 and 2.0 QC scripts. Then you can run the script of each analyses indepentently by respecting the order in each given analyses. Note that you must first run the single-cell ChIPseq analyse before running the script '5.0_ComparisonChIPseq.R'.

2.0 scChIP-seq
Run the scripts of H3K27me3 and H3K4me3 in any order.

3.0 ChIPreChIP (a.k.a. Sequential ChIP-seq)
Run these scripts to retrieve bivalent genes as well as bivalent pathways.

4.0 bulk ChIPseq
These scripts are mainly to produce snapshots of specific genes from the bigwigs.

Additional Files:
input/bulk_ChIPseq/MM468/chromatin_indexing_qc.csv -> ratios of Chromatin Indexing IP / input
input/scChIPseq/MM468/Raw_Counts/ -> raw counts used to calculate FrIP for scChIP
input/scChIPseq/MM468/BigWigs/ clusters C2 / C3 / C4 to produce Extened Figure 5h
