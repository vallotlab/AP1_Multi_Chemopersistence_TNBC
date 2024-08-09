# AP-1 drives a recurrent drug persister state in triple negative breast cancer 

Scripts to reproduce the analysis in the paper **'AP-1 drives a recurrent drug persister state in triple negative breast cancer'** by _Baudre et al_.

## 0.0 Setup

In order to re-run the analysis from the paper you must first download this repository. Then, at the base of the repository, download the processed data (e.g. count matrices) from GSEXXXXX / EGAXXXX, and place it in the **_Input_** folder. In the **_Input_** folder, each kind of data should be placed in the appropriate directory. The following hiearchy and naming should be kept:

```
.
├── Annotations
│   ├── EnsDB.All.bed.df
│   ├── H12CORE_meme_format.meme
│   ├── Peaks_Common
│   │   ├── Peaks_MACS2_Common_FOSL1.bed
│   │   ├── Peaks_MACS2_Common_H3K27ac.bed
│   │   └── Peaks_MACS2_Common_H3K27ac_&_FOSL1.bed
│   ├── Peaks_Individual
│   │   ├── Peaks_MACS2_FOSL1_5FU_VF2_peaks.narrowPeak
│   │   ├── Peaks_MACS2_FOSL1_GOF_VF2_peaks.narrowPeak
│   │   ├── Peaks_MACS2_FOSL1_GOFctrl_VF2_peaks.narrowPeak
│   │   ├── Peaks_MACS2_FOSL1_VPR_VF2_peaks.narrowPeak
│   │   ├── Peaks_MACS2_FOSL1_VPRctrl_VF2_peaks.narrowPeak
│   │   ├── Peaks_MACS2_FOSL1_WT_VF2_peaks.narrowPeak
│   │   ├── Peaks_MACS2_H3K27ac_5FU_VF2_peaks.narrowPeak
│   │   ├── Peaks_MACS2_H3K27ac_GOF_VF2_peaks.narrowPeak
│   │   └── Peaks_MACS2_H3K27ac_WT_VF2_peaks.narrowPeak
│   └── TFs_network_CollecTRI.csv
├── Input
│   └── hg38
│       ├── Cell_lines
│       │   ├── Bulk_CutTag
│       │   │   ├── BAM
│       │   │   │   ├── MM468_ATCC_hu_5FU_d35_AM_m10y22_H3K27ac.sorted.bam
│       │   │   │   ├── MM468_ATCC_hu_5FU_d35_AM_m10y22_H3K27ac.sorted.bam.bai
│       │   │   │   ├── MM468_ATCC_hu_WT_AM_m10y22_H3K27ac.sorted.bam
│       │   │   │   ├── MM468_ATCC_hu_WT_AM_m10y22_H3K27ac.sorted.bam.bai
│       │   │   │   ├── MM468_NB_hu_FOSL1_GOF_AM_m10y22_HA.sorted.bam
│       │   │   │   ├── MM468_NB_hu_FOSL1_GOF_AM_m10y22_HA.sorted.bam.bai
│       │   │   │   ├── MM468_NB_hu_FOSL1_GOF_n1_AM_m04y23_H3K27ac.sorted.bam
│       │   │   │   ├── MM468_NB_hu_FOSL1_GOF_n1_AM_m04y23_H3K27ac.sorted.bam.bai
│       │   │   │   ├── MM468_NB_hu_FOSL1_GOF_n2_AM_m04y23_H3K27ac.sorted.bam
│       │   │   │   ├── MM468_NB_hu_FOSL1_GOF_n2_AM_m04y23_H3K27ac.sorted.bam.bai
│       │   │   │   ├── MM468_NB_hu_FOSL1_GOF_n2_AM_m04y23_HA.sorted.bam
│       │   │   │   ├── MM468_NB_hu_FOSL1_GOF_n2_AM_m04y23_HA.sorted.bam.bai
│       │   │   │   ├── MM468_NB_hu_FOSL1_GOF_n3_AM_m04y23_H3K27ac.sorted.bam
│       │   │   │   ├── MM468_NB_hu_FOSL1_GOF_n3_AM_m04y23_H3K27ac.sorted.bam.bai
│       │   │   │   ├── MM468_NB_hu_FOSL1_GOF_n3_AM_m04y23_HA.sorted.bam
│       │   │   │   ├── MM468_NB_hu_FOSL1_GOF_n3_AM_m04y23_HA.sorted.bam.bai
│       │   │   │   ├── MM468_NB_hu_VPRFOSL1g3_GOF_n1_AM_m04y23_FOSL1.sorted.bam
│       │   │   │   ├── MM468_NB_hu_VPRFOSL1g3_GOF_n1_AM_m04y23_FOSL1.sorted.bam.bai
│       │   │   │   ├── MM468_NB_hu_VPRctrl_GOF_n1_AM_m04y23_FOSL1.sorted.bam
│       │   │   │   ├── MM468_NB_hu_VPRctrl_GOF_n1_AM_m04y23_FOSL1.sorted.bam.bai
│       │   │   │   ├── MM468_NB_hu_WT_5FU_d35_n1_AM_m04y23_FOSL1.sorted.bam
│       │   │   │   ├── MM468_NB_hu_WT_5FU_d35_n1_AM_m04y23_FOSL1.sorted.bam.bai
│       │   │   │   ├── MM468_NB_hu_WT_5FU_d35_n2_AM_m04y23_FOSL1.sorted.bam
│       │   │   │   ├── MM468_NB_hu_WT_5FU_d35_n2_AM_m04y23_FOSL1.sorted.bam.bai
│       │   │   │   ├── MM468_NB_hu_WT_5FU_d35_n2_AM_m04y23_H3K27ac.sorted.bam
│       │   │   │   ├── MM468_NB_hu_WT_5FU_d35_n2_AM_m04y23_H3K27ac.sorted.bam.bai
│       │   │   │   ├── MM468_NB_hu_WT_5FU_d35_n3_AM_m04y23_FOSL1.sorted.bam
│       │   │   │   ├── MM468_NB_hu_WT_5FU_d35_n3_AM_m04y23_FOSL1.sorted.bam.bai
│       │   │   │   ├── MM468_NB_hu_WT_5FU_d35_n3_AM_m04y23_H3K27ac.sorted.bam
│       │   │   │   ├── MM468_NB_hu_WT_5FU_d35_n3_AM_m04y23_H3K27ac.sorted.bam.bai
│       │   │   │   ├── MM468_NB_hu_WT_d0_n1_AM_m04y23_FOSL1.sorted.bam
│       │   │   │   ├── MM468_NB_hu_WT_d0_n1_AM_m04y23_FOSL1.sorted.bam.bai
│       │   │   │   ├── MM468_NB_hu_WT_d0_n2_AM_m04y23_FOSL1.sorted.bam
│       │   │   │   ├── MM468_NB_hu_WT_d0_n2_AM_m04y23_FOSL1.sorted.bam.bai
│       │   │   │   ├── MM468_NB_hu_WT_d0_n2_AM_m04y23_H3K27ac.sorted.bam
│       │   │   │   ├── MM468_NB_hu_WT_d0_n2_AM_m04y23_H3K27ac.sorted.bam.bai
│       │   │   │   ├── MM468_NB_hu_WT_d0_n3_AM_m04y23_FOSL1.sorted.bam
│       │   │   │   ├── MM468_NB_hu_WT_d0_n3_AM_m04y23_FOSL1.sorted.bam.bai
│       │   │   │   ├── MM468_NB_hu_WT_d0_n3_AM_m04y23_H3K27ac.sorted.bam
│       │   │   │   ├── MM468_NB_hu_WT_d0_n3_AM_m04y23_H3K27ac.sorted.bam.bai
│       │   │   │   ├── MM468_NB_hu_pcr275_GOF_AM_m10y22_HA.sorted.bam
│       │   │   │   ├── MM468_NB_hu_pcr275_GOF_AM_m10y22_HA.sorted.bam.bai
│       │   │   │   ├── MM468_NB_hu_pcr276_GOF_n2_AM_m04y23_HA.sorted.bam
│       │   │   │   ├── MM468_NB_hu_pcr276_GOF_n2_AM_m04y23_HA.sorted.bam.bai
│       │   │   │   ├── MM468_NB_hu_pcr276_GOF_n3_AM_m04y23_HA.sorted.bam
│       │   │   │   └── MM468_NB_hu_pcr276_GOF_n3_AM_m04y23_HA.sorted.bam.bai
│       │   │   └── Bigwigs
│       │   │       ├── MM468_ATCC_hu_5FU_d35_AM_m10y22_H3K27ac.bw
│       │   │       ├── MM468_ATCC_hu_WT_AM_m10y22_H3K27ac.bw
│       │   │       ├── MM468_NB_hu_FOSL1_GOF_AM_m10y22_HA.bw
│       │   │       ├── MM468_NB_hu_FOSL1_GOF_n1_AM_m04y23_H3K27ac.bw
│       │   │       ├── MM468_NB_hu_FOSL1_GOF_n2_AM_m04y23_H3K27ac.bw
│       │   │       ├── MM468_NB_hu_FOSL1_GOF_n2_AM_m04y23_HA.bw
│       │   │       ├── MM468_NB_hu_FOSL1_GOF_n3_AM_m04y23_H3K27ac.bw
│       │   │       ├── MM468_NB_hu_FOSL1_GOF_n3_AM_m04y23_HA.bw
│       │   │       ├── MM468_NB_hu_VPRFOSL1g3_GOF_n1_AM_m04y23_FOSL1.bw
│       │   │       ├── MM468_NB_hu_VPRctrl_GOF_n1_AM_m04y23_FOSL1.bw
│       │   │       ├── MM468_NB_hu_WT_5FU_d35_n1_AM_m04y23_FOSL1.bw
│       │   │       ├── MM468_NB_hu_WT_5FU_d35_n2_AM_m04y23_FOSL1.bw
│       │   │       ├── MM468_NB_hu_WT_5FU_d35_n2_AM_m04y23_H3K27ac.bw
│       │   │       ├── MM468_NB_hu_WT_5FU_d35_n3_AM_m04y23_FOSL1.bw
│       │   │       ├── MM468_NB_hu_WT_5FU_d35_n3_AM_m04y23_H3K27ac.bw
│       │   │       ├── MM468_NB_hu_WT_d0_n1_AM_m04y23_FOSL1.bw
│       │   │       ├── MM468_NB_hu_WT_d0_n2_AM_m04y23_FOSL1.bw
│       │   │       ├── MM468_NB_hu_WT_d0_n2_AM_m04y23_H3K27ac.bw
│       │   │       ├── MM468_NB_hu_WT_d0_n3_AM_m04y23_FOSL1.bw
│       │   │       ├── MM468_NB_hu_WT_d0_n3_AM_m04y23_H3K27ac.bw
│       │   │       ├── MM468_NB_hu_pcr275_GOF_AM_m10y22_HA.bw
│       │   │       ├── MM468_NB_hu_pcr276_GOF_n2_AM_m04y23_HA.bw
│       │   │       └── MM468_NB_hu_pcr276_GOF_n3_AM_m04y23_HA.bw
│       │   └── scRNAseq
│       │       ├── BT20_chemonaive_barcodes.tsv.gz
│       │       ├── BT20_chemonaive_features.tsv.gz
│       │       ├── BT20_chemonaive_matrix.mtx.gz
│       │       ├── BT20_persister_barcodes.tsv.gz
│       │       ├── BT20_persister_features.tsv.gz
│       │       ├── BT20_persister_matrix.mtx.gz
│       │       ├── HCC38_chemonaive_barcodes.tsv.gz
│       │       ├── HCC38_chemonaive_features.tsv.gz
│       │       ├── HCC38_chemonaive_matrix.mtx.gz
│       │       ├── HCC38_persister_barcodes.tsv.gz
│       │       ├── HCC38_persister_features.tsv.gz
│       │       ├── HCC38_persister_matrix.mtx.gz
│       │       ├── MM468_5FU1_day33_barcodes.tsv.gz
│       │       ├── MM468_5FU1_day33_features.tsv.gz
│       │       ├── MM468_5FU1_day33_matrix.mtx.gz
│       │       ├── MM468_5FU2_day67_barcodes.tsv.gz
│       │       ├── MM468_5FU2_day67_features.tsv.gz
│       │       ├── MM468_5FU2_day67_matrix.mtx.gz
│       │       ├── MM468_5FU3_day50_barcodes.tsv.gz
│       │       ├── MM468_5FU3_day50_features.tsv.gz
│       │       ├── MM468_5FU3_day50_matrix.mtx.gz
│       │       ├── MM468_5FU3_day77_barcodes.tsv.gz
│       │       ├── MM468_5FU3_day77_features.tsv.gz
│       │       ├── MM468_5FU3_day77_matrix.mtx.gz
│       │       ├── MM468_GOF
│       │       │   ├── MM468_GOF_assignment_confidence_table.csv
│       │       │   ├── barcodes.tsv.gz
│       │       │   ├── features.tsv.gz
│       │       │   └── matrix.mtx.gz
│       │       ├── MM468_chemonaive_barcodes.tsv.gz
│       │       ├── MM468_chemonaive_features.tsv.gz
│       │       └── MM468_chemonaive_matrix.mtx.gz
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
│               │   ├── HBCx33_Capecitabine_rawmat.csv
│               │   ├── HBCx33_Cisplatin_rawmat.csv
│               │   └── HBCx33_chemonaive_rawmat.csv
│               └── scRNAseq
│                   ├── HBCx172_Capecitabine
│                   │   ├── barcodes.tsv.gz
│                   │   ├── features.tsv.gz
│                   │   └── matrix.mtx.gz
│                   ├── HBCx172_chemonaive
│                   │   ├── barcodes.tsv.gz
│                   │   ├── features.tsv.gz
│                   │   └── matrix.mtx.gz
│                   ├── HBCx218_AC
│                   │   ├── features.tsv.gz
│                   │   ├── matrix.mtx.gz
│                   │   └── web_summary.html
│                   ├── HBCx218_chemonaive
│                   │   ├── barcodes.tsv.gz
│                   │   ├── features.tsv.gz
│                   │   └── matrix.mtx.gz
│                   ├── HBCx221_Capecitabine
│                   │   ├── barcodes.tsv.gz
│                   │   ├── features.tsv.gz
│                   │   └── matrix.mtx.gz
│                   ├── HBCx221_chemonaive
│                   │   ├── barcodes.tsv.gz
│                   │   ├── features.tsv.gz
│                   │   └── matrix.mtx.gz
│                   ├── HBCx39_Capecitabine
│                   │   ├── barcodes.tsv.gz
│                   │   ├── features.tsv.gz
│                   │   └── matrix.mtx.gz
│                   ├── HBCx39_chemonaive
│                   │   ├── barcodes.tsv.gz
│                   │   ├── features.tsv.gz
│                   │   └── matrix.mtx.gz
│                   ├── HBCx95_Capecitabine
│                   │   ├── barcodes.tsv.gz
│                   │   ├── features.tsv.gz
│                   │   └── matrix.mtx.gz
│                   └── HBCx95_chemonaive
│                       ├── barcodes.tsv.gz
│                       ├── features.tsv.gz
│                       └── matrix.mtx.gz
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
├── README.md
├── Scripts
│   ├── global_var.R
│   └── hg38
│       ├── Cell_lines
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
All R packages used in the analysis are loaded by sourcing **_global_var.R_**. If you don't have them already installed on your laptop, please install them before running the analysis.

## 1.0 - PDX RNA Analysis

For this Analysis, please start by running each of the 3 scripts **_QCs_&_Objects_XXX.Rmd_**, and then run **_Common_Analysis.Rmd_**. Scripts in this section mainly refers to what was used to produce _Fig.1_, _Fig.2_, _Supp.Fig.1_ & _Supp.Fig.2_ in the paper. 

## 2.0 - Cell_lines scRNAseq Analysis 
For this Analysis, please start by running the script **_QCs_&_Objects.Rmd_**, and then run **_Analysis.Rmd_**. Scripts in this section mainly refers to what was used to produce _Fig.3_, _Fig.5_, _Supp.Fig.3_ & _Supp.Fig.5_ in the paper. 

## 3.0 - Cell_lines Bulk CUT&Tag Analysis 
For this Analysis, please start by running the script **_QCs_&_Objects.Rmd_**, and then run **_Analysis.Rmd_**. Scripts in this section mainly refers to what was used to produce _Fig.3_, _Fig.5_, _Supp.Fig.3_ & _Supp.Fig.5_ in the paper.

