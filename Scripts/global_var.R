
###### LIBRARIES ######

# Packages
library(ChromSCape)          
library(devtools)            
library(DropletUtils)        
library(irlba)               
library(corrplot)           
library(R.utils)            
library(scater)             
library(Rtsne)               
library(ccRemover)           
library(viridis)             
library(colorRamps)          
library(RColorBrewer)       
library(edgeR)                
library(gplots)              
library(ggplot2)             
library(RColorBrewer)       
library(genefilter)          
library(xtable)              
library(WriteXLS)             
library(data.table)           
library(stringr)             
library(limma)                
library(edgeR)               
library(monocle3)           
library(Seurat)              
library(dendextend)          
library(ape)                
library(ConsensusClusterPlus) 
library(Matrix)              
library(genefilter)          
library(ggpubr)              
library(eulerr)              
library(GenomicRanges)      
library(kableExtra)           
library(colorspace)          
library(forcats)             
library(dplyr)
library(tidyr)
library(scTools)
library(EnsDb.Hsapiens.v86)
library(Signac)
library(decoupleR)
library(factoextra)
library(ggdendro)
library(msigdbr)
library(fgsea)
library(Nebulosa)


###### VARIABLES ######

set.seed(1)
PDX_scRNA_samples = c("HBCx172_Capecitabin", "HBCx172_chemonaive", "HBCx218_AC","HBCx218_chemonaive", "HBCx221_Capecitabin",
                      "HBCx221_chemonaive","HBCx39_Capecitabin","HBCx39_chemonaive","HBCx95_Capecitabin","HBCx95_chemonaive")

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

AnnColor_Model = c("HBCx172" = "#e4be3f",
            "HBCx218" = "#191970",
            "HBCx221" = "tan4",
            "HBCx10" = "#7D26CD", 
            "HBCx14" = "#CD00CD", 
            "HBCx33" = "#00B2EE",
            "HBCx39" = "#FF7F00", 
            "HBCx95" = "#CD3700")

AnnColor_State = c("chemonaive" = "grey10",
            "Persister" = "chartreuse4", 
            "Relapse" = "grey40")

AnnShape_treatment = c("AC" = 21,
                   "Capecitabin" = 22, 
                   "Carboplatin" = 24, 
                   "chemonaive" = 3, 
                   "Cisplatin" = 23)

AnnShape_Cluster = c("1" = 6, 
                     "2" = 13)

TFs_top10_Allsamples = c("IRF1", "STAT1", "AP1", "JUN", "NFKB1", "NFKB","RELA","SPI1","ETS1", "CTNNB1")

List_genes_boxplot = c("KRT14", "CDH2", "JUNB", "E2F7", "ABCA1", "IDI1", "LDLR", "SQLE","IRF1", "STAT1", 
                       "JUN","FOS","NFKB2","RELA","JUND","FOSL1","FOSL2","FOSB","JUNB","NFKB1","REL",
                       "RELB","IRF3","IRF5","IRF7","IRF9","STAT2","STAT3","STAT5A","STAT5B")
