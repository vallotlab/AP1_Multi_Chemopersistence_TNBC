
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
library(memes)
library(universalmotif)


###### VARIABLES ######

set.seed(1)
PDX_scRNA_samples = c("HBCx172_Capecitabine", "HBCx172_chemonaive", "HBCx218_AC","HBCx218_chemonaive", "HBCx221_Capecitabine",
                      "HBCx221_chemonaive","HBCx39_Capecitabine","HBCx39_chemonaive","HBCx95_Capecitabine","HBCx95_chemonaive")

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
                   "Capecitabine" = 22, 
                   "Carboplatin" = 24, 
                   "chemonaive" = 3, 
                   "Cisplatin" = 23)

AnnShape_Cluster = c("1" = 6, 
                     "2" = 13)

TFs_top10_Allsamples = c("IRF1", "STAT1", "AP1", "JUN", "NFKB1", "NFKB","RELA","SPI1","ETS1", "CTNNB1")

List_genes_boxplot = c("KRT14", "CDH2", "JUNB", "E2F7", "ABCA1", "IDI1", "LDLR", "SQLE","IRF1", "STAT1", 
                       "JUN","FOS","NFKB2","RELA","JUND","FOSL1","FOSL2","FOSB","JUNB","NFKB1","REL",
                       "RELB","IRF3","IRF5","IRF7","IRF9","STAT2","STAT3","STAT5A","STAT5B")

Cell_lines_samples = c("MM468_chemonaive", "MM468_5FU1_day33","MM468_5FU2_day67","MM468_5FU3_day50","MM468_5FU3_day77","BT20_chemonaive", "BT20_persister","HCC38_chemonaive", "HCC38_persister")

List_MM468_persister = c("MM468_5FU1_day33","MM468_5FU2_day67","MM468_5FU3_day50","MM468_5FU3_day77")


###### FUNCTIONS ######

Associate_Annot = function(masterPeak, Annotations){
  Peak$peak_loca = paste0(Peak@seqnames, "_",Peak@ranges)
  genes.coords.all = read.table(file.path(Annotations, "EnsDB.All.bed.df"), header = T)
  genes.coords.all$seqnames = paste0("chr", genes.coords.all$seqnames)
  genes.coords.all = plyranges::as_granges(genes.coords.all)
  
  Peak = plyranges::join_overlap_inner(Peak, genes.coords.all)
  Peak = as.data.frame(Peak)
  Order_loca = c("Promoter_1K","Promoter_2K","Promoter_3K","Promoter_5K","TSS_1K","TSS_2K","TSS_5K","TSS_10K","Enhancers","Genes","Intergenic")
  Peak$seq_coord_system = factor(Peak$seq_coord_system, levels = Order_loca)
  Peak = Peak %>% dplyr::group_by(peak_loca) %>% dplyr::arrange(seq_coord_system, .by_group = T) %>% dplyr::slice_head(n = 1)
  Peak = plyranges::as_granges(Peak)
  
  genes.coords.all_onlygenes = genes.coords.all %>% dplyr::filter(seq_coord_system == "Genes")
  Intergenic_matdiff = Peak %>% dplyr::filter(seq_coord_system == "Intergenic") %>% plyranges::as_granges()
  Intergenic_matdiff$Nearest_gene = genes.coords.all_onlygenes$gene_name[nearest(Intergenic_matdiff, genes.coords.all_onlygenes)]
  Intergenic_matdiff = as.data.frame(Intergenic_matdiff)
  gene_MatDiff_fullannot = Peak %>% dplyr::filter(seq_coord_system != "Intergenic")
  gene_MatDiff_fullannot$Nearest_gene = gene_MatDiff_fullannot$gene_name
  gene_MatDiff_fullannot = as.data.frame(gene_MatDiff_fullannot)
  gene_MatDiff_fullannot = rbind(gene_MatDiff_fullannot, Intergenic_matdiff)
  gene_MatDiff_fullannot = plyranges::as_granges(gene_MatDiff_fullannot)
  return(gene_MatDiff_fullannot)
}

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}
