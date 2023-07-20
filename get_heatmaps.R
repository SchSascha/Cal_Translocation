# Description: generate heatmaps for genes of interest
#
# author Sascha Schäuble
# date of creation: Mon Jun 12 15:49:58 2023
# license: MIT

library("tidyverse")
library("magrittr")
library("ggpubr")
# library("ggsci")
library("ComplexHeatmap")
library("circlize")
library("org.Hs.eg.db")

PROJ_PATH <- "./" # set, e.g. getwd()
DAT_PATH <- paste0(PROJ_PATH, "dat/")
RES_PATH <- paste0(PROJ_PATH, "res/")

DATE_STR <- format(Sys.time(), "%y%m%d")

FN0 <- "Cal_zinc.tsv"
FN.cal.1 <- "dat_expr_1.tsv"
FN.cal.2 <- "dat_expr_2.tsv"
FN.cal.3 <- "dat_expr_3.tsv"
FN.hsa.1 <- "dat_expr_iecs_1.tsv"
FN.hsa.2 <- "dat_expr_iecs_2.tsv"

# ================================================================ #


# convenience function for ensembl id mapping
get_hsa_map <- function(from = "ENSEMBL", to = "SYMBOL") {
  geneID.map <- AnnotationDbi::select(org.Hs.eg.db,
                keys = keys(org.Hs.eg.db),
                columns = c(from,to),
                keytype = 'ENTREZID') %>% dplyr::select(-ENTREZID) %>% 
  drop_na() %>% distinct()
  
  return(geneID.map)
  }



#### Data wrangling ################################################
goi <- read_tsv(paste0(DAT_PATH, FN0))

## Calb data
# timecourse
dat.expr.1 <- read_tsv(paste0(DAT_PATH, FN.cal.1))
dat.expr.1.wide <- dat.expr.1 %>% dplyr::select(-c(ORF,DESeq2_adj_pval)) %>% 
  pivot_wider(names_from = comparison, values_from = log2_fc_mrn) %>% 
  column_to_rownames("Gene_name") %>% 
  dplyr::select(c(4,3,5,1,2)) # order the time course
dat.expr.1.wide %>% colnames()
dat.expr.1.wide %<>% dplyr::rename("0.75h"  = "Calbicans_air_CC_45_VS_air_Ca_45",
                           "3h"   = "Calbicans_air_CC_3_VS_air_Ca_3",
                           "6h"   = "Calbicans_air_CC_6_VS_air_Ca_6",
                           "12h"  = "Calbicans_air_CC_12_VS_air_Ca_12",
                           "24h"  = "Calbicans_air_CC_24_VS_air_Ca_24"
                           )
dat.expr.1.mat <-  dat.expr.1.wide %>% as.matrix()
dat.expr.1.mat %>% summary()

dat.expr.1.wide.p <- dat.expr.1 %>% dplyr::select(-c(ORF, log2_fc_mrn)) %>% 
  pivot_wider(names_from = comparison, values_from = DESeq2_adj_pval) %>% 
  column_to_rownames("Gene_name") %>% 
  dplyr::select(c(4,3,5,1,2)) # order the time course
dat.expr.1.wide.p %<>%
  mutate(Calbicans_air_CC_45_VS_air_Ca_45 = if_else(Calbicans_air_CC_45_VS_air_Ca_45 <= 0.05,1,0),
         Calbicans_air_CC_3_VS_air_Ca_3 = if_else(Calbicans_air_CC_3_VS_air_Ca_3 <= 0.05,1,0),
         Calbicans_air_CC_6_VS_air_Ca_6 = if_else(Calbicans_air_CC_6_VS_air_Ca_6 <= 0.05,1,0),
         Calbicans_air_CC_12_VS_air_Ca_12 = if_else(Calbicans_air_CC_12_VS_air_Ca_12 <= 0.05,1,0),
         Calbicans_air_CC_24_VS_air_Ca_24 = if_else(Calbicans_air_CC_24_VS_air_Ca_24 <= 0.05,1,0)
  )


## data - Ca vs pre
dat.expr.2 <- read_tsv(paste0(DAT_PATH, FN.cal.2))
dat.expr.2.wide <- dat.expr.2 %>% dplyr::select(-c(ORF, DESeq2_adj_pval)) %>% 
  pivot_wider(names_from = comparison, values_from = log2_fc_mrn) %>% 
  column_to_rownames("Gene_name") %>% 
  dplyr::select(c(4,3,5,1,2)) # order the time course
dat.expr.2.wide %>% colnames()
dat.expr.2.wide %<>% dplyr::rename("0.75h"  = "Calbicans_air_Ca_45_VS_air_Ca_0",
                            "3h"   = "Calbicans_air_Ca_3_VS_air_Ca_0",
                            "6h"   = "Calbicans_air_Ca_6_VS_air_Ca_0",
                            "12h"  = "Calbicans_air_Ca_12_VS_air_Ca_0",
                            "24h"  = "Calbicans_air_Ca_24_VS_air_Ca_0"
)
dat.expr.2.mat <-  dat.expr.2.wide %>% as.matrix()
dat.expr.2.mat %>% summary()

dat.expr.2.wide.p <- dat.expr.2 %>% dplyr::select(-c(ORF, log2_fc_mrn)) %>% 
  pivot_wider(names_from = comparison, values_from = DESeq2_adj_pval) %>% 
  column_to_rownames("Gene_name") %>% 
  dplyr::select(c(4,3,5,1,2)) # order the time course
dat.expr.2.wide.p %<>%
  mutate(Calbicans_air_Ca_45_VS_air_Ca_0 = if_else(Calbicans_air_Ca_45_VS_air_Ca_0 <= 0.05,1,0),
         Calbicans_air_Ca_3_VS_air_Ca_0 = if_else(Calbicans_air_Ca_3_VS_air_Ca_0 <= 0.05,1,0),
         Calbicans_air_Ca_6_VS_air_Ca_0 = if_else(Calbicans_air_Ca_6_VS_air_Ca_0 <= 0.05,1,0),
         Calbicans_air_Ca_12_VS_air_Ca_0 = if_else(Calbicans_air_Ca_12_VS_air_Ca_0 <= 0.05,1,0),
         Calbicans_air_Ca_24_VS_air_Ca_0 = if_else(Calbicans_air_Ca_24_VS_air_Ca_0 <= 0.05,1,0)
  )

## data - CC vs pre
dat.expr.3 <- read_tsv(paste0(DAT_PATH, FN.cal.3))
dat.expr.3.wide <- dat.expr.3 %>% dplyr::select(-c(ORF, DESeq2_adj_pval)) %>% 
  pivot_wider(names_from = comparison, values_from = log2_fc_mrn) %>% 
  column_to_rownames("Gene_name") %>% 
  dplyr::select(c(4,3,5,1,2)) # order the time course
dat.expr.3.wide %>% colnames()
dat.expr.3.wide %<>% dplyr::rename("0.75h"  = "Calbicans_air_CC_45_VS_air_Ca_0",
                            "3h"   = "Calbicans_air_CC_3_VS_air_Ca_0",
                            "6h"   = "Calbicans_air_CC_6_VS_air_Ca_0",
                            "12h"  = "Calbicans_air_CC_12_VS_air_Ca_0",
                            "24h"  = "Calbicans_air_CC_24_VS_air_Ca_0"
)
dat.expr.3.mat <-  dat.expr.3.wide %>% as.matrix()
dat.expr.3.mat %>% summary()

dat.expr.3.wide.p <- dat.expr.3 %>% dplyr::select(-c(ORF, log2_fc_mrn)) %>% 
  pivot_wider(names_from = comparison, values_from = DESeq2_adj_pval) %>% 
  column_to_rownames("Gene_name") %>% 
  dplyr::select(c(4,3,5,1,2)) # order the time course
dat.expr.3.wide.p %<>%
  mutate(Calbicans_air_CC_45_VS_air_Ca_0 = if_else(Calbicans_air_CC_45_VS_air_Ca_0 <= 0.05,1,0),
         Calbicans_air_CC_3_VS_air_Ca_0 = if_else(Calbicans_air_CC_3_VS_air_Ca_0 <= 0.05,1,0),
         Calbicans_air_CC_6_VS_air_Ca_0 = if_else(Calbicans_air_CC_6_VS_air_Ca_0 <= 0.05,1,0),
         Calbicans_air_CC_12_VS_air_Ca_0 = if_else(Calbicans_air_CC_12_VS_air_Ca_0 <= 0.05,1,0),
         Calbicans_air_CC_24_VS_air_Ca_0 = if_else(Calbicans_air_CC_24_VS_air_Ca_0 <= 0.05,1,0)
  )

### IEC data
geneID.map <- get_hsa_map()

## main heatmap
dat.expr.iecs.1 <- read_tsv(paste0(DAT_PATH, FN.hsa.1))
dat.expr.iecs.1.wide <- dat.expr.iecs.1 %>% dplyr::select(-c(id, DESeq2_adj_pval)) %>% 
  pivot_wider(names_from = comparison, values_from = log2_fc_mrn) %>% 
  column_to_rownames("SYMBOL")
dat.expr.iecs.1.wide %>% colnames()
dat.expr.iecs.1.wide %<>% dplyr::rename("ece1 vs. bwp17" = "Hsapiens_Ca_ece1_12h_VS_Ca_bwp17_12h",
                            "efg1cph1 vs. SC5314" = "Hsapiens_Ca_efg1cph1_12h_VS_Ca_SC5314_12h")
dat.expr.iecs.1.wide.p <- dat.expr.iecs.1 %>% dplyr::select(-c(id, log2_fc_mrn)) %>% 
  pivot_wider(names_from = comparison, values_from = DESeq2_adj_pval) %>% 
  column_to_rownames("SYMBOL")
dat.expr.iecs.1.wide.p %<>%
  mutate(Hsapiens_Ca_ece1_12h_VS_Ca_bwp17_12h = if_else(Hsapiens_Ca_ece1_12h_VS_Ca_bwp17_12h <= 0.05,1,0),
         Hsapiens_Ca_efg1cph1_12h_VS_Ca_SC5314_12h = if_else(Hsapiens_Ca_efg1cph1_12h_VS_Ca_SC5314_12h <= 0.05,1,0))

dat.expr.iecs.1.mat <-  dat.expr.iecs.1.wide %>% as.matrix()
dat.expr.iecs.1.mat %>% summary()


## supplementary heatmap
dat.expr.iecs.2 <- read_tsv(paste0(DAT_PATH, FN.hsa.2))
dat.expr.iecs.2.wide <- dat.expr.iecs.2 %>% dplyr::select(-c(id, DESeq2_adj_pval)) %>% 
  pivot_wider(names_from = comparison, values_from = log2_fc_mrn) %>% 
  column_to_rownames("SYMBOL")
dat.expr.iecs.2.wide %>% colnames()
dat.expr.iecs.2.wide %<>% dplyr::select(c(1,2,4,3)) 
dat.expr.iecs.2.wide %<>% dplyr::rename("bwp17"    = "Hsapiens_Ca_bwp17_12h_VS_EC_only_12h",
                                        "ece1"     = "Hsapiens_Ca_ece1_12h_VS_EC_only_12h",
                                        "efg1cph1" = "Hsapiens_Ca_efg1cph1_12h_VS_EC_only_12h",
                                        "SC5314"   = "Hsapiens_Ca_SC5314_12h_VS_EC_only_12h"
                                        )
dat.expr.iecs.2.wide.p <- dat.expr.iecs.2 %>% dplyr::select(-c(id, log2_fc_mrn)) %>% 
  pivot_wider(names_from = comparison, values_from = DESeq2_adj_pval) %>% 
  column_to_rownames("SYMBOL")
dat.expr.iecs.2.wide.p %<>%
  mutate(Hsapiens_Ca_bwp17_12h_VS_EC_only_12h = if_else(Hsapiens_Ca_bwp17_12h_VS_EC_only_12h <= 0.05,1,0),
         Hsapiens_Ca_ece1_12h_VS_EC_only_12h = if_else(Hsapiens_Ca_ece1_12h_VS_EC_only_12h <= 0.05,1,0),
         Hsapiens_Ca_efg1cph1_12h_VS_EC_only_12h = if_else(Hsapiens_Ca_efg1cph1_12h_VS_EC_only_12h <= 0.05,1,0),
         Hsapiens_Ca_SC5314_12h_VS_EC_only_12h = if_else(Hsapiens_Ca_SC5314_12h_VS_EC_only_12h <= 0.05,1,0)
         )
dat.expr.iecs.2.wide.p %>% colnames()
dat.expr.iecs.2.wide.p %<>% dplyr::select(c(1,2,4,3)) 

dat.expr.iecs.2.mat <-  dat.expr.iecs.2.wide %>% as.matrix()
dat.expr.iecs.2.mat %>% summary()

# ================================================================ #


#### heatmaps Ca time course #######################################
## CC vs Ca for CC vs Ca per time point
p.Ca.1 <- ComplexHeatmap::Heatmap(dat.expr.1.mat,
                        col = colorRamp2(breaks = c(-5, 0, 2), 
                                         hcl.colors(3, "Blue-Red")),
                        heatmap_legend_param = list(title = expression(paste("Log"[2], "(fc)"))),
                        cluster_columns = F,
                        cluster_rows = F,
                        column_names_rot = 45,
                        width = unit(5*2, "char"),
                        height = unit(6*2, "char"),
                        row_names_gp = gpar(fontface = "italic"),
                        column_title = "infection vs. medium",
                        cell_fun = function(j, i, x, y, width, height, fill) {
                          if(dat.expr.1.wide.p[i,j]==1)
                            grid.points(pch = 8, x, y, size = unit(1.5, "mm"),
                                        gp = gpar(col = "black"))
                        })
p.Ca.1 %>% ggexport(filename = paste0(RES_PATH, "Ca_zinc_CC_vs_Ca_", DATE_STR, ".pdf"), 
                    height = 3, width = 4) 
tiff(paste0(RES_PATH, "Ca_zinc_CC_vs_Ca_", DATE_STR, ".tiff"), units="in", width=4, height=3, res=300,
     compression = "lzw")
draw(p.Ca.1)
dev.off()

## Ca vs 0h tp (preculture)
p.Ca.2 <- ComplexHeatmap::Heatmap(dat.expr.2.mat,
                                  col = colorRamp2(breaks = c(-3, 0, 11), 
                                                   hcl.colors(3, "Blue-Red")),
                                  heatmap_legend_param = list(title = expression(paste("Log"[2], "(fc)"))),
                                  cluster_columns = F,
                                  cluster_rows = F,
                                  column_names_rot = 45,
                                  width = unit(5*2.5, "char"),
                                  height = unit(6*2.5, "char"),
                                  row_names_gp = gpar(fontface = "italic"),
                                  column_title = "medium vs. preculture",
                                  cell_fun = function(j, i, x, y, width, height, fill) {
                                    if(dat.expr.2.wide.p[i,j]==1)
                                      grid.points(pch = 8, x, y, size = unit(1.5, "mm"),
                                                  gp = gpar(col = "black"))
                                  })
p.Ca.2 %>% ggexport(filename = paste0(RES_PATH, "Ca_zinc_Ca_vs_pre_", DATE_STR, ".pdf"), 
                    height = 4, width = 4) 
tiff(paste0(RES_PATH, "Ca_zinc_Ca_vs_pre_", DATE_STR, ".tiff"), units="in", width=4, height=4, res=300,
     compression = "lzw")
draw(p.Ca.2)
dev.off()
## Cc vs 0h tp (preculture)
p.Ca.3 <- ComplexHeatmap::Heatmap(dat.expr.3.mat,
                                  col = colorRamp2(breaks = c(-3, 0, 11), 
                                                   hcl.colors(3, "Blue-Red")),
                                  heatmap_legend_param = list(title = expression(paste("Log"[2], "(fc)"))),
                                  cluster_columns = F,
                                  cluster_rows = F,
                                  show_heatmap_legend  =F,
                                  column_names_rot = 45, 
                                  width = unit(5*2.5, "char"),
                                  height = unit(6*2.5, "char"),
                                  row_names_gp = gpar(fontface = "italic"),
                                  column_title = "infection vs. preculture",
                                  cell_fun = function(j, i, x, y, width, height, fill) {
                                    if(dat.expr.3.wide.p[i,j]==1)
                                      grid.points(pch = 8, x, y, size = unit(1.5, "mm"),
                                                  gp = gpar(col = "black"))
                                  })

p.Ca.3 %>% ggexport(filename = paste0(RES_PATH, "Ca_zinc_CC_vs_pre_", DATE_STR, ".pdf"), 
                    height = 4, width = 4) 
tiff(paste0(RES_PATH, "Ca_zinc_CC_vs_pre_", DATE_STR, ".tiff"), units="in", width=4, height=4, res=300,
     compression = "lzw")
draw(p.Ca.3)
dev.off()
ggexport(plotlist = p.Ca.3 + p.Ca.2, 
         filename = paste0(RES_PATH, "Ca_zinc_CCCa_vs_pre_", DATE_STR, ".pdf"),
         height = 4, width = 7) 
tiff(paste0(RES_PATH, "Ca_zinc_CCCa_vs_pre_", DATE_STR, ".tiff"), 
     units="in", width=7, height=4, res=300,
     compression = "lzw")
draw(p.Ca.3 + p.Ca.2)
dev.off()
# ================================================================ #


#### heatmaps IECs mutant vs wt ####################################
col_labels = c("*ece1*Δ/Δ vs. WT (BWP17)", "*efg1*ΔΔ/*cph1*ΔΔ vs. WT (SC5314)")
p.IEC.1 <- ComplexHeatmap::Heatmap(dat.expr.iecs.1.mat,
                                  col = colorRamp2(breaks = c(-3, 0, 3), 
                                                   hcl.colors(3, "Blue-Red")),
                                  heatmap_legend_param = list(title = expression(paste("Log"[2], "(fc)"))),
                                  cluster_columns = F,
                                  width = unit(2*1.5, "char"),
                                  height = unit((dat.expr.iecs.1.mat %>% dim())[1]*1, "char"),
                                  cell_fun = function(j, i, x, y, width, height, fill) {
                                    if(dat.expr.iecs.1.wide.p[i,j]==1)
                                      grid.points(pch = 8, x, y, size = unit(0.4, "char"),
                                                  gp = gpar(col = "black"))
                                  },
                                  # column_labels = c("ece1Δ/Δ vs. WT (BWP17)", "efg1ΔΔ/cph1ΔΔ vs. WT (SC5314)")
                                  column_labels = gt_render(col_labels),
                                  row_names_gp = gpar(fontface = "italic")
                                  )
# p.IEC.1 %>% ggexport(filename = paste0(RES_PATH, "IEC_DEGs_1.pdf"), height = 18, width = 5)
cairo_pdf(filename = paste0(RES_PATH, "IEC_DEGs_1_logfc1_", DATE_STR, ".pdf"), height = 9, width = 4)
draw(p.IEC.1)
dev.off()
tiff(filename = paste0(RES_PATH, "IEC_DEGs_1_logfc1_", DATE_STR, ".tiff"), units="cm", height=23, width=10, res=300,
     compression = "lzw")
draw(p.IEC.1)
dev.off()

col_labels = c("WT (BWP17)", "*ece1*Δ/Δ", "WT (SC5314)", "*efg1*ΔΔ/*cph1*ΔΔ")
p.IEC.2 <- ComplexHeatmap::Heatmap(dat.expr.iecs.2.mat,
                                   col = colorRamp2(breaks = c(-3, 0, 8), 
                                                    hcl.colors(3, "Blue-Red")),
                                   heatmap_legend_param = list(title = expression(paste("Log"[2], "(fc)"))),
                                   cluster_columns = F,
                                   width = unit(4*1.5, "char"),
                                   height = unit((dat.expr.iecs.2.mat %>% dim())[1]*1, "char"),
                                   cell_fun = function(j, i, x, y, width, height, fill) {
                                     if(dat.expr.iecs.2.wide.p[i,j]==1)
                                       grid.points(pch = 8, x, y, size = unit(0.4, "char"),
                                                   gp = gpar(col = "black"))
                                   },
                                   column_labels = gt_render(col_labels),
                                   row_names_gp = gpar(fontface = "italic")
)
cairo_pdf(filename = paste0(RES_PATH, "IEC_DEGs_2_logfc1_", DATE_STR, ".pdf"), height = 8, width = 4)
draw(p.IEC.2)
dev.off()
tiff(filename = paste0(RES_PATH, "IEC_DEGs_2_logfc1_", DATE_STR, ".tiff"), units="cm", height=20, width=11, res=300,
     compression = "lzw")
draw(p.IEC.2)
dev.off()

# ================================================================ #

