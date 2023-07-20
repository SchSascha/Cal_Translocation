# Description: vizualise enrichment analysis results
#
# author Sascha Schäuble
# date of creation: Wed Jun 14 11:45:34 2023
# license: MIT

library("tidyverse")
library("magrittr")
library("ggpubr")
library("data.table")
library("rvest")

PROJ_PATH <- "./" # set, e.g. getwd()
DAT_PATH <- paste0(PROJ_PATH, "dat/")
RES_PATH <- paste0(PROJ_PATH, "res/")

DATE_STR <- format(Sys.time(), "%y%m%d")

FN0 <- "DEGs.rdata"
FN1 <- "Calb_ORA.rds"
FN2 <- "IECs_ORA_tc_KEGG.rds"
FN3 <- "IECs_ORA_wtVSmut_KEGG.rds"
FN4 <- "KEGG_pathways2beremoved.txt"
# CONFIGs
SAVE_OUT <- T
# ================================================================ #


#### functions #####################################################

#
#### Description: add a column with meanLog2Fc to each enriched category
#
# IN:
#    variables
# OUT:
#   variables
# PRE:
#   df.enrich has the genes as comma separated string
#
add_meanFcPerCategory <- function(df.enrich, df.comps, inputGeneCol = "result_gene_list") {
  df.enrich %<>% mutate(meanLog2Fc = NA) 
  for (i in 1:dim(df.enrich)[1]) {
    dummy_genes <- df.enrich[[inputGeneCol]][i] %>% 
      str_split_1(pattern = ",") %>% 
      str_remove("_[A-Z]$")
    df.enrich$meanLog2Fc[i] <- df.comps %>% 
      filter(id %in% dummy_genes) %>% pull(log2_fc_mrn) %>% mean()
  }
  return(df.enrich)
}


#
#### Description: given a list of GO terms return revigo output
#
# IN:
#    goTerms : character vector
# OUT:
#   df with revigo terms
# PRE:
#   none
#
getRevigoByAPIcall <- function(goTerms, cutoff=0.7, rmObs = T) {
  rmObsStr <- if_else(rmObs, "true", "false")
  write_tsv(x = goTerms, col_names = F,file = "tmpGOterms.lst")
  system(paste("bash runRevigoAPI.sh tmpGOterms.lst", cutoff, rmObsStr, " > tmpRvgRes.html"))
  res.revigo <- rvest::html_table(x= read_html("tmpRvgRes.html"))
  Sys.sleep(0.5)   
  system("rm tmpGOterms.lst tmpRvgRes.html")
  stopifnot(length(res.revigo) == 1)
  
  return(res.revigo[[1]]) # results provided as tibble in a list
}
# ================================================================ #

#
#### Description: given a list of GO terms return revigo output
#
# IN:
#    goTerms : character vector
# OUT:
#   df with revigo terms
# PRE:
#   none
#
getReducedGOterms <- function(revigoDf) {
  return(revigoDf %>% filter(Representative == "null") %>% 
           dplyr::pull("Term ID"))
  
}
# ================================================================ #

# ================================================================ #

#### Data wrangling ################################################
cond_vec <- c("0.75h", "3h", "6h", "12h", "24h")

## load DEGs
load(paste0(DAT_PATH, FN0))

## load ORA results
dat.Calb              <- read_rds(paste0(DAT_PATH, FN1))
dat.IECs.KEGG         <- read_rds(paste0(DAT_PATH, FN2))
dat.IECs.wtVsMut.KEGG <- read_rds(paste0(DAT_PATH, FN3))

### filter and extend ORA results for plotting

## Calb
dat.Calb[["0.75h"]] %<>% mutate(geneRatio = result_count / bgd_count) 
dat.Calb[["0.75h"]] <- add_meanFcPerCategory(df.enrich = dat.Calb[["0.75h"]], 
                                             df.comps = dat.DEGs.Calb$x45)
dat.Calb[["3h"]] %<>% mutate(geneRatio = result_count / bgd_count)
dat.Calb[["3h"]] <- add_meanFcPerCategory(df.enrich = dat.Calb[["3h"]], 
                                          df.comps = dat.DEGs.Calb$x3)
dat.Calb[["6h"]] %<>% mutate(geneRatio = result_count / bgd_count)
dat.Calb[["6h"]] <- add_meanFcPerCategory(df.enrich = dat.Calb[["6h"]], 
                                          df.comps = dat.DEGs.Calb$x6)
dat.Calb[["12h"]] %<>% mutate(geneRatio = result_count / bgd_count)
dat.Calb[["12h"]] <- add_meanFcPerCategory(df.enrich = dat.Calb[["12h"]], 
                                           df.comps = dat.DEGs.Calb$x12)
dat.Calb[["24h"]] %<>% mutate(geneRatio = result_count / bgd_count)
dat.Calb[["24h"]] <- add_meanFcPerCategory(df.enrich = dat.Calb[["24h"]], 
                                           df.comps = dat.DEGs.Calb$x24)

# discard unwanted background count
for (e in 1:length(dat.Calb)) {
  dat.Calb[[e]] %<>% filter(bgd_count <= 1000 & bgd_count >= 3) 
}

# cleaning redundant terms with revigo
#  note that we KEEP obsolete terms
for (e in 1:length(dat.Calb)) {
  tmp_ids <- getReducedGOterms(revigoDf = getRevigoByAPIcall(
    goTerms = dat.Calb[[e]] %>% dplyr::select(id), cutoff = 0.5, rmObs = T))
  dat.Calb[[e]] %<>% filter(id %in% tmp_ids)
}

# clean terms with exactly the same genes
# keep the ones with less bg gene count
dat.Calb[["0.75h"]]$result_gene_list %>% sort()
dat.Calb[["0.75h"]]$name
dat.Calb[["0.75h"]]$result_gene_list %>% duplicated()
dat.Calb[["0.75h"]] %<>% filter( !(id %in% c("GO:0036171"))) # NOTE: not exact genes, but a better subclass is enriched


dat.Calb[["3h"]]$result_gene_list %>% sort() %>% duplicated()
dat.Calb[["3h"]]$name

dat.Calb[["6h"]] %>% 
  dplyr::arrange(result_gene_list) %>% pull %>% 
  duplicated()
dat.Calb[["6h"]]$name
dat.Calb[["6h"]] %<>% filter( !(id %in% c(
  "GO:0015942", "GO:0006629", "GO:0006694", "GO:1902652", "GO:0008202", "GO:1901615"
  )))

dat.Calb[["12h"]] %>% 
  dplyr::arrange(result_gene_list) %>% pull %>% 
  duplicated()

dat.Calb[["24h"]] %>% 
  dplyr::arrange(result_gene_list) %>% pull %>% 
  duplicated()
dat.Calb[["24h"]] %<>% filter( !(id %in% c(
  "GO:0140115"
)))

dat.Calb.df <- data.table::rbindlist(dat.Calb, idcol = "tp")
stopifnot((dat.Calb.df$result_gene_list %>% duplicated() %>% sum) == 0)
dat.Calb.df$tp %>% class()
dat.Calb.df$tp %>% table()
dat.Calb.df$tp %<>% factor(level = cond_vec)

## IECs KEGG 
# which unneccessary KEGG pathways to remove
keggPWYs2bRemoved <- paste0("KEGG:",
                            read_tsv(paste0(DAT_PATH, FN4), col_names = F) %>% 
                              pull()
)

## timecourse
# filter by background
for (t in 1:length(dat.IECs.KEGG)) {
  dat.IECs.KEGG[[t]] %<>% 
  filter(term_size <= 1000 & term_size >= 3)
}
# add mean fold change of relevant genes
dat.IECs.KEGG[["0.75h"]] <- add_meanFcPerCategory(df.enrich = dat.IECs.KEGG[["0.75h"]], 
                                                  df.comps = dat.DEGs.Hsap$x45, 
                                                  inputGeneCol = "intersection")
dat.IECs.KEGG[["3h"]] <- add_meanFcPerCategory(df.enrich = dat.IECs.KEGG[["3h"]], 
                                               df.comps = dat.DEGs.Hsap$x3, 
                                               inputGeneCol = "intersection")
dat.IECs.KEGG[["6h"]] <- add_meanFcPerCategory(df.enrich = dat.IECs.KEGG[["6h"]], 
                                               df.comps = dat.DEGs.Hsap$x6, 
                                               inputGeneCol = "intersection")
dat.IECs.KEGG[["12h"]] <- add_meanFcPerCategory(df.enrich = dat.IECs.KEGG[["12h"]], 
                                                df.comps = dat.DEGs.Hsap$x12, 
                                                inputGeneCol = "intersection")
dat.IECs.KEGG[["24h"]] <- add_meanFcPerCategory(df.enrich = dat.IECs.KEGG[["24h"]], 
                                                df.comps = dat.DEGs.Hsap$x24, 
                                                inputGeneCol = "intersection")

dat.IECs.KEGG.df <- data.table::rbindlist(dat.IECs.KEGG, idcol = "tp")

# filter unneccessary KEGG pathways
dat.IECs.KEGG.df %<>% filter(!(term_id %in% keggPWYs2bRemoved))
# add geneRatio column
dat.IECs.KEGG.df %<>% mutate(geneRatio = intersection_size / term_size)
dat.IECs.KEGG.df$tp %<>% factor(level = cond_vec)


## IECs KEGG wt vs mut data
# filter by background
for (t in 1:length(dat.IECs.wtVsMut.KEGG)) {
  dat.IECs.wtVsMut.KEGG[[t]] %<>% 
    filter(term_size <= 1000 & term_size >= 3)
}
# add mean fc column
dat.IECs.wtVsMut.KEGG$Hsa_Ca_bwp17_12h_VS_EC_only_12h <- add_meanFcPerCategory(
  df.enrich = dat.IECs.wtVsMut.KEGG$Hsa_Ca_bwp17_12h_VS_EC_only_12h, 
  df.comps = dat.DEGs.Hsap$bwp17, inputGeneCol = "intersection")
dat.IECs.wtVsMut.KEGG$Hsa_Ca_ece1_12h_VS_EC_only_12h <- add_meanFcPerCategory(
  df.enrich = dat.IECs.wtVsMut.KEGG$Hsa_Ca_ece1_12h_VS_EC_only_12h, 
  df.comps = dat.DEGs.Hsap$ece1, inputGeneCol = "intersection")
dat.IECs.wtVsMut.KEGG$Hsa_Ca_SC5314_12h_VS_EC_only_12h <- add_meanFcPerCategory(
  df.enrich = dat.IECs.wtVsMut.KEGG$Hsa_Ca_SC5314_12h_VS_EC_only_12h, 
  df.comps = dat.DEGs.Hsap$sc5314, inputGeneCol = "intersection")
dat.IECs.wtVsMut.KEGG$Hsa_Ca_efg1cph1_12h_VS_EC_only_12h <- add_meanFcPerCategory(
  df.enrich = dat.IECs.wtVsMut.KEGG$Hsa_Ca_efg1cph1_12h_VS_EC_only_12h, 
  df.comps = dat.DEGs.Hsap$efg1cph1, inputGeneCol = "intersection")

dat.IECs.wtVsMut.KEGG.df <- data.table::rbindlist(dat.IECs.wtVsMut.KEGG, idcol = "comparison")

# filter unneccessary KEGG pathways
dat.IECs.wtVsMut.KEGG.df %<>% filter(!(term_id %in% keggPWYs2bRemoved))
# add geneRatio column
dat.IECs.wtVsMut.KEGG.df %<>% mutate(geneRatio = intersection_size / term_size)
dat.IECs.wtVsMut.KEGG.df$comparison %<>% as_factor()

# ================================================================ #



#### vizualise GO ##################################################
p1 <- ggdotchart(dat.Calb.df %>% dplyr::rename("gene ratio" = "geneRatio"),
           y = "tp",
           x = "name",
           color = "meanLog2Fc",
           size = "gene ratio",
           xlab = "",
           ylab = "",
           rotate = T,
           legend.title = expression(paste("mean(log"[2], "(fc))")),
           ggtheme = theme_bw()
           ) +
  scale_color_gradient2(low  = hcl.colors(3, palette = "Blue-Red")[1],
                        mid  = hcl.colors(3, palette = "Blue-Red")[2],
                        high = hcl.colors(3, palette = "Blue-Red")[3]
                        )

ggexport(p1, filename = paste0(RES_PATH, "Calb_tc_GO_", DATE_STR, ".pdf"), height = 6, width = 10)
tiff(filename = paste0(RES_PATH, "Calb_tc_GO_", DATE_STR, ".tiff"), units="cm", 
     height=15, width=23, res=300,
     compression = "lzw")
print(p1)
dev.off()

## IECs KEGG
p2 <- ggdotchart(dat.IECs.KEGG.df %>% dplyr::rename("gene ratio" = "geneRatio"), 
                 y = "tp",
                 x = "term_name",
                 color = "meanLog2Fc",
                 size = "gene ratio",
                 xlab = "",
                 ylab = "",
                 rotate = T,
                 legend.title = expression(paste("mean(log"[2], "(fc))")),
                 ggtheme = theme_bw()
) +
scale_color_gradient2(low  = hcl.colors(3, palette = "Blue-Red")[1],
                      mid  = hcl.colors(3, palette = "Blue-Red")[2],
                      high = hcl.colors(3, palette = "Blue-Red")[3]
)


ggexport(p2, filename = paste0(RES_PATH, "IECs_tc_KEGG_", DATE_STR, ".pdf"), height = 5, width = 6.5)
tiff(filename = paste0(RES_PATH, "IECs_tc_KEGG_", DATE_STR, ".tiff"), units="cm", 
     height=12, width=16, res=300,
     compression = "lzw")
print(p2)
dev.off()

# wt vs mutant
p3 <- ggdotchart(dat.IECs.wtVsMut.KEGG.df %>% dplyr::rename("gene ratio" = "geneRatio"), 
                 y = "comparison",
                 x = "term_name",
                 color = "meanLog2Fc",
                 size = "gene ratio",
                 xlab = "",
                 ylab = "",
                 rotate = T,
                 legend.title = expression(paste("mean(log"[2], "(fc))")),
                 ggtheme = theme_bw(),
) %>%  ggpar(x.text.angle = 45) +
scale_color_gradient2(low  = hcl.colors(3, palette = "Blue-Red")[1],
                      mid  = hcl.colors(3, palette = "Blue-Red")[2],
                      high = hcl.colors(3, palette = "Blue-Red")[3]
)

ggexport(p3, filename = paste0(RES_PATH, "IECs_tc_KEGG_wt_vs_mut_", DATE_STR, ".pdf"), height = 5, width = 6.5)
tiff(filename = paste0(RES_PATH, "IECs_tc_KEGG_wt_vs_mut_", DATE_STR, ".tiff"), units="cm", 
     height=13, width=16, res=300,
     compression = "lzw")
print(p3)
dev.off()
# ================================================================ #
