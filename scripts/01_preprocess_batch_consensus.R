
# ============================================================
# Script 01 - Preprocessing, batch correction, consensus clustering
# (formerly figure1.R). Run from the repository root.
#
# Inputs : data/GSE60427_RAW.rds, data/GSE130823_RAW.rds,
#          data/GSE60427.rds, data/GSE130823.rds
# Outputs: data/expr_combat_intergration.rds,
#          data/expr_combat_intergration_df.rds,
#          data/sample_info_intergration.rds,
#          consensus_out_intergration/clusters.csv,
#          results/figs/batch_before.png, results/figs/batch_after.png
# ============================================================

##################################data preprocess and intergration

# setwd("D:\\document\\20260425_pujiang_project\\submit\\PLOS One\\revision_v4")

library(data.table)
library(dplyr)
library(tidyr)
library(GEOquery)
library('factoextra')
library("FactoMineR")

#
rename_probe <- function(rds_raw,rds_series){
  tmp_exp=readRDS(rds_raw)
  series_matrix=readRDS(rds_series)
  rds_feature = fData(series_matrix)
  tmp_exp %>% 
    as.data.frame() %>% 
    mutate(ID=rownames(.)) %>% 
    left_join(.,
              rds_feature %>% dplyr::select(c(ID,GENE_SYMBOL))
    ) %>% 
    dplyr::select(-c('ID')) %>% 
    filter(GENE_SYMBOL!="") %>%
    group_by(
      GENE_SYMBOL
    ) %>%
    summarise(
      across(where(is.numeric),median,na.rm=T)
    ) -> exp_filter
  return(exp_filter)
}
GSE60427 <- rename_probe('data/GSE60427_RAW.rds','data/GSE60427.rds') #32 x 29548
GSE130823 <- rename_probe('data/GSE130823_RAW.rds','data/GSE130823.rds') #94 x 31893

GSE60427 %>% 
  inner_join(.,
             GSE130823
  ) -> expr_mat # 127 samples, 29388 genes
gse=c('GSE60427','GSE130823')
lapply(
  gse,
  FUN=function(x){
    pData(readRDS(paste0('data/',x,'.rds'))) %>% 
      dplyr::select(c('title','geo_accession','source_name_ch1')) %>% 
      mutate(GSE=x)
  }
) -> res_list
rbindlist(res_list) -> res_df


#PCA分析
expr_mat %>% dplyr::select(-c('GENE_SYMBOL')) -> df
res.pca <- prcomp(t(df))
fviz_pca_ind(res.pca,
             geom = c("point")
)+
  geom_point(aes(
    #shape = factor(res_df$GSE), 
    colour = factor(res_df$GSE)
  ),
  size=5,alpha=0.7)+
  scale_color_manual(
    values=c(
      'GSE130823'='#956EB5',
      'GSE60427'='#6A9750'
      
    )
  )+
  theme_classic(base_size = 20) +
  theme(
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.8
    )
  )+
  theme(legend.position = 'none',
        text = element_text(family = "Arial")) -> p1
ggsave('results/figs/batch_before.png',p1,width = 6,height = 6)

#批次矫正
library(sva)
batch <- res_df$GSE
expr_combat <- ComBat(dat = as.matrix(df), batch = batch, par.prior = TRUE, prior.plots = FALSE)
res.pca <- prcomp(t(expr_combat))
fviz_pca_ind(res.pca,
             geom = c("point")
)+
  geom_point(aes(
    #shape = factor(res_df$GSE), 
    colour = factor(res_df$GSE)
  ),
  size=5,alpha=0.7)+
  scale_color_manual(
    values=c(
      'GSE130823'='#956EB5',
      'GSE60427'='#6A9750'
      
    )
  )+
  theme_classic(base_size = 20) +
  theme(
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.8
    )
  )+
  theme(legend.position = 'none',
        text = element_text(family = "Arial")) -> p1
ggsave('results/figs/batch_after.png',p1,width = 6,height = 6)


expr_combat_df=as.data.frame(expr_combat)
rownames(expr_combat_df) = expr_mat$GENE_SYMBOL
saveRDS(expr_combat,'data/expr_combat_intergration.rds')
saveRDS(expr_combat_df,'data/expr_combat_intergration_df.rds')
saveRDS(res_df,'data/sample_info_intergration.rds')

##################################consensus analysis
library(ConsensusClusterPlus)
#GSE130823,GSE60427
res_df=readRDS('data/sample_info_intergration.rds')
res_df %>% 
  filter(GSE %in% c('GSE130823','GSE60427')) %>%
  group_by(source_name_ch1) %>% 
  summarise(N=n())
#GSE130823,Gastritis
#GSE60427, remove normal
c(res_df %>% 
    filter(GSE == 'GSE130823' & 
             source_name_ch1 =='Gastritis') %>% 
    pull(geo_accession),
  res_df %>% 
    filter(GSE == 'GSE60427') %>% 
    filter(!grepl('normal',source_name_ch1)) %>% 
    pull(geo_accession)
) -> select_gsm #71

expr_combat_df=readRDS('data/expr_combat_intergration_df.rds')
expr_combat_df[,select_gsm] -> df
vars <- apply(df, 1, sd, na.rm=TRUE)
top_genes <- names(sort(vars, decreasing=TRUE))[1:2000]
mat_sub <- df[top_genes, ]

# optional:  z-score normalization
mat_z <- t(scale(t(mat_sub)))  
out_dir <- "consensus_out_intergration"
dir.create(out_dir, showWarnings = FALSE)
# run consensus clustering 
results = ConsensusClusterPlus(
  mat_z,
  maxK = 5,                
  reps = 2000,             
  pItem = 0.8,             
  pFeature = 1,            
  clusterAlg = "hc",       # "hc" (hierarchical) 
  distance = "pearson",    # "pearson", "spearman", "euclidean"
  seed = 1234,
  plot = "pdf",            
  title = out_dir
)
dev.off()

# select k = 2 for downstream analysis
k <- 2
clust_assign <- results[[k]]$consensusClass  # names = sample names
table(clust_assign)
# extract the consensus matrix results
cons_mat <- results[[k]]$consensusMatrix
pheatmap::pheatmap(cons_mat)
# # PCA
# library(ggplot2)
# library(ggrepel)
# pca <- prcomp(t(mat_z), scale=FALSE)
# df_pca <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], cluster = factor(clust_assign)) %>% mutate(id=rownames(.))
# ggplot(df_pca, aes(PC1, PC2, color=cluster)) + geom_point(size=3)+
#   ggrepel::geom_text_repel(aes(label=id))
# 
as.data.frame(
  clust_assign
) %>% write.csv(.,'consensus_out_intergration/clusters.csv')
clust_assign<-fread('consensus_out_intergration/clusters.csv')
# df_pca <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], cluster = factor(clust_assign$clust_assign)) %>% mutate(id=rownames(.))
# ggplot(df_pca, aes(PC1, PC2, color=cluster)) + geom_point(size=3)+
#   ggrepel::geom_text_repel(aes(label=id))

res_df %>% 
  inner_join(.,
             clust_assign,
             by=c('geo_accession'='V1')
             ) -> cluster_res
#dplyr::select(c('geo_accession','GSE','clust_assign')) -> cluster_res
