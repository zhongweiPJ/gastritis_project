
# ============================================================
# Script 03 - ORA and GSEA enrichment (formerly enrichment_ORA_GSEA.R)
# Run from the repository root; also sourced by Script 02.
# Input: results/csv/deg_list.csv (written by Script 02).
# ============================================================

library(clusterProfiler)
library(org.Hs.eg.db)
library(data.table)
library(dplyr)

library(msigdbr)
c2=readRDS('data/c2.rds')
c5=readRDS('data/c5.rds')
hallmark=readRDS('data/hallmark.rds')

pathway2gene=readRDS('data/pathway2gene.rds')
pathway2name=readRDS('data/pathway2name.rds')
pathway2gene_symbol=readRDS('data/pathway2gene_symbol.rds')

fread('results/csv/deg_list.csv') %>% 
  setnames(c('gene','log2fc','pvalue')) %>% 
  filter(pvalue<0.05) %>% 
  filter(abs(log2fc)>0.5)  %>% 
  distinct(.,gene,.keep_all = T) %>% 
  arrange(desc(log2fc))-> tmp

lapply(
  c('BP','CC','MF','KEGG'),
  FUN = function(term){
    if(term != 'KEGG'){
      ego <- enrichGO(
        gene = tmp %>% filter(log2fc>0) %>% pull(gene),
        OrgDb = 'org.Hs.eg.db',
        keyType = 'SYMBOL',
        ont = term,
        pvalueCutoff = 1
      )
      saveRDS(object = ego,file = paste0('results/rds/',term,'_ego_up.rds'))
    }else{
      ego <- enricher(
        gene = tmp %>% filter(log2fc>0) %>% pull(gene),
        TERM2GENE = pathway2gene_symbol,
        TERM2NAME = pathway2name
      )
      saveRDS(object = ego,file = paste0('results/rds/',term,'_ego_up.rds'))
    }
    return(ego@result %>% mutate(term_type=term))
  }
) -> go_up_list
rbindlist(go_up_list) %>% write.csv(.,'results/csv/go_kegg_up_results.csv')

#down
lapply(
  c('BP','CC','MF','KEGG'),
  FUN = function(term){
    if(term != 'KEGG'){
      ego <- enrichGO(
        gene = tmp %>% filter(log2fc<0) %>% pull(gene),
        OrgDb = 'org.Hs.eg.db',
        keyType = 'SYMBOL',
        ont = term,
        pvalueCutoff = 1
      )
      saveRDS(object = ego,file = paste0('results/rds/',term,'_ego_down.rds'))
    }else{
      ego <- enricher(
        gene = tmp %>% filter(log2fc<0) %>% pull(gene),
        TERM2GENE = pathway2gene_symbol,
        TERM2NAME = pathway2name
      )
      saveRDS(object = ego,file = paste0('results/rds/',term,'_ego_down.rds'))
    }
    return(ego@result %>% mutate(term_type=term))
  }
) -> go_down_list
rbindlist(go_down_list) %>% write.csv(.,'results/csv/go_kegg_down_results.csv')

#overall
lapply(
  c('BP','CC','MF','KEGG'),
  FUN = function(term){
    if(term != 'KEGG'){
      ego <- enrichGO(
        gene = tmp %>% pull(gene),
        OrgDb = 'org.Hs.eg.db',
        keyType = 'SYMBOL',
        ont = term,
        pvalueCutoff = 1
      )
      saveRDS(object = ego,file = paste0('results/rds/',term,'_ego_overall.rds'))
    }else{
      ego <- enricher(
        gene = tmp %>% pull(gene),
        TERM2GENE = pathway2gene_symbol,
        TERM2NAME = pathway2name
      )
      saveRDS(object = ego,file = paste0('results/rds/',term,'_ego_overall.rds'))
    }
    return(ego@result %>% mutate(term_type=term))
  }
) -> go_overall_list
rbindlist(go_overall_list) %>% write.csv(.,'results/csv/go_kegg_overall_results.csv')


#尝试使用clusterprofile， org.Hs.eg.db中的GO数据进行GSEA
tmp %>% pull(log2fc) -> deg_genes
tmp %>% pull(gene) -> names(deg_genes)
go_df=readRDS('data/go_df.rds')
lapply(
  c('BP','CC','MF','KEGG'),
  FUN = function(term){
    if(term != 'KEGG'){
      gsea <- GSEA(
        geneList = deg_genes,
        TERM2GENE = switch(EXPR = term,
                           BP=go_df %>% filter(ONTOLOGY =='BP') %>% dplyr::select(c(GOID,SYMBOL)),
                           CC=go_df %>% filter(ONTOLOGY =='CC') %>% dplyr::select(c(GOID,SYMBOL)),
                           MF=go_df %>% filter(ONTOLOGY =='MF') %>% dplyr::select(c(GOID,SYMBOL)),
        ),
        TERM2NAME = switch(
          EXPR = term,
          BP=go_df %>% filter(ONTOLOGY =='BP') %>% dplyr::select(c(GOID,TERM)),
          CC=go_df %>% filter(ONTOLOGY =='CC') %>% dplyr::select(c(GOID,TERM)),
          MF=go_df %>% filter(ONTOLOGY =='MF') %>% dplyr::select(c(GOID,TERM)),
        ),
        maxGSSize = 1001,
        minGSSize = 2,
        pvalueCutoff = 1
      )
      saveRDS(object = gsea,file = paste0('results/rds/',term,'_gsea_go_org.hs.eg.db.rds'))
    }else{
      gsea <- GSEA(
        geneList = deg_genes,
        TERM2GENE = pathway2gene_symbol,
        TERM2NAME = pathway2name,
        maxGSSize = 1001,
        minGSSize = 2,
        pvalueCutoff = 1
      )
      saveRDS(object = gsea,file = paste0('results/rds/',term,'_gsea_kegg_org.hs.eg.db.rds'))
    }
    return(gsea@result %>% mutate(Id=row_number(),term_type=term))
  }
)-> gsea_list
rbindlist(gsea_list) %>% write.csv(.,'results/csv/gsea_GO_KEGG_org.hs.eg.db_results.csv')


###########################################GSEA
#C2/C5/hallmark/kegg
tmp %>% pull(log2fc) -> deg_genes
tmp %>% pull(gene) -> names(deg_genes)

lapply(
  c('c2','c5','hallmark','KEGG'),
  FUN = function(term){
    if(term != 'KEGG'){
      gsea <- GSEA(
        geneList = deg_genes,
        TERM2GENE = switch(EXPR = term,
                           c2=c2 %>% dplyr::select(c(gs_name,gene_symbol)),
                           c5=c5 %>% dplyr::select(c(gs_name,gene_symbol)),
                           hallmark=hallmark %>% dplyr::select(c(gs_name,gene_symbol)),
        ),
        maxGSSize = 1001,
        minGSSize = 2,
        pvalueCutoff = 1
      )
      saveRDS(object = gsea,file = paste0('results/rds/',term,'_gsea_c2c5hallmark.rds'))
    }else{
      gsea <- GSEA(
        geneList = deg_genes,
        TERM2GENE = pathway2gene_symbol,
        TERM2NAME = pathway2name,
        maxGSSize = 1001,
        minGSSize = 2,
        pvalueCutoff = 1
      )
      saveRDS(object = gsea,file = paste0('results/rds/',term,'_gsea_kegg.rds'))
    }
    return(gsea@result %>% mutate(Id=row_number(),term_type=term))
  }
) -> gsea_list

rbindlist(gsea_list) %>% write.csv(.,'results/csv/gsea_c2_c5_hall_kegg_results.csv')
