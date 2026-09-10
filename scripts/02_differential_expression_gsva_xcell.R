
# ============================================================
# Script 02 - Differential expression, GSVA, and xCell
# (formerly figure2.R). Run from the repository root.
# Uses objects produced by Script 01 (expr_mat, res_df, sample_info,
# clust_assign) and sources Script 03 for ORA/GSEA enrichment.
# ============================================================

# setwd('D:\\document\\20260425_pujiang_project\\submit\\PLOS One\\revision_v4')
##################################comparative analysis
library(limma)
#1: A
#2: B
expr_mat %>% dplyr::select(-c('GENE_SYMBOL')) %>% as.data.frame()-> df
rownames(df) =expr_mat$GENE_SYMBOL

clust_assign<-fread('consensus_out_intergration/clusters.csv')
sample_info=data.frame(
  clust_assign=clust_assign$clust_assign
) %>% mutate(
  Group = ifelse(clust_assign==1,'subtypeA','subtypeB')
)
rownames(sample_info)=clust_assign$V1
design <- model.matrix(
  ~ 0 + sample_info$Group
)
colnames(design) = c('subtypeA','subtypeB')

fit <- lmFit(df, design)
cm <- makeContrasts(subtypeBVSsubtypeA=subtypeB - subtypeA, levels=c('subtypeA','subtypeB'))
fit2 <- contrasts.fit(fit, cm)
fit2 <- eBayes(fit2)
res <- topTable(fit2, number=Inf)
head(res)
saveRDS(res,'results/rds/degs.rds')
res %>% write.csv(.,'results/csv/degs.csv')
res %>% 
  mutate(Group=ifelse(adj.P.Val<0.05,
                      ifelse(
                        logFC > 0.5,'Up',
                        ifelse(logFC < -0.5,'Down','Non')
                      ),'Non'
  )) %>% 
  ggplot(aes(x=logFC,y=-log10(adj.P.Val)))+
  geom_point(aes(color=Group))+
  scale_color_manual(
    values=c(
      'Down'='#557EB6',
      'Up'='#EA6361',
      'Non'='grey'
    )
  )+
  geom_hline(yintercept = c(-log10(0.05)),lty=2)+
  geom_vline(xintercept = c(-0.5,0.5),lty=2)+
  theme_classic(base_size = 20)+
  theme(legend.position = 'none',
        text = element_text(family = "Arial"))+
  labs(x='Log2 Fold Change',y='-Log10 adjust P value') -> p
ggsave('results/figs/degs_volcano.png',p,width = 6,height = 6)
res %>% 
  mutate(Group=ifelse(adj.P.Val<0.05,
                      ifelse(
                        logFC > 0.5,'Up',
                        ifelse(logFC < -0.5,'Down','Non')
                      ),'Non'
  )) %>% group_by(Group) %>% summarise(N=n())
#Down   1421
#Up     1532
res %>%   filter(adj.P.Val<0.05) %>% 
  filter(abs(logFC)>0.5) %>% 
  mutate(id=rownames(.)) -> degs 
nrow(degs)


####heatmap
library(pheatmap)
res %>% 
  filter(adj.P.Val<0.05) %>% 
  filter(abs(logFC) >1) %>% 
  mutate(id=rownames(.)) %>% 
  pull(id) -> x1
sample_info %>% arrange(desc(Group)) %>% 
  mutate(id=rownames(.)) %>% 
  left_join(.,
            res_df,
            by=c('id'='geo_accession')
  ) %>% 
  dplyr::select(c(2,3,6)) -> tmp
tmp %>% dplyr::select(c('Group','GSE')) %>% setnames(c('Subtype','Dataset')) %>% 
  as.data.frame() -> annotation_col
rownames(annotation_col) = tmp %>% pull(id)
pheatmap(
  scale(df[x1,tmp$id]),
  annotation_col = annotation_col,
  cluster_cols = F,
  scale = 'row',
  show_rownames = F,
  show_colnames = F,
  annotation_colors = list(
    Subtype=c(subtypeA='#557EB6',subtypeB='#EA6361'),
    Dataset=c(GSE130823='#956EB5',GSE60427='#6A9750')
  ),
  treeheight_row = 0,
  fontsize = 15,
  filename = 'results/figs//degs_heatmap.png'
) 
dev.off()


res %>% mutate(gene=rownames(.)) %>% 
  dplyr::select(c(gene,logFC,P.Value)) %>% 
  write.csv(.,'results/csv/deg_list.csv',row.names = F)
source('scripts/03_enrichment_ORA_GSEA.R')

go_bp_up=readRDS('results/rds/BP_ego_up.rds')
go_bp_up@result %>% head(5) %>% 
  mutate(id=forcats::fct_reorder(
    stringr::str_wrap(stringr::str_to_title(Description),width = 25),
    dplyr::desc(pvalue)
  )
  ) %>% 
  ggplot(aes(x=-log10(p.adjust),y=id))+
  geom_bar(stat='identity',
           width=0.5,
           fill='#EA6361',
           color='black')+
  theme_classic(base_size = 15)+
  theme(axis.text.y = element_text(hjust = 0))+
  labs(x='-Log10 adjust P Value',y='') -> p1
ggsave('results/figs/go_bp_up_barplot.png',p1,width = 5,height = 3)
kegg_up = readRDS('results/rds/KEGG_ego_up.rds')
kegg_up@result %>% 
  filter(ID %in% c('hsa04110','hsa04115','hsa04060','hsa04657','hsa04979')) %>% 
  mutate(id=forcats::fct_reorder(
    stringr::str_wrap(stringr::str_to_title(Description),width = 30),
    dplyr::desc(pvalue)
  )
  ) %>% 
  ggplot(aes(x=-log10(p.adjust),y=id))+
  geom_bar(stat='identity',
           width=0.5,
           fill='#EA6361',
           color='black')+
  theme_classic(base_size = 15)+
  theme(axis.text.y = element_text(hjust = 0))+
  labs(x='-Log10 adjust P Value',y='')-> p1
ggsave('results/figs/go_kegg_up_barplot.png',p1,width = 5,height = 3)


go_bp_down=readRDS('results/rds//BP_ego_down.rds')
go_bp_down@result %>% head(5) %>% 
  mutate(id=forcats::fct_reorder(
    stringr::str_wrap(stringr::str_to_title(Description),width = 30),
    dplyr::desc(pvalue)
  )
  ) %>% 
  ggplot(aes(x=-log10(p.adjust),y=id))+
  geom_bar(stat='identity',
           width=0.5,
           fill='#557EB6',
           color='black')+
  theme_classic(base_size = 15)+
  theme(axis.text.y = element_text(hjust = 0))+
  labs(x='-Log10 adjust P Value',y='') ->p1
ggsave('results/figs/go_bp_down_barplot.png',p1,width = 5,height = 3)

kegg_down = readRDS('results/rds/KEGG_ego_down.rds')
kegg_down@result %>% 
  filter(ID %in% c('hsa04971','hsa04081','hsa04974','hsa04972','hsa04082')) %>% 
  mutate(id=forcats::fct_reorder(
    stringr::str_wrap(stringr::str_to_title(Description),width = 30),
    dplyr::desc(pvalue)
  )
  ) %>% 
  ggplot(aes(x=-log10(p.adjust),y=id))+
  geom_bar(stat='identity',
           width=0.5,
           fill='#557EB6',
           color='black')+
  theme_classic(base_size = 15)+
  theme(axis.text.y = element_text(hjust = 0))+
  labs(x='-Log10 adjust P Value',y='') -> p1
ggsave('results/figs/go_kegg_down_barplot.png',p1,width = 5,height = 3)

##################################GSVA enrichment
m=load('data/dat.RData')
lapply(
  c('REACTOME_ION_CHANNEL_TRANSPORT','REACTOME_POTASSIUM_CHANNELS','REACTOME_CELL_CYCLE'),
  FUN=function(x){
    c2 %>% filter(gs_name==x) %>% pull(gene_symbol)
  }
) -> c2_gs
names(c2_gs)=c('REACTOME_ION_CHANNEL_TRANSPORT','REACTOME_POTASSIUM_CHANNELS','REACTOME_CELL_CYCLE')
lapply(
  c('GOBP_GASTRIC_ACID_SECRETION','GOBP_INFLAMMATORY_RESPONSE','GOBP_POTASSIUM_ION_TRANSPORT',
    'GOMF_POTASSIUM_CHANNEL_ACTIVITY'),
  FUN=function(x){
    c5 %>% filter(gs_name==x) %>% pull(gene_symbol)
  }
) -> c5_gs
names(c5_gs)=c('GOBP_GASTRIC_ACID_SECRETION','GOBP_INFLAMMATORY_RESPONSE','GOBP_POTASSIUM_ION_TRANSPORT',
               'GOMF_POTASSIUM_CHANNEL_ACTIVITY')
lapply(
  c('HALLMARK_E2F_TARGETS ','HALLMARK_INTERFERON_GAMMA_RESPONSE'),
  FUN=function(x){
    hallmark %>% filter(gs_name==x) %>% pull(gene_symbol)
  }
) -> hallmark_gs
names(hallmark_gs)=c('HALLMARK_E2F_TARGETS ','HALLMARK_INTERFERON_GAMMA_RESPONSE')
m=load('data/kegg_hsa.RData')
lapply(
  c('hsa04110','hsa04657','hsa04668','hsa04971','hsa00190'),
  FUN=function(x){
    pathway2genesymbol %>% filter(from == x) %>% pull(SYMBOL)
  }
) -> kegg_gs
names(kegg_gs)=c('Cell cycle','IL-17 signaling pathway','TNF signaling pathway','Gastric acid secretion','Oxidative phosphorylation')
select_gs = c(c2_gs,c5_gs,hallmark_gs,kegg_gs)
library(GSVA)
gsva_param=gsvaParam(
  as.matrix(df), 
  select_gs, 
  verbose=TRUE
)
gsva_score <- gsva(gsva_param)
write.csv(gsva_score,'results/csv/gsva_score_figure2.csv')
sample_info %>% arrange(desc(Group)) %>% 
  dplyr::select(-c(1)) %>% setnames('Subtype')-> annotation_col
pheatmap(gsva_score[,rownames(annotation_col)],
         annotation_col = annotation_col,
         show_colnames = F,
         gaps_col = c(39),
         cutree_rows  = 2,
         treeheight_row = 0,
         annotation_colors = list(
           Subtype=c(subtypeA='#557EB6',subtypeB='#EA6361')
         ),
         border_color = "NA",
         cluster_cols = F,
         filename = 'results/figs/subtypeA_subtypeB_gsva_heatmap.png',
         height = 3,
         width = 8,
         scale = 'row'
)
dev.off()

##################################xCell 分析
library(xCell)
xcell_res = xCellAnalysis(df,rnaseq = F,scale = F)
apply(
  xcell_res,
  1,
  FUN = function(x){
    wilcox.test(as.numeric(x)~sample_info$Group)
  }
) -> tmp #67

tmp[
  unlist(lapply(
    tmp,
    FUN=function(x){
      x$p.value<0.05
    }
  ))
] -> tmpx #30


xcell_res[names(tmpx),] %>% t() %>% 
  as.data.frame() %>% 
  mutate(id=rownames(.)) %>% 
  left_join(.,
            sample_info %>% mutate(id=rownames(.)) %>% dplyr::select(-c('clust_assign'))
  )-> xcell_res_select
select_types=c(
  'CD4+ naive T-cells',
  'CD8+ T-cells',
  'Endothelial cells',#（组织稳态）
  'NKT',
  'Th2 cells',#（免疫调控）
  'Neutrophils',#（炎症）
  'Macrophages M1',#（修复/重塑）
  'pro B-cells'#（体液免疫）
)
lapply(
  select_types,
  FUN=function(ccl){
    xcell_res_select %>% 
      ggplot(aes(x=Group,y=.data[[ccl]]))+
      geom_boxplot(aes(fill=Group),width=0.4)+
      scale_fill_manual(
        values = c(
          'subtypeA'='#557EB6',
          'subtypeB'='#EA6361'
        ))+
      ggpubr::stat_compare_means(method='wilcox',label = 'p.format',size=5)+
      geom_jitter(width=0.2,size=2)+
      labs(x='',y='',title=sprintf('%s',ccl))+
      theme_classic(base_size = 15)+
      theme(legend.position = 'none',
            text = element_text(family = "Arial"),
            plot.title = element_text(size = 15),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank())
  }
) -> ps
library(patchwork)
#wrap_plots(ps, nrow = 2, ncol = 5)
ggsave("results/figs/xcell_boxplots.png",
       wrap_plots(ps, nrow = 2, ncol = 4),
       width = 14,
       height = 7)
xcell_res_select %>% 
  write.csv(.,'results/csv/xcell_results.csv')

