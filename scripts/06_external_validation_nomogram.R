
# ============================================================
# Script 06 - External validation, nomogram, ROC (formerly figure3D.R)
# Run from the repository root.
# Inputs: results/rds/candidate_signatures.rds,
#         results/rds/expr_combat_intergration_df.rds,
#         results/rds/sample_info_intergration.rds,
#         consensus_out_intergration/clusters.csv, data/step04_ML.RData
# ============================================================

library(ggplot2)
library(ggpubr)


clust_assign<-fread('consensus_out_intergration/clusters.csv')

candidate_signatures<- readRDS('results/rds/candidate_signatures.rds')
library(GSVA)
expr_combat_df=readRDS('results/rds/expr_combat_intergration_df.rds')
gsva_param=ssgseaParam(
  as.matrix(expr_combat_df), 
  candidate_signatures, 
  #candidate_signatures_overall,
  normalize = TRUE,
  verbose=TRUE
)
ssgsea_res <- gsva(gsva_param)
sample_info = readRDS('results/rds/sample_info_intergration.rds')


as.data.frame(t(ssgsea_res)) %>% 
  mutate(id=rownames(.)) %>% 
  left_join(.,
            sample_info,
            by=c('id'='geo_accession')
  ) %>% 
  filter(GSE=='GSE60427') %>% 
  left_join(.,
            clust_assign,
            by=c('id'='V1')
  ) %>% 
  filter(!is.na(clust_assign)) %>% 
  mutate(
    Group = ifelse(clust_assign==1,'subtypeA','subtypeB')
  ) -> df
df %>% ggplot(aes(x=Group,y=gastritis_subtypeB))+
  geom_boxplot(aes(fill=Group),width=0.5)+
  geom_jitter(width=0.3)+
  scale_fill_manual(
    values = c(
      'subtypeA'='#EA6361',
      'subtypeB'='#557EB6'
    )
  )+
  ggpubr::stat_compare_means(method='t.test',label = "p.format")+
  theme_classic(base_size = 20)+
  labs(x='',y='subtype B \n signature score')+
  theme(legend.position = 'none') -> p1
df %>% ggplot(aes(x=Group,y=gastritis_subtypeA))+
  geom_boxplot(aes(fill=Group),width=0.5)+
  geom_jitter(width=0.3)+
  scale_fill_manual(
    values = c(
      'subtypeA'='#EA6361',
      'subtypeB'='#557EB6'
    )
  )+
  ggpubr::stat_compare_means(method='t.test',label = "p.format")+
  theme_classic(base_size = 20)+
  labs(x='',y='subtype A \n signature score')+
  theme(legend.position = 'none') -> p2
ggsave('figures/GSE60427_boxplot1.png',p1,width = 5,height = 5)
ggsave('figures/GSE60427_boxplot2.png',p2,width = 5,height = 5)
#和疾病特征进行比较
# 指定两两比较
my_comparisons <- list(
  c("mild", "severe"),
  c("severe", "IM"),
  c("mild", "IM")
)
df %>% 
  #filter(grepl('DR',source_name_ch1)) %>% 
  #mutate(Grade=factor(rep(rep(c('mild','severe','IM'),each=4),2),levels = c('mild','severe','IM'))) %>% 
  filter(grepl('DR',source_name_ch1)) %>% 
  mutate(Grade=factor(rep(c('mild','severe','IM'),each=4),levels = c('mild','severe','IM'))) %>% 
  ggplot(
    aes(x=Grade,y=gastritis_subtypeB)
  )+
  geom_boxplot(aes(fill=Grade),width=0.5)+
  geom_jitter(width=0.3)+
  stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox",
    label = "p.format"
  ) +
  scale_fill_manual(values = c(IM="#D95F5F", severe="#4C78A8", mild="#59A96A"))+
  theme_classic(base_size = 20)+
  labs(x='',y='subtype B \n signature score')+
  theme(legend.position = 'none') -> p1
ggsave('results/figs/GSE60427_boxplot1.png',p1,width = 5,height = 5)
ggsave('results/figs/GSE60427_boxplot2.png',p1,width = 5,height = 5)
#GSE130823结果可视化
as.data.frame(t(ssgsea_res)) %>% 
  mutate(id=rownames(.)) %>% 
  left_join(.,
            sample_info,
            by=c('id'='geo_accession')
  ) %>% 
  filter(GSE=='GSE130823') %>%
  left_join(.,
            clust_assign,
            by=c('id'='V1')
  ) %>% 
  #filter(!is.na(clust_assign)) %>% 
  mutate(
    Group = ifelse(clust_assign==1,'subtypeA','subtypeB')
  ) -> df
df %>% 
  filter(!is.na(Group)) %>% 
  ggplot(aes(x=Group,y=gastritis_subtypeB))+
  geom_boxplot(aes(fill=Group),width=0.5)+
  geom_jitter(width=0.3)+
  scale_fill_manual(
    values = c(
      'subtypeA'='#EA6361',
      'subtypeB'='#557EB6'
    )
  )+
  ggpubr::stat_compare_means(method='wilcox',label = "p.format")+
  theme_classic(base_size = 20)+
  labs(x='',y='subtype B \n signature score')+
  theme(legend.position = 'none') -> p1
ggsave('results/figs/GSE130823_boxplot1.png',p1,width = 5,height = 5)
df %>% 
  filter(!is.na(Group)) %>% 
  ggplot(aes(x=Group,y=gastritis_subtypeA))+
  geom_boxplot(aes(fill=Group),width=0.5)+
  geom_jitter(width=0.3)+
  scale_fill_manual(
    values = c(
      'subtypeA'='#EA6361',
      'subtypeB'='#557EB6'
    )
  )+
  ggpubr::stat_compare_means(method='t.test',label = "p.format")+
  theme_classic(base_size = 20)+
  labs(x='',y='subtype A \n signature score')+
  theme(legend.position = 'none') -> p1
ggsave('results/figs/GSE130823_boxplot2.png',p1,width = 5,height = 5)
#同样地，看一看疾病分级的关系
#HGN, LGN, GC, gastritis
my_comparisons <- list(
  c("Gastritis", "LGN"),
  c("Gastritis", "HGN"),
  c("Gastritis", "GC")
)
df %>% 
  mutate(
    Grade = factor(case_when(
      source_name_ch1 =='GastricHighGradeIntraepithelialNeoplasia' ~ 'HGN',
      source_name_ch1 =='GastricLowGradeIntraepithelialNeoplasia' ~ 'LGN',
      source_name_ch1 =='IntestinalGastricCancer' ~ 'GC',
      .default='Gastritis'
    ),levels=c('Gastritis','LGN','HGN','GC') )) %>% 
  ggplot(
    aes(x=Grade,y=gastritis_subtypeB)
  )+
  geom_boxplot(aes(fill=Grade),width=0.5)+
  geom_jitter(width=0.3)+
  stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox",
    label = "p.format"
  ) +
  scale_fill_manual(values = c(
    "Gastritis" = "#557EB6",  # 蓝（最轻）
    "LGN"       = "#7AA6D1",  # 浅蓝（过渡）
    "HGN"       = "#F39A8F",  # 浅红（进展）
    "GC"        = "#EA6361"   # 红（最重）
  ))+
  theme_classic(base_size = 20)+
  labs(x='',y='subtype B \n signature score')+
  theme(legend.position = 'none') -> p1
ggsave('results/figs/GSE130823_boxplot3.png',p1,width = 5,height = 5)


df %>% 
  mutate(
    Grade = factor(case_when(
      source_name_ch1 =='GastricHighGradeIntraepithelialNeoplasia' ~ 'HGN',
      source_name_ch1 =='GastricLowGradeIntraepithelialNeoplasia' ~ 'LGN',
      source_name_ch1 =='IntestinalGastricCancer' ~ 'GC',
      .default='Gastritis'
    ),levels=c('Gastritis','LGN','HGN','GC') )) %>% 
  ggplot(
    aes(x=Grade,y=gastritis_subtypeA)
  )+
  geom_boxplot(aes(fill=Grade),width=0.5)+
  geom_jitter(width=0.3)+
  stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox",
    label = "p.format"
  ) +
  scale_fill_manual(values = c(
    "Gastritis" = "#557EB6",  # 蓝（最轻）
    "LGN"       = "#7AA6D1",  # 浅蓝（过渡）
    "HGN"       = "#F39A8F",  # 浅红（进展）
    "GC"        = "#EA6361"   # 红（最重）
  ))+
  theme_classic(base_size = 20)+
  labs(x='',y='subtype A \n signature score')+
  theme(legend.position = 'none') -> p1
ggsave('results/figs/GSE130823_boxplot4.png',p1,width = 5,height = 5)

#external dataset的应用：GSE60662，GSE55696
#GSE60662
ssgsea_res=fread('results/csv/ssgsea_result.csv',index = 'id')
ssgsea_res %>% 
  #mutate(id=rownames(.)) %>% 
  left_join(.,
            sample_info,
            by=c('id'='geo_accession')
  ) %>% write.csv('results/csv/ssgsea_result_sampleid.csv')
#as.data.frame(t(ssgsea_res)) %>% 
  ssgsea_res %>% 
  #mutate(id=rownames(.)) %>% 
  left_join(.,
            sample_info,
            by=c('id'='geo_accession')
  ) %>% 
  filter(GSE=='GSE60662') -> df
my_comparisons <- list(
  c("mild", "severe"),
  c("severe", "IM"),
  c("mild", "IM")
)
df %>% 
  filter(!grepl('control',source_name_ch1)) %>% 
  mutate(
    Grade = factor(case_when(
      grepl('metaplasia',source_name_ch1) ~ 'IM',
      grepl('mild',source_name_ch1) ~ 'mild',
      grepl('severe',source_name_ch1) ~ 'severe',
    ),levels=c('mild','severe','IM') )) %>% 
  ggplot(
    aes(x=Grade,y=gastritis_subtypeB)
  )+
  geom_boxplot(aes(fill=Grade),width=0.5,outlier.shape = NA,)+
  geom_jitter(width=0.3)+
  stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox",
    label = "p.format"
  ) +
  scale_fill_manual(values = c(IM="#D95F5F", severe="#4C78A8", mild="#59A96A"))+
  theme_classic(base_size = 20)+
  labs(x='',y='subtype B \n signature score')+
  theme(legend.position = 'none') -> p1
ggsave('results/figs//GSE60662_boxplot1.png',p1,width = 5,height = 5)
df %>% 
  filter(!grepl('control',source_name_ch1)) %>% 
  mutate(
    Grade = factor(case_when(
      grepl('metaplasia',source_name_ch1) ~ 'IM',
      grepl('mild',source_name_ch1) ~ 'mild',
      grepl('severe',source_name_ch1) ~ 'severe',
    ),levels=c('mild','severe','IM') )) %>% 
  ggplot(
    aes(x=Grade,y=gastritis_subtypeA)
  )+
  geom_boxplot(aes(fill=Grade),width=0.5,outlier.shape = NA)+
  geom_jitter(width=0.3)+
  stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox.test",#
    label = "p.format"
  ) +
  scale_fill_manual(values = c(IM="#D95F5F", severe="#4C78A8", mild="#59A96A"))+
  theme_classic(base_size = 20)+
  labs(x='',y='subtype A \n signature score')+
  theme(legend.position = 'none') ->p1
ggsave('results/figs/GSE60662_boxplot2.png',p1,width = 5,height = 5)
#是不是可以使用RF模型来个预测呢？

#GSE55696
as.data.frame(t(ssgsea_res)) %>% 
  mutate(id=rownames(.)) %>% 
  left_join(.,
            sample_info,
            by=c('id'='geo_accession')
  ) %>% 
  filter(GSE=='GSE55696') -> df
my_comparisons <- list(
  c("Gastritis", "LGN"),
  c("Gastritis", "HGN"),
  c("Gastritis", "EGC")
)
df %>% 
  mutate(
    Grade = factor(case_when(
      source_name_ch1 =='HGD' ~ 'HGN',
      source_name_ch1 =='LGD' ~ 'LGN',
      source_name_ch1 =='EGC' ~ 'EGC',
      .default='Gastritis'
    ),levels=c('Gastritis','LGN','HGN','EGC') ))  %>% 
  #filter(Grade!='GC') %>% 
  ggplot(
    aes(x=Grade,y=gastritis_subtypeB)
  )+
  geom_boxplot(aes(fill=Grade),width=0.5,outlier.shape = NA)+
  geom_jitter(width=0.3)+
  stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox",#wilcox.test,t.test
    label = "p.format"
  ) +
  scale_fill_manual(values = c(
    "Gastritis" = "#557EB6",  # 蓝（最轻）
    "LGN"       = "#7AA6D1",  # 浅蓝（过渡）
    "HGN"       = "#F39A8F",  # 浅红（进展）
    "EGC"        = "#EA6361"   # 红（最重）
  ))+
  theme_classic(base_size = 20)+
  labs(x='',y='subtype B \n signature score')+
  theme(legend.position = 'none') -> p1
ggsave('results/figs/GSE55696_boxplot1.png',p1,width = 5,height = 5)
df %>% 
  mutate(
    Grade = factor(case_when(
      source_name_ch1 =='HGD' ~ 'HGN',
      source_name_ch1 =='LGD' ~ 'LGN',
      source_name_ch1 =='EGC' ~ 'EGC',
      .default='Gastritis'
    ),levels=c('Gastritis','LGN','HGN','EGC') ))  %>% 
  #filter(Grade!='GC') %>% 
  ggplot(
    aes(x=Grade,y=gastritis_subtypeA)
  )+
  geom_boxplot(aes(fill=Grade),width=0.5,outlier.shape = NA)+
  geom_jitter(width=0.3)+
  stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox.test",#wilcox.test,t.test
    label = "p.format"
  ) +
  scale_fill_manual(values = c(
    "Gastritis" = "#557EB6",  # 蓝（最轻）
    "LGN"       = "#7AA6D1",  # 浅蓝（过渡）
    "HGN"       = "#F39A8F",  # 浅红（进展）
    "EGC"        = "#EA6361"   # 红（最重）
  ))+
  theme_classic(base_size = 20)+
  labs(x='',y='subtype A \n signature score')+
  theme(legend.position = 'none') -> p2
ggsave('results/figs//GSE55696_boxplot2.png',p2,width = 5,height = 5)

#是否可以考虑基于基因建立一个nanogram模型呢？
model_glmnet
m=load('data/step04_ML.RData')
#3个模型共有的genes
intersect(
  lasso_genes,
  intersect(
    svm_top,rf_top
  )
) #7个genes

#"AGXT2L1"  "C11orf86" "IGF1R"    "FAR2"     "MBIP"     "CLIC6"    "ZCWPW1" 
#用这些基因建立一个lasso模型，然后nomogram可视化即可
library(rms)
#nomogram(lasso_model)
expr_mat_matrix=as.data.frame(t(expr_mat_df))
expr_mat_matrix$Group = ifelse(sample_info$clust_assign==2,1,0)#subtypeB=1, subtypeA=0
fit <- lrm(
  Group ~ AGXT2L1 + C11orf86 + IGF1R + FAR2 + MBIP + CLIC6 + ZCWPW1,
  data = expr_mat_matrix,
  x = TRUE,
  y = TRUE,
  penalty = 1 #这个参数很有用
)
fit
dd <- datadist(expr_mat_matrix)
options(datadist = "dd")
nom <- nomogram(
  fit,
  fun = plogis,
  fun.at = c(0.1, 0.3, 0.5, 0.7, 0.9),
  funlabel = "Risk probability"
)
plot(nom)
plot(calibrate(fit))
plot(fit)
# pdf("nomogram_7gene_model.pdf", width = 12, height = 8)
# plot(
#   nom,
#   xfrac = 0.2,
#   cex.axis = 1.5,
#   cex.var = 1.5,
#   cex.main = 1.5,
#   lmgp = 0.25
# )
# dev.off()
png(
  "nomogram_7gene_model.png",
  width = 4200,
  height = 3000,
  res = 300
)
plot(
  nom,
  xfrac = 0.2,
  cex.axis = 2,
  cex.var = 2,
  cex.main = 2,
  lmgp = 0.25
)
dev.off()
png(
  "calibration_curve.png",
  width = 4000,
  height = 4000,
  res = 800
)
plot(calibrate(fit))
dev.off()


library(pROC)
pred_prob <- predict(fit, type = "fitted")
roc_obj <- roc(expr_mat_matrix$Group, pred_prob)
png(
  "figures/logistic_roc.png",
  width = 5000,
  height = 4000,
  res = 800
)
plot(
  roc_obj,
  print.auc = TRUE,
  main = "ROC curve of logistic model"
)
dev.off()

##"AGXT2L1"  "C11orf86" "IGF1R"    "FAR2"     "MBIP"     "CLIC6"    "ZCWPW1"
genes
library(pROC)
roc_AGXT2L1 <- roc(expr_mat_matrix$Group, expr_mat_matrix$AGXT2L1)
roc_C11orf86 <- roc(expr_mat_matrix$Group, expr_mat_matrix$C11orf86)
roc_IGF1R <- roc(expr_mat_matrix$Group, expr_mat_matrix$IGF1R)
roc_FAR2 <- roc(expr_mat_matrix$Group, expr_mat_matrix$FAR2)
roc_MBIP <- roc(expr_mat_matrix$Group, expr_mat_matrix$MBIP)
roc_CLIC6 <- roc(expr_mat_matrix$Group, expr_mat_matrix$CLIC6)
roc_ZCWPW1 <- roc(expr_mat_matrix$Group, expr_mat_matrix$ZCWPW1)
auc(roc_lasso)

png(
  "gene_roc.png",
  width = 5000,
  height = 4000,
  res = 800
)
plot(roc_AGXT2L1,col = "#4E79A7",
     xlim = c(1,0),
     ylim = c(0, 1),
     legacy.axes = FALSE,
     main = "") 
lines(roc_C11orf86, col = "#F28E2B")
lines(roc_IGF1R, col = "#59A14F")
lines(roc_FAR2, col = "#E15759")
lines(roc_MBIP, col = "#B07AA1")
lines(roc_CLIC6, col = "#76B7B2")
lines(roc_ZCWPW1, col = "#EDC948")
auc(roc_AGXT2L1) #0.985
auc(roc_C11orf86) #0.979
auc(roc_IGF1R) #0.986
auc(roc_FAR2) #0.965
auc(roc_MBIP) #0.971
auc(roc_CLIC6) #0.992
auc(roc_ZCWPW1) #0.968
dev.off()
gene_colors <- c(
  "AGXT2L1"  = "#4E79A7",  # muted blue
  "C11orf86" = "#F28E2B",  # orange
  "IGF1R"    = "#59A14F",  # green
  "FAR2"     = "#E15759",  # red
  "MBIP"     = "#B07AA1",  # purple
  "CLIC6"    = "#76B7B2",  # cyan
  "ZCWPW1"   = "#EDC948"   # yellow
)

expr_mat_matrix[,genes]  %>% t() %>% as.data.frame() ->tmp
sample_info %>%
  arrange(desc(Group)) %>%
  dplyr::select(c('Group')) %>% as.data.frame()-> annotation_col
rownames(annotation_col) = sample_info %>% arrange(desc(Group)) %>% pull(1)
pheatmap::pheatmap(
  tmp[,rownames(annotation_col)],
  annotation_col = annotation_col,
  cluster_cols = F,
  treeheight_row = 0,
  show_colnames = F,
  scale = 'row',
  annotation_colors = list(
    Group=c(subtypeA='#557EB6',subtypeB='#EA6361')
  ),
  gaps_col = c(37),
  cellwidth = 4,
  cellheight = 30,
  filename = 'figures/7genes.png'
)
dev.off()


