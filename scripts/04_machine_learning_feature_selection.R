
# ============================================================
# Script 04 - Machine-learning feature prioritization
# (LASSO / SVM / random forest; formerly figure3.R)
# Run from the repository root.
# Inputs: data/expr_combat_intergration_df.rds,
#         data/sample_info_intergration.rds,
#         consensus_out_intergration/clusters.csv,
#         results/csv/deg_list.csv
# ============================================================

#figure3
library(caret)
library(pROC)
##############
expr_mat <- readRDS('data/expr_combat_intergration_df.rds')
sample_info  <- readRDS('data/sample_info_intergration.rds')
#GSE130823, GSE60427
c(sample_info %>% 
    filter(GSE == 'GSE130823' & 
             source_name_ch1 =='Gastritis') %>% 
    pull(geo_accession),
  sample_info %>% 
    filter(GSE == 'GSE60427') %>% 
    filter(!grepl('normal',source_name_ch1)) %>% 
    pull(geo_accession)
) -> select_gsm #71

sample_info<-fread('consensus_out_intergration/clusters.csv') %>% 
  mutate(
    Group = ifelse(clust_assign==1,'subtypeA','subtypeB')
  ) 
expr_mat_df = expr_mat[,sample_info$V1] #gene(29388) x sample(71)
expr_mat_df = expr_mat_df %>%mutate(id=rownames(.)) %>%  filter(!grepl('LOC',id))
#high variable genes
var_per_gene <- apply(expr_mat_df, 1, var)
length(var_per_gene[var_per_gene==0])  # check zero-variance genes
topN <- 3000
top_genes <- names(sort(var_per_gene, decreasing = TRUE))[1:topN]
ml_data <- expr_mat_df[top_genes, , drop=FALSE]
expr <- as.matrix(ml_data)
label <- factor(sample_info$Group)
# invert matrix to: sample × gene
data <- t(expr)
# data split: 7:3
set.seed(123)
trainIndex <- createDataPartition(label, p = 0.7, list = FALSE)
train_data <- data[trainIndex, ]
test_data  <- data[-trainIndex, ]
train_label <- label[trainIndex]
test_label  <- label[-trainIndex]

train_ctrl <- trainControl(method = "repeatedcv",
                           number = 10,
                           repeats = 5,
                           classProbs = TRUE,
                           summaryFunction = twoClassSummary,
                           savePredictions = TRUE,
                           verboseIter = FALSE)
train_label <- factor(train_label, levels = c("subtypeA","subtypeB"))
test_label  <- factor(test_label,  levels = c("subtypeA","subtypeB"))

# ------------- model training --------------
# 1) Logistic regression with LASSO (glmnet)
set.seed(123)
grid_glmnet <- expand.grid(alpha = 1,    # alpha = 1，LASSO；alpha = 0，Ridge；alpha = 0.5，elastic。
                           lambda = 10^seq(-4, 1, length = 50))
model_glmnet <- train(x = train_data, 
                      y = train_label,
                      method = "glmnet",
                      metric = "ROC",
                      trControl = train_ctrl,
                      tuneGrid = grid_glmnet)

print(model_glmnet)
# best model
best_glmnet <- model_glmnet$bestTune
# LASSO non-zero gene effect
best_lambda <- model_glmnet$bestTune$lambda
coef_lasso <- coef(model_glmnet$finalModel, s = best_lambda)
lasso_genes <- rownames(coef_lasso)[as.numeric(coef_lasso) != 0]
lasso_genes <- setdiff(lasso_genes, "(Intercept)")
length(lasso_genes)
lasso_genes


#2） SVM
set.seed(123)
svm_grid <- expand.grid(
  sigma = c(0.001, 0.005, 0.01, 0.05, 0.1),
  C = c(0.25, 0.5, 1, 2, 4, 8)
)
svm_model <- train(
  x = train_data,
  y = train_label,
  method = "svmRadial",
  metric = "ROC",
  trControl = train_ctrl,
  tuneGrid = svm_grid
)
svm_model$bestTune

# SVM importance
svm_imp <- varImp(svm_model, scale = FALSE)
# 如果没有 Overall，就用 A/B 两列的平均值作为综合重要性
if (!"Overall" %in% colnames(svm_imp$importance)) {
  svm_imp$importance$Overall <- rowMeans(svm_imp$importance)
}
svm_ranked <- rownames(svm_imp$importance)[
  order(svm_imp$importance$Overall, decreasing = TRUE)
]
svm_top <- svm_ranked[1:50]
length(svm_top)
svm_top

# 3) Random Forest
set.seed(123)
grid_rf <- expand.grid(mtry = c(5, 10, 30, 50, 100))  # mtry 
grid_rf <- expand.grid(mtry = c(5, 10, 30, 50, 100))
model_rf <- train(x = train_data, 
                  y = train_label,
                  method = "rf",
                  metric = "ROC",
                  trControl = train_ctrl,
                  tuneGrid = grid_rf,
                  importance = TRUE,
                  ntree = 500)   # 指定树的数量
model_rf$bestTune

# 提取RF重要基因
rf_imp <- varImp(model_rf, scale = FALSE)
if (!"Overall" %in% colnames(rf_imp$importance)) {
  rf_imp$importance$Overall <- rowMeans(rf_imp$importance)
}
rf_ranked <- rownames(rf_imp$importance)[
  order(rf_imp$importance$Overall, decreasing = TRUE)
]

rf_top <- rf_ranked[1:50]

length(rf_top)
rf_top

save(
  model_glmnet,
  svm_model,
  model_rf,
  file = 'data/figure3_ML.RData'
)

#########visualization importance
svm_imp$importance %>% head()
rf_imp$importance %>% head()
coef_lasso

data.frame(gene=coef_lasso[lasso_genes,]) %>% 
  mutate(
    lasso_importance=ifelse(
      (gene-min(gene)) / (max(gene)-min(gene))==0,0.1,(gene-min(gene)) / (max(gene)-min(gene))
    )
  ) %>% mutate(
    id=rownames(.)
  ) %>% dplyr::select(c(id,lasso_importance))-> df1
svm_imp$importance %>% 
  mutate(
    id=rownames(.)
  ) %>% 
  filter(id %in% svm_top) %>% 
  mutate(
    svm_importance=ifelse(
      (Overall-min(Overall)) / (max(Overall)-min(Overall))==0,0.1,(Overall-min(Overall)) / (max(Overall)-min(Overall))
    )) %>% 
  dplyr::select(c(id,svm_importance))  -> df2

rf_imp$importance %>% 
  mutate(
    id=rownames(.)
  ) %>% 
  filter(id %in% rf_top) %>% 
  mutate(
    rf_importance=ifelse(
      (Overall-min(Overall)) / (max(Overall)-min(Overall))==0,0.1,(Overall-min(Overall)) / (max(Overall)-min(Overall))
    )) %>% 
  dplyr::select(c(id,rf_importance))  -> df3

df1 %>% 
  full_join(.,
            df2
  ) %>% 
  full_join(.,
            df3
  ) %>% mutate_all(~replace(., is.na(.), 0))-> df
# library(pheatmap)
# df0=df %>% dplyr::select(-c(id))
# rownames(df0)=df %>% pull(id)
# pheatmap(t(df0),
#          color = colorRampPalette(c("white", "darkred"))(100),
#          cluster_rows = F,
#          cluster_cols = F,
#          cellheight  = 35,
#          cellwidth = 8,
#          filename = 'test.png')
# dev.off()

dfx=df[,-c(1)]
rownames(dfx)=df$id
mat <- dfx %>% mutate(id=rownames(.)) %>% filter(!grepl('LOC',id)) %>% dplyr::select(-c('id'))
library(circlize)
library(ComplexHeatmap)
library(RColorBrewer)
circos.clear()
col_fun <- colorRamp2(
  c(0, 0.5, 1),
  c("white", "#F4A6A6", "#B2182B")
)

pdf("results/figs/circular_importance_heatmap_with_labels.pdf", width = 5, height = 5)
circos.par(
  start.degree = 90,
  gap.degree = 1,
  track.margin = c(0.01, 0.01),
  points.overflow.warning = FALSE
)
circos.heatmap(
  mat,
  col = col_fun,
  cluster = FALSE,
  show.sector.labels = TRUE,
  rownames.side = "outside",
  rownames.cex = 0.8,
  cell.border = "grey80",
  track.height = 0.25
)
grid.text(
  "Inner: LASSO\nMiddle: SVM\nOuter: RF",
  x = unit(0.5, "npc"),
  y = unit(0.5, "npc"),
  gp = gpar(fontsize = 11, fontface = "bold")
)
dev.off()
circos.clear()
#######################################candidate signatures
ml_features = unique(c(
  lasso_genes,
  svm_top,
  rf_top
))

fread('results/csv/deg_list.csv') %>%
  filter(gene %in% ml_features) %>% 
  filter(!grepl('LOC',gene)) %>% 
  filter(logFC >0) %>% pull(gene) -> signature_subtypeB

fread('results/csv/deg_list.csv') %>% 
  filter(gene %in% ml_features) %>% 
  filter(!grepl('LOC',gene)) %>% 
  filter(logFC <0) %>% pull(gene)-> signature_subtypeA

candidate_signatures = list('gastritis_subtypeB'=signature_subtypeB,
                            'gastritis_subtypeA'=signature_subtypeA)
candidate_signatures_overall = list('gastritis_signature'=c(signature_subtypeA,signature_subtypeB))
saveRDS(candidate_signatures,'results/rds/candidate_signatures.rds')
saveRDS(candidate_signatures_overall,'results/rds/candidate_signatures_overall.rds')

