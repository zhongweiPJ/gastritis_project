
#figure3 
# ============================================================
# Script 05 - Core seven-gene signature (formerly figure3ML.R)
# Run from the repository root.
# Input: data/step04_ML.RData (ML workspace)
# Derives the intersection of LASSO, SVM, and RF feature sets (7 genes).
# ============================================================

m=load(file = 'data/step04_ML.RData')
####################model_glmnet
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
####################SVM
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
##################Random Forest
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

#3个模型共有的genes
intersect(
  lasso_genes,
  intersect(
    svm_top,rf_top
  )
) #7个genes
