############################################################
# xCell statistical analysis for SCI manuscript revision
#
# Input:
#   xcell_results.csv
#
# Required columns:
#   - id
#   - Group
#   - xCell-derived features
#
# Group:
#   - SubtypeA
#   - SubtypeB
#
# Output:
#   - Median and IQR for each subtype
#   - Wilcoxon rank-sum test P value
#   - Benjamini-Hochberg adjusted P value
#   - Rank-biserial correlation effect size
#   - 95% bootstrap confidence interval
#   - Direction of difference
#
############################################################


############################
# 1. Load packages
############################

library(tidyverse)
library(rstatix)
library(openxlsx)


############################
# 2. Input / output paths
############################

input_file <- "xcell_results.csv"

output_dir <- "xCell_statistics"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}


############################
# 3. Read xCell results
############################

xcell <- read.csv(
  input_file,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


############################
# 4. Inspect data
############################

cat("========================================\n")
cat("xCell statistical analysis\n")
cat("========================================\n\n")

cat("Number of samples:", nrow(xcell), "\n")

cat("Columns:\n")
print(colnames(xcell))


############################
# 5. Check Group column
############################

if (!"Group" %in% colnames(xcell)) {
  stop("Column 'Group' was not found in the xCell input file.")
}

cat("\nGroup distribution:\n")
print(table(xcell$Group))


############################
# 6. Standardize group names
############################

xcell$Group <- trimws(as.character(xcell$Group))

# Convert possible naming variations
xcell$Group[xcell$Group %in% c(
  "Subtype A",
  "Subtype_A",
  "subtypeA",
  "A"
)] <- "SubtypeA"

xcell$Group[xcell$Group %in% c(
  "Subtype B",
  "Subtype_B",
  "subtypeB",
  "B"
)] <- "SubtypeB"


############################
# 7. Keep only Subtype A/B
############################

xcell <- xcell %>%
  filter(Group %in% c("SubtypeA", "SubtypeB"))

xcell$Group <- factor(
  xcell$Group,
  levels = c("SubtypeA", "SubtypeB")
)

cat("\nFinal sample numbers:\n")
print(table(xcell$Group))


############################
# 8. Identify xCell features
############################

# Metadata columns that should NOT be tested
metadata_columns <- c(
  "id",
  "ID",
  "sample",
  "Sample",
  "Group",
  "group",
  "Unnamed: 0"
)

cell_features <- setdiff(
  colnames(xcell),
  metadata_columns
)

cat("\nNumber of xCell-derived features:",
    length(cell_features), "\n")

cat("\nFeatures tested:\n")
print(cell_features)


############################
# 9. Convert to numeric
############################

for (feature in cell_features) {
  
  xcell[[feature]] <- as.numeric(
    xcell[[feature]]
  )
  
}


############################################################
# 10. Function for rank-biserial correlation
############################################################

rank_biserial <- function(x, y) {
  
  x <- x[!is.na(x)]
  y <- y[!is.na(y)]
  
  n_x <- length(x)
  n_y <- length(y)
  
  if (n_x == 0 || n_y == 0) {
    return(NA_real_)
  }
  
  # Mann-Whitney U
  test <- wilcox.test(
    x,
    y,
    exact = FALSE
  )
  
  U <- as.numeric(test$statistic)
  
  # Rank-biserial correlation
  #
  # Positive:
  #   Subtype A tends to have higher scores
  #
  # Negative:
  #   Subtype B tends to have higher scores
  
  rbc <- (
    2 * U / (n_x * n_y)
  ) - 1
  
  return(rbc)
}


############################################################
# 11. Bootstrap 95% CI for rank-biserial correlation
############################################################

bootstrap_rbc <- function(
    x,
    y,
    n_boot = 2000,
    conf_level = 0.95
) {
  
  x <- x[!is.na(x)]
  y <- y[!is.na(y)]
  
  n_x <- length(x)
  n_y <- length(y)
  
  if (n_x == 0 || n_y == 0) {
    return(
      c(
        lower = NA_real_,
        upper = NA_real_
      )
    )
  }
  
  boot_values <- numeric(n_boot)
  
  for (i in seq_len(n_boot)) {
    
    x_boot <- sample(
      x,
      size = n_x,
      replace = TRUE
    )
    
    y_boot <- sample(
      y,
      size = n_y,
      replace = TRUE
    )
    
    boot_values[i] <- rank_biserial(
      x_boot,
      y_boot
    )
  }
  
  alpha <- 1 - conf_level
  
  lower <- quantile(
    boot_values,
    probs = alpha / 2,
    na.rm = TRUE
  )
  
  upper <- quantile(
    boot_values,
    probs = 1 - alpha / 2,
    na.rm = TRUE
  )
  
  return(
    c(
      lower = as.numeric(lower),
      upper = as.numeric(upper)
    )
  )
}


############################################################
# 12. Statistical analysis
############################################################

results <- list()

set.seed(20260829)

for (feature in cell_features) {
  
  cat(
    "Analyzing:",
    feature,
    "\n"
  )
  
  ########################################
  # Extract data
  ########################################
  
  dat <- xcell %>%
    select(
      Group,
      all_of(feature)
    ) %>%
    rename(
      Score = all_of(feature)
    ) %>%
    filter(
      !is.na(Score)
    )
  
  
  ########################################
  # Separate groups
  ########################################
  
  group_A <- dat %>%
    filter(Group == "SubtypeA") %>%
    pull(Score)
  
  group_B <- dat %>%
    filter(Group == "SubtypeB") %>%
    pull(Score)
  
  
  ########################################
  # Sample size
  ########################################
  
  n_A <- length(group_A)
  n_B <- length(group_B)
  
  
  ########################################
  # Median and IQR
  ########################################
  
  median_A <- median(
    group_A,
    na.rm = TRUE
  )
  
  median_B <- median(
    group_B,
    na.rm = TRUE
  )
  
  IQR_A <- IQR(
    group_A,
    na.rm = TRUE
  )
  
  IQR_B <- IQR(
    group_B,
    na.rm = TRUE
  )
  
  
  ########################################
  # Wilcoxon rank-sum test
  ########################################
  
  wilcox_result <- wilcox.test(
    group_A,
    group_B,
    alternative = "two.sided",
    exact = FALSE
  )
  
  raw_p <- wilcox_result$p.value
  
  
  ########################################
  # Effect size
  ########################################
  
  effect_size <- rank_biserial(
    group_A,
    group_B
  )
  
  
  ########################################
  # 95% CI
  ########################################
  
  ci <- bootstrap_rbc(
    group_A,
    group_B,
    n_boot = 2000,
    conf_level = 0.95
  )
  
  
  ########################################
  # Direction
  ########################################
  
  if (effect_size > 0) {
    
    direction <- "Subtype A higher"
    
  } else if (effect_size < 0) {
    
    direction <- "Subtype B higher"
    
  } else {
    
    direction <- "No difference"
    
  }
  
  
  ########################################
  # Save result
  ########################################
  
  results[[feature]] <- data.frame(
    
    `Cell population` = feature,
    
    `N_Subtype_A` = n_A,
    
    `N_Subtype_B` = n_B,
    
    `Median_Subtype_A` = median_A,
    
    `IQR_Subtype_A` = IQR_A,
    
    `Median_Subtype_B` = median_B,
    
    `IQR_Subtype_B` = IQR_B,
    
    `Effect_size_rank_biserial_r` =
      effect_size,
    
    `CI95_low` =
      ci["lower"],
    
    `CI95_high` =
      ci["upper"],
    
    `Raw_P` =
      raw_p,
    
    check.names = FALSE
    
  )
}


############################################################
# 13. Combine results
############################################################

results_df <- bind_rows(
  results
)


############################################################
# 14. Benjamini-Hochberg correction
############################################################

results_df$BH_adjusted_P <- p.adjust(
  results_df$Raw_P,
  method = "BH"
)


############################################################
# 15. Statistical significance
############################################################

results_df$FDR_significant <-
  results_df$BH_adjusted_P < 0.05


############################################################
# 16. Add direction
############################################################

results_df$Direction <-
  ifelse(
    results_df$Effect_size_rank_biserial_r > 0,
    "Subtype A higher",
    ifelse(
      results_df$Effect_size_rank_biserial_r < 0,
      "Subtype B higher",
      "No difference"
    )
  )


############################################################
# 17. Sort by adjusted P value
############################################################

results_df <- results_df %>%
  arrange(
    BH_adjusted_P
  )


############################################################
# 18. Create compact supplementary table
############################################################

supplementary_table <- results_df %>%
  select(
    `Cell population`,
    N_Subtype_A,
    N_Subtype_B,
    
    Median_Subtype_A,
    IQR_Subtype_A,
    
    Median_Subtype_B,
    IQR_Subtype_B,
    
    Effect_size_rank_biserial_r,
    CI95_low,
    CI95_high,
    
    Raw_P,
    BH_adjusted_P,
    
    Direction
  )


############################################################
# 19. Create summary table
############################################################

summary_table <- data.frame(
  
  Metric = c(
    
    "Number of samples",
    "Subtype A samples",
    "Subtype B samples",
    "Number of xCell-derived features tested",
    "FDR-significant features",
    "FDR-significant features higher in subtype A",
    "FDR-significant features higher in subtype B"
    
  ),
  
  Value = c(
    
    nrow(xcell),
    
    sum(
      xcell$Group == "SubtypeA"
    ),
    
    sum(
      xcell$Group == "SubtypeB"
    ),
    
    length(cell_features),
    
    sum(
      results_df$BH_adjusted_P < 0.05
    ),
    
    sum(
      results_df$BH_adjusted_P < 0.05 &
        results_df$Effect_size_rank_biserial_r > 0
    ),
    
    sum(
      results_df$BH_adjusted_P < 0.05 &
        results_df$Effect_size_rank_biserial_r < 0
    )
    
  )
  
)


############################################################
# 20. Statistical methods table
############################################################

methods_table <- data.frame(
  
  Item = c(
    
    "Statistical test",
    "Multiple-testing correction",
    "Effect size",
    "Confidence interval",
    "Bootstrap resampling",
    "Direction convention",
    "Significance threshold"
    
  ),
  
  Specification = c(
    
    "Wilcoxon rank-sum test",
    
    "Benjamini-Hochberg correction",
    
    "Rank-biserial correlation",
    
    "95% confidence interval",
    
    "2,000 bootstrap resamples",
    
    "Positive = higher in subtype A; negative = higher in subtype B",
    
    "BH-adjusted P < 0.05"
    
  )
  
)


############################################################
# 21. Write CSV files
############################################################

write.csv(
  
  results_df,
  
  file = file.path(
    output_dir,
    "xCell_statistics_full.csv"
  ),
  
  row.names = FALSE
  
)


write.csv(
  
  supplementary_table,
  
  file = file.path(
    output_dir,
    "Supplementary_xCell_statistics.csv"
  ),
  
  row.names = FALSE
  
)


write.csv(
  
  summary_table,
  
  file = file.path(
    output_dir,
    "xCell_statistics_summary.csv"
  ),
  
  row.names = FALSE
  
)


############################################################
# 22. Write Excel workbook
############################################################

excel_file <- file.path(
  
  output_dir,
  
  "xCell_statistics_Wilcoxon_BH_effect_CI.xlsx"
  
)


write.xlsx(
  
  list(
    
    xCell_statistics =
      supplementary_table,
    
    Summary =
      summary_table,
    
    Methods =
      methods_table
    
  ),
  
  file = excel_file,
  
  overwrite = TRUE
  
)


############################################################
# 23. Print summary
############################################################

cat("\n\n")
cat("========================================\n")
cat("Analysis completed\n")
cat("========================================\n\n")

print(summary_table)

cat("\n\nTop statistically significant features:\n\n")

print(
  results_df %>%
    select(
      `Cell population`,
      Effect_size_rank_biserial_r,
      CI95_low,
      CI95_high,
      Raw_P,
      BH_adjusted_P,
      Direction
    ) %>%
    head(15)
)


cat("\n\nOutput files:\n")

cat(
  file.path(
    output_dir,
    "xCell_statistics_full.csv"
  ),
  "\n"
)

cat(
  file.path(
    output_dir,
    "Supplementary_xCell_statistics.csv"
  ),
  "\n"
)

cat(
  file.path(
    output_dir,
    "xCell_statistics_summary.csv"
  ),
  "\n"
)

cat(
  excel_file,
  "\n"
)

############################################################
# END
############################################################