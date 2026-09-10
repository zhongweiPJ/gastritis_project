############################################################
# Q3: ssGSEA pairwise group comparisons
# Input:
#   ssgsea_result_sampleid(2).csv
#
# Columns required:
#   GSE
#   Group
#   gastritis_subtypeA
#   gastritis_subtypeB
#
# Analysis:
#   Within each GSE and each ssGSEA signature:
#   - all pairwise comparisons among Groups
#   - two-sided Wilcoxon rank-sum test
#   - rank-biserial correlation effect size
#   - 95% bootstrap CI
#   - BH correction across the six pairwise comparisons
#
############################################################

library(tidyverse)
library(openxlsx)

input_file <- "ssgsea_result_sampleid(2).csv"
output_dir <- "Q3_ssGSEA_statistics"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

dat <- read.csv(
  input_file,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

signatures <- c(
  "gastritis_subtypeA",
  "gastritis_subtypeB"
)

############################################################
# Rank-biserial correlation
#
# Positive = Group 1 higher
# Negative = Group 2 higher
############################################################

rank_biserial <- function(x, y) {

  x <- x[!is.na(x)]
  y <- y[!is.na(y)]

  n1 <- length(x)
  n2 <- length(y)

  wt <- wilcox.test(
    x, y,
    alternative = "two.sided",
    exact = FALSE
  )

  U <- as.numeric(wt$statistic)

  2 * U / (n1 * n2) - 1
}

############################################################
# Bootstrap 95% CI
############################################################

bootstrap_rbc_ci <- function(
    x,
    y,
    B = 3000
) {

  x <- x[!is.na(x)]
  y <- y[!is.na(y)]

  n1 <- length(x)
  n2 <- length(y)

  boot <- numeric(B)

  for (i in seq_len(B)) {

    xb <- sample(
      x,
      size = n1,
      replace = TRUE
    )

    yb <- sample(
      y,
      size = n2,
      replace = TRUE
    )

    boot[i] <- rank_biserial(
      xb,
      yb
    )
  }

  quantile(
    boot,
    probs = c(0.025, 0.975),
    na.rm = TRUE
  )
}

############################################################
# All pairwise comparisons
############################################################

results <- list()

set.seed(12345)

for (gse in unique(dat$GSE)) {

  dat_gse <- dat %>%
    filter(GSE == gse)

  group_levels <- sort(
    unique(dat_gse$Group)
  )

  pair_list <- combn(
    group_levels,
    2,
    simplify = FALSE
  )

  for (signature in signatures) {

    family_results <- list()

    for (pair in pair_list) {

      group1 <- pair[1]
      group2 <- pair[2]

      x <- dat_gse %>%
        filter(Group == group1) %>%
        pull(all_of(signature)) %>%
        as.numeric() %>%
        na.omit()

      y <- dat_gse %>%
        filter(Group == group2) %>%
        pull(all_of(signature)) %>%
        as.numeric() %>%
        na.omit()

      ######################################################
      # Wilcoxon rank-sum test
      ######################################################

      wt <- wilcox.test(
        x,
        y,
        alternative = "two.sided",
        exact = FALSE
      )

      ######################################################
      # Effect size
      ######################################################

      effect <- rank_biserial(
        x,
        y
      )

      ######################################################
      # 95% CI
      ######################################################

      ci <- bootstrap_rbc_ci(
        x,
        y,
        B = 3000
      )

      family_results[[paste(
        group1,
        group2,
        sep = "_vs_"
      )]] <- data.frame(

        GSE = gse,

        Signature = signature,

        Group_1 = group1,

        Group_2 = group2,

        N_Group1 = length(x),

        N_Group2 = length(y),

        Median_Group1 = median(x),

        IQR_Group1 = IQR(x),

        Median_Group2 = median(y),

        IQR_Group2 = IQR(y),

        Effect_size_RBC = effect,

        CI95_low = as.numeric(ci[1]),

        CI95_high = as.numeric(ci[2]),

        Raw_P = wt$p.value,

        stringsAsFactors = FALSE

      )
    }

    ########################################################
    # BH correction within GSE × signature
    ########################################################

    family_results <- bind_rows(
      family_results
    )

    family_results$BH_FDR <- p.adjust(
      family_results$Raw_P,
      method = "BH"
    )

    family_results$FDR_significant <-
      family_results$BH_FDR < 0.05

    family_results$Direction <- ifelse(
      family_results$Effect_size_RBC > 0,
      paste(
        family_results$Group_1,
        "higher"
      ),
      paste(
        family_results$Group_2,
        "higher"
      )
    )

    results[[paste(
      gse,
      signature,
      sep = "_"
    )]] <- family_results
  }
}

############################################################
# Combine results
############################################################

results_df <- bind_rows(results) %>%
  arrange(
    GSE,
    Signature,
    BH_FDR
  )

############################################################
# Summary
############################################################

summary_df <- results_df %>%
  group_by(
    GSE,
    Signature
  ) %>%
  summarise(
    Comparisons = n(),
    FDR_significant =
      sum(FDR_significant),
    .groups = "drop"
  )

############################################################
# Save
############################################################

write.csv(
  results_df,
  file.path(
    output_dir,
    "ssGSEA_Figure3D_pairwise_Wilcoxon_BH.csv"
  ),
  row.names = FALSE
)

write.csv(
  summary_df,
  file.path(
    output_dir,
    "ssGSEA_statistics_summary.csv"
  ),
  row.names = FALSE
)

write.xlsx(
  list(
    Pairwise_statistics = results_df,
    Summary = summary_df
  ),
  file.path(
    output_dir,
    "ssGSEA_Figure3D_pairwise_statistics.xlsx"
  ),
  overwrite = TRUE
)

print(results_df)
print(summary_df)

############################################################
# END
############################################################
