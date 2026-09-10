############################################################
# Q5_subtype_stability_sensitivity_FINAL_v3.R
#
# Complete, robust Q5 stability/sensitivity analysis.
#
# INPUT:
#   expr_combat_intergration_df.rds
#   sample_info_intergration.rds
#
# RUN:
#   source("Q5_subtype_stability_sensitivity.R")
#
# IMPORTANT:
#   Dataset-stratified ConsensusClusterPlus is run with plot="none".
#   This avoids the "clusterTrackingPlot ... end fraction clustered"
#   plotting error that can occur with small/specific cohort sizes.
#   Clustering itself is unchanged: hc + Pearson, 2000 reps,
#   pItem=0.8, pFeature=1, seed=1234.
############################################################

options(stringsAsFactors = FALSE)

############################################################
# 0. Packages
############################################################

cran_pkgs <- c(
  "data.table",
  "dplyr",
  "tidyr",
  "cluster",
  "mclust",
  "openxlsx"
)

bioc_pkgs <- "ConsensusClusterPlus"

missing_cran <- cran_pkgs[
  !vapply(
    cran_pkgs,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_cran) > 0) {
  stop(
    paste0(
      "Missing CRAN package(s): ",
      paste(missing_cran, collapse = ", "),
      ". Install with install.packages()."
    )
  )
}

missing_bioc <- bioc_pkgs[
  !vapply(
    bioc_pkgs,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_bioc) > 0) {
  stop(
    "Missing Bioconductor package: ConsensusClusterPlus"
  )
}

library(data.table)
library(dplyr)
library(tidyr)
library(cluster)
library(mclust)
library(openxlsx)
library(ConsensusClusterPlus)

############################################################
# 1. Parameters
############################################################

DISCOVERY_GSES <- c(
  "GSE130823",
  "GSE60427"
)

TOP_N_GENES <- 2000

CC_MAX_K <- 5
CC_K <- 2
CC_REPS <- 2000
CC_PITEM <- 0.8
CC_PFEATURE <- 1
CC_ALG <- "hc"
CC_DISTANCE <- "pearson"
CC_SEED <- 1234

OUTDIR <- "Q5_results"

dir.create(
  OUTDIR,
  recursive = TRUE,
  showWarnings = FALSE
)

############################################################
# 2. Find input files
############################################################

find_input <- function(filename) {

  candidates <- c(
    filename,
    file.path("data", filename),
    file.path("datasets", filename),
    file.path("input", filename)
  )

  found <- candidates[
    file.exists(candidates)
  ]

  if (length(found) == 0) {
    stop(
      paste0(
        "Cannot find: ",
        filename,
        "\nCurrent working directory: ",
        getwd()
      )
    )
  }

  normalizePath(
    found[1]
  )
}

expr_file <- find_input(
  "expr_combat_intergration_df.rds"
)

info_file <- find_input(
  "sample_info_intergration.rds"
)

############################################################
# 3. Utilities
############################################################

safe_numeric_cluster <- function(x) {

  nm <- names(x)

  if (is.factor(x)) {
    x <- as.character(x)
  }

  if (is.character(x)) {
    x <- suppressWarnings(
      as.numeric(x)
    )
  }

  x <- as.numeric(x)

  if (!is.null(nm)) {
    names(x) <- nm
  }

  x
}


safe_character_cluster <- function(x) {

  nm <- names(x)

  if (is.factor(x)) {
    x <- as.character(x)
  }

  x <- as.character(x)

  if (!is.null(nm)) {
    names(x) <- nm
  }

  x
}


select_top_variable_genes <- function(
    expr,
    top_n = 2000
) {

  expr <- as.matrix(expr)

  storage.mode(expr) <- "numeric"

  gene_sd <- apply(
    expr,
    1,
    sd,
    na.rm = TRUE
  )

  gene_sd <- gene_sd[
    is.finite(gene_sd)
  ]

  if (length(gene_sd) == 0) {
    stop(
      "No genes with finite SD."
    )
  }

  top_n <- min(
    top_n,
    length(gene_sd)
  )

  names(
    sort(
      gene_sd,
      decreasing = TRUE
    )
  )[seq_len(top_n)]
}


zscore_by_gene <- function(expr) {

  expr <- as.matrix(expr)

  storage.mode(expr) <- "numeric"

  z <- t(
    scale(
      t(expr),
      center = TRUE,
      scale = TRUE
    )
  )

  keep <- apply(
    z,
    1,
    function(x) {
      all(is.finite(x))
    }
  )

  z[
    keep,
    ,
    drop = FALSE
  ]
}


make_pearson_distance <- function(expr) {

  cor_mat <- cor(
    expr,
    method = "pearson",
    use = "pairwise.complete.obs"
  )

  cor_mat[
    !is.finite(cor_mat)
  ] <- 0

  cor_mat[
    cor_mat > 1
  ] <- 1

  cor_mat[
    cor_mat < -1
  ] <- -1

  diag(
    cor_mat
  ) <- 1

  as.dist(
    1 - cor_mat
  )
}

############################################################
# 4. ConsensusClusterPlus wrapper
############################################################
#
# plot_mode:
#   "pdf"  = ordinary plotting
#   "none" = NO plotting; recommended for dataset-stratified
#            sensitivity analysis.
#
# The clustering parameters are otherwise identical.
############################################################

run_cc <- function(
    expr,
    output_dir,
    maxK = CC_MAX_K,
    reps = CC_REPS,
    seed = CC_SEED,
    plot_mode = "pdf"
) {

  expr <- as.matrix(expr)

  storage.mode(expr) <- "numeric"

  genes <- select_top_variable_genes(
    expr,
    TOP_N_GENES
  )

  mat <- expr[
    genes,
    ,
    drop = FALSE
  ]

  mat_z <- zscore_by_gene(
    mat
  )

  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  cc <- ConsensusClusterPlus(
    mat_z,
    maxK = maxK,
    reps = reps,
    pItem = CC_PITEM,
    pFeature = CC_PFEATURE,
    clusterAlg = CC_ALG,
    distance = CC_DISTANCE,
    seed = seed,
    plot = plot_mode,
    title = output_dir
  )

  list(
    cc = cc,
    mat_z = mat_z,
    genes = rownames(mat_z)
  )
}

############################################################
# 5. Consensus stability
############################################################

calculate_consensus_stability <- function(
    consensus_matrix,
    cluster_labels,
    sample_ids
) {

  consensus_matrix <- as.matrix(
    consensus_matrix
  )

  storage.mode(
    consensus_matrix
  ) <- "numeric"

  cluster_labels <-
    safe_character_cluster(
      cluster_labels
    )

  sample_ids <- as.character(
    sample_ids
  )

  if (
    nrow(consensus_matrix) !=
      length(sample_ids) ||
    ncol(consensus_matrix) !=
      length(sample_ids)
  ) {

    stop(
      "Consensus matrix dimensions do not match sample IDs."
    )
  }

  rownames(
    consensus_matrix
  ) <- sample_ids

  colnames(
    consensus_matrix
  ) <- sample_ids

  cluster_labels <-
    cluster_labels[
      sample_ids
    ]

  clusters <- sort(
    unique(
      cluster_labels
    )
  )

  if (
    length(clusters) != 2
  ) {
    stop(
      paste0(
        "Expected two clusters; found ",
        length(clusters)
      )
    )
  }

  within_values <- c()

  within_by_cluster <- list()

  for (cl in clusters) {

    idx <- which(
      cluster_labels == cl
    )

    if (length(idx) > 1) {

      sub <- consensus_matrix[
        idx,
        idx,
        drop = FALSE
      ]

      vals <- sub[
        upper.tri(sub)
      ]

      within_values <- c(
        within_values,
        vals
      )

      within_by_cluster[[cl]] <-
        vals
    }
  }

  idx1 <- which(
    cluster_labels == clusters[1]
  )

  idx2 <- which(
    cluster_labels == clusters[2]
  )

  between_values <-
    as.vector(
      consensus_matrix[
        idx1,
        idx2,
        drop = FALSE
      ]
    )

  summary <- data.frame(

    N =
      length(cluster_labels),

    Cluster1 =
      clusters[1],

    Cluster1_N =
      length(idx1),

    Cluster2 =
      clusters[2],

    Cluster2_N =
      length(idx2),

    Mean_within_consensus =
      mean(
        within_values,
        na.rm = TRUE
      ),

    Median_within_consensus =
      median(
        within_values,
        na.rm = TRUE
      ),

    Mean_between_consensus =
      mean(
        between_values,
        na.rm = TRUE
      ),

    Median_between_consensus =
      median(
        between_values,
        na.rm = TRUE
      )

  )

  list(
    summary = summary,
    within_values = within_values,
    between_values = between_values,
    matrix = consensus_matrix
  )
}

############################################################
# 6. Robust silhouette
############################################################

calculate_silhouette_robust <- function(
    expr,
    cluster_labels
) {

  expr <- as.matrix(expr)

  storage.mode(expr) <- "numeric"

  sample_ids <- colnames(
    expr
  )

  if (
    is.null(sample_ids)
  ) {
    stop(
      "Expression matrix has no sample IDs."
    )
  }

  cluster_labels <-
    safe_character_cluster(
      cluster_labels
    )

  common <- intersect(
    sample_ids,
    names(cluster_labels)
  )

  if (
    length(common) < 3
  ) {
    stop(
      "Too few common samples."
    )
  }

  expr <- expr[
    ,
    common,
    drop = FALSE
  ]

  cluster_labels <-
    cluster_labels[
      common
    ]

  keep <-
    !is.na(
      cluster_labels
    )

  expr <- expr[
    ,
    keep,
    drop = FALSE
  ]

  cluster_labels <-
    cluster_labels[
      keep
    ]

  sample_ids <-
    colnames(expr)

  cluster_levels <- sort(
    unique(
      cluster_labels
    )
  )

  if (
    length(cluster_levels) != 2
  ) {
    stop(
      "Silhouette requires exactly two clusters."
    )
  }

  cluster_numeric <- match(
    cluster_labels,
    cluster_levels
  )

  cluster_numeric <- as.integer(
    cluster_numeric
  )

  names(
    cluster_numeric
  ) <- NULL

  d <- make_pearson_distance(
    expr
  )

  sil <- cluster::silhouette(
    x = cluster_numeric,
    dist = d
  )

  sil_matrix <- as.matrix(
    sil
  )

  if (
    nrow(sil_matrix) !=
      length(sample_ids)
  ) {

    stop(
      "Silhouette rows do not match sample IDs."
    )
  }

  # NEVER use rownames(sil).
  sil_df <- data.frame(

    Sample =
      sample_ids,

    Cluster_numeric =
      as.integer(
        sil_matrix[
          ,
          "cluster"
        ]
      ),

    Neighbour_cluster =
      as.integer(
        sil_matrix[
          ,
          "neighbor"
        ]
      ),

    Silhouette =
      as.numeric(
        sil_matrix[
          ,
          "sil_width"
        ]
      ),

    stringsAsFactors = FALSE
  )

  sil_df$Cluster <-
    ifelse(
      sil_df$Cluster_numeric == 1,
      "subtypeA",
      "subtypeB"
    )

  sil_df$Neighbour_cluster_name <-
    ifelse(
      sil_df$Neighbour_cluster == 1,
      "subtypeA",
      "subtypeB"
    )

  summary_df <-
    sil_df %>%
    group_by(
      Cluster
    ) %>%
    summarise(

      N = n(),

      Mean_silhouette =
        mean(
          Silhouette,
          na.rm = TRUE
        ),

      Median_silhouette =
        median(
          Silhouette,
          na.rm = TRUE
        ),

      SD_silhouette =
        sd(
          Silhouette,
          na.rm = TRUE
        ),

      .groups = "drop"

    )

  list(
    object = sil,
    sample = sil_df,
    summary = summary_df,
    mean =
      mean(
        sil_df$Silhouette,
        na.rm = TRUE
      ),
    median =
      median(
        sil_df$Silhouette,
        na.rm = TRUE
      )
  )
}

############################################################
# 7. Load data
############################################################

cat(
  "\nLoading data...\n"
)

expr_all <-
  readRDS(
    expr_file
  )

sample_info_all <-
  readRDS(
    info_file
  )

expr_all <-
  as.matrix(
    expr_all
  )

storage.mode(
  expr_all
) <- "numeric"

sample_info_all <-
  as.data.frame(
    sample_info_all
  )

cat(
  "Expression:",
  nrow(expr_all),
  "genes x",
  ncol(expr_all),
  "samples\n"
)

############################################################
# 8. Identify sample ID
############################################################

id_candidates <- c(
  "geo_accession",
  "GSM",
  "gsm",
  "sample",
  "sample_id",
  "Sample",
  "V1"
)

id_found <- id_candidates[
  id_candidates %in%
    colnames(
      sample_info_all
    )
]

if (
  length(id_found) == 0
) {
  stop(
    "No sample ID column found."
  )
}

id_col <- id_found[1]

sample_info_all$.sample_id <-
  as.character(
    sample_info_all[[id_col]]
  )

if (
  !"GSE" %in%
  colnames(
    sample_info_all
  )
) {
  stop(
    "Metadata lacks GSE column."
  )
}

############################################################
# 9. Select discovery samples
############################################################

info_discovery <-
  sample_info_all %>%
  filter(
    GSE %in%
      DISCOVERY_GSES
  )

if (
  "source_name_ch1" %in%
  colnames(
    info_discovery
  )
) {

  select_gsm <- c(

    info_discovery %>%
      filter(
        GSE == "GSE130823",
        source_name_ch1 ==
          "Gastritis"
      ) %>%
      pull(
        .sample_id
      ),

    info_discovery %>%
      filter(
        GSE == "GSE60427",
        !grepl(
          "normal",
          source_name_ch1,
          ignore.case = TRUE
        )
      ) %>%
      pull(
        .sample_id
      )

  )

} else {

  select_gsm <-
    info_discovery$.sample_id

}

select_gsm <- unique(
  select_gsm[
    select_gsm %in%
      colnames(
        expr_all
      )
  ]
)

if (
  length(select_gsm) == 0
) {
  stop(
    "No discovery samples found."
  )
}

info_discovery <-
  info_discovery %>%
  filter(
    .sample_id %in%
      select_gsm
  )

info_discovery <-
  info_discovery[
    match(
      select_gsm,
      info_discovery$.sample_id
    ),
    ,
    drop = FALSE
  ]

expr_discovery <-
  expr_all[
    ,
    info_discovery$.sample_id,
    drop = FALSE
  ]

############################################################
# 10. Integrated consensus clustering
############################################################

cat(
  "\n============================================\n",
  "Integrated discovery clustering\n",
  "============================================\n"
)

integrated <-
  run_cc(
    expr_discovery,
    file.path(
      OUTDIR,
      "integrated_consensus"
    ),
    maxK = CC_MAX_K,
    reps = CC_REPS,
    seed = CC_SEED,
    plot_mode = "pdf"
  )

############################################################
# 11. Extract K=2
############################################################

consensus_matrix <-
  integrated$cc[[CC_K]]$consensusMatrix

cluster_numeric <-
  safe_numeric_cluster(
    integrated$cc[[CC_K]]$consensusClass
  )

if (
  is.null(
    names(cluster_numeric)
  )
) {
  names(cluster_numeric) <-
    colnames(
      integrated$mat_z
    )
}

sample_ids <-
  names(cluster_numeric)

if (
  nrow(consensus_matrix) !=
    length(sample_ids) ||
  ncol(consensus_matrix) !=
    length(sample_ids)
) {

  stop(
    "Consensus matrix dimensions do not match cluster assignments."
  )
}

rownames(
  consensus_matrix
) <- sample_ids

colnames(
  consensus_matrix
) <- sample_ids

############################################################
# 12. Subtype labels
############################################################

cluster_levels <- sort(
  unique(
    cluster_numeric
  )
)

if (
  length(cluster_levels) != 2
) {
  stop(
    "Integrated clustering did not produce exactly two clusters."
  )
}

subtype_labels <-
  ifelse(
    cluster_numeric ==
      cluster_levels[1],
    "subtypeA",
    "subtypeB"
  )

names(
  subtype_labels
) <- sample_ids

info_discovery$Molecular_subtype <-
  subtype_labels[
    info_discovery$.sample_id
  ]

############################################################
# 13. Dataset/subtype association
############################################################

dataset_table <-
  table(
    info_discovery$GSE,
    info_discovery$Molecular_subtype
  )

dataset_proportion <-
  prop.table(
    dataset_table,
    margin = 1
  )

dataset_chisq <-
  suppressWarnings(
    chisq.test(
      dataset_table
    )
  )

dataset_association <-
  data.frame(

    Test =
      "Chi-square test",

    Chi_square =
      unname(
        dataset_chisq$statistic
      ),

    df =
      unname(
        dataset_chisq$parameter
      ),

    P_value =
      dataset_chisq$p.value

  )

############################################################
# 14. Sex
############################################################

sex_result <- data.frame()

sex_candidates <- c(
  "sex",
  "Sex",
  "gender",
  "Gender"
)

sex_found <- sex_candidates[
  sex_candidates %in%
    colnames(
      info_discovery
    )
]

if (
  length(sex_found) > 0
) {

  sex_var <- sex_found[1]

  tmp <-
    info_discovery %>%
    filter(
      !is.na(
        .data[[sex_var]]
      ),
      !is.na(
        Molecular_subtype
      )
    )

  tab <-
    table(
      tmp[[sex_var]],
      tmp$Molecular_subtype
    )

  if (
    nrow(tab) >= 2 &&
    ncol(tab) == 2
  ) {

    chi <-
      suppressWarnings(
        chisq.test(tab)
      )

    if (
      any(
        chi$expected < 5
      )
    ) {

      test <-
        fisher.test(tab)

      method <-
        "Fisher's exact test"

    } else {

      test <- chi

      method <-
        "Chi-square test"

    }

    sex_result <-
      data.frame(

        Variable =
          sex_var,

        Test =
          method,

        N =
          nrow(tmp),

        P_value =
          test$p.value

      )
  }
}

############################################################
# 15. Age
############################################################

age_result <- data.frame()

age_candidates <- colnames(
  info_discovery
)[
  grepl(
    "age",
    colnames(
      info_discovery
    ),
    ignore.case = TRUE
  )
]

age_candidates <-
  setdiff(
    age_candidates,
    ".sample_id"
  )

if (
  length(age_candidates) > 0
) {

  age_var <-
    age_candidates[1]

  age_x <-
    suppressWarnings(
      as.numeric(
        info_discovery[[age_var]]
      )
    )

  age_y <-
    info_discovery$Molecular_subtype

  keep <-
    is.finite(age_x) &
    !is.na(age_y)

  if (
    sum(keep) >= 3 &&
    length(
      unique(
        age_y[keep]
      )
    ) == 2
  ) {

    wt <-
      wilcox.test(
        age_x[keep] ~
          age_y[keep],
        exact = FALSE
      )

    age_result <-
      data.frame(

        Variable =
          age_var,

        Test =
          "Wilcoxon rank-sum test",

        N =
          sum(keep),

        Median_subtypeA =
          median(
            age_x[
              keep &
                age_y ==
                "subtypeA"
            ],
            na.rm = TRUE
          ),

        IQR_subtypeA =
          IQR(
            age_x[
              keep &
                age_y ==
                "subtypeA"
            ],
            na.rm = TRUE
          ),

        Median_subtypeB =
          median(
            age_x[
              keep &
                age_y ==
                "subtypeB"
            ],
            na.rm = TRUE
          ),

        IQR_subtypeB =
          IQR(
            age_x[
              keep &
                age_y ==
                "subtypeB"
            ],
            na.rm = TRUE
          ),

        P_value =
          wt$p.value

      )
  }
}

############################################################
# 16. Pathological variables
############################################################

path_keywords <- c(
  "inflamm",
  "atrophy",
  "metaplas",
  "intestinal",
  "h.?pylori",
  "helicobacter",
  "pylori"
)

path_candidates <- colnames(
  info_discovery
)[
  vapply(
    colnames(
      info_discovery
    ),
    function(x) {

      grepl(
        paste(
          path_keywords,
          collapse = "|"
        ),
        x,
        ignore.case = TRUE
      )

    },
    logical(1)
  )
]

path_candidates <-
  setdiff(
    path_candidates,
    c(
      ".sample_id",
      "GSE",
      "Molecular_subtype"
    )
  )

path_results <- list()

for (
  v in path_candidates
) {

  tmp <-
    info_discovery %>%
    filter(
      !is.na(
        .data[[v]]
      ),
      !is.na(
        Molecular_subtype
      )
    )

  if (
    nrow(tmp) < 3
  ) next

  tab <-
    table(
      tmp[[v]],
      tmp$Molecular_subtype
    )

  if (
    nrow(tab) < 2 ||
    ncol(tab) < 2
  ) next

  chi <-
    suppressWarnings(
      chisq.test(tab)
    )

  if (
    any(
      chi$expected < 5
    )
  ) {

    test <-
      fisher.test(tab)

    method <-
      "Fisher's exact test"

  } else {

    test <- chi

    method <-
      "Chi-square test"

  }

  path_results[[v]] <-
    data.frame(

      Variable =
        v,

      Test =
        method,

      N =
        nrow(tmp),

      P_value =
        test$p.value

    )
}

path_results_df <-
  if (
    length(path_results) > 0
  ) {
    bind_rows(
      path_results
    )
  } else {
    data.frame()
  }

############################################################
# 17. Consensus stability
############################################################

stability <-
  calculate_consensus_stability(
    consensus_matrix,
    subtype_labels,
    sample_ids
  )

stability_summary <-
  stability$summary

############################################################
# 18. Silhouette
############################################################

sil_result <-
  calculate_silhouette_robust(
    integrated$mat_z,
    cluster_numeric
  )

stability_summary$Mean_silhouette <-
  sil_result$mean

stability_summary$Median_silhouette <-
  sil_result$median

############################################################
# 19. Dataset-stratified clustering
############################################################
#
# IMPORTANT:
#   plot="none" prevents ConsensusClusterPlus plotting code
#   from failing for a smaller cohort. It does NOT change
#   clustering or consensus calculations.
############################################################

############################################################
# 19. Dataset-stratified clustering
############################################################
#
# IMPORTANT:
# We do NOT call ConsensusClusterPlus for the dataset-stratified
# sensitivity analysis.
#
# Reason:
# Some versions of ConsensusClusterPlus invoke internal plotting
# functions even when plot="none", producing:
#
#   end fraction
#   clustered
#   Error in 1:nrow(m): argument of length zero
#
# The sensitivity analysis below reproduces the essential
# consensus-clustering logic directly:
#
#   1. top variable genes
#   2. gene-wise Z-score
#   3. repeated subsampling of samples (pItem = 0.8)
#   4. hierarchical clustering using Pearson distance
#   5. accumulate co-clustering / co-occurrence counts
#   6. divide co-clustering by co-occurrence
#   7. obtain the final K=2 partition by hierarchical clustering
#
# This avoids the plotting bug while preserving the statistical
# resampling logic needed for the sensitivity analysis.
############################################################

run_custom_consensus <- function(
    expr,
    reps = 2000,
    pItem = 0.8,
    pFeature = 1,
    K = 2,
    seed = 1234
) {

  expr <- as.matrix(expr)

  storage.mode(expr) <- "numeric"

  sample_ids <- colnames(expr)

  if (
    is.null(sample_ids)
  ) {
    stop(
      "Expression matrix must have sample IDs."
    )
  }

  n_samples <- ncol(expr)

  if (
    n_samples < 4
  ) {
    stop(
      "At least 4 samples are required."
    )
  }

  ##########################################################
  # Select top variable genes ONCE within the dataset
  ##########################################################

  genes <- select_top_variable_genes(
    expr,
    TOP_N_GENES
  )

  mat <- expr[
    genes,
    ,
    drop = FALSE
  ]

  mat_z <- zscore_by_gene(
    mat
  )

  n_genes <- nrow(mat_z)

  ##########################################################
  # pFeature = 1 means all selected genes are used in every
  # resampling iteration.
  ##########################################################

  n_item <- max(
    2,
    floor(
      n_samples *
        pItem
    )
  )

  co_cluster <- matrix(
    0,
    n_samples,
    n_samples,
    dimnames = list(
      sample_ids,
      sample_ids
    )
  )

  co_occurrence <- matrix(
    0,
    n_samples,
    n_samples,
    dimnames = list(
      sample_ids,
      sample_ids
    )
  )

  ##########################################################
  # Repeated consensus resampling
  ##########################################################

  set.seed(
    seed
  )

  for (
    r in seq_len(reps)
  ) {

    idx <- sort(
      sample(
        seq_len(n_samples),
        size = n_item,
        replace = FALSE
      )
    )

    sub <- mat_z[
      ,
      idx,
      drop = FALSE
    ]

    d <- make_pearson_distance(
      sub
    )

    hc <- hclust(
      d,
      method = "complete"
    )

    cl <- cutree(
      hc,
      k = K
    )

    ########################################################
    # Update co-occurrence and co-clustering matrices
    ########################################################

    present <- matrix(
      0,
      n_samples,
      n_samples
    )

    present[
      idx,
      idx
    ] <- 1

    co_occurrence <-
      co_occurrence +
      present

    same_cluster <- outer(
      cl,
      cl,
      FUN = "=="
    )

    co_cluster[
      idx,
      idx
    ] <-
      co_cluster[
        idx,
        idx
      ] +
      same_cluster

  }

  ##########################################################
  # Consensus matrix
  ##########################################################

  consensus_matrix <- matrix(
    0,
    n_samples,
    n_samples,
    dimnames = list(
      sample_ids,
      sample_ids
    )
  )

  valid <- co_occurrence > 0

  consensus_matrix[
    valid
  ] <-
    co_cluster[
      valid
    ] /
    co_occurrence[
      valid
    ]

  diag(
    consensus_matrix
  ) <- 1

  ##########################################################
  # Final clustering of consensus matrix
  #
  # Distance = 1 - consensus.
  ##########################################################

  consensus_dist <- as.dist(
    1 -
      consensus_matrix
  )

  hc_final <- hclust(
    consensus_dist,
    method = "complete"
  )

  final_cluster <- cutree(
    hc_final,
    k = K
  )

  names(
    final_cluster
  ) <-
    sample_ids

  ##########################################################
  # Return object
  ##########################################################

  list(

    consensusMatrix =
      consensus_matrix,

    consensusClass =
      final_cluster,

    mat_z =
      mat_z,

    genes =
      rownames(mat_z),

    reps =
      reps,

    pItem =
      pItem,

    pFeature =
      pFeature,

    clusterAlg =
      "hierarchical clustering (complete linkage)",

    distance =
      "Pearson",

    seed =
      seed

  )
}


run_dataset_clustering <- function(
    dataset_name
) {

  cat(
    "\nDataset-stratified clustering:",
    dataset_name,
    "\n"
  )

  samples <-
    info_discovery$.sample_id[
      info_discovery$GSE ==
        dataset_name
    ]

  expr <-
    expr_discovery[
      ,
      samples,
      drop = FALSE
    ]

  cat(
    "  samples:",
    length(samples),
    "\n"
  )

  if (
    length(samples) < 4
  ) {

    stop(
      paste0(
        dataset_name,
        " has fewer than 4 samples."
      )
    )
  }

  cat(
    "  Running custom consensus clustering:",
    CC_REPS,
    "repetitions\n"
  )

  result <-
    run_custom_consensus(

      expr =
        expr,

      reps =
        CC_REPS,

      pItem =
        CC_PITEM,

      pFeature =
        CC_PFEATURE,

      K =
        2,

      seed =
        CC_SEED

    )

  cat(
    "  Finished:",
    dataset_name,
    "\n"
  )

  result
}

res_130823 <-
  run_dataset_clustering(
    "GSE130823"
  )

res_60427 <-
  run_dataset_clustering(
    "GSE60427"
  )


############################################################
# 20. Dataset-stratified ARI
############################################################

calculate_dataset_ari <- function(
    result,
    dataset_name
) {

  target_cluster <-
    safe_numeric_cluster(
      result$consensusClass
    )

  if (
    is.null(
      names(target_cluster)
    )
  ) {
    names(target_cluster) <-
      colnames(
        result$mat_z
      )
  }

  levels_target <- sort(
    unique(
      target_cluster
    )
  )

  target_subtype <-
    ifelse(
      target_cluster ==
        levels_target[1],
      "subtypeA",
      "subtypeB"
    )

  names(
    target_subtype
  ) <- names(
    target_cluster
  )

  original <-
    subtype_labels[
      names(target_subtype)
    ]

  keep <-
    !is.na(original) &
    !is.na(target_subtype)

  original <-
    original[keep]

  target_subtype <-
    target_subtype[keep]

  tab <-
    table(
      target_subtype,
      original
    )

  if (
    nrow(tab) == 2 &&
    ncol(tab) == 2
  ) {

    same <-
      tab[1,1] +
      tab[2,2]

    flip <-
      tab[1,2] +
      tab[2,1]

    if (
      flip > same
    ) {

      target_subtype <-
        ifelse(
          target_subtype ==
            rownames(tab)[1],
          "subtypeB",
          "subtypeA"
        )
    }
  }

  ari <-
    adjustedRandIndex(
      original,
      target_subtype
    )

  concordance <-
    mean(
      original ==
        target_subtype
    )

  data.frame(

    Dataset =
      dataset_name,

    N =
      length(original),

    ARI =
      ari,

    Concordance =
      concordance

  )
}

dataset_ari <-
  bind_rows(

    calculate_dataset_ari(
      res_130823,
      "GSE130823"
    ),

    calculate_dataset_ari(
      res_60427,
      "GSE60427"
    )

  )

############################################################
# 21. Leave-one-dataset-out analysis
############################################################

derive_centroids <- function(
    expr,
    labels
) {

  labels <-
    safe_character_cluster(
      labels
    )

  clusters <-
    sort(
      unique(
        labels
      )
    )

  centroids <-
    sapply(
      clusters,
      function(cl) {

        rowMeans(
          expr[
            ,
            labels == cl,
            drop = FALSE
          ],
          na.rm = TRUE
        )

      }
    )

  centroids <-
    as.matrix(
      centroids
    )

  colnames(
    centroids
  ) <-
    clusters

  centroids
}


assign_to_centroid <- function(
    expr,
    centroids
) {

  expr <-
    as.matrix(
      expr
    )

  pred <-
    sapply(
      colnames(expr),
      function(s) {

        x <-
          expr[
            ,
            s
          ]

        cors <-
          apply(
            centroids,
            2,
            function(y) {

              suppressWarnings(
                cor(
                  x,
                  y,
                  method = "pearson",
                  use = "pairwise.complete.obs"
                )
              )

            }
          )

        if (
          all(
            !is.finite(cors)
          )
        ) {
          return(NA_character_)
        }

        names(
          which.max(
            cors
          )
        )

      }
    )

  names(pred) <-
    colnames(expr)

  pred
}


run_LOO <- function(
    reference_gse,
    heldout_gse
) {

  cat(
    "\nLOO:",
    reference_gse,
    " -> ",
    heldout_gse,
    "\n"
  )

  ref_samples <-
    info_discovery$.sample_id[
      info_discovery$GSE ==
        reference_gse
    ]

  held_samples <-
    info_discovery$.sample_id[
      info_discovery$GSE ==
        heldout_gse
    ]

  ref_expr <-
    expr_discovery[
      ,
      ref_samples,
      drop = FALSE
    ]

  held_expr <-
    expr_discovery[
      ,
      held_samples,
      drop = FALSE
    ]

  genes <-
    select_top_variable_genes(
      ref_expr,
      TOP_N_GENES
    )

  ref_mat <-
    ref_expr[
      genes,
      ,
      drop = FALSE
    ]

  ref_mean <-
    rowMeans(
      ref_mat,
      na.rm = TRUE
    )

  ref_sd <-
    apply(
      ref_mat,
      1,
      sd,
      na.rm = TRUE
    )

  ref_sd[
    !is.finite(ref_sd) |
      ref_sd == 0
  ] <- 1

  ref_z <-
    sweep(
      ref_mat,
      1,
      ref_mean,
      "-"
    )

  ref_z <-
    sweep(
      ref_z,
      1,
      ref_sd,
      "/"
    )

  keep_gene <-
    apply(
      ref_z,
      1,
      function(x) {
        all(
          is.finite(x)
        )
      }
    )

  ref_z <-
    ref_z[
      keep_gene,
      ,
      drop = FALSE
    ]

  genes <-
    rownames(ref_z)

  outdir <-
    file.path(
      OUTDIR,
      "leave_one_dataset_out",
      paste0(
        reference_gse,
        "_reference"
      )
    )

  dir.create(
    outdir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  # Use the same custom consensus-resampling implementation as the
  # dataset-stratified sensitivity analysis. This deliberately avoids
  # ConsensusClusterPlus plotting internals, which can fail in some
  # package versions even when plot="none".
  loo_ref <-
    run_custom_consensus(

      expr =
        ref_z,

      reps =
        CC_REPS,

      pItem =
        CC_PITEM,

      pFeature =
        CC_PFEATURE,

      K =
        2,

      seed =
        CC_SEED

    )

  ref_cluster <-
    safe_numeric_cluster(
      loo_ref$consensusClass
    )

  if (
    is.null(
      names(ref_cluster)
    )
  ) {
    names(ref_cluster) <-
      colnames(ref_z)
  }

  ref_original <-
    subtype_labels[
      names(ref_cluster)
    ]

  ref_levels <-
    sort(
      unique(
        ref_cluster
      )
    )

  ref_new <-
    ifelse(
      ref_cluster ==
        ref_levels[1],
      "subtypeA",
      "subtypeB"
    )

  names(
    ref_new
  ) <- names(
    ref_cluster
  )

  tab <-
    table(
      ref_new,
      ref_original
    )

  if (
    nrow(tab) == 2 &&
    ncol(tab) == 2
  ) {

    same <-
      tab[1,1] +
      tab[2,2]

    flip <-
      tab[1,2] +
      tab[2,1]

    if (
      flip > same
    ) {

      ref_new <-
        ifelse(
          ref_new ==
            rownames(tab)[1],
          "subtypeB",
          "subtypeA"
        )
    }
  }

  centroids <-
    derive_centroids(
      ref_z,
      ref_new
    )

  held_mat <-
    held_expr[
      genes,
      ,
      drop = FALSE
    ]

  held_z <-
    sweep(
      held_mat,
      1,
      ref_mean[
        genes
      ],
      "-"
    )

  held_z <-
    sweep(
      held_z,
      1,
      ref_sd[
        genes
      ],
      "/"
    )

  predicted <-
    assign_to_centroid(
      held_z,
      centroids
    )

  original <-
    subtype_labels[
      names(predicted)
    ]

  common <-
    intersect(
      names(original),
      names(predicted)
    )

  keep <-
    !is.na(
      original[
        common
      ]
    ) &
    !is.na(
      predicted[
        common
      ]
    )

  common <-
    common[
      keep
    ]

  ari <-
    adjustedRandIndex(
      original[
        common
      ],
      predicted[
        common
      ]
    )

  concordance <-
    mean(
      original[
        common
      ] ==
        predicted[
          common
        ]
    )

  data.frame(

    Reference_dataset =
      reference_gse,

    Heldout_dataset =
      heldout_gse,

    N_heldout =
      length(common),

    ARI =
      ari,

    Concordance =
      concordance

  )
}

loo_results <-
  bind_rows(

    run_LOO(
      "GSE130823",
      "GSE60427"
    ),

    run_LOO(
      "GSE60427",
      "GSE130823"
    )

  )

############################################################
# 22. ConsensusClusterPlus parameters
############################################################

cc_parameters <-
  data.frame(

    Parameter = c(

      "Discovery cohorts",
      "Number of discovery samples",
      "Gene selection",
      "Number of genes",
      "Transformation",
      "Maximum K",
      "Selected K",
      "Repetitions",
      "pItem",
      "pFeature",
      "Clustering algorithm",
      "Distance metric",
      "Random seed",
      "Dataset-stratified / LOO sensitivity method"

    ),

    Value = c(

      paste(
        DISCOVERY_GSES,
        collapse = "; "
      ),

      ncol(
        expr_discovery
      ),

      "Top genes ranked by SD",
      nrow(integrated$mat_z),
      "Gene-wise Z-score",
      CC_MAX_K,
      CC_K,
      CC_REPS,
      CC_PITEM,
      CC_PFEATURE,
      CC_ALG,
      CC_DISTANCE,
      CC_SEED,
      "Custom consensus resampling; hierarchical complete-linkage clustering of consensus matrix"

    )

  )

############################################################
# 23. Final assignment
############################################################

assignment_table <-
  info_discovery %>%
  select(
    .sample_id,
    GSE,
    Molecular_subtype,
    everything()
  )

############################################################
# 24. Write CSVs
############################################################

write.csv(
  as.data.frame(dataset_table),
  file.path(
    OUTDIR,
    "dataset_subtype_counts.csv"
  ),
  row.names = FALSE
)

write.csv(
  as.data.frame(dataset_proportion),
  file.path(
    OUTDIR,
    "dataset_subtype_proportions.csv"
  ),
  row.names = FALSE
)

write.csv(
  dataset_association,
  file.path(
    OUTDIR,
    "dataset_subtype_association.csv"
  ),
  row.names = FALSE
)

write.csv(
  sex_result,
  file.path(
    OUTDIR,
    "sex_subtype_association.csv"
  ),
  row.names = FALSE
)

write.csv(
  age_result,
  file.path(
    OUTDIR,
    "age_subtype_association.csv"
  ),
  row.names = FALSE
)

write.csv(
  path_results_df,
  file.path(
    OUTDIR,
    "pathological_subtype_associations.csv"
  ),
  row.names = FALSE
)

write.csv(
  stability_summary,
  file.path(
    OUTDIR,
    "cluster_stability_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  sil_result$sample,
  file.path(
    OUTDIR,
    "silhouette_by_sample.csv"
  ),
  row.names = FALSE
)

write.csv(
  sil_result$summary,
  file.path(
    OUTDIR,
    "silhouette_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  dataset_ari,
  file.path(
    OUTDIR,
    "dataset_stratified_ARI.csv"
  ),
  row.names = FALSE
)

write.csv(
  loo_results,
  file.path(
    OUTDIR,
    "leave_one_dataset_out_results.csv"
  ),
  row.names = FALSE
)

write.csv(
  cc_parameters,
  file.path(
    OUTDIR,
    "ConsensusClusterPlus_parameters.csv"
  ),
  row.names = FALSE
)

write.csv(
  assignment_table,
  file.path(
    OUTDIR,
    "final_subtype_assignments.csv"
  ),
  row.names = FALSE
)

############################################################
# 25. Excel
############################################################

wb <- createWorkbook()

tables <- list(

  Dataset_subtype_counts =
    as.data.frame(dataset_table),

  Dataset_proportions =
    as.data.frame(dataset_proportion),

  Dataset_association =
    dataset_association,

  Sex_association =
    sex_result,

  Age_association =
    age_result,

  Pathological_association =
    path_results_df,

  Cluster_stability =
    stability_summary,

  Silhouette_by_sample =
    sil_result$sample,

  Silhouette_summary =
    sil_result$summary,

  Dataset_stratified_ARI =
    dataset_ari,

  Leave_one_dataset_out =
    loo_results,

  ConsensusClusterPlus_parameters =
    cc_parameters,

  Subtype_assignments =
    assignment_table

)

for (
  nm in names(tables)
) {

  addWorksheet(
    wb,
    nm
  )

  writeData(
    wb,
    nm,
    tables[[nm]]
  )

}

saveWorkbook(
  wb,
  file.path(
    OUTDIR,
    "Q5_stability_sensitivity_results.xlsx"
  ),
  overwrite = TRUE
)

############################################################
# 26. Save objects
############################################################

saveRDS(

  list(

    integrated =
      integrated,

    consensus_matrix =
      consensus_matrix,

    cluster_numeric =
      cluster_numeric,

    subtype_labels =
      subtype_labels,

    stability =
      stability,

    silhouette =
      sil_result,

    dataset_stratified_ARI =
      dataset_ari,

    leave_one_dataset_out =
      loo_results,

    metadata =
      info_discovery

  ),

  file.path(
    OUTDIR,
    "Q5_stability_sensitivity_objects.rds"
  )

)

############################################################
# 27. Silhouette PDF
############################################################

pdf(
  file.path(
    OUTDIR,
    "silhouette_plot.pdf"
  ),
  width = 8,
  height = 6
)

plot(
  sil_result$object,
  border = NA,
  main =
    "Silhouette plot for integrated K = 2 clustering"
)

dev.off()

############################################################
# 28. Final console report
############################################################

cat(
  "\n\n========================================================\n"
)

cat(
  "Q5 ANALYSIS COMPLETED SUCCESSFULLY\n"
)

cat(
  "========================================================\n\n"
)

cat(
  "Discovery samples:\n"
)

print(
  table(
    info_discovery$GSE
  )
)

cat(
  "\nSubtype counts:\n"
)

print(
  table(
    info_discovery$Molecular_subtype
  )
)

cat(
  "\nDataset x subtype:\n"
)

print(
  dataset_table
)

cat(
  "\nDataset association:\n"
)

print(
  dataset_association
)

cat(
  "\nCluster stability:\n"
)

print(
  stability_summary
)

cat(
  "\nDataset-stratified ARI:\n"
)

print(
  dataset_ari
)

cat(
  "\nLeave-one-dataset-out:\n"
)

print(
  loo_results
)

cat(
  "\nOutput directory:\n"
)

cat(
  normalizePath(
    OUTDIR
  ),
  "\n"
)

cat(
  "\nExcel workbook:\n"
)

cat(
  normalizePath(
    file.path(
      OUTDIR,
      "Q5_stability_sensitivity_results.xlsx"
    )
  ),
  "\n"
)

cat(
  "\n========================================================\n"
)

cat(
  "All reported results were calculated from the supplied RDS files.\n"
)

cat(
  "No numerical results were simulated.\n"
)

cat(
  "========================================================\n"
)
