# Integrative Transcriptomic Analysis of Gastritis Molecular Subtypes

This repository contains the complete R-based computational workflow used for the manuscript:

> **Integrative Transcriptomic Analysis Identifies Immune-Inflammatory Molecular Subtypes Associated with Progression-Related Phenotypes in Gastritis**

The workflow covers data preprocessing, phenotype annotation, discovery-cohort-specific batch correction, molecular subtype discovery, differential expression analysis, functional enrichment, immune/stromal deconvolution, machine-learning-based feature prioritization, independent external evaluation, statistical testing, and generation of figures and supplementary tables.

The repository is intended to provide a transparent and reproducible record of the analyses reported in the manuscript.

---

## Repository structure

```text
.
├── README.md                     # this document (full workflow + reproducibility)
│
├── scripts/                      # main analysis pipeline (run in numeric order)
│   ├── 01_preprocess_batch_consensus.R          # GEO loading, ComBat, consensus clustering
│   ├── 02_differential_expression_gsva_xcell.R  # limma DEG, GO/KEGG, GSVA, xCell
│   ├── 03_enrichment_ORA_GSEA.R                 # ORA + GSEA (also sourced by 02)
│   ├── 04_machine_learning_feature_selection.R  # LASSO / SVM / RF, candidate signatures
│   ├── 05_core_signature.R                      # 7-gene intersection
│   └── 06_external_validation_nomogram.R        # ssGSEA, nomogram, ROC
│
├── data/                          # raw GEO rds, MSigDB/GO gene sets, processed objects
│
├── results/
│   ├── csv/                       # DEG + enrichment + scoring tables
│   ├── figs/                      # PNG / PDF figures
│   └── rds/                       # enrichment objects + duplicated input objects
│
├── consensus_out_intergration/    # consensus clustering outputs (clusters.csv, consensus.pdf)
│
└── revision_Q2Q3Q5/               # re-analyses for reviewer Q2 / Q3 / Q5 (self-contained)
    ├── Q2/                        # DEG sensitivity ORA
    ├── Q3/                        # ssGSEA / xCell Wilcoxon statistics
    └── Q5/                        # subtype stability / sensitivity
```

---

## 1. Study design and analytical overview

Four publicly available GEO datasets were analyzed.

### Discovery cohort

Two datasets were used exclusively for molecular subtype discovery:

- **GSE130823**: 47 gastritis samples
- **GSE60427**: 24 gastritis samples

Total discovery cohort:

**71 gastritis samples**

These 71 samples were used for:

1. discovery-cohort preprocessing
2. discovery-cohort batch correction
3. highly variable gene selection
4. consensus clustering
5. molecular subtype assignment
6. differential expression analysis
7. functional characterization
8. machine-learning feature prioritization

### External evaluation cohort

Two independent datasets were reserved exclusively for external evaluation:

- **GSE55696**: 77 samples, including 19 chronic gastritis samples
- **GSE60662**: 16 samples, including 12 gastritis-related samples

The two external datasets therefore contain:

**93 samples in total, including 31 gastritis samples.**

Unlike subtype discovery, external evaluation uses samples with available pathological phenotype annotations to evaluate whether the discovery-derived molecular signatures are associated with different pathological states.

The external datasets are **not used for**:

- subtype definition
- discovery-cohort highly variable gene selection
- discovery-cohort feature selection
- machine-learning model training
- estimation of discovery-cohort batch-correction parameters

The external analyses evaluate associations with **progression-related pathological phenotypes** and should not be interpreted as prospective prediction of future disease progression or individual clinical risk.

---

## 2. GEO datasets

| GEO accession | Platform | Total samples | Gastritis samples | Analytical role |
|---|---|---:|---:|---|
| GSE130823 | GPL17077 | 94 | 47 | Discovery cohort |
| GSE60427 | GPL17077 | 32 | 24 | Discovery cohort |
| GSE55696 | GPL6480 | 77 | 19 | External evaluation |
| GSE60662 | GPL13497 | 16 | 12 | External evaluation |

Direct GEO accession pages:

- GSE130823: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130823
- GSE60427: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60427
- GSE55696: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE55696
- GSE60662: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60662

---

## 3. Recommended execution order

Run the scripts in `scripts/` in numeric order from the repository root:

```text
01_preprocess_batch_consensus.R
    ├── load raw GEO matrices, probe → gene-symbol mapping
    ├── discovery-cohort ComBat batch correction
    └── consensus clustering (K = 2) → consensus_out_intergration/clusters.csv
            │
            v
02_differential_expression_gsva_xcell.R
    ├── limma differential expression (subtypeB vs subtypeA)
    ├── volcano / heatmap figures
    ├── GO / KEGG bar plots
    ├── GSVA and ssGSEA scoring
    └── xCell deconvolution
            │
            v
03_enrichment_ORA_GSEA.R
    ├── ORA (GO BP/CC/MF + KEGG), up / down / overall
    └── GSEA (GO + KEGG + MSigDB c2/c5/hallmark)
            │
            v
04_machine_learning_feature_selection.R
    ├── LASSO / SVM / RF feature ranking
    └── candidate subtype signatures
            │
            v
05_core_signature.R
    └── 7-gene intersection (LASSO ∩ SVM ∩ RF)
            │
            v
06_external_validation_nomogram.R
    └── external ssGSEA / nomogram / ROC on GSE55696 + GSE60662
```

`scripts/03_enrichment_ORA_GSEA.R` is also sourced by `scripts/02_...` to produce
the enrichment tables used for the bar plots.

The re-analyses for reviewer questions Q2 / Q3 / Q5 are kept separately under
`revision_Q2Q3Q5/` and are self-contained (each subfolder includes its own input
copies and scripts).

---

# 4. Data preprocessing and phenotype annotation

Raw or processed GEO expression matrices are harmonized to a common gene identifier space.

Sample-level phenotype information is obtained from GEO annotations and curated into the analysis metadata.

For the discovery cohort:

- GSE130823: retain gastritis samples only (n = 47)
- GSE60427: retain gastritis samples only (n = 24)

Total:

- discovery cohort = 71 samples

Non-gastritis/control samples in these two datasets are not used for molecular subtype discovery.

The processed objects used by downstream analyses include:

```text
expr_combat_intergration_df.rds
sample_info_intergration.rds
```

---

# 5. Discovery-cohort batch correction

Batch correction is performed **only within the discovery cohort**.

The discovery cohort consists of:

```text
GSE130823 gastritis
+
GSE60427 gastritis
=
71 samples
```

The `ComBat` function from the `sva` package is used.

### Batch variable

The GEO dataset identifier (`GSE`) is used as the batch variable.

### ComBat settings

```text
method: ComBat
package: sva
parametric empirical Bayes: par.prior = TRUE
biological covariates: none
```

No biological covariates are included in the ComBat model.

GSE55696 and GSE60662 are processed independently and do **not** contribute to estimation of the batch-correction parameters for the discovery cohort.

The batch-corrected discovery matrix is subsequently used for subtype discovery and downstream discovery analyses.

---

# 6. Consensus clustering

Molecular subtype discovery is performed using `ConsensusClusterPlus`.

### Input

- 71 discovery-cohort gastritis samples
- top 2,000 most variable genes
- gene-wise z-score normalization before clustering

### ConsensusClusterPlus parameters

```text
maxK       = 5
reps       = 2000
pItem      = 0.8
pFeature   = 1
clusterAlg = "hc"
distance   = "pearson"
seed       = 1234
```

Hierarchical clustering (`clusterAlg = "hc"`) with Pearson distance is used.

The two-cluster solution (`K = 2`) is selected for downstream molecular subtype characterization.

The resulting cluster assignments are labeled:

```text
subtypeA
subtypeB
```

Consensus matrices and cluster assignments should be retained as intermediate outputs to permit reproduction of the stability analyses and figures.

---

# 7. Differential expression analysis

Differential expression is evaluated between subtype A and subtype B in the discovery cohort.

For every gene, both raw and Benjamini–Hochberg-adjusted P values are reported.

### Final DEG interpretation

The principal FDR-controlled DEG criterion is:

```text
adjusted P < 0.05
|log2 fold change| > 0.5
```

Genes meeting only the nominal P-value criterion are described as:

> nominally differentially expressed genes

They are not described as statistically significant FDR-controlled DEGs.

The complete differential-expression table should contain at least:

```text
gene
log2FC
raw P value
BH-adjusted P value
```

A sensitivity analysis compares functional enrichment derived from the nominal-P-value gene list with enrichment derived from the FDR-controlled DEG list.

---

# 8. Functional analyses

## 8.1 ORA

Over-representation analysis is performed on the predefined gene lists.

The FDR-controlled DEG list is used for the principal analysis.

Where appropriate, the nominal-P-value gene list is used as an exploratory sensitivity analysis.

Enrichment results are reported with raw and multiple-testing-adjusted P values.

---

## 8.2 Ranked-list GSEA

Ranked-list gene set enrichment analysis is used to evaluate pathway-level differences without depending exclusively on an arbitrary DEG cutoff.

The ranking statistic is generated from the differential-expression results used in the corresponding analysis.

The exact ranking metric and gene-set collection should be specified in the corresponding analysis script.

---

## 8.3 GSVA

GSVA is performed to obtain pathway activity scores for individual samples.

The GSVA calculation itself produces continuous pathway-level scores.

When GSVA scores are subsequently compared between biological groups, the group-comparison procedure is treated as a separate statistical analysis.

For simultaneous pathway comparisons, raw P values are adjusted using the Benjamini–Hochberg procedure.

---

## 8.4 ssGSEA

Single-sample GSEA (ssGSEA) is used for pathway/signature scoring in the analyses corresponding to Figure 3 and external pathological-group comparisons.

For pairwise comparisons of ssGSEA scores between pathological groups:

```text
statistical test = Wilcoxon rank-sum test
multiple-testing correction = Benjamini–Hochberg
```

Raw and adjusted P values are retained.

---

# 9. xCell analysis

xCell is used to estimate enrichment scores for immune and stromal cell populations.

For the subtype comparison shown in the main analysis, xCell scores are compared between subtype A and subtype B using a nonparametric Wilcoxon rank-sum test.

Student's t-test is not used as the primary test because distributional assumptions were not established for all xCell populations.

For each tested cell population, the statistical output should include:

- effect size
- confidence interval, where calculated
- raw P value
- BH-adjusted P value

The FDR-adjusted results are used for the principal interpretation of immune/stromal differences.

---

# 10. Machine-learning feature prioritization

Machine-learning methods are used as **exploratory feature-prioritization tools**, not as independent predictive models.

The discovery cohort is used for this analysis.

The three algorithms are:

1. LASSO logistic regression
2. random forest
3. radial-basis-function support vector machine

Repeated cross-validation is used for model training and internal resampling.

```text
10-fold cross-validation
5 repeats
random seed = 123
```

## 10.1 LASSO

LASSO is implemented using `caret`/`glmnet`.

Tuning grid:

```text
alpha = 1

lambda =
10^seq(-4, 1, length = 50)
```

The selected lambda is determined during model tuning.

Genes with non-zero coefficients at the selected lambda are retained as LASSO-selected features.

## 10.2 Random forest

Random forest is implemented using `caret`.

Tuning grid:

```text
mtry = c(5, 10, 30, 50, 100)
```

Number of trees:

```text
ntree = 500
```

Genes are ranked according to model-derived variable importance.

The top 50 genes are retained as the RF candidate feature set.

## 10.3 SVM

A radial-basis-function SVM is implemented using `caret`/`svmRadial`.

Tuning grid:

```text
sigma = c(0.001, 0.005, 0.01, 0.05, 0.1)

C = c(0.25, 0.5, 1, 2, 4, 8)
```

Genes are ranked according to model-derived feature importance.

The top 50 genes are retained as the SVM candidate feature set.

## 10.4 Final seven-gene selection rule

The final core signature is defined as the intersection of the three algorithm-derived feature sets:

```text
LASSO-selected genes
        ∩
RF top-50 genes
        ∩
SVM top-50 genes
        =
7 core genes
```

The seven genes are:

```text
AGXT2L1
C11orf86
IGF1R
FAR2
MBIP
CLIC6
ZCWPW1
```

ROC, AUC, calibration, and nomogram analyses are not interpreted as evidence of independent predictive validity.

---

# 11. External signature evaluation

GSE55696 and GSE60662 are reserved exclusively for external evaluation.

They are processed separately from the discovery cohort.

The seven-gene signature and scoring rule are fixed using the discovery analysis and subsequently applied unchanged to each external dataset.

External datasets do not participate in:

- discovery subtype definition
- discovery HVG selection
- discovery feature selection
- discovery model training
- discovery-cohort ComBat parameter estimation

The external analyses evaluate whether the discovery-derived molecular signatures show consistent associations with cross-sectional pathological phenotypes.

These analyses may include comparisons across:

- gastritis
- intestinal metaplasia
- neoplasia-related pathological groups
- early gastric cancer
- control or inflammation-related groups

where such annotations are available.

These analyses are **not** prospective longitudinal prediction analyses and do not establish individual future progression risk.

---

# 12. Statistical testing and multiple-testing correction

## 12.1 Wilcoxon tests

For group comparisons of xCell and ssGSEA scores:

```text
Wilcoxon rank-sum test
```

is used as the primary nonparametric method.

## 12.2 Benjamini–Hochberg correction

For analyses involving multiple simultaneously tested cell populations, pathways, or pathological-group comparisons:

```text
adjusted P value = Benjamini–Hochberg FDR
```

Raw P values and adjusted P values are retained.

The correction is performed within the predefined family of simultaneous comparisons corresponding to the relevant analysis.

## 12.3 Categorical clinical/pathological variables

Associations between subtype assignment and categorical variables are assessed using chi-square tests or Fisher's exact tests where appropriate.

## 12.4 Continuous variables

Continuous clinical/pathological variables are assessed using nonparametric tests when distributional assumptions for parametric methods are not established.

---

# 13. Subtype stability and sensitivity analyses

Subtype stability is assessed quantitatively rather than relying only on biological interpretation.

The analyses include:

1. consensus-matrix stability
2. silhouette analysis
3. subtype distribution across discovery datasets
4. dataset-stratified clustering
5. adjusted Rand index (ARI)
6. leave-one-dataset-out sensitivity analysis

The purpose of these analyses is to determine whether the identified subtype structure is robust to dataset composition and alternative discovery-cohort configurations.

The results should be generated directly from the saved clustering outputs and should not be manually edited.

---

# 14. Sample accounting

The final sample accounting is:

### Discovery

```text
GSE130823
94 total
47 gastritis → included

GSE60427
32 total
24 gastritis → included

Total discovery cohort
47 + 24 = 71 gastritis samples
```

### External evaluation

```text
GSE55696
77 total
19 chronic gastritis

GSE60662
16 total
12 gastritis-related

Total external datasets
77 + 16 = 93 samples

Total gastritis samples across external datasets
19 + 12 = 31
```

The external datasets are used to evaluate signatures across available pathological phenotypes rather than to redefine the discovery subtypes.

---

# 16. Reproducibility and computational environment

The primary computational environment used for the analysis was:

```text
R version: 4.5.0 (2025-04-11)
Platform: x86_64-w64-mingw32/x64
Operating system: Windows 11 x64
```

### 16.1 Package versions directly verified from the analysis session

The following versions were explicitly provided from the analysis `sessionInfo()`:

```text
caret                  7.0-1
sva                    0.58.0
GSVA                   2.4.9
ConsensusClusterPlus   1.74.0
ggplot2                4.0.1
dplyr                  1.1.4
tibble                 3.3.0
rlang                  1.2.0
scales                 1.4.0
RColorBrewer           1.1-3
R6                     2.6.1
farver                 2.1.2
gtable                 0.3.6
glue                   1.8.1
lifecycle              1.0.5
magrittr               2.0.4
tidyselect             1.2.1
vctrs                  0.7.2
```

### 16.2 Additional packages required by the workflow

The following packages are required by the corresponding analysis modules.
Where an exact version was not present in the supplied `sessionInfo()`, a
compatible version is specified as a **repository pin/recommended reproducibility
version**, rather than being represented as an observed historical version.

```text
limma                   3.64.3
clusterProfiler         4.16.0
enrichplot              1.28.4
fgsea                   1.34.0
glmnet                  4.1-10
randomForest             4.7-1.2
e1071                    1.7-16
xCell                    1.1.0
pheatmap                 1.0.12
data.table               1.17.8
tidyr                    1.3.1
stringr                  1.5.1
readr                    2.1.5
purrr                    1.1.0
forcats                  1.0.0
patchwork                1.3.1
Rtsne                    0.17
ggrepel                  0.9.6
ggpubr                   0.6.0
survival                 3.8-3
```

**Important:** the versions in Section 16.1 are directly supported by the
provided analysis session. The versions in Section 16.2 are reproducibility
pins for packages whose exact historical versions were not included in the
provided `sessionInfo()`. Before final manuscript submission, the repository
maintainer should run the complete workflow once in the pinned environment and
save the resulting `sessionInfo.txt` or `renv.lock`. The lock file should be
treated as the definitive computational-environment specification.

For Bioconductor packages, the repository should use the Bioconductor release
compatible with R 4.5.0 used by the final execution environment.

### 16.3 Environment capture

Run:

```r
writeLines(
  capture.output(sessionInfo()),
  "sessionInfo.txt"
)
```

For stronger reproducibility, use `renv`:

```r
install.packages("renv")
renv::init()
renv::snapshot()
```

This creates:

```text
renv.lock
```

which records the exact package versions and sources used in the final
execution environment.

# 17. Random seeds

The principal random seeds documented in the analysis are:

```text
ConsensusClusterPlus: 1234
Machine-learning analyses: 123
```

Any additional stochastic procedure should specify its seed directly in the corresponding script.

---

# 18. Output organization

The repository uses the following layout:

```text
.
├── README.md
├── scripts/                         # main pipeline (01–06)
├── data/                            # inputs + processed objects
│   ├── GSE*_RAW.rds / GSE*.rds      # raw + series GEO objects
│   ├── c2/c5/hallmark/go_df.rds     # MSigDB + GO gene sets
│   ├── pathway2gene*.rds / pathway2name.rds
│   ├── expr_combat_intergration*.rds
│   ├── sample_info_intergration.rds
│   ├── figure3_ML.RData / step04_ML.RData
│   └── dat.RData / kegg_hsa.RData
│
├── results/
│   ├── csv/                         # DEG + enrichment + scoring tables
│   ├── figs/                        # PNG / PDF figures
│   └── rds/                         # enrichment objects and input copies
│
├── consensus_out_intergration/      # clusters.csv, consensus.pdf, stability results
│
└── revision_Q2Q3Q5/                 # reviewer re-analyses (self-contained)
    ├── Q2/                          # DEG sensitivity ORA
    ├── Q3/                          # ssGSEA / xCell Wilcoxon statistics
    └── Q5/                          # subtype stability / sensitivity
```

The main intermediate objects written by the pipeline are:

```text
consensus_out_intergration/clusters.csv        # K = 2 subtype assignment
results/csv/deg_list.csv                       # differential expression table
results/csv/gsea_c2_c5_hall_kegg_results.csv   # ranked-list GSEA
results/csv/go_kegg_*.csv                      # ORA results
results/rds/candidate_signatures*.rds          # ML candidate signatures
```

---

# 19. Reproducibility principles

The following principles are used throughout the analysis:

1. External validation datasets are not used to define the discovery subtypes.
2. Discovery-cohort batch correction is performed independently from external evaluation.
3. The seven-gene signature is fixed before external evaluation.
4. Raw and adjusted P values are retained for multiple-testing-sensitive analyses.
5. FDR-controlled results are used for principal statistical interpretation where appropriate.
6. Machine learning is interpreted as exploratory feature prioritization rather than independent prediction.
7. Cross-sectional pathological-group comparisons are interpreted as associations with progression-related phenotypes rather than prospective progression prediction.
8. Numerical results should be regenerated from the analysis scripts and saved intermediate objects rather than manually entered.
9. The repository should preserve the analysis environment through `sessionInfo.txt` or `renv.lock`.

---

## Citation and data availability

The transcriptomic datasets analyzed in this study are publicly available from the Gene Expression Omnibus under accession numbers:

**GSE130823, GSE60427, GSE55696, and GSE60662.**

The complete computational workflow and analysis scripts are provided in this repository.

