Q2 sensitivity analysis package

Input:
- degs.csv supplied by the author

BH/FDR recalculation:
- Genes tested: 29388
- Nominal DEG (P < 0.05 and |log2FC| > 0.5): 2984
- Nominal, B higher: 1554
- Nominal, A higher: 1430
- FDR DEG (BH-adjusted P < 0.05 and |log2FC| > 0.5): 2953
- FDR, B higher: 1532
- FDR, A higher: 1421

Core genes:
- AGXT2L1
- C11orf86
- IGF1R
- FAR2
- MBIP
- CLIC6
- ZCWPW1

The R script performs:
1. GO-BP ORA
2. KEGG ORA
3. MSigDB Hallmark ORA
4. Nominal vs FDR pathway comparison
5. Concordant / nominal-only / FDR-only classification
6. Supplementary figure generation

Important:
- Ranked-list GSEA should remain the primary pathway analysis.
- ORA is used here as the requested sensitivity/complementary analysis.
- Do not claim nominal-only pathways as statistically significant.
