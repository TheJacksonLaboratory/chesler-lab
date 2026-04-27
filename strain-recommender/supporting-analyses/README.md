# Guide to `supporting-analyses` files

Below is a guide to each file in the `supporting-analyses` folder.

| File | Description |
|------|-------------|
| `retrospective_function_check_exclusions.R` | Removes rows with missing disease names or study IDs, non-exact strain matches, duplicate entries, and self-comparisons |
| `retrospective_function_get_sensitivity.R` | Computes per-disease sensitivity as the proportion of correct predictions |
| `retrospective_function_process_data.R` | Carries out sensitivity calculations using `check_exclusions()` and `get_sensitivity()` |
| `retrospective_manuscript_mESC_strain_ranks_across_diseases_final.Rmd` | Analysis of strain rankings with respect to both extremeness (vulnerable or resistant) and vulnerability |
| `retrospective_manuscript_other_tissues_analyses.Rmd` | Analysis of other tissue types; finds mean sensitivity and 95% CI for each tissue |
| `retrospective_manuscript_processing_other_tissue_types.Rmd` | Preprocesses other tissue result data for `retrospective_manuscript_other_tissues_analyses.Rmd` |
| `retrospective_manuscript_results_disease-sensitivities.Rmd` | Reads in the final filtered result set and calculates disease sensitivity distribution across all diseases and the rare disease subset |
| `retrospective_manuscript_sensitivity-analysis_fpr-threshold.Rmd` | Assesses mean sensitivity when setting FPR cutoffs in 0.05 increments |
| `retrospective_manuscript_sensitivity-analysis_impact-of-fixed-factors-on-performance.Rmd` | Runs ANOVAs testing whether gene set size, species, or tier impacts correctness in the global validation |
| `retrospective_manuscript_sensitivity-analysis_non-redundant_comparisons.Rmd` | Assesses sensitivity on non-redundant comparisons introduced by the ontology mapping |
| `CTD/retrospective_manuscript_CTD_analysis.Rmd` | Multivariate linear regression analyses to assess the effect of mouse strain on EKG and echocardiogram phenotypes relevant to Marfan syndrome and cardiovascular health, using scaled variables, ANOVA-based F-statistics, and BH-corrected p-values to evaluate both overall and variable-specific strain effects |
| `CTD/retrospective_manuscript_CTD_marfan_gs_SR.Rmd` | Code was used to run Marfan Syndrome (MFS) gene sets through Strain Recommender to predict strain vulnerability |
| `AMD/Camd_recommender_manuscript_2026-02-23.Rmd` | Cross-Species Gene Expression Analysis for AMD Model Strain Recommendation |
