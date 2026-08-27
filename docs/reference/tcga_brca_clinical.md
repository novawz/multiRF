# TCGA BRCA Clinical Data

Clinical metadata matched to the bundled TCGA BRCA omics cohort.

## Format

A data frame with 674 rows (one per sample) and 9 columns:

- sampleID:

  TCGA sample barcode matching the row names of the `tcga_brca` omics
  blocks.

- sample_type:

  Sample type: `"Primary Tumor"`, `"Solid Tissue Normal"`, or
  `"Metastatic"`.

- BRCA_Subtype_PAM50:

  PAM50 intrinsic subtype call (`"Basal"`, `"Her2"`, `"LumA"`, `"LumB"`,
  `"Normal"`, or `NA`).

- OS:

  Overall survival event indicator (1 = death, 0 = censored).

- OS.time:

  Overall survival time in days.

- PFI:

  Progression-free interval event indicator (1 = event, 0 = censored).

- PFI.time:

  Progression-free interval time in days.

- age_at_initial_pathologic_diagnosis:

  Age in years at initial pathologic diagnosis.

- pathologic_stage:

  AJCC pathologic tumor stage (e.g. `"Stage IIA"`).

## Source

<https://www.cancer.gov/tcga>
