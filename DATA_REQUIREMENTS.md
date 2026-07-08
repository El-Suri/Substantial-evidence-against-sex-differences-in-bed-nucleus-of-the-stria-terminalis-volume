# Data requirements

This code-only release deliberately excludes data. To rerun the manuscript analyses, create the following files locally.

## Public data after HCP-mapped IDs are approved

These files should use HCP-approved study-specific IDs, not original HCP subject IDs.

### `data/open_derived/bnst_open_derived.csv`

Required columns:

- `study_id`: HCP-approved mapped study ID.
- `Gender`: sex/gender variable used in the manuscript analysis (`F`, `M`).
- `Age_group`: broad HCP open-access age group, if released.
- `TBV_mm3`: total brain volume in cubic millimetres.
- `BNST_left_mm3`: left BNST mask volume in cubic millimetres.
- `BNST_right_mm3`: right BNST mask volume in cubic millimetres.
- `BNST_bilateral_mm3`: left plus right BNST volume in cubic millimetres.
- `segmentation_file`: relative path to the segmentation mask.

### `data/segmentations/<study_id>.nii.gz`

Final BNST manual segmentation masks, one file per participant, named with the mapped `study_id`.

Mask labels:

- `0`: background.
- `1`: right BNST.
- `2`: left BNST.

### `data/prior_sources/`

Literature-derived prior inputs used by `scripts/02_create_prior_table.R`. The script contains the summary-statistic derivations directly. Participant-level Allen/Gorski and Chung tables, and either a private Guma participant-level HCP file or a cached Guma adjusted-effect CSV, are required for a full rerun. The Guma et al. participant-level HCP file should not be redistributed.

## Local restricted data for exact reruns

Create this file locally and keep it out of the public repository:

```text
data/local_inputs/restricted_covariates.csv
```

Required columns:

- `study_id`: same mapped ID used in the public open-derived table.
- `Age_in_Yrs`: exact age in years. HCP Restricted Access variable.
- `Family_ID`: HCP family identifier. HCP Restricted Access variable.
- `ZygositySR`: twin/zygosity variable, if using the script to derive `monozygotic`.
- `monozygotic`: optional derived indicator where `1` marks MZ twins and `0` marks everyone else.

The scripts can use either `ZygositySR` or `monozygotic`.

## Local HCP-ID mapping

If recomputing TBV or rebuilding tables before the HCP-approved public mapping is available, keep any HCP-ID mapping file outside version control. For local FAST recomputation, `scripts/00_calculate_tbv_from_fast_pves.R` can use:

```text
data/local_inputs/study_id_mapping.csv
```

with columns:

- `study_id`: mapped study ID, or a private local ID.
- `hcp_subject_id`: original HCP subject ID used only for local file lookup.

Do not upload this mapping file unless HCP explicitly approves the mapped-ID release mechanism.
