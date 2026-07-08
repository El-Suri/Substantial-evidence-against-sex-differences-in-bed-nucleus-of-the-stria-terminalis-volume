# Human BNST sex-difference analysis code release

This is a code-only staging folder for the manuscript "Substantial evidence against sex differences in bed nucleus of the stria terminalis volume in the human brain".

No participant-level HCP data, HCP-derived segmentations, HCP subject IDs, derived individual-level volume tables, or manuscript result tables are included in this pre-acceptance release. Data files will be added only after HCP-approved study-specific IDs are available.

## Included

- `scripts/`: R scripts for BNST mask volume extraction, FAST-derived TBV calculation, open-derived table assembly, prior construction, bilateral model fitting, prior-width sensitivity checks, unilateral sensitivity checks, and release checks.
- `DATA_REQUIREMENTS.md`: the files required to rerun the analysis once HCP-approved study IDs and local HCP Restricted Access data are available.
- `DATA_USE.md`: data-use notes for the future data release.
- `TBV_PROVENANCE.md`: description of the FSL FAST partial-volume method used to calculate total brain volume.
- `data/`: placeholder folders only.

## Not Included

- HCP raw anatomical images.
- FAST partial-volume images.
- BNST segmentation masks.
- HCP subject IDs.
- HCP Restricted Access variables, including exact age, family IDs, twin/zygosity information, handedness, race/ethnicity, height, weight, BMI, SSAGA, endocrine, or menstrual variables.
- Individual-level derived BNST/TBV tables.
- Manuscript result tables.
- Saved model objects (`.rds`, `.RDS`, `.RData`).

## Expected Future Data Layout

After publication acceptance and HCP approval for mapped study IDs, the repository can be completed with:

```text
data/open_derived/bnst_open_derived.csv
data/segmentations/<study_id>.nii.gz
data/prior_sources/
```

The public join key will be `study_id`, not the original HCP subject ID.

Exact model reruns also require a local restricted covariate file that will require access to the HCP Young Adults restricted dataset:

```text
study_id,Age_in_Yrs,Family_ID,ZygositySR,monozygotic
```

## Running

Install the R packages used by the scripts:

```r
install.packages(c("brms", "posterior", "polspline", "metafor", "RNifti", "bayesplot"))
```

