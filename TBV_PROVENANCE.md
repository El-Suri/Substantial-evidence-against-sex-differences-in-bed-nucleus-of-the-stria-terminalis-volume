# Total brain volume provenance

The manuscript states that total brain volume (TBV) was calculated from FSL FAST tissue segmentation outputs. For each participant, grey- and white-matter partial-volume images were used to estimate tissue volumes, and TBV was defined as:

```text
TBV = grey matter tissue volume + white matter tissue volume
```

For a FAST partial-volume image, tissue volume is calculated as the sum of the voxelwise partial-volume fractions multiplied by voxel volume. This is equivalent to multiplying the mean partial-volume value by the corresponding image volume, as described in the FSL FAST tissue-volume quantification documentation.

In FAST's standard three-tissue output naming:

- `*_pve_1.nii.gz` is the grey-matter partial-volume image.
- `*_pve_2.nii.gz` is the white-matter partial-volume image.

The code-only release does not include raw HCP anatomical images, FAST partial-volume images, HCP subject IDs, or individual-level TBV values. The manuscript analysis used precomputed TBV values derived from this FAST method.

To recompute TBV from local FAST outputs, place or point to a directory containing the relevant `*_pve_1.nii.gz` and `*_pve_2.nii.gz` files and run:

```bash
Rscript scripts/00_calculate_tbv_from_fast_pves.R --subjects=data/local_inputs/study_id_mapping.csv --fast_root=/path/to/local/FAST_outputs
```

This writes `results/tbv_from_fast_recomputed.csv`. The downstream open-derived table can then be rebuilt with:

```bash
Rscript scripts/00_prepare_open_derived_data.R --tbv=results/tbv_from_fast_recomputed.csv
```

The optional `study_id_mapping.csv` file should contain `study_id` and, for private local use only, `hcp_subject_id`. The script uses `hcp_subject_id` to locate local FAST files and writes only `study_id` to the output.
