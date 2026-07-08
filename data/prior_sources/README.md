# Prior-source data placeholder

No prior-source data are included in this code-only release.

`scripts/02_create_prior_table.R` documents how the literature-informed prior was derived. To rerun it locally, add the following files when source terms allow:

- `Allen_Gorski.csv`: participant-level table with `BNST.Volume`, `Sex`, `Brain.Weight`, and `Age`.
- `Chung_etal.csv`: participant-level table with `BNST.Volume`, `Sex`, `Brain.Weight`, and `Age`.
- `Guma_2024_adjusted_effect.csv`: study-level adjusted effect table with `n_male`, `n_female`, `effect_size`, and `effect_se`, if the participant-level HCP RDS cannot be used.

Alternatively, rerun the Guma re-analysis from a private local copy of the participant-level HCP-derived file:

```bash
Rscript scripts/02_create_prior_table.R --guma_participant_file=/path/to/df_HCP_volumes_clean.RDS
```

Do not upload Guma participant-level HCP data unless HCP explicitly permits it.
