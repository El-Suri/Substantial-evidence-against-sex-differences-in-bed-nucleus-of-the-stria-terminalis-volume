# Data-use and sharing notes

This is a code-only pre-acceptance release. It intentionally contains no HCP participant-level data, HCP subject IDs, HCP-derived segmentation masks, derived individual-level volume tables, or manuscript result tables.

The exact manuscript model uses HCP Restricted Access variables, including exact age, family membership, and twin/zygosity information. These variables must not be redistributed publicly.

When the manuscript is accepted, request the appropriate HCP-approved mechanism for sharing study-specific mapped IDs. The future public data release should use those mapped IDs as `study_id`, not original HCP subject IDs.

Do not publicly upload:

- filled restricted covariate files;
- any file linking mapped study IDs to original HCP subject IDs, unless HCP explicitly authorises that route;
- raw HCP anatomical images;
- FAST partial-volume images;
- broad HCP phenotype files containing unrelated variables.

Public uses of future HCP-derived files should acknowledge HCP according to the relevant HCP data-use terms:

- https://www.humanconnectome.org/study/hcp-young-adult/document/wu-minn-hcp-consortium-open-access-data-use-terms
- https://www.humanconnectome.org/study/hcp-young-adult/document/wu-minn-hcp-consortium-restricted-data-use-terms
