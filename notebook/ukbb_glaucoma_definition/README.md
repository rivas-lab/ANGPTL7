## Glaucoma definition in UK Biobank

We performed this analysis using 3 scripts.

1. `glaucoma_def_extract_data.R`: this script extracts the relevant columns from the raw table data in `/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/9796/9797/download/ukb9797.tab` and save to `ukb9797.glaucoma.tsv` (this is not part of the GitHub repo because it has individual-level data).

2. `glaucoma_def_check.R`: this short script checks if the extracted data matches with our current phenotype definition (stored in the phe file).

3. `glaucoma_def_upset_plot.R`: this script generates UpSet plots.

We used the following data fields:

- self-reported
    - UKB Field: 20002
    - UKB coding: '1277'
- ICD-10
    - UKB Fields: 41202, 42104, 40001, 41078, 41143, 41142, 41104, 41079, 41105, 40002, 40006, 41204
    - ICD-10 codes: H40, H41, H400, H401, H402, H403, H404, H405, H406, H408, H409, H42, H420, H428, Q150
