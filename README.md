# Adherence and persistence to novel glucose-lowering medications in persons with type 2 diabetes mellitus undergoing routine care
DOI: https://doi.org/10.1016/j.diabres.2024.111745

## Abstract
**Aims:** To assess adherence and persistence to sodium-glucose cotransporter-2 inhibitors (SGLT2i), glucagon-like peptide-1 receptor agonists (GLP1-RA), and dipeptidyl peptidase-4 inhibitors (DPP4i) in routine care.

**Methods:** Using retrospective healthcare data from the Stockholm region, Sweden, we evaluated new-users of these agents during 2015-2020. We investigated adherence (≥80 % of days covered by an active supply), persistence (no treatment gap ≥ 60 days), and predictors for non-adherence and non-persistence.

**Results:** We identified 24,470 new-users of SGLT2i (10,743), GLP1-RA (10,315), and/or DPP4i (9,488). Over 2.8 years median follow-up, the proportion demonstrating adherence was higher for SGLT2i (57 %) than DPP4i (53 %, comparison p < 0.001), and for GLP1-RA than DPP4i (54 % vs 53 %, p < 0.001). Similarly, persistence was higher for both SGLT2i and GLP-RA than DPP4i (respectively, 50 % vs 44 %, p < 0.001; 49 % vs 44 %, p < 0.001). Overall adherence was better among users who were older, had a history of high blood pressure, used more non-diabetic medications, had lower Hba1c, had better kidney function, and had completed secondary schooling or university. Women had worse adherence to SGLT2i and GLP1-RA than DPP4i.

**Conclusions:** We report adherence and persistence to SGLT2i, GLP1-RA and DPP4i in routine care, and identify prognostic factors that could inform implementation interventions to improve uptake of these important therapies.

## Codes
- 1-cohort_derivation-alternative-20230202.R: Derivation of cohort and flowchart. Alternative because the structuring differed from a previous script to get the flowchart per drug of interest.
- 2-covariate_derivation-20221118.R: Derivation of covariates for all individuals.
- 3-outcomes-20221125.R: Derivation of outcomes for all individuals, specifically persistence and adherence (see below).
- 4-analysis-20221129.Rmd: Final analyses.
- adhper.R: Script for calculating adherence and persistence. See below.

## Determining adherence and persistence
Adherence and persistence were calculated,  taking into account stockpiling :pill:. This repository also contains the separate [code](./adhper.R) for this.

The script creates two functions:

```
adh(df, drug_df, yr)
per(df, drug_df, yr, gp)
```

The following data-specific variable names are used in the code and mean the following:
- lopnr           =   participant ID
- drug            =   drug of interest
- index_dt        =   participant's cohort entry date
- censor_dt_cod   =   participant's censor date
- edatum          =   dispensation date
- antal           =   amount of packages dispensed
- antnum          =   amount of units (e.g., pills) per package

The  diagrams shown below give an explanation of how the code works, with on the left a text-based explanation and on the right a visual example.

### Adherence
![Adherence function](https://github.com/user-attachments/assets/d91fb74e-530d-4d65-b140-d44d2603e3eb)

### Persistence
![Persistence function](https://github.com/user-attachments/assets/fbbe4b1b-ca3c-4e41-8d9f-f614cb97186b)

