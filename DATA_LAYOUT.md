# Input data

Participant-level data are not distributed in this repository. Data access is
subject to the study's consent and sharing conditions; consult the accompanying
paper's data-availability statement and corresponding author. The code does not
download participant data.

Each MAT file contains one finite numeric variable. Connectomes are square
node × node × participant arrays; behaviors are participant × 1 vectors.
Subjects are aligned by row/index, not matched by identifiers.

```text
data/
├── behav/
├── connectomes/
├── connectomes_altparcel/
│   ├── no_cerebellum/
│   └── with_cerebellum/
├── covars.mat
└── covars_indicator.mat
```

## Covariates and ordering

`covars_indicator.mat` is 67 × 5:
`[study_2, study_3, sex, age, mean_FD]`.
`covars.mat` is 67 × 4:
`[study_id, sex, age, mean_FD]`.
Study labels define groups only; regression uses indicators selected for the
studies represented. Rows 1–19, 20–42, and 43–67 correspond to the 5-HT, LAM,
and LPM studies, respectively.

## Internal validation

Connectomes have 67 participants and 416 nodes:
`LSD.mat`, `LSD_difference.mat`, and `LSD_GSR_difference.mat` for the GSR
control. Only the connectomes needed by the chosen preset are required.

Behavior files are `LSD_<scale>.mat` and, where available,
`LSD_<scale>_difference.mat`, for BDE, GDE, OBN, AED, VRS, and MEQ30.
LSD MEQ30 vectors contain 48 values corresponding to rows 20–67 of the LSD
connectomes; other LSD vectors contain 67 values.

Difference analyses prefer the difference behavior file and otherwise use the
standard vector. In the study data this fallback is intentional for scales
whose placebo ratings are all zero, including BDE and AED: the vectors are
arithmetically identical. It is not a missing-data procedure.

## Cross-drug LOSO

Provide `psilocybin_difference.mat` and `mescaline_difference.mat` (25
participants each), and `amphetamine_difference.mat` and
`MDMA_difference.mat` (23 participants each), under connectomes/.
For GDE, MEQ30, and VRS, provide the corresponding
`<drug>_<scale>_difference.mat` behavior vector, or its equivalent standard
`<drug>_<scale>.mat` vector. Each individual-drug behavior vector has 25 or
23 rows, including MEQ30.

The distributed `metadata/cross_drug_row_map.csv` defines the correspondence
to LSD rows. Pooled results concatenate amphetamine then MDMA, and psilocybin
then mescaline. Pooled `all_amphs_difference.mat` (46 sessions) and
`other_psych_difference.mat` (50 sessions), and their behavior equivalents,
are optional consistency checks; the engine derives pooling from the
single-drug inputs. Standard-connectome overrides use the analogous filenames
without `_difference`.

## Alternative parcellations and metadata

Alternative-parcellation folders contain LSD and four individual-drug difference
connectomes with suffix `_altparcel_no_cerebellum` (432 nodes) or
`_altparcel_with_cerebellum` (439 nodes), before `.mat`.
For example: `LSD_difference_altparcel_no_cerebellum.mat`.
The same behavior and covariate files are used. Override the connectome directory
explicitly if these files are stored elsewhere.

The repository includes the 416-node network assignment and descriptive parcel
metadata, plus structure membership for the 439-node atlas. Sources are listed
in [metadata/atlas/SOURCES.md](metadata/atlas/SOURCES.md). Atlas images are not
distributed.
