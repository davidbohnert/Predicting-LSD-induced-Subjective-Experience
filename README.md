# Neural Signatures of Lysergic Acid Diethylamide–Induced Subjective Experience Identified via Connectome-Based Predictive Modeling

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22692846.svg)](https://doi.org/10.5281/zenodo.22692846)

Analysis and visualization code accompanying the revised manuscript by Bohnert
and colleagues. This repository includes internal validation, subject-wise
cross-drug generalization, and supplementary sensitivity analyses. Changes from
the earlier release are described in [CHANGELOG.md](CHANGELOG.md).

## Requirements and data

Requires **MATLAB R2022b or newer**, with the Statistics and Machine Learning
and Parallel Computing Toolboxes; tested with R2026a. R2022b is the minimum
for the `Processes` pool profile used here. Earlier compatible syntax includes
`arguments` blocks and name-value arguments. R2022b has not been tested directly.

Use Python 3.12 in a virtual environment for figure and table generation.
Run these commands from the repository root:

```bash
python3.12 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
```

On Windows, use `py -3.12 -m venv .venv` and `.venv\Scripts\python.exe`
in place of `.venv/bin/python`. Python commands below assume the repository
root is the working directory.

Participant-level data are not included. See [DATA_LAYOUT.md](DATA_LAYOUT.md)
for processed inputs, subject ordering, and data access information.

## Usage

```matlab
repo = '/path/to/Predicting-LSD-induced-Subjective-Experience';
addpath(genpath(repo));
cfg = paper_analysis_config('/path/to/data', '/path/to/outputs');
cfg.runtime.python_executable = fullfile(repo, '.venv', 'bin', 'python');
% Windows: fullfile(repo, '.venv', 'Scripts', 'python.exe')
internal_validation_main(cfg, "primary");
```

MATLAB invokes Python itself when rendering and verifying network figures.
Set `cfg.runtime.python_executable` to the absolute path of the environment
where you installed the requirements; activating it in a separate terminal
is not sufficient.

The primary specification uses positive edges from LSD–placebo difference
connectomes, adjusted for study, sex, age, and mean framewise displacement.
MEQ30 uses LSD rows 20–67. Study indicators are selected for the represented
studies; MEQ30 therefore uses study 2 as its reference. Feature selection uses
Spearman correlation for BDE, GDE, and AED and Pearson correlation for OBN,
VRS, and MEQ30. Prediction accuracy is evaluated with Pearson correlation.
Cross-drug LOSO excludes the matching participant's LSD observation before
feature selection and fitting.

Main inferential analyses default to 10,000 iterations; partition and
learning-curve analyses use 1,000 per fit. For internal CV (including partition,
learning-curve, and study-held-out tests) and cross-drug LOSO, the iteration
count **includes the observed ordering**: 10,000 means one observed fit and
9,999 shuffled fits; 1,000 means one observed and 999 shuffled fits.
Freedman–Lane network-interaction and overlap tests instead use 10,000 null
permutations **in addition to** the observed fit, with the usual +1 correction
in their empirical p-values. These conventions preserve the paper calculations.

Five process workers are requested.
Run one analysis at a time. Each pipeline records settings and saves progress;
incompatible checkpoints require a fresh output directory.

| Paper output | Analysis and instructions |
|---|---|
| Table 1; Tables S5–S9, S11–S12, S14 | [Internal validation](pipelines/internal_validation/README.md) |
| Table 2; Table S10 | [Cross-drug LOSO](pipelines/cross_drug_loso/README.md) |
| Figure S1; Table S17 | [Partition sensitivity and stability](pipelines/partition_sensitivity/README.md) |
| Figure S2 | [Learning curves](pipelines/learning_curves/README.md) |
| Table S13 | [Parcellation sensitivity](pipelines/parcellation_sensitivity/README.md) |
| Figures 1–3; Table S15 | [Network interactions](pipelines/network_interactions/README.md) |
| Table S16 | [Network overlap](pipelines/network_overlap/README.md) |

Tables S1–S4 and preprocessing/QC analyses are outside the CPM code's scope.
This is a code release; aggregate results and participant data are not bundled.
See [release_assets/README.md](release_assets/README.md) for the optional
aggregate-output assembler.

## Repository structure

```text
config/          Shared paths, analysis presets, and reproducibility settings
pipelines/       Internal validation, cross-drug LOSO, and sensitivity analyses
scripts/         Shared CPM functions, utilities, and optional release assembler
metadata/        Atlas descriptions and cross-drug row mapping
tests/           Synthetic MATLAB and Python checks; reference implementations
release_assets/  Instructions for optionally packaging aggregate results
```

## Checks

Small deterministic checks use synthetic data and do not run the paper analyses:

```matlab
run_all_tests
```

```bash
.venv/bin/python -m unittest discover -s tests -p 'test_*.py'
```

For a short serial input check, use a separate temporary output directory:

```matlab
internal_validation_main(cfg, "primary", no_iterations=3, workers=0, ...
    use_parallel=false, bootstrap_iterations=100, output_path=tempname);
```

## Citation and contributions

Please cite the accompanying manuscript and the software version used.
Machine-readable software citation information is in [CITATION.cff](CITATION.cff).
Version 2.0.0 accompanies the revised manuscript submission and is archived on
Zenodo: https://doi.org/10.5281/zenodo.22692846. The corresponding source is
available from the
[GitHub release](https://github.com/davidbohnert/Predicting-LSD-induced-Subjective-Experience/releases/tag/v2.0.0).
Earlier code remains available in Git history and the original `v1.0.0` tag.

The original CPM pipeline was developed collaboratively by Olivia M. F. Rapp
and David Bohnert under the supervision of Mihai Avram. David Bohnert developed
and maintained the revised analysis and visualization code accompanying this
manuscript revision.

The implementation was adapted for these datasets and partly guided by Boyle,
R., & Weng, Y. (2025), *Studying the Connectome at a Large Scale*, in R. Whelan
& H. Lemaître (Eds.), *Methods for Analyzing Large Neuroimaging Datasets*,
pp. 365–394. https://doi.org/10.1007/978-1-0716-4260-3_15.

Code is distributed under the [MIT License](LICENSE). Atlas sources are
documented in [metadata/atlas/SOURCES.md](metadata/atlas/SOURCES.md).
