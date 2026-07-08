# Data folder

ADNI data are not redistributed here.

Place the local input files in `data/raw/` before running the scripts:

- `Data.csv`: the baseline analytic file used for the main cross-sectional cohort.
- `Data_all.csv`: the companion file containing absolute choroid plexus volume and imaging/source metadata.
- `manual_chp_qc_integrated_80.csv` (optional): the 80-participant manually corrected ChP subset used for segmentation-sensitivity checks.

The scripts write cleaned analysis files to `data/clean/` and timestamped outputs to `results/`.
