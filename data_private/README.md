# Private data folder

This folder is intentionally empty in the public release.

Authorized users may place local patient-level data here to rerun selected
scripts. Do not commit or upload files from this folder.

Expected local filenames for the main analysis scripts:

- `mimic队列插补后.csv`: MIMIC-IV analysis dataset after local extraction and
  allowed local imputation.
- `临沂队列插补后.xlsx`: local single-center exploratory cohort dataset after
  approved local preprocessing.

These files contain or derive from patient-level clinical data and are not part
of the public repository.
