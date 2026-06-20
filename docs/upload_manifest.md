# Upload manifest

## Included

- Public R analysis scripts.
- MIMIC-IV cohort extraction SQL.
- Aggregate result CSV files.
- Compressed PNG figures.
- Documentation for data access, variable definitions, and workflow.

## Excluded

- MIMIC-derived patient-level CSV files.
- Linyi raw or imputed patient-level files.
- Names, hospital numbers, and other patient identifiers.
- MICE RDS objects and RData workspaces.
- Cover letters, author submission forms, proof replies, similarity checks, and
  manuscript DOCX/PDF files.
- Large journal-formatted TIFF figures.

## Final check before GitHub upload

Before pushing publicly, run:

```powershell
Get-ChildItem -Recurse github_release_pni_scmdeath
rg -n -i "姓名|住院号|Cover_Letter|Similarity|authors_proof|\\.RData|\\.rds|mids|医院数据|医院脓毒症" github_release_pni_scmdeath
```

The second command should return no actual included private files. It may still
find warning text in documentation or ignored private-data filename examples.

