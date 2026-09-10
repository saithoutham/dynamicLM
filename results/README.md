# Result files

I keep the computed tables in this directory so that the reports can be checked
without rerunning the full simulation study. `trace.csv` is the index: it links
reported values to the script, function, seed, runtime, and commit that produced
them. Files ending in `seed_manifest.csv` contain the corresponding random
seeds.

## Why some failure files contain only a header

These five files intentionally have column names but no data rows:

- `simulation_failures.csv`
- `simulation_auc_coverage_failures.csv`
- `timescale_diagnostic_failures.csv`
- `informative_entry_failures.csv`
- `informative_auc_coverage_failures.csv`

They are not unfinished placeholders. Each analysis writes the same failure-log
schema whether or not a fit fails. A header with no rows records that the run
captured no failures. I kept these files because deleting an empty failure log
would make “no failures” harder to distinguish from “failure logging was never
run.”

## Large raw tables

The replicate-level files `simulation_raw.csv` and
`informative_entry_raw.csv` are not tracked because of their size. They can be
rebuilt from the committed seed manifests with:

```sh
REGENERATE_CORES=8 Rscript analysis/regenerate_raw.R --write
```

The byte-comparison record made before their removal is in
[`check/raw_regeneration_verification.csv`](check/raw_regeneration_verification.csv).
The smaller summaries and all seed manifests remain in Git.

## Reading the directory

- `*_summary.csv` files contain the report-level summaries.
- `*_raw.csv` files contain replicate- or test-level records when their size is
  reasonable for Git.
- `*_failures.csv` files retain captured error records, including the valid
  zero-row case described above.
- `check/` contains package-check and clean-clone records.
- `session.txt` records the R session used for the analyses.
