# Study metadata for step 13 (DESeq2)

Write one CSV file for each study. The file name, without `.csv`, is the study name used in all
later steps. Files whose names start with `example_` are templates, and step 13 skips them.

| Column | Required | Meaning |
|---|---|---|
| `sample` | yes | Sample ID, the name of the Salmon folder `results/11_salmon/quants/<sample>/` |
| `condition` | yes | `C` = control, `T` = treated (other labels: set `REF_LEVEL` / `TREAT_LEVEL`) |
| `time` | no | Time point or sampling stage. Each level also gets its own T vs C model |
| `genotype` | no | Host genotype. Enters as main effect and `genotype:condition` |
| `cultivar` | no | Cultivar. Enters as main effect and `cultivar:condition` |
| `batch` | no | Technical batch (main effect only) |

An optional column is used only when it has at least two levels in the study. Factor levels are
sorted alphabetically, so the first level is the reference.

Samples whose `condition` is neither the control nor the treated level are ignored. Several
BioProjects that form a single experiment can share one metadata file; add a `batch` column if
needed.

Examples:

- `example_condition.csv`: control vs treated only
- `example_condition_time.csv`: control vs treated at two time points
- `example_condition_genotype_batch.csv`: two genotypes and two batches
