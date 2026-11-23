# RTD-Expression-Analysis

A pipeline for public RNA-seq data from many independent studies. It builds a
**Reference Transcript Dataset (RTD)** from the data themselves, then quantifies gene expression
on it. Next it runs a differential expression analysis for each study and combines the studies in a
**gene-wise multilevel meta-analysis**. Finally it tests the resulting gene sets for GO enrichment.

The pipeline is generic and works for any organism with a reference genome and annotation.
Each step is a PBS job script that you can also run with plain `bash`.

```
raw reads
  └─ 01 FastQC/MultiQC
  └─ 02 fastp trimming
       └─ 03 contaminant screening ───────────────────────────┐
            └─ 05 STAR 2-pass  (◄ 04 STAR index)              │
                 └─ 06 Scallop + StringTie                    │
                      └─ 07 RTDmaker ─► RTD + tx2gene         │
                           ├─ 08 BUSCO (transcripts)          │
                           ├─ 09 TransDecoder (+ 10 BUSCO)    │
                           │    └─ 12 InterProScan (GO) ──────┼──────────┐
                           └─ 11 Salmon ◄─────────────────────┘          │
                                └─ 13 DESeq2 per study                   │
                                     └─ 14 meta-analysis                 │
                                          └─ 15 GO enrichment ◄──────────┘
```

## Contents

| Path | What it is |
|---|---|
| `script/NN_*.pbs` | The 15 pipeline steps |
| `script/common.sh` | Helpers that every step sources: environment set-up, input checks, sample discovery |
| `script/R/` | R code for steps 13–15 |
| `config/trimming_config.txt` | fastp parameters for each BioProject (template) |
| `config/example_contaminants.tsv` | Contaminant genomes for the screening step (template) |
| `config/example_study_groups.tsv` | Example grouping of studies for the meta-analysis |
| `metadata/` | One CSV for each study for DESeq2 (format: [metadata/README.md](metadata/README.md)) |
| `paper_config/` | Inputs and settings of the published analysis ([paper_config/README.md](paper_config/README.md)) |
| `logs/` | Scheduler logs |

Raw reads, reference files and results are not tracked by Git (see `.gitignore`).

## Running a step

Run every step **from the repository root**. All paths are relative to the root, and the defaults
chain the steps together, so you normally only set the inputs of the first steps.

```bash
# With a PBS scheduler (PBS Pro / OpenPBS / Torque): set parameters as environment variables with -v
qsub -v GENOME_FA=reference/genome.fa,REF_GTF=reference/annotation.gtf script/04_build_star_index.pbs

# Without a scheduler, or interactively: use the command-line options
bash script/04_build_star_index.pbs -f reference/genome.fa -g reference/annotation.gtf -t 8

# Help for any step
bash script/04_build_star_index.pbs -h
```

`qsub` does not pass command-line arguments to a job script, so on a cluster always use `-v`.
Each option of a step has an environment variable of the same meaning; `-h` lists both.

**Software.** Each step checks that its tools are in `PATH` and stops with a clear message if they
are not. To load them inside the job, set either variable:

- `MODULES="star/2.7.8a samtools/1.13"` to load environment modules
- `CONDA_ENV=myenv` to activate a conda environment

For example: `qsub -v MODULES="salmon/1.9.0" script/11_salmon_quantification.pbs`.

**Resources.** The `#PBS -l select=...:ncpus=...:mem=...` and `walltime` lines in each script are
starting values. Change them to fit your data and your cluster. Add `#PBS -q <queue>`
if your site requires a queue. Thread counts default to the CPUs the scheduler grants (`NCPUS`).

**Re-runs.** Steps 03, 05, 06 and 11 skip samples that already have an output; `-f` forces a re-run.

## Input layout

```
raw_reads/
  paired/<BIOPROJECT>/<SAMPLE>_1.fastq.gz  <SAMPLE>_2.fastq.gz
  single/<BIOPROJECT>/<SAMPLE>.fastq.gz
reference/
  genome.fa            # genome FASTA
  annotation.gtf       # reference annotation (GTF)
config/trimming_config.txt
config/contaminants.tsv  # contaminant genomes for step 03
metadata/<STUDY>.csv   # one file per study, sample IDs = <SAMPLE> above
```

`<SAMPLE>` is normally the SRA run accession. It must be unique across all studies, because it
names the Salmon quantification folder and links reads to metadata.

## Steps

| # | Script | Does | Main output (under `results/`) |
|---|---|---|---|
| 01 | `01_fastqc_multiqc.pbs` | FastQC and MultiQC on raw reads | `01_fastqc_raw/` |
| 02 | `02_trimming_qc.pbs` | fastp trimming with per-BioProject parameters, then FastQC/MultiQC | `02_trimming/trimmed/{paired,single}/<BP>/` |
| 03 | `03_contamination_screening.pbs` | Remove reads of contaminant organisms (HISAT2 to contaminants, re-check on the host) | `03_screening/clean/{paired,single}/<BP>/`, `03_screening/summary.tsv` |
| 04 | `04_build_star_index.pbs` | STAR genome index (`OVERHANG` = read length − 1) | `04_star_index/` |
| 05 | `05_star_2pass_alignment.pbs` | Per-sample STAR 2-pass alignment (paired and single-end) | `05_star/pass1/*SJ.out.tab`, `05_star/pass2/*.bam` |
| 06 | `06_assembly.pbs` | Reference-guided assembly of each BAM with Scallop and/or StringTie | `06_assembly/{scallop,stringtie}/*.gtf` |
| 07 | `07_rtdmaker.pbs` | Merge the assemblies into the RTD with RTDmaker ShortReads; build `tx2gene.tsv` | `07_rtdmaker/<OUTNAME>_RTDmaker_output/`, `07_rtdmaker/tx2gene.tsv` |
| 08 | `08_busco_transcriptome.pbs` | BUSCO on RTD transcripts (`LINEAGE` required) | `08_busco_rtd/` |
| 09 | `09_transdecoder.pbs` | ORF and peptide prediction (TransDecoder LongOrfs + Predict) | `09_transdecoder/*.pep` |
| 10 | `10_busco_peptides.pbs` | BUSCO on predicted peptides (`LINEAGE` required) | `10_busco_peptides/` |
| 11 | `11_salmon_quantification.pbs` | Salmon index on the RTD and quantification of each sample | `11_salmon/quants/<SAMPLE>/quant.sf` |
| 12 | `12_interproscan.pbs` | InterProScan annotation with GO terms (stop codons removed first) | `12_interproscan/*.tsv` |
| 13 | `13_deseq2.pbs` | DESeq2 for each study, one study per metadata file | `13_deseq2/<STUDY>/`, `13_deseq2/contrasts.tsv` |
| 14 | `14_meta_analysis.pbs` | Gene-wise multilevel random-effects meta-analysis | `14_meta_analysis/<GROUP>/` |
| 15 | `15_go_enrichment.pbs` | GO over-representation analysis (BP/MF/CC) against the tested genes | `15_go_enrichment/` |

### Step details

**02 – trimming.** `config/trimming_config.txt` holds one line per BioProject:
`BIOPROJECT TRIM_FRONT1 TRIM_FRONT2 MIN_PHRED N_BASE_LIMIT`. fastp also runs with
`--unqualified_percent_limit 50` and, for paired-end reads, `--detect_adapter_for_pe`. Choose the 5′
trimming from the FastQC per-base sequence content of step 01.

**03 – contaminant screening.** Reads are aligned with HISAT2 to each genome in `CONTAMINANTS`
(`name<TAB>source`, where source is an NCBI assembly accession or a FASTA path). Reads with a primary
alignment to any contaminant are re-aligned to the host genome (`HOST_FA`). Those that do not align to
the host are removed; reads that align to both are kept. For paired-end data the whole pair is
removed. Output files keep the step-02 names, so later steps read `results/03_screening/clean` by
default. To skip screening, set `READS_DIR=results/02_trimming/trimmed` in steps 05, 07 and 11.

**05 – alignment.** Each sample gets two passes. Pass 1 finds splice junctions. Pass 2 re-aligns with
that sample's junctions, using `--outSAMstrandField intronMotif` and the `NH HI AS nM XS` attributes
that the assemblers need. Pass-1 `SJ.out.tab` files are also inputs to RTDmaker.

**07 – RTD.** The step collects the assemblies, junction files and screened reads into flat input
folders (RTDmaker's expected layout). It then runs RTDmaker ShortReads with
`--SJ-reads 5 1 --tpm 0.5 1 --fragment-len 0.70 --antisense-len 0.50`; the environment variables
`SJ_READS`, `TPM`, `FRAGMENT_LEN` and `ANTISENSE_LEN` change these thresholds. `OUTNAME` and `PREFIX`
(default `RTD`) set the file names and the prefix of new gene and transcript IDs. Downstream steps use:

- `<OUTNAME>.fixed.fa`: transcript sequences, input to BUSCO and TransDecoder
- `<OUTNAME>.salmon_unique.fa`: transcripts with unique sequences, the Salmon index
- `tx2gene.tsv`: `transcript_id → gene_id` from `<OUTNAME>.fixed.gtf`

If you change `OUTNAME`, pass the new file names to steps 08, 09 and 11, and to every step that
follows them.

**13 – differential expression.** Each `metadata/<STUDY>.csv` is analysed on its own. Counts are
imported with `tximport` (`countsFromAbundance = "lengthScaledTPM"`) and summarised to genes. The
model is

```
~ [batch] + [time] + [genotype] + [cultivar] + condition + [time:condition] + [genotype:condition] + [cultivar:condition]
```

Each optional factor enters only if it has at least 2 levels. The step exports these results:

- **Main contrast** `condition T vs C`, at the reference level of the other factors.
- **Interaction coefficients**.
- **One separate model per time point**, if the study has a `time` column:
  `~ [batch] + [genotype] + [cultivar] + condition`, fitted on the samples of that time point.

Log2 fold changes are **maximum-likelihood estimates, not shrunken**. Their standard errors
(`lfcSE`) are the within-study variances that the meta-analysis needs. Significance tables
(`DE_sig_*`) use padj < 0.05 and |log2FC| ≥ 1.

`contrasts.tsv` lists the tables that enter the meta-analysis: one per time point for studies with
`time`, otherwise the main contrast. To use other levels than `C`/`T`, set `REF_LEVEL` and
`TREAT_LEVEL`.

**14 – meta-analysis.** For each gene: `yi = log2FC`, `vi = lfcSE²` and

```r
metafor::rma.mv(yi, vi, random = ~ 1 | study/contrast, method = "REML")
```

Contrasts (for example, time points) are nested within studies. A gene is tested only when it has at
least 3 contrasts from at least 2 studies (`MIN_CONTRASTS`, `MIN_STUDIES`). P-values are corrected by
Benjamini–Hochberg over the genes tested. Significant genes have FDR < 0.05 and |estimate| ≥ 1
(`FDR_CUT`, `LFC_CUT`). The output also has between-study and within-study heterogeneity
(`tau2_study`, `tau2_contrast`) and Cochran's Q.

With `STUDY_GROUPS=<file>` (a `study<TAB>group` table; see `config/example_study_groups.tsv`), each
group is meta-analysed separately. Use one group for each biological question that can be compared.
Do not pool unrelated treatments in one model in order to compare them; compare the results of
separate groups instead.

**15 – GO enrichment.** GO terms come from the InterProScan output. They pass from each peptide to
its transcript, then to its gene. The test is an over-representation test (`clusterProfiler::enricher`,
BH correction, terms with 06–500 genes, q < 0.2). It runs separately for BP, MF and CC. The universe
is the **genes actually tested** in the analysis that produced the set, limited to annotated genes;
sets with fewer than 5 annotated genes are not tested. By default the sets are the UP and DOWN genes
of each group from step 14. For other sets (a single study, or the overlap between groups with the
intersection of their tested genes as universe), list them in your own `SETS` table with the columns
`set genes_file universe_file`.

## Software

Versions used for the reference run:

| Tool | Version | | Tool | Version |
|---|---|---|---|---|
| FastQC | 0.11.9 | | RTDmaker | 1.0.0 |
| MultiQC | 1.11 | | BUSCO | 6.0.0 |
| fastp | 1.0.1 | | TransDecoder | 5.5.0 |
| STAR | 2.7.8a | | Salmon | 1.9.0 |
| Scallop | 0.10.5 | | InterProScan | 5.68-100.0 |
| StringTie | 2.2.1 | | R (step 13, cluster) | 4.4.1 |
| samtools | 1.13 | | HISAT2 | 2.2.1 |
| seqtk | 1.2-r95 | | NCBI datasets | 16.6.0 |

R packages: `tximport`, `DESeq2`, `ggplot2` (step 13); `metafor` (step 14); `clusterProfiler`,
`GO.db`, `AnnotationDbi` (step 15). The reference run used R 4.3.1 / Bioconductor 3.18 for steps 14–15
(metafor 4.8.0, clusterProfiler 4.10.1, GO.db 3.18.0, AnnotationDbi 1.64.1). Steps 13–15 write `sessionInfo.txt` with the exact package
versions of each run. RTDmaker is installed from <https://github.com/anonconda/RTDmaker>;
`RTDmaker.py` must be in `PATH`, or its folder set in `RTDMAKER_PREFIX`.

## License

See [LICENSE](LICENSE).
