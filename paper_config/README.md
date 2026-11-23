# Configuration of the published analysis

This folder holds the study-specific inputs and settings used for the associated article. The
pipeline code in `script/` is generic; these files make it reproduce the published run.

| File | Content |
|---|---|
| `trimming_config.txt` | fastp parameters (the same for all BioProjects, as in the fastp logs of the original run) |
| `contaminants.tsv` | Genomes screened for contaminant reads in step 03 (phiX174, *S. cerevisiae*, *E. coli*) |
| `metadata/<STUDY>.csv` | DESeq2 sample tables, one per study (step 13) |
| `study_groups.tsv` | Studies meta-analysed together in step 14 |
| `gene2go.tsv` | Gene–GO table used for enrichment (63,270 associations, 16,858 genes, 4,326 terms; from the InterProScan annotation) |

## Reference files

| Role | File |
|---|---|
| Host genome (`GENOME_FA`, `HOST_FA`) | *Citrus sinensis* DVS_A1.0, NCBI RefSeq assembly GCF_022201045.2 (genomic FASTA) |
| STAR index annotation (step 04) | `GCF_022201045.2_DVS_A1.0_genomic.gtf` (RefSeq GTF of the same assembly) |
| StringTie guide (step 06) | `GCF_022201045.2_refseq_nuclear.stringtie.gtf` (the RefSeq GTF restricted to nuclear sequences) |
| RTDmaker references (step 07) | both GTFs above |

## Run settings

Run from the repository root, after placing the raw reads under `raw_reads/` and the files above
under `reference/` (`genome.fa`, `annotation.gtf` = RefSeq GTF, `annotation_nuclear.gtf` = nuclear GTF):

```bash
qsub script/01_fastqc_multiqc.pbs
qsub -v CONFIG=paper_config/trimming_config.txt                          script/02_trimming_qc.pbs
qsub -v CONTAMINANTS=paper_config/contaminants.tsv                         script/03_contamination_screening.pbs
qsub -v OVERHANG=99                                                         script/04_build_star_index.pbs
qsub                                                                        script/05_star_2pass_alignment.pbs
qsub -v REF_GTF=reference/annotation_nuclear.gtf                           script/06_assembly.pbs
qsub -v REF_GTF="reference/annotation.gtf reference/annotation_nuclear.gtf",OUTNAME=Citrus_RTD,PREFIX=CitrusRTD script/07_rtdmaker.pbs
R=results/07_rtdmaker/Citrus_RTD_RTDmaker_output
qsub -v LINEAGE=eudicots_odb10,INPUT=$R/Citrus_RTD.fixed.fa              script/08_busco_transcriptome.pbs
qsub -v INPUT_FA=$R/Citrus_RTD.fixed.fa                                    script/09_transdecoder.pbs
qsub -v LINEAGE=eudicots_odb10,INPUT=results/09_transdecoder/Citrus_RTD.fixed.fa.transdecoder.pep script/10_busco_peptides.pbs
qsub -v RTD_FA=$R/Citrus_RTD.salmon_unique.fa                             script/11_salmon_quantification.pbs
qsub -v INPUT_PEP=results/09_transdecoder/Citrus_RTD.fixed.fa.transdecoder.pep script/12_interproscan.pbs
qsub -v METADATA_DIR=paper_config/metadata                                 script/13_deseq2.pbs
qsub -v STUDY_GROUPS=paper_config/study_groups.tsv                         script/14_meta_analysis.pbs
qsub -v GENE2GO=paper_config/gene2go.tsv                                  script/15_go_enrichment.pbs
```

Each `qsub` must wait for the previous step (for example with `qsub -W depend=afterok:<jobid>`).

Notes:

- The two BioProjects PRJNA739184 and PRJNA739186 form one experiment and are analysed as the
  single study `PRJNA739184_739186_MERGED` (one metadata file).
- PRJNA838230 is the only study of its group. It has no meta-analysis, so it is not listed in
  `study_groups.tsv`; its DESeq2 results (`results/13_deseq2/PRJNA838230/`) are used directly.
- Step 15 runs by default on the UP/DOWN genes of each meta-analysis group. The other published
  gene sets are tested with the same script through a custom `SETS` table: viral UP/DOWN genes
  (universe: genes with a non-missing adjusted p-value in PRJNA838230), commonly regulated genes
  of each pair of classes (universe: genes tested in both classes) and the exploratory bacterial
  set (nominal p < 0.05 and |estimate| >= 1 in `results/14_meta_analysis/BACTERIA/meta_results.tsv`).
- All other thresholds are the pipeline defaults described in the main README.
