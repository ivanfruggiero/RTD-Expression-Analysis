# RNA-seq Pipeline Scripts

This folder contains clean, reproducible PBS scripts for RNA-seq data analysis.

## Recommended repository structure

- `github/raw_reads/`
  - `paired/`
    - `<BIOPROJECT>/`
      - `sample_1.fastq.gz`
      - `sample_2.fastq.gz`
  - `single/`
    - `<BIOPROJECT>/`
      - `sample.fastq.gz`
- `github/results/`
  - analysis outputs
- `github/logs/`
  - execution logs
- `github/reference/`
  - reference files for STAR and other steps
- `github/script/`
  - executable PBS scripts
- `github/trimming_config.txt`
  - trimming parameter configuration


## Running scripts with parameters

The scripts are designed to receive command-line options, so you can submit jobs flexibly:

```bash
qsub script/02_trimming_qc.pbs -i ../raw_reads -o ../results/trimming_qc -c trimming_config.txt -t 10
```

### Supported options for `02_trimming_qc.pbs`

- `-i RAWDIR` : base folder of raw reads (default `../raw_reads`)
- `-o BASEDIR` : base folder for results (default `../results/trimming_qc`)
- `-c CONFIG` : trimming configuration file (default `trimming_config.txt`)
- `-t THREADS` : total number of threads to use (default `10`)
- `-h` : show this help message

## STAR index

The script `script/03_build_star_index.pbs` builds a STAR index reproducibly using generic reference files.

### Example execution

```bash
qsub script/03_build_star_index.pbs -r ../reference -f ../reference/genome.fa -g ../reference/annotations.gtf -o ../results/star_index -t 12 -s 99
```

### Supported options for `03_build_star_index.pbs`

- `-r REFDIR` : reference file directory (default `../reference`)
- `-f FASTA` : genome FASTA file (default `../reference/genome.fa`)
- `-g GTF` : annotation GTF file (default `../reference/annotations.gtf`)
- `-o OUTDIR` : STAR index output directory (default `../results/star_index`)
- `-t THREADS` : STAR thread count (default `12`)
- `-s OVERHANG` : `sjdbOverhang` value (default `99`)
- `-h` : show this help message

## STAR 2-pass alignment

The script `script/04_star_2pass_alignment.pbs` aligns trimmed files produced by `02_trimming_qc.pbs` in two passes.

### Example execution

```bash
qsub script/04_star_2pass_alignment.pbs -i ../results/trimming_qc/trimmed -r ../results/star_index -o ../results/star_alignment -t 12
```

### Supported options for `04_star_2pass_alignment.pbs`

- `-i INPUTDIR` : base folder for trimmed reads (default `../results/trimming_qc/trimmed`)
- `-r STAR_INDEX` : STAR index directory (default `../results/star_index`)
- `-o OUTDIR` : base folder for outputs (default `../results/star_alignment`)
- `-t THREADS` : STAR thread count (default `12`)
- `-s SJDB` : `sjdbOverhang` value (default `99`)
- `-m RAM` : `limitBAMsortRAM` in bytes (default `110000000000`)
- `-h` : show this help message

### Expected output

- `../results/star_alignment/pass1/` : STAR first-pass output
- `../results/star_alignment/pass2/` : aligned BAMs from the second pass
- `../results/star_alignment/logs/` : per-sample logs

## Assembly with Scallop and StringTie

The script `script/05_assembly.pbs` assembles transcripts from STAR pass2 BAMs. It supports Scallop, StringTie, or both in a single job.

### Example execution

```bash
qsub script/05_assembly.pbs -b ../results/star_alignment/pass2 -o ../results/assembly -a both -g ../reference/annotations.gtf -t 12
```

### Supported options for `05_assembly.pbs`

- `-b BAMDIR` : directory of STAR pass2 BAM files (default `../results/star_alignment/pass2`)
- `-o OUTDIR` : base folder for assembly results (default `../results/assembly`)
- `-a ASSEMBLERS` : `scallop`, `stringtie`, or `both` (default `both`)
- `-g GTF` : reference GTF for StringTie (default `../reference/annotations.gtf`)
- `-t THREADS` : threads for StringTie / samtools (default `8`)
- `-f` : force rerun even if output already exists
- `-h` : show this help message

### Expected output

- `../results/assembly/scallop/` : Scallop output GTFs
- `../results/assembly/stringtie/` : StringTie output GTFs
- `../results/assembly/logs/` : job logs

## RTDmaker

The script `script/06_rtdmaker.pbs` builds a Reference Transcriptome Database (RTD) using assembly GTFs, STAR SJ files, and cleaned FASTQ.

### Example execution

```bash
qsub script/06_rtdmaker.pbs -a ../results/assembly -s ../results/star_alignment/pass1 -f ../results/trimming_qc/trimmed -r ../reference -g ../reference/genome.fa -t ../reference/annotations.gtf -o ../results/rtdmaker -n MyRTD -p MyPrefix -m 8
```

### Supported options for `06_rtdmaker.pbs`

- `-a ASSEMBLY_DIR` : directory containing assembly GTFs (default `../results/assembly`)
- `-l ASSEMBLY_LIST` : file listing paths to assembly GTFs (default: built automatically from `ASSEMBLY_DIR`)
- `-s SJ_DIR` : directory containing `SJ.out.tab` files (default `../results/star_alignment/pass1`)
- `-f FASTQ_DIR` : directory of cleaned FASTQ files (default `../results/trimming_qc/trimmed`)
- `-r REFDIR` : reference folder (default `../reference`)
- `-g GENOME_FA` : genome FASTA file (default `../reference/genome.fa`)
- `-t REF_GTF` : reference GTF file for RTDmaker (default `../reference/annotations.gtf`)
- `-o OUTDIR` : RTDmaker output directory (default `../results/rtdmaker`)
- `-n OUTNAME` : RTD output name (default `Citrus_RTD`)
- `-p PREFIX` : RTDmaker prefix (default `CitrusRTD`)
- `-m RAM` : RAM in GB for RTDmaker (default `8`)
- `-k KEEP` : keep level (`intermediary` by default)
- `-h` : show this help message

### Expected output

- `../results/rtdmaker/` : full RTDmaker output
- `../results/rtdmaker/input/assemblies/` : linked assembly GTFs
- `../results/rtdmaker/reference_gtf/` : reference GTF for RTDmaker
- `../results/rtdmaker/logs/` : job logs

## BUSCO transcriptome

The script `script/07_busco_transcriptome.pbs` assesses transcriptome completeness with BUSCO.

### Example execution

```bash
qsub script/07_busco_transcriptome.pbs -i ../results/rtdmaker/Citrus_RTD.fa -d ../reference/eudicots_odb10 -o ../results/busco -n Citrus_RTD_busco -c 12
```

### Supported options for `07_busco_transcriptome.pbs`

- `-i INPUT_FA` : input transcriptome FASTA (default `../results/rtdmaker/${OUTNAME}.fa`)
- `-d DBDIR` : BUSCO lineage database directory (default `../reference/eudicots_odb10`)
- `-m MODE` : BUSCO mode (default `transcriptome`)
- `-o OUTDIR` : BUSCO output directory (default `../results/busco`)
- `-n OUTNAME` : BUSCO run name (default `busco_transcriptome`)
- `-c THREADS` : BUSCO thread count (default `12`)
- `-h` : show this help message

### Expected output

- `../results/busco/<OUTNAME>/` : full BUSCO output

## TransDecoder

The script `script/08_transdecoder.pbs` predicts ORFs (Open Reading Frames) in the transcriptome using TransDecoder.

### Example execution

```bash
qsub script/08_transdecoder.pbs -i ../results/rtdmaker/Citrus_RTD.fa -o ../results/transdecoder
```

### Supported options for `08_transdecoder.pbs`

- `-i INPUT_FA` : input transcriptome FASTA (default `../results/rtdmaker/Citrus_RTD.fa`)
- `-o OUTDIR` : TransDecoder output directory (default `../results/transdecoder`)
- `-h` : show this help message

### Expected output

- `../results/transdecoder/` : full TransDecoder output (.pep, .cds, .gff3, etc.)

## BUSCO peptides

The script `script/09_busco_peptides.pbs` assesses the completeness of peptides predicted by TransDecoder using BUSCO.

### Example execution

```bash
qsub script/09_busco_peptides.pbs -i ../results/transdecoder/Citrus_RTD.fa.transdecoder.pep -d ../reference/eudicots_odb10 -o ../results/busco_peptides -n Citrus_pep_busco -c 12
```

### Supported options for `09_busco_peptides.pbs`

- `-i INPUT_PEP` : input peptide FASTA (default `../results/transdecoder/Citrus_RTD.fa.transdecoder.pep`)
- `-d DBDIR` : BUSCO lineage database directory (default `../reference/eudicots_odb10`)
- `-m MODE` : BUSCO mode (default `protein`)
- `-o OUTDIR` : BUSCO output directory (default `../results/busco_peptides`)
- `-n OUTNAME` : BUSCO run name (default `busco_peptides`)
- `-c THREADS` : BUSCO thread count (default `12`)
- `-h` : show this help message

### Expected output

- `../results/busco_peptides/<OUTNAME>/` : full BUSCO peptide output

## Salmon quantification

The script `script/10_salmon_quantification.pbs` builds the Salmon index and quantifies expression for each sample.

### Example execution

```bash
qsub script/10_salmon_quantification.pbs -t ../results/rtdmaker/Citrus_RTD.fa -s salmon_samples.tsv -i ../results/salmon/index -o ../results/salmon/quants -p 12
```

### Supported options for `10_salmon_quantification.pbs`

- `-t TRANSCRIPTOME_FA` : transcriptome FASTA for indexing (default `../results/rtdmaker/Citrus_RTD.fa`)
- `-s SAMPLES_TSV` : TSV file with sample info: sample_id\tr1\tr2 (default `salmon_samples.tsv`)
- `-i INDEX_DIR` : Salmon index output directory (default `../results/salmon/index`)
- `-o QUANTDIR` : quantification output directory (default `../results/salmon/quants`)
- `-m MAP_FILE` : TSV file mapping sample to quant directory (default `../results/salmon/sample_to_quantdir.tsv`)
- `-p THREADS` : Salmon threads (default `12`)
- `-f` : force rerun even if output exists
- `-h` : show this help message

### Expected output

- `../results/salmon/index/` : Salmon index
- `../results/salmon/quants/<sample>/` : per-sample quantifications (`quant.sf`)
- `../results/salmon/sample_to_quantdir.tsv` : sample-to-quant directory mapping

### Samples TSV format

The file `salmon_samples.tsv` should have a header and TSV rows:

```tsv
sample_id	r1	r2
sample1	/path/to/sample1_R1.fastq.gz	/path/to/sample1_R2.fastq.gz
sample2	/path/to/sample2_R1.fastq.gz	/path/to/sample2_R2.fastq.gz
```

## InterProScan annotation

The script `script/11_interproscan.pbs` annotates predicted peptides with InterProScan using multiple databases.

### Example execution

```bash
qsub script/11_interproscan.pbs -i ../results/transdecoder/Citrus_RTD.fa.transdecoder.pep -o ../results/interproscan -c 12
```

### Supported options for `11_interproscan.pbs`

- `-i INPUT_PEP` : input peptide FASTA (default `../results/transdecoder/Citrus_RTD.fa.transdecoder.pep`)
- `-o OUTDIR` : InterProScan output directory (default `../results/interproscan`)
- `-t TMPDIR` : temporary directory (default `../results/interproscan/tmp`)
- `-c THREADS` : CPU threads (default `12`)
- `-f FORMATS` : output formats: TSV,XML,GFF3 (default `TSV,XML,GFF3`)
- `-h` : show this help message

### Expected output

- `../results/interproscan/` : InterProScan output in TSV, XML, GFF3 formats

## DESeq2 analysis

The script `script/12_deseq_analysis.pbs` performs differential expression analysis using DESeq2 on Salmon data.

### Example execution

```bash
qsub script/12_deseq_analysis.pbs -m ../metadata -q ../results/salmon/quants -t ../results/rtdmaker/tx2gene.tsv -o ../results/deseq
```

### Supported options for `12_deseq_analysis.pbs`

- `-m METADATA_DIR` : directory with metadata CSV files (default `../metadata`)
- `-q QUANTS_DIR` : directory with Salmon quantifications (default `../results/salmon/quants`)
- `-t TX2GENE` : tx2gene TSV file (default `../results/rtdmaker/tx2gene.tsv`)
- `-o OUTDIR` : DESeq2 output directory (default `../results/deseq`)
- `-r RSCRIPT` : path to Rscript (default `Rscript`)
- `-h` : show this help message

### Metadata CSV format

CSV files in `metadata/` should have columns such as:

```csv
sample,condition,time,genotype,batch
sample1,C,0,WT,1
sample2,T,0,WT,1
sample3,C,24,WT,2
sample4,T,24,WT,2
```

### Expected output

- `../results/deseq/<dataset>/results/` : DESeq2 results for each contrast (CSV, RDS)
- `../results/deseq/<dataset>/figures/` : MA plots for each contrast

## Trimming configuration example

The `trimming_config.txt` file should contain lines with the following fields separated by spaces:

```text
BIOPROJECT TRIM1 TRIM2 PHRED NLIMIT
```

Example:

```text
PRJ001 10 10 20 5
PRJ002 5 5 20 5
```

Lines starting with `#` are ignored.

## Software requirements

The scripts require these tools to be available in `PATH`:

- `fastp`
- `fastqc`
- `multiqc`

If using a module system or Conda environment, activate it before submitting jobs.

