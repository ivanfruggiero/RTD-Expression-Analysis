#!/usr/bin/env Rscript
# Gene-level differential expression for ONE study (one metadata file) with DESeq2.
#
# Usage: Rscript deseq2_per_study.R <metadata.csv> <quants_dir> <tx2gene.tsv> <out_dir> [ref_level] [treat_level]
#
# metadata.csv: one row per sample. Required columns: sample, condition.
#   Optional columns: time, genotype, cultivar, batch (used only when they have >= 2 levels).
#   'sample' must match the Salmon quantification folders (<quants_dir>/<sample>/quant.sf).
# tx2gene.tsv : transcript_id <TAB> gene_id, no header (step 07).
#
# Model (full data):  ~ [batch] + [time] + [genotype] + [cultivar] + condition
#                       + [time:condition] + [genotype:condition] + [cultivar:condition]
# Exported results (log2 fold changes are maximum-likelihood estimates, NOT shrunken, so that
# log2FoldChange and lfcSE can be used directly as effect size and standard error in step 14):
#   DE_all_main_<T>_vs_<C>.csv           treatment vs control at the reference level of the other factors
#   INTERACTION_<coef>.csv               interaction coefficients, if any
#   DE_all_<T>_vs_<C>_time_<level>.csv   if 'time' exists: separate model per time point,
#                                        ~ [batch] + [genotype] + [cultivar] + condition
#   DE_sig_*.csv                         subsets with padj < 0.05 and |log2FC| >= 1
#   counts_raw_global.csv, vst_matrix_global.csv, dds.rds, figures/PCA_global.png, figures/MA_*.png
#   contrasts.tsv                        index of the treatment-vs-control tables (input of step 14)

suppressPackageStartupMessages({
  library(tximport); library(DESeq2); library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) stop("usage: deseq2_per_study.R <metadata.csv> <quants_dir> <tx2gene.tsv> <out_dir> [ref] [treat]")
META <- args[1]; QUANTS <- args[2]; TX2GENE <- args[3]; OUTROOT <- args[4]
REF <- if (length(args) >= 5) args[5] else "C"
TRT <- if (length(args) >= 6) args[6] else "T"
TAG <- sub("\\.[^.]+$", "", basename(META))

PADJ_CUT <- 0.05; LFC_CUT <- 1
RESDIR <- file.path(OUTROOT, "results"); FIGDIR <- file.path(OUTROOT, "figures")
dir.create(RESDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGDIR, recursive = TRUE, showWarnings = FALSE)
lab <- paste0(TRT, "_vs_", REF)

message("[INFO] study=", TAG, " metadata=", META, " reference=", REF, " treatment=", TRT)
message("[INFO] R ", R.version.string, " | DESeq2 ", packageVersion("DESeq2"), " | tximport ", packageVersion("tximport"))

# ---------------- metadata
meta <- read.csv(META, stringsAsFactors = FALSE, check.names = FALSE)
if (!all(c("sample", "condition") %in% names(meta))) stop("metadata must contain 'sample' and 'condition'")
meta <- meta[!is.na(meta$sample) & meta$sample != "", , drop = FALSE]
meta$sample <- as.character(meta$sample)
if (anyDuplicated(meta$sample)) stop("duplicated sample IDs in ", META)
if (!all(c(REF, TRT) %in% meta$condition))
  stop("condition must contain both '", REF, "' and '", TRT, "' (found: ", paste(unique(meta$condition), collapse = ", "), ")")
meta <- meta[meta$condition %in% c(REF, TRT), , drop = FALSE]
meta$condition <- relevel(factor(meta$condition), ref = REF)
rownames(meta) <- meta$sample

# Optional factors: kept only if they have >= 2 levels.
opt_factors <- c("batch", "time", "genotype", "cultivar")
for (f in opt_factors) {
  if (f %in% names(meta)) {
    if (length(unique(meta[[f]])) >= 2) {
      meta[[f]] <- factor(meta[[f]])  # default (alphabetical) level order: first level = reference
    } else {
      message("[INFO] '", f, "' has a single level: ignored"); meta[[f]] <- NULL
    }
  }
}
has <- function(f) f %in% names(meta)

# ---------------- import
qf <- file.path(QUANTS, meta$sample, "quant.sf")
if (!all(file.exists(qf))) stop("missing quant.sf for: ", paste(meta$sample[!file.exists(qf)], collapse = ", "))
names(qf) <- meta$sample
tx2gene <- read.delim(TX2GENE, header = FALSE, stringsAsFactors = FALSE)[, 1:2]
colnames(tx2gene) <- c("transcript", "gene")
txi <- tximport(qf, type = "salmon", tx2gene = tx2gene,
                countsFromAbundance = "lengthScaledTPM", ignoreTxVersion = FALSE, ignoreAfterBar = TRUE)

# ---------------- full model
main_terms <- c(if (has("batch")) "batch", if (has("time")) "time",
                if (has("genotype")) "genotype", if (has("cultivar")) "cultivar", "condition")
int_factors <- intersect(c("time", "genotype", "cultivar"), names(meta))
int_terms <- if (length(int_factors)) paste0(int_factors, ":condition") else character(0)
design_full <- as.formula(paste("~", paste(c(main_terms, int_terms), collapse = " + ")))
message("[INFO] design = ", deparse(design_full))

dds <- DESeqDataSetFromTximport(txi, colData = meta, design = design_full)
dds <- DESeq(dds)
saveRDS(dds, file.path(RESDIR, "dds.rds"))

vsd <- if (nrow(dds) >= 1000) vst(dds, blind = TRUE) else varianceStabilizingTransformation(dds, blind = TRUE)
ig <- intersect(c("condition", "time", "genotype", "cultivar", "batch"), names(colData(vsd)))
ggsave(file.path(FIGDIR, "PCA_global.png"), plotPCA(vsd, intgroup = ig) + ggtitle(paste0("PCA (", TAG, ")")),
       width = 7, height = 5, dpi = 220)
write.csv(as.data.frame(counts(dds)), file.path(RESDIR, "counts_raw_global.csv"))
write.csv(as.data.frame(assay(vsd)), file.path(RESDIR, "vst_matrix_global.csv"))

plot_ma <- function(df, title, file) {
  df$sig <- ifelse(!is.na(df$padj) & df$padj < PADJ_CUT, "sig", "ns")
  p <- ggplot(df, aes(log10(baseMean + 1), log2FoldChange, alpha = sig)) + geom_point(size = 0.6) +
    geom_hline(yintercept = c(-LFC_CUT, LFC_CUT), linetype = "dashed") +
    scale_alpha_manual(values = c(sig = 1, ns = 0.2)) + theme_bw() +
    labs(title = title, x = "log10(baseMean + 1)", y = "log2 fold change")
  ggsave(file, p, width = 6, height = 4, dpi = 220)
}
write_res <- function(res, name, title) {
  df <- as.data.frame(res); df$gene_id <- rownames(df)
  write.csv(df, file.path(RESDIR, paste0("DE_all_", name, ".csv")), row.names = FALSE)
  sig <- df[!is.na(df$padj) & df$padj < PADJ_CUT & abs(df$log2FoldChange) >= LFC_CUT, ]
  write.csv(sig[order(sig$padj), ], file.path(RESDIR, paste0("DE_sig_", name, ".csv")), row.names = FALSE)
  plot_ma(df, title, file.path(FIGDIR, paste0("MA_", name, ".png")))
  message("[OK] ", name, ": ", sum(!is.na(df$padj)), " genes tested, ", nrow(sig), " DE (padj<", PADJ_CUT, ", |LFC|>=", LFC_CUT, ")")
  paste0("DE_all_", name, ".csv")
}

index <- data.frame(study = character(), contrast = character(), file = character())

# 1) treatment vs control at the reference level of the other factors
f_main <- write_res(results(dds, contrast = c("condition", TRT, REF)), paste0("main_", lab), paste(TAG, lab))

# 2) interaction coefficients
for (nm in grep(paste0("\\.condition", TRT, "$"), resultsNames(dds), value = TRUE)) {
  df <- as.data.frame(results(dds, name = nm)); df$gene_id <- rownames(df)
  write.csv(df, file.path(RESDIR, paste0("INTERACTION_", gsub("[^A-Za-z0-9_]+", "_", nm), ".csv")), row.names = FALSE)
}

# 3) per time point: separate model on the subset
if (has("time")) {
  for (tl in levels(meta$time)) {
    sel <- meta$time == tl
    mt <- droplevels(meta[sel, , drop = FALSE])
    if (length(unique(mt$condition)) < 2) { message("[WARN] time ", tl, ": only one condition, skipped"); next }
    terms <- c(Filter(function(f) has(f) && length(unique(mt[[f]])) >= 2, c("batch", "genotype", "cultivar")), "condition")
    txi_t <- list(counts = txi$counts[, sel, drop = FALSE], abundance = txi$abundance[, sel, drop = FALSE],
                  length = txi$length[, sel, drop = FALSE], countsFromAbundance = txi$countsFromAbundance)
    dt <- DESeqDataSetFromTximport(txi_t, colData = mt, design = as.formula(paste("~", paste(terms, collapse = " + "))))
    dt <- DESeq(dt)
    nm <- paste0(lab, "_time_", gsub("[^A-Za-z0-9_.-]+", "_", tl))
    f <- write_res(results(dt, contrast = c("condition", TRT, REF)), nm, paste(TAG, lab, "time", tl))
    index <- rbind(index, data.frame(study = TAG, contrast = tl, file = file.path(RESDIR, f)))
  }
} else {
  index <- rbind(index, data.frame(study = TAG, contrast = "main", file = file.path(RESDIR, f_main)))
}

# Contrasts entering the meta-analysis: one per time point if the study has a time factor,
# otherwise the single main contrast.
write.table(index, file.path(OUTROOT, "contrasts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(OUTROOT, "sessionInfo.txt"))
message("[OK] ", TAG, " completed")
