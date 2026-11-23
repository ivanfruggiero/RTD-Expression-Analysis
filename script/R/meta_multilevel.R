#!/usr/bin/env Rscript
# Gene-wise multilevel random-effects meta-analysis of DESeq2 log2 fold changes (metafor::rma.mv).
#
# Usage: Rscript meta_multilevel.R <contrasts.tsv> <out_dir> [groups.tsv] [threads]
#
# contrasts.tsv : study <TAB> contrast <TAB> file   (written by step 13; one row per
#                 treatment-vs-control table; each file needs gene_id, log2FoldChange, lfcSE, padj)
# groups.tsv    : optional, study <TAB> group. Studies are meta-analysed separately within each
#                 group; studies not listed are excluded. Without it all studies form one group "all".
#                 Use one group per comparable biological question: pooling unrelated treatments in a
#                 single model is not a valid basis for comparing them.
#
# For each gene and group:
#   yi = log2FoldChange, vi = lfcSE^2
#   rma.mv(yi, vi, random = ~ 1 | study/contrast, method = "REML")
#   i.e. contrasts (e.g. time points) nested within studies.
# A gene is tested only with >= MIN_CONTRASTS contrasts from >= MIN_STUDIES studies (env vars,
# defaults 3 and 2). FDR = Benjamini-Hochberg over the tested genes of the group.
#
# Output per group (out_dir/<group>/):
#   master_table.tsv          input effect sizes (long format)
#   meta_results.tsv          all genes: estimate, se, z, p, CI, tau2_study, tau2_contrast, QE, FDR, ...
#   meta_sig.tsv              FDR < FDR_CUT and |estimate| >= LFC_CUT (defaults 0.05 and 1)
#   genes_UP.txt, genes_DOWN.txt, genes_tested.txt   (gene lists for step 15)

suppressPackageStartupMessages({ library(metafor); library(parallel) })

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: meta_multilevel.R <contrasts.tsv> <out_dir> [groups.tsv] [threads]")
IDX <- args[1]; OUT <- args[2]
GROUPS <- if (length(args) >= 3 && nzchar(args[3])) args[3] else NA
THREADS <- if (length(args) >= 4) as.integer(args[4]) else 1L
env_num <- function(x, d) { v <- Sys.getenv(x, ""); if (nzchar(v)) as.numeric(v) else d }
MIN_CONTRASTS <- env_num("MIN_CONTRASTS", 3); MIN_STUDIES <- env_num("MIN_STUDIES", 2)
FDR_CUT <- env_num("FDR_CUT", 0.05); LFC_CUT <- env_num("LFC_CUT", 1)

message("[INFO] metafor ", packageVersion("metafor"), " | min contrasts ", MIN_CONTRASTS,
        ", min studies ", MIN_STUDIES, ", FDR < ", FDR_CUT, ", |estimate| >= ", LFC_CUT, ", threads ", THREADS)

idx <- read.delim(IDX, stringsAsFactors = FALSE)
stopifnot(all(c("study", "contrast", "file") %in% names(idx)))
if (anyDuplicated(idx[, c("study", "contrast")])) stop("duplicated study/contrast pairs in ", IDX)

if (is.na(GROUPS)) {
  idx$group <- "all"
} else {
  g <- read.delim(GROUPS, stringsAsFactors = FALSE)
  stopifnot(all(c("study", "group") %in% names(g)))
  miss <- setdiff(idx$study, g$study)
  if (length(miss)) message("[WARN] studies without a group, excluded: ", paste(miss, collapse = ", "))
  idx <- merge(idx, g[, c("study", "group")], by = "study")
}

fit_gene <- function(d) {
  base <- data.frame(n_contrasts = nrow(d), n_studies = length(unique(d$study)),
                     estimate = NA_real_, se = NA_real_, zval = NA_real_, pval = NA_real_,
                     ci_lb = NA_real_, ci_ub = NA_real_, tau2_study = NA_real_, tau2_contrast = NA_real_,
                     QE = NA_real_, QEp = NA_real_, converged = FALSE)
  if (nrow(d) < MIN_CONTRASTS || base$n_studies < MIN_STUDIES) return(base)
  fit <- tryCatch(rma.mv(yi = yi, V = vi, random = ~ 1 | study / contrast, method = "REML", data = d),
                  error = function(e) NULL)
  if (is.null(fit)) return(base)
  base$estimate <- as.numeric(fit$b[1]); base$se <- fit$se; base$zval <- fit$zval; base$pval <- fit$pval
  base$ci_lb <- fit$ci.lb; base$ci_ub <- fit$ci.ub
  base$tau2_study <- fit$sigma2[1]; base$tau2_contrast <- if (length(fit$sigma2) >= 2) fit$sigma2[2] else NA
  base$QE <- fit$QE; base$QEp <- fit$QEp; base$converged <- TRUE
  base
}

for (grp in sort(unique(idx$group))) {
  sub <- idx[idx$group == grp, ]
  gdir <- file.path(OUT, gsub("[^A-Za-z0-9_.-]+", "_", grp)); dir.create(gdir, recursive = TRUE, showWarnings = FALSE)
  message("[INFO] group '", grp, "': ", length(unique(sub$study)), " studies, ", nrow(sub), " contrasts")

  master <- do.call(rbind, lapply(seq_len(nrow(sub)), function(i) {
    x <- read.csv(sub$file[i], stringsAsFactors = FALSE)
    x <- x[, c("gene_id", "log2FoldChange", "lfcSE", "padj")]
    x$study <- sub$study[i]; x$contrast <- sub$contrast[i]; x
  }))
  master <- master[!is.na(master$gene_id) & is.finite(master$log2FoldChange) &
                     is.finite(master$lfcSE) & master$lfcSE > 0, ]
  write.table(master, file.path(gdir, "master_table.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  master$yi <- master$log2FoldChange; master$vi <- master$lfcSE^2
  master$study <- factor(master$study)
  master$contrast <- factor(paste(master$study, master$contrast, sep = "__"))
  by_gene <- split(master, master$gene_id)
  message("[INFO] fitting ", length(by_gene), " genes")
  res <- mclapply(by_gene, fit_gene, mc.cores = THREADS)
  res <- do.call(rbind, res); res <- cbind(gene_id = names(by_gene), res)

  res$FDR <- p.adjust(res$pval, method = "BH")
  res$direction <- ifelse(is.na(res$estimate), NA, ifelse(res$estimate > 0, "UP", ifelse(res$estimate < 0, "DOWN", "ZERO")))
  res <- res[order(res$FDR, res$pval, na.last = TRUE), ]
  write.table(res, file.path(gdir, "meta_results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  sig <- res[!is.na(res$FDR) & res$FDR < FDR_CUT & abs(res$estimate) >= LFC_CUT, ]
  write.table(sig, file.path(gdir, "meta_sig.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  writeLines(sig$gene_id[sig$estimate > 0], file.path(gdir, "genes_UP.txt"))
  writeLines(sig$gene_id[sig$estimate < 0], file.path(gdir, "genes_DOWN.txt"))
  writeLines(res$gene_id[res$converged & !is.na(res$pval)], file.path(gdir, "genes_tested.txt"))
  message("[OK] group '", grp, "': ", sum(res$converged), " genes tested, ", nrow(sig), " significant (",
          sum(sig$estimate > 0), " up, ", sum(sig$estimate < 0), " down)")
}
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
