#!/usr/bin/env Rscript
# GO over-representation analysis (ORA) of gene sets against the genes actually tested.
#
# Usage: Rscript go_enrichment.R <interproscan.tsv> <tx2gene.tsv> <sets.tsv> <out_dir>
#
# interproscan.tsv : InterProScan TSV output (step 12), run with -goterms.
# tx2gene.tsv      : transcript_id <TAB> gene_id (step 07).
# sets.tsv         : set <TAB> genes_file <TAB> universe_file   (header required)
#                    genes_file    = one gene_id per line (the set to test)
#                    universe_file = one gene_id per line (background: the genes that were tested
#                                    in the analysis that produced the set)
#
# Annotation: GO terms of each peptide are transferred to its transcript (TransDecoder suffix
# ".pN" removed) and then to its gene; the gene-to-GO table is the union over the gene's peptides.
# Alternatively, set GENE2GO to a precomputed gene <TAB> GO table (no header) to use it directly.
# Test: clusterProfiler::enricher per ontology (BP, MF, CC), restricted to annotated genes of the
# universe; BH adjustment; terms with 5-500 annotated universe genes; sets with < 5 annotated genes
# are not tested. Thresholds can be changed with MIN_GS, MAX_GS, P_CUT, Q_CUT (env vars).
#
# Output: gene2go.tsv, enrichment_summary.tsv (one row per set x ontology), enrichment_terms.tsv
#         (significant terms), dotplots/<set>_<ONT>.png (top 20 terms)

suppressPackageStartupMessages({
  library(clusterProfiler); library(GO.db); library(AnnotationDbi); library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) stop("usage: go_enrichment.R <interproscan.tsv> <tx2gene.tsv> <sets.tsv> <out_dir>")
IPS <- args[1]; TX2GENE <- args[2]; SETS <- args[3]; OUT <- args[4]
env_num <- function(x, d) { v <- Sys.getenv(x, ""); if (nzchar(v)) as.numeric(v) else d }
MIN_GS <- env_num("MIN_GS", 5); MAX_GS <- env_num("MAX_GS", 500)
P_CUT <- env_num("P_CUT", 0.05); Q_CUT <- env_num("Q_CUT", 0.2)
dir.create(file.path(OUT, "dotplots"), recursive = TRUE, showWarnings = FALSE)
message("[INFO] clusterProfiler ", packageVersion("clusterProfiler"), " | GO.db ", packageVersion("GO.db"))

# ---------------- gene -> GO
GENE2GO <- Sys.getenv("GENE2GO", "")
if (nzchar(GENE2GO)) {
  # Precomputed gene <TAB> GO table (no header); InterProScan and tx2gene are then not used.
  g2g <- read.delim(GENE2GO, header = FALSE, stringsAsFactors = FALSE)[, 1:2]
  names(g2g) <- c("gene_id", "GO")
  g2g <- unique(g2g[grepl("^GO:[0-9]{7}$", g2g$GO), ])
  message("[INFO] gene2go read from ", GENE2GO)
} else {
  ips <- read.delim(IPS, header = FALSE, sep = "\t", quote = "", fill = TRUE, col.names = paste0("V", 1:15),
                    stringsAsFactors = FALSE)
  ips <- ips[!is.na(ips$V14) & ips$V14 != "" & ips$V14 != "-", c("V1", "V14")]
  go_list <- strsplit(ips$V14, "|", fixed = TRUE)
  p2g <- data.frame(protein = rep(ips$V1, lengths(go_list)), GO = sub("\\(.*\\)$", "", unlist(go_list)))
  p2g <- p2g[grepl("^GO:[0-9]{7}$", p2g$GO), ]
  p2g$transcript <- sub("\\.p[0-9]+$", "", p2g$protein)
  tx <- read.delim(TX2GENE, header = FALSE, stringsAsFactors = FALSE)[, 1:2]; names(tx) <- c("transcript", "gene_id")
  g2g <- unique(merge(p2g, tx, by = "transcript")[, c("gene_id", "GO")])
}
if (!nrow(g2g)) stop("no gene-GO pair: check the annotation inputs (InterProScan/tx2gene or GENE2GO)")
write.table(g2g, file.path(OUT, "gene2go.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
message("[INFO] gene2go: ", nrow(g2g), " pairs, ", length(unique(g2g$gene_id)), " annotated genes")

gometa <- suppressMessages(AnnotationDbi::select(GO.db, keys = unique(g2g$GO), columns = c("TERM", "ONTOLOGY"),
                                                 keytype = "GOID"))
gometa <- gometa[!is.na(gometa$ONTOLOGY), ]; names(gometa)[1] <- "GO"

# ---------------- ORA
read_ids <- function(f) { x <- trimws(readLines(f, warn = FALSE)); unique(x[x != ""]) }
sets <- read.delim(SETS, stringsAsFactors = FALSE)
stopifnot(all(c("set", "genes_file", "universe_file") %in% names(sets)))

summary_rows <- list(); term_rows <- list()
for (i in seq_len(nrow(sets))) {
  genes <- read_ids(sets$genes_file[i]); univ <- read_ids(sets$universe_file[i])
  for (ont in c("BP", "MF", "CC")) {
    om <- gometa[gometa$ONTOLOGY == ont, ]
    t2g <- g2g[g2g$GO %in% om$GO, c("GO", "gene_id")]
    annot_u <- intersect(univ, t2g$gene_id); annot_g <- intersect(genes, annot_u)
    row <- data.frame(set = sets$set[i], ontology = ont, input = length(genes),
                      outside_universe = sum(!genes %in% univ), input_annotated = length(annot_g),
                      universe = length(univ), universe_annotated = length(annot_u), significant_terms = NA,
                      status = "")
    if (length(annot_g) < 5) {
      row$status <- "not tested (< 5 annotated genes)"
    } else {
      ego <- enricher(annot_g, universe = annot_u, TERM2GENE = t2g, TERM2NAME = unique(om[, c("GO", "TERM")]),
                      pAdjustMethod = "BH", pvalueCutoff = P_CUT, qvalueCutoff = Q_CUT,
                      minGSSize = MIN_GS, maxGSSize = MAX_GS)
      r <- if (is.null(ego)) data.frame() else as.data.frame(ego)
      row$significant_terms <- nrow(r); row$status <- "tested"
      if (nrow(r)) {
        term_rows[[length(term_rows) + 1]] <- cbind(set = sets$set[i], ontology = ont, r)
        top <- head(r[order(r$p.adjust), ], 20)
        top$ratio <- sapply(strsplit(top$GeneRatio, "/"), function(z) as.numeric(z[1]) / as.numeric(z[2]))
        p <- ggplot(top, aes(ratio, reorder(Description, ratio), size = Count, colour = p.adjust)) +
          geom_point() + scale_colour_continuous(trans = "reverse") + theme_bw() +
          labs(x = "Gene ratio", y = NULL, title = paste(sets$set[i], ont))
        ggsave(file.path(OUT, "dotplots", paste0(gsub("[^A-Za-z0-9_.-]+", "_", sets$set[i]), "_", ont, ".png")),
               p, width = 8, height = 2 + 0.3 * nrow(top), dpi = 200, limitsize = FALSE)
      }
    }
    summary_rows[[length(summary_rows) + 1]] <- row
    message("[OK] ", sets$set[i], " ", ont, ": ", row$status,
            if (!is.na(row$significant_terms)) paste0(", ", row$significant_terms, " significant terms") else "")
  }
}
write.table(do.call(rbind, summary_rows), file.path(OUT, "enrichment_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
terms <- if (length(term_rows)) do.call(rbind, term_rows) else data.frame()
write.table(terms, file.path(OUT, "enrichment_terms.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
