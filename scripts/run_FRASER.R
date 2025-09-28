#!/usr/bin/env Rscript
## Parse args 
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript run_FRASER.R <sample_table.tsv> <out_dir> [workers=16] [q=2] [filter=FALSE] [minExp=5]")
}
sample_table <- args[[1]]
out_dir      <- args[[2]]
workers      <- if (length(args) >= 3) as.integer(args[[3]]) else 16L
q_dim        <- if (length(args) >= 4) as.integer(args[[4]]) else 2L
do_filter    <- if (length(args) >= 5) as.logical(args[[5]]) else FALSE
min_exp      <- if (length(args) >= 6) as.integer(args[[6]]) else 5L

## Require packages
needed <- c(
  "FRASER", "BiocParallel", "data.table", "S4Vectors",
  "TxDb.Hsapiens.UCSC.hg38.knownGene", "org.Hs.eg.db", "GenomicRanges"
)
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Missing packages (please install them yourself, script will not install): ",
       paste(missing, collapse=", "))
}
library(FRASER)
library(BiocParallel)
library(data.table)
library(S4Vectors)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(GenomicRanges)

## --- I/O & parallel ---
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
register(MulticoreParam(workers = workers))
message("Workers: ", workers, " | q(jaccard): ", q_dim,
        " | filter: ", do_filter, " | minExpressionInOneSample: ", min_exp)

## Read sample table
tbl <- data.table::fread(sample_table)

req <- c("sampleID","bamFile","pairedEnd","strandSpecific")
if (!all(req %in% names(tbl))) {
  stop("Your sample table must contain columns: ", paste(req, collapse=", "))
}

## Types & paths
tbl[, bamFile := normalizePath(bamFile, mustWork = TRUE)]
tbl[, pairedEnd := as.logical(pairedEnd)]
tbl[, strandSpecific := as.integer(strandSpecific)]   # 0 (unstranded), 1 (forward), 2 (reverse)

## Ensure unique sampleID for FRASER
if (any(duplicated(tbl$sampleID))) {
  if (!"replicate" %in% names(tbl)) {
    tbl[, replicate := paste0("R", seq_len(.N)), by = sampleID]
  }
  tbl[, sampleID := make.unique(paste0(sampleID, "_", replicate))]
  warning("Duplicate sampleIDs found; created unique internal IDs by appending replicate.")
}

## Warn for missing BAM indexes
has_bai <- function(b) file.exists(paste0(b, ".bai")) || file.exists(sub("\\.bam$", ".bai", b))
no_idx <- tbl[!vapply(bamFile, has_bai, logical(1))]
if (nrow(no_idx)) warning("Missing .bai for:\n  ", paste(no_idx$bamFile, collapse="\n  "))

## Keep optional columns if present
cols_for_fraser <- c("sampleID","bamFile","pairedEnd","strandSpecific",
                     intersect(c("group","replicate"), names(tbl)))

## colData as S4Vectors::DataFrame
cd <- S4Vectors::DataFrame(tbl[, ..cols_for_fraser])
data.table::fwrite(as.data.table(cd), file.path(out_dir, "colData_used.tsv"), sep = "\t")

## Build dataset
fds <- FraserDataSet(colData = cd, workingDir = out_dir)
name(fds) <- "FRASER_run"

## Counting from BAMs
fds <- countRNAData(
  fds,
  NcpuPerSample = 10,  # total parallelism via BiocParallel
  filter = do_filter
)

# compute stats
fds <- calculatePSIValues(fds)

## Annotate junction ranges (hg38 knownGene + org.Hs.eg.db)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
orgDb <- org.Hs.eg.db
fds  <- annotateRangesWithTxDb(fds, txdb = txdb, orgDb = orgDb)

## Save annotated junction ranges to TSV (coords + mcols)
rr <- rowRanges(fds, type = "j")
rr_dt <- data.table(
  seqnames = as.character(GenomicRanges::seqnames(rr)),
  start    = GenomicRanges::start(rr),
  end      = GenomicRanges::end(rr),
  strand   = as.character(GenomicRanges::strand(rr))
)
mcols_dt <- as.data.table(S4Vectors::mcols(rr))
data.table::fwrite(cbind(rr_dt, mcols_dt),
                   file.path(out_dir, "junction_rowRanges_annotated.tsv.gz"),
                   sep = "\t")

## Fit FRASER (Jaccard)
fds <- FRASER(fds)

## Results (junction × sample)
res <- results(fds, all=TRUE, padjCutoff = 1, deltaPsiCutoff = 0)
data.table::fwrite(as.data.table(res),
                   file.path(out_dir, "fraser_results_jaccard.tsv.gz"),
                   sep = "\t")

## Save fitted object
dir.create(file.path(out_dir, "savedObjects"), showWarnings = FALSE, recursive = TRUE)
saveRDS(fds, file.path(out_dir, "savedObjects", "fds_fitted.rds"))

message("Done.\n  Results: ", file.path(out_dir, "fraser_results_jaccard.tsv.gz"),
        "\n  ColData: ", file.path(out_dir, "colData_used.tsv"),
        "\n  Junctions: ", file.path(out_dir, "junction_rowRanges_annotated.tsv.gz"),
        if ("geneSymbol" %in% colnames(res))
          paste0("\n  Gene counts: ", file.path(out_dir, "per_sample_gene_outlier_counts.tsv"))
        else "")