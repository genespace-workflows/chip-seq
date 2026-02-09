#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(DESeq2)
})

opt_list <- list(
  make_option("--samples", type="character", help="samples.tsv with columns: sample, group, counts"),
  make_option("--out_results", type="character", default="deseq2_results.tsv"),
  make_option("--out_norm", type="character", default="normalized_counts.tsv")
)
opt <- parse_args(OptionParser(option_list = opt_list))

meta <- fread(opt$samples)
stopifnot(all(c("sample","group","counts") %in% colnames(meta)))

# Read counts: each file has columns: chr start end peak_id count
# We'll use peak_id as row index (NR from counting step).
count_list <- lapply(meta$counts, function(fn){
  dt <- fread(fn, header=FALSE)
  setnames(dt, c("chr","start","end","peak_id","count"))
  dt[, .(peak_id, count)]
})

# Merge into matrix by peak_id
mat <- Reduce(function(x,y) merge(x,y,by="peak_id",all=TRUE), count_list)
mat[is.na(mat)] <- 0

# Assign column names as samples
for(i in seq_along(meta$sample)){
  setnames(mat, old = names(mat)[i+1], new = meta$sample[i])
}

count_mat <- as.matrix(mat[, -1, with=FALSE])
rownames(count_mat) <- mat$peak_id

coldata <- data.frame(row.names = meta$sample, group = factor(meta$group))

dds <- DESeqDataSetFromMatrix(countData = count_mat, colData = coldata, design = ~ group)
dds <- DESeq(dds)

# Default contrast: second level vs first level
# Ensure deterministic: reorder levels if you want: coldata$group <- relevel(coldata$group, "control")
res <- results(dds)

res_dt <- as.data.table(res, keep.rownames = "peak_id")
res_dt <- res_dt[order(pvalue)]
fwrite(res_dt, opt$out_results, sep="\t")

norm_counts <- counts(dds, normalized=TRUE)
norm_dt <- as.data.table(norm_counts, keep.rownames = "peak_id")
fwrite(norm_dt, opt$out_norm, sep="\t")
