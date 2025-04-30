#!/usr/bin/env Rscript
library(data.table)
library(Matrix)

# Retrieve command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: matrix_to_counts.R <bam_dir> <out_file>")
}
bamdir <- args[1]
out_file <- args[2]

# Build the path to the STARsolo output matrix folder.
matrix_dir <- file.path(bamdir, "Solo.out", "Gene", "raw", "")

# Read in the sparse matrix
f <- file(paste0(matrix_dir, "umiDedup-1MM_All.mtx"), "r")
mat <- as.data.frame(as.matrix(readMM(f)))
close(f)

# Read in feature and barcode files
feature.names <- fread(paste0(matrix_dir, "features.tsv"), header = FALSE)
barcode.names <- fread(paste0(matrix_dir, "barcodes.tsv"), header = FALSE, stringsAsFactors = FALSE, data.table = FALSE)

# Assign row and column names to the matrix
colnames(mat) <- barcode.names$V1
rownames(mat) <- feature.names$V1

# Write the count matrix to a file
fwrite(mat, file = out_file, sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

