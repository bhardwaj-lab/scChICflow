#!/usr/bin/env Rscript

library(magrittr)
library(ggplot2)

Args <- commandArgs(trailingOnly = T)
bcfile <- Args[1]
countdir <- Args[2]
outpdf <- Args[3]

idx <- read.delim(bcfile, header = FALSE)
files <- list.files(countdir, pattern = "*\\.per_barcode.tsv", full.names = TRUE)
names <- list.files(countdir, pattern = "*\\.per_barcode.tsv") %>% gsub(".per_barcode.tsv", "", .)
ncounts <- lapply(files, read.delim, header = FALSE, sep = "\t")
names(ncounts) <- names

ncounts %<>% lapply(function(x) {
  x$idx <- match(x$V1, idx$V1)
  colnames(x) <- c("barcode", "total_count", "idx")
  return(x) })

pdf(outpdf, width = 8, height = 5)
mapply(function(x, y) {
  g <- expand.grid(x = LETTERS[1:16], y = 1:24)
  g$n <- rep(0:15, 24)
  g$idx <- 24*g$n + g$y
  g$count <- x[match(g$idx, x$idx), "total_count"]
  g$y %<>% as.factor()
  g$x %<>% as.factor()
  g$x <- with(g, factor(x, levels = rev(levels(x))))
  ggplot(g, aes(y, x, fill = log10(count))) + 
    geom_tile() + ggtitle(label = y) +
    scale_fill_viridis_c() + 
    theme_bw(base_size = 16)
  
}, ncounts, names(ncounts), SIMPLIFY = FALSE)
dev.off()