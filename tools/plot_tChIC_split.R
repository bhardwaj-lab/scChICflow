#!/usr/bin/env Rscript

library(magrittr)
library(ggplot2)

Args <- commandArgs(trailingOnly = T)
outfile <- Args[1]
indir <- Args[2]

f <- list.files(indir, pattern = "_split_tChIC.log", full.names=TRUE)

lapply(f, function(x) read.delim(x, sep=" ", header=F, col.names = c("sample", "count"))) %>%
  Reduce(rbind, .) %>% dplyr::mutate(sample=gsub(":", "", sample)) %>%
  tidyr::separate(sample, c("sample", "protocol"), "_") %>%
  ggplot(., aes(sample, count, fill=protocol)) +
  geom_bar(position="dodge", stat = "identity") +
  scale_y_continuous(breaks = seq(0, 100000000, 5000000)) +
  theme_gray(base_size = 14) + coord_flip() -> p

ggsave(plot = p, file=outfile, width=10, height = 6)
