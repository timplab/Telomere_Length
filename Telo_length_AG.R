#!/usr/bin/env Rscript

# laod libraries
if(!require(optparse)){
  install.packages("optparse")
  library(optparse)
}

# load functions
source("utils.R")


option_list = list(
  make_option(c("-c", "--coverage"), type="numeric", default=NULL, 
              help="minimum coverage value (default is 20X)"),
  make_option(c("-o", "--out"), type="character", default=NULL, 
              help="Output file"), 
  make_option(c("-s", "--starts"), type="character", default=NULL, 
              help="path to telomere starts file"),
  make_option(c("-i", "--data"), type="character", default=NULL, 
              help="path to telomere data bed file")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


if (is.null(opt$starts)){
  print_help(opt_parser)
  stop("provide path to telomere starts file.\n", call.=FALSE)
}

if (is.null(opt$out)){
  print_help(opt_parser)
  stop("Specify an output directory.\n", call.=FALSE)
}
if (is.null(opt$data)){
  print_help(opt_parser)
  stop("provide path to bed file.\n", call.=FALSE)
}

if (is.null(opt$coverage)){
  opt$coverage = 20
}
# load data
telo_starts <- read_delim(opt$starts, delim= "\t", col_names =  F) %>%
  dplyr::rename("TEL" = X1, "telo_start" = X2) %>%
  mutate(telo_start = as.numeric(telo_start)) %>%
  dplyr::select(c(TEL, telo_start))

#load data (bedfile) 
pulse <- read_delim(opt$data, delim= "\t", col_names =  F) %>%
  dplyr::rename("chrom" = X1, "start" = X2, "end" = X3, "qname" = X4, "score" = X5, "strand" = X6) %>%
  mutate(score = as.numeric(score)) %>%
  mutate(start = as.numeric(start)) %>%
  mutate(end = as.numeric(end)) %>%
  dplyr::select(c(chrom, start, end, qname, score, strand))
pulse <- pulse[!duplicated(pulse$qname), ]

# run lengths function
dat.out <- get_lengths(opt$coverage)

# save data 

write.table(dat.out,
            file = opt$out, 
            sep = "\t", 
            quote = F, 
            row.names = F,
            col.names = T)
