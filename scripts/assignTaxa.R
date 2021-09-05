#!/usr/bin/env Rscript

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

input <- snakemake@input[[1]]
output <- snakemake@output[[1]]
refFasta <- snakemake@input[[2]]
minBoot <- snakemake@params$minboot
outputBootstraps <- snakemake@params$outputBootstraps
threads <- snakemake@threads

library(dada2)
set.seed(100)
training_fasta <- system.file("extdata", refFasta, package="dada2")
taxonomy <- assignTaxonomy(seqs, training_fasta, minBoot=minBoot,
                           outputBootstraps=outputBootstraps,
                           multithread=threads, verbose=TRUE)
rownames(taxonomy) <- names(rownames(taxonomy))
write.table(taxonomy, output, sep="\t", quote=FALSE)