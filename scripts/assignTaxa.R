#!/usr/bin/env Rscript

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

seqs <- snakemake@input$seqs
output <- snakemake@output[[1]]
training_fasta <- snakemake@input$training_fasta
minBoot <- snakemake@params$minboot
outputBootstraps <- snakemake@params$outputBootstraps
threads <- snakemake@threads

library(dada2)
set.seed(100)
refFasta <- training_fasta #system.file("extdata", training_fasta, package="dada2")
taxonomy <- assignTaxonomy(seqs=seqs, refFasta=refFasta, minBoot=minBoot,
                           outputBootstraps=outputBootstraps,
                           multithread=threads, verbose=TRUE)
rownames(taxonomy) <- names(rownames(taxonomy))
write.table(taxonomy, output, sep="\t", quote=FALSE)
