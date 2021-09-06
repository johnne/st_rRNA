#!/usr/bin/env Rscript

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

seqs <- snakemake@input$seqs
taxdf <- snakemake@output$taxdf
bootdf <- snakemake@output$bootdf
training_fasta <- snakemake@input$training_fasta
minBoot <- snakemake@params$minboot
outputBootstraps <- snakemake@params$outputBootstraps
threads <- snakemake@threads

library(dada2)
set.seed(100)
refFasta <- training_fasta
taxonomy <- assignTaxonomy(seqs=seqs, refFasta=refFasta, minBoot=minBoot,
                           outputBootstraps=outputBootstraps,
                           multithread=threads, verbose=TRUE)
df <- as.data.frame(taxonomy)
rownames(df) <- names(dimnames(taxonomy$tax)[[1]])
taxdf <- df[, grepl("tax.", colnames(df))]
bootdf <- df[, grepl("boot", colnames(df))]
write.table(taxdf, taxdf, sep="\t", quote=FALSE)
write.table(bootdf, bootdf, sep="\t", quote=FALSE)
print("DONE")