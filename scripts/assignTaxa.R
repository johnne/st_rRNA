#!/usr/bin/env Rscript

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

seqs <- snakemake@input$seqs
taxdf <- snakemake@output$taxdf
bootdf <- snakemake@output$bootdf
refFasta <- snakemake@input$refFasta
minBoot <- snakemake@params$minboot
outputBootstraps <- snakemake@params$outputBootstraps
threads <- snakemake@threads

library(dada2)
set.seed(100)
taxonomy <- assignTaxonomy(seqs=seqs, refFasta=refFasta, minBoot=minBoot,
                           outputBootstraps=outputBootstraps,
                           multithread=threads, verbose=TRUE)
df <- as.data.frame(taxonomy)
rownames(df) <- names(dimnames(taxonomy$tax)[[1]])
taxdf <- df[, grepl("tax.", colnames(df))]
colnames(taxdf) <- gsub("tax.", "", colnames(taxdf))
bootdf <- df[, grepl("boot", colnames(df))]
colnames(bootdf) <- gsub("boot.", "", colnames(bootdf))
write.table(taxdf, taxdf, sep="\t", quote=FALSE)
write.table(bootdf, bootdf, sep="\t", quote=FALSE)
print("DONE")