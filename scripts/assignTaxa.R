#!/usr/bin/env Rscript

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

seqs <- snakemake@input$seqs
taxdf_file <- snakemake@output$taxdf
bootdf_file <- snakemake@output$bootdf
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
write.table(taxdf, taxdf_file, sep="\t", quote=FALSE)
write.table(bootdf, bootdf_file, sep="\t", quote=FALSE)
print("DONE")