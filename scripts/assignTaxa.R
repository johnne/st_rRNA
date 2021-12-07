#!/usr/bin/env Rscript

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

subunit <- snakemake@wildcards$subunit
seqs <- snakemake@input$seqs
taxdf_file <- snakemake@output$taxdf
bootdf_file <- snakemake@output$bootdf
refFasta <- snakemake@input$refFasta
spFasta <- snakemake@input$spFasta
minBoot <- snakemake@params$minBoot
taxLevels <- snakemake@params$taxLevels
tryRC <- snakemake@params$tryRC
outputBootstraps <- snakemake@params$outputBootstraps
threads <- snakemake@threads

library(dada2)
library(ShortRead)
reads <- readFastq(seqs)
if (length(reads) == 0) {
    writeLines("", taxdf_file)
    writeLines("", bootdf_file)
    quit(save = "no", status = 0)
}
set.seed(100)
taxonomy <- assignTaxonomy(seqs=seqs, refFasta=refFasta, minBoot=minBoot, taxLevels=taxLevels,
                           outputBootstraps=outputBootstraps, tryRC=tryRC,
                           multithread=threads, verbose=TRUE)
df <- as.data.frame(taxonomy)
# Remove sequences with numbered suffix (resulting from duplicates)
df <- df[grep("[0-9]", rownames(df), invert=TRUE), ]
sequences <- as.character(dimnames(taxonomy$tax)[[1]])
seqids <- names(dimnames(taxonomy$tax)[[1]])
seqidmap <- as.data.frame(seqids, row.names=sequences)
# Add species
df <- addSpecies(df, spFasta, allowMultiple=TRUE)
colnames(df) <- gsub("Species", "tax.species", colnames(df))
rownames(df) <- seqidmap[rownames(df),]
taxdf <- df[, grepl("tax.", colnames(df))]
colnames(taxdf) <- gsub("tax.", "", colnames(taxdf))
rownames(taxdf) <- gsub(" .+", "", rownames(taxdf))
bootdf <- df[, grepl("boot", colnames(df))]
colnames(bootdf) <- gsub("boot.", "", colnames(bootdf))
rownames(bootdf) <- gsub(" .+", "", rownames(bootdf))
write.table(taxdf, taxdf_file, sep="\t", quote=FALSE)
write.table(bootdf, bootdf_file, sep="\t", quote=FALSE)
print("DONE")