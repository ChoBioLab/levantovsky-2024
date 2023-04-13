#!/usr/bin/env Rscript


library(Seurat)
library(SeuratObject)
library(rtracklayer)
library(GenomeInfoDb)
library(Signac)



# read in the integrated object
int_obj = readRDS("/mnt/cho_lab/disk2/jiayuzh/projects/perianal-cd/data/2023-01-19_MultiomeReRun_Annotated.rds")
meta_dt = int_obj@meta.data


DefaultAssay(int_obj) = "ATAC"

# call peaks for the subsetted objects
# group by parameter

variables = levels(as.factor(meta_dt$Bin))

for (i in 1:length(variables)){ 
  peaks <- CallPeaks(int_obj, assay = "ATAC", group.by = "Bin", idents = variables[i]) # more accurate peak calling using MACS2
  peaks <- keepStandardChromosomes(
    peaks,
    pruning.mode = "coarse"
  ) # remove peaks on nonstandard chromosomes and in genomic blacklist regions
  peaks <- subsetByOverlaps(
    x = peaks,
    ranges = blacklist_hg38_unified,
    invert = TRUE
  )
rtracklayer::export.bed(peaks, paste0(variables[i], "_peaks.bed"))
}
