library(Seurat)
library(dplyr)

# read in the integrated object
int_obj = readRDS("/mnt/cho_lab/disk2/jiayuzh/projects/perianal-cd/data/2023-01-19_MultiomeReRun_Annotated.rds")
sample_info = read.csv(file = "/mnt/cho_lab/disk2/ctastad/projects/perianal-cd/analysis/coreSC/output/output_2023-01-16_12.27.51/samples.csv")
meta_dt = int_obj@meta.data

# extract sample_id 
samples = levels(as.factor(meta_dt$object))

# set the variable of interests

meta_dt$combined_bin = meta_dt$Bin
meta_dt$combined_bin[which(meta_dt$combined_bin == "Fibroblasts")] = "Stromal"
meta_dt$combined_bin = as.factor(meta_dt$combined_bin)
variables = levels(as.factor(meta_dt$Bin)) # example with Fibroblasts vs. Myeloid
path_out = "/mnt/cho_lab/disk1/jiayuzh/tmp/"



for(j in 1:length(variables)){
    # make a dataframe to be input as pseudobulking commands
    df = data.frame(matrix(ncol = 6, nrow = 0))
    colnames(df) = c("file", "bam", "sam_body", "sam_header", "filtered_sam", "filtered_bam")
    non_df = data.frame(matrix(ncol = 6, nrow = 0))
    colnames(non_df) = c("file", "bam", "sam_body", "sam_header", "filtered_sam", "filtered_bam")
    non_var = paste0("non", variables[j])

    for (i in 1:length(samples)){
        barcodes = meta_dt[which(meta_dt$object == samples[i] & meta_dt$Bin == variables[j]), "gex_barcode"]
        barcodes = paste("CB:Z:",barcodes, sep = "")
        file_name = paste0(path_out, paste(samples[i], variables[j], sep = "_"), ".txt")
        write.table(barcodes, file = file_name, quote = F, row.names = F, col.names = F)
        bam_loc = paste(sample_info[sample_info$name == samples[i], "dir"], "atac_possorted_bam.bam", sep = "/")
        sam_body = paste0(path_out, samples[i], "_", variables[j], "_filtered_SAM_body")
        sam_header = paste0(path_out, samples[i], "_", variables[j], "_SAM_header")
        filtered_sam = paste0(path_out, samples[i], "_", variables[j], "_filtered.sam")
        filtered_bam = paste0(path_out, samples[i], "_", variables[j], "_filtered.bam")
        row = c(file_name, bam_loc, sam_body, sam_header, filtered_sam, filtered_bam)
        df[nrow(df)+1,] = row

        barcodes = meta_dt[which(meta_dt$object == samples[i] & meta_dt$Bin != variables[j]), "gex_barcode"]
        barcodes = paste("CB:Z:",barcodes, sep = "")
        file_name = paste0(path_out, paste(samples[i], non_var, sep = "_"), ".txt")
        write.table(barcodes, file = file_name, quote = F, row.names = F, col.names = F)
        bam_loc = paste(sample_info[sample_info$name == samples[i], "dir"], "atac_possorted_bam.bam", sep = "/")
        sam_body = paste0(path_out, samples[i], "_", non_var, "_filtered_SAM_body")
        sam_header = paste0(path_out, samples[i], "_", non_var, "_SAM_header")
        filtered_sam = paste0(path_out, samples[i], "_", non_var, "_filtered.sam")
        filtered_bam = paste0(path_out, samples[i], "_", non_var, "_filtered.bam")
        row = c(file_name, bam_loc, sam_body, sam_header, filtered_sam, filtered_bam)
        non_df[nrow(non_df)+1,] = row
    }
    write.table(df, file = paste0(path_out, variables[j], "_bam_loc.txt"), quote = F, row.names = F, col.names = F)
    write.table(non_df, file = paste0(path_out, non_var, "_bam_loc.txt"), quote = F, row.names = F, col.names = F)
}

# barcodes for Stromal separately

    df = data.frame(matrix(ncol = 6, nrow = 0))
    colnames(df) = c("file", "bam", "sam_body", "sam_header", "filtered_sam", "filtered_bam")
    non_df = data.frame(matrix(ncol = 6, nrow = 0))
    colnames(non_df) = c("file", "bam", "sam_body", "sam_header", "filtered_sam", "filtered_bam")
    non_var = paste0("non", "Stromal")

    for (i in 1:length(samples)){
        barcodes = meta_dt[which(meta_dt$object == samples[i] & meta_dt$combined_bin == "Stromal"), "gex_barcode"]
        barcodes = paste("CB:Z:",barcodes, sep = "")
        file_name = paste0(path_out, paste(samples[i], "Stromal", sep = "_"), ".txt")
        write.table(barcodes, file = file_name, quote = F, row.names = F, col.names = F)
        bam_loc = paste(sample_info[sample_info$name == samples[i], "dir"], "atac_possorted_bam.bam", sep = "/")
        sam_body = paste0(path_out, samples[i], "_", "Stromal", "_filtered_SAM_body")
        sam_header = paste0(path_out, samples[i], "_", "Stromal", "_SAM_header")
        filtered_sam = paste0(path_out, samples[i], "_", "Stromal", "_filtered.sam")
        filtered_bam = paste0(path_out, samples[i], "_", "Stromal", "_filtered.bam")
        row = c(file_name, bam_loc, sam_body, sam_header, filtered_sam, filtered_bam)
        df[nrow(df)+1,] = row

        barcodes = meta_dt[which(meta_dt$object == samples[i] & meta_dt$combined_bin != "Stromal"), "gex_barcode"]
        barcodes = paste("CB:Z:",barcodes, sep = "")
        file_name = paste0(path_out, paste(samples[i], non_var, sep = "_"), ".txt")
        write.table(barcodes, file = file_name, quote = F, row.names = F, col.names = F)
        bam_loc = paste(sample_info[sample_info$name == samples[i], "dir"], "atac_possorted_bam.bam", sep = "/")
        sam_body = paste0(path_out, samples[i], "_", non_var, "_filtered_SAM_body")
        sam_header = paste0(path_out, samples[i], "_", non_var, "_SAM_header")
        filtered_sam = paste0(path_out, samples[i], "_", non_var, "_filtered.sam")
        filtered_bam = paste0(path_out, samples[i], "_", non_var, "_filtered.bam")
        row = c(file_name, bam_loc, sam_body, sam_header, filtered_sam, filtered_bam)
        non_df[nrow(non_df)+1,] = row
    }
    write.table(df, file = paste0(path_out, "Stromal", "_bam_loc.txt"), quote = F, row.names = F, col.names = F)
    write.table(non_df, file = paste0(path_out, non_var, "_bam_loc.txt"), quote = F, row.names = F, col.names = F)
