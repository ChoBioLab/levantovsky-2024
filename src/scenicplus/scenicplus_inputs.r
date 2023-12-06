library(remotes)
library(BiocManager)

remotes::install_github("aertslab/cisTopic")
BiocManager::install(c('Rsubread', 'umap', 'Rtsne', 'ComplexHeatmap', 'fastcluster', 'data.table', 'rGREAT', 'ChIPseeker', 'TxDb.Hsapiens.UCSC.hg19.knownGene', 'org.Hs.eg.db'))
library(cisTopic)
library(Rsubread)
library(umap)
library(Rtsne)
library(ComplexHeatmap)
library(fastcluster)
library(data.table)
library(rGREAT)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(Seurat)
library(Signac)
library(SeuratDisk)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(ggplot2)
library(future)
library(chromVAR)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(harmony)
library(limma)
library(DropletUtils)

plan(
  multicore,
  workers = 20
) # parallelization
options(future.globals.maxSize = 20000 * 1024^2 * 1000)

# read in the object
obj <- readRDS("../data/2023-01-19_MultiomeReRun_Annotated.rds")
atac_assay <- obj@assays$ATAC@counts
meta <- obj@meta.data
# write into tsv and csv file as inputs to python
write.csv(meta, file = "scenicplus/resources/PerianalCD_MultiomeReRun_metadata.csv")
write.table(atac_assay, file = "scenicplus/resources/PerianalCD_MultiomeReRun_ATAC.tsv", sep = "\t")

# RNA assay write 
count_mt <- obj@assays$RNA@counts

DropletUtils::write10xCounts(path = "scenicplus/resources/PerianalCD_MultiomeReRun_RNA/", count_mt, barcodes = colnames(count_mt), 
                            gene.id = rownames(count_mt), gene.symbol = rownames(count_mt), 
                            gene.type = "Gene Expression", type = "sparse", version = "3")

# subset to test run with AA3 Ileum
sub <- subset(obj, subset = Libraries == "AA3_Ileum")
atac_assay <- sub@assays$ATAC@counts
meta <- sub@meta.data
rna_assay <- sub@assays$RNA@counts
write.csv(x = meta, file = "scenicplus/resources/AA3_Ileum_metadata.csv")
write.table(x = atac_assay, file = "scenicplus/resources/AA3_Ileum_ATAC.tsv", sep = "\t")
DropletUtils::write10xCounts(path = "scenicplus/resources/AA3_Ileum_RNA/", rna_assay, barcodes = colnames(rna_assay), 
                            gene.id = rownames(rna_assay), gene.symbol = rownames(rna_assay), 
                            gene.type = "Gene Expression", type = "sparse", version = "3")

# initialize and setup the cisTopic object
new_rownames <- sub("-", ":", rownames(atac_assay))
rownames(atac_assay) <- new_rownames
cisTopicObject <- createcisTopicObject(atac_assay, project = "PerianalCD_multi")

# add cell metadata

cisTopicObject <- addCellMetadata(cisTopicObject, cell.data = meta)

# Building the models
cisTopicObject <- runWarpLDAModels(cisTopicObject, topic=c(2, 10, 20, 30, 40), seed=12345, nCores=12, tmp = "scenicplus/tmp", addModels=FALSE, returnType = "selectedModel")

# select the best model (topic)
#par(mfrow=c(3,3))
#cisTopicObject <- selectModel(cisTopicObject, type='maximum')
#cisTopicObject <- selectModel(cisTopicObject, type='perplexity')
#cisTopicObject <- selectModel(cisTopicObject, type='derivative')

# interpreting the models
# 1. identification of cell states using the cell-cistopic distribution

cisTopicObject <- runUmap(cisTopicObject, target='cell')
par(mfrow=c(1,2))
plotFeatures(cisTopicObject, method='Umap', target='cell', topic_contr=NULL, colorBy=c('Bin', 'ClusterAnnotation'), cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE, col.low='darkgreen', col.mid='yellow', col.high='brown1', intervals=20)

par(mfrow=c(2,5))
plotFeatures(cisTopicObject, method='Umap', target='cell', topic_contr='Probability', colorBy=NULL, cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE)

cellTopicHeatmap(cisTopicObject, method='Probability', colorBy=c('Bin', 'ClusterAnnotation'))


# Export data to python
library(arrow)
path <- "scenicplus/resources/"
#cisTopicObject <- readRDS(PATH_TO_YOUR_CISTOPICOBJECT)
modelMat <- modelMatSelection(cisTopicObject, 'cell', 'Probability')
modelMat <- as.data.frame(modelMat)
write_feather(modelMat, sink=paste0(path, 'model_to_pycisTopic/cell_topic.feather'))
modelMat <- modelMatSelection(cisTopicObject, 'region', 'Probability', all.regions=TRUE)
modelMat <- as.data.frame(modelMat)
# Function to check and convert non-character columns to character
convert_to_character <- function(data_frame) {
  for (i in 1:ncol(data_frame)) {
    if (!is.character(data_frame[, i])) {
      data_frame[, i] <- as.character(data_frame[, i])
    }
  }
  return(data_frame)
}

# Apply the function to your data frame
modelMat <- convert_to_character(modelMat)

write_feather(modelMat, sink=paste0(path, 'model_to_pycisTopic/topic_region.feather'))
# You can also save the count.matrix
ct <- cisTopicObject@count.matrix
ct <- as.data.frame(ct)
write_feather(ct, sink=paste0(path, 'model_to_pycisTopic/count_matrix.feather'))

# test run with first 100 rows
atac_100 <- atac_assay[1:100, ]
cisTopicObject <- createcisTopicObject(atac_100, project = "test")

new_rownames <- sub("-", ":", rownames(atac_100))
rownames(atac_100) <- new_rownames
# get the row names of atac_assay and split the strings by "-", extract the second element and see if any of them cannot be turned into numeric
row_names <- rownames(atac_assay)
row_names1 <- strsplit(row_names, "-")
row_names <- unlist(row_names)

row_coordinates <- sapply(strsplit(rownames(atac_assay), "-"), `[`, 2)
row_coordinates <- as.numeric(row_coordinates)



