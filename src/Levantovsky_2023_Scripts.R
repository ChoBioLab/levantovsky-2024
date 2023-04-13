#code to generate figures in Levantovsky 2023

#libraries
library(dplyr)
library(DESeq2)
library(ggplot2)
library(ggpubr)
library(ggrepel)
options(ggrepel.max.overlaps = Inf) #https://ggrepel.slowkow.com/articles/examples.html
library(EnhancedVolcano)
library(Seurat)
library(SCpubr)
library(patchwork)
library(RColorBrewer)
library(Signac)
library(BSgenome.Hsapiens.UCSC.hg38)



##########################
####### RDS Files ########
##########################
#bulk RNA seq
dds.CD <- readRDS('/data0/RDS_Files/Med_Submission/MSCCR_dds-CD_AllBiopsies.rds') #Figures 1A-D
dds.ileum <- readRDS('/data0/RDS_Files/Med_Submission/MSCCR_dds_Ileum_Inf-Non.rds') #Figure S1B
dds.rectum <- readRDS('/data0/RDS_Files/Med_Submission/MSCCR_dds_Rectum_Inf-Non.rds') #Figure S1B

#scRNAseq
combined.cd <- readRDS('/data0/RDS_Files/Med_Submission/scRNA_combinedCD.RDS') #Figures S1D, S1E, S1F, 1H, 1I
perianal.cd <- readRDS('/data0/RDS_Files/Med_Submission/scRNA_perianalCD.RDS') #Figures 1F, 1G, 
fistula.clusters <- readRDS('/data0/RDS_Files/Med_Submission/scRNA_fistula_n-2.rds') #n = 2 patients
ms.int <- readRDS('/data0/RDS_Files/Med_Submission/scRNA_myeloid-stromal_integrated.RDS') #myeloid-stromal subclustering n = 14 patients, annotated
ileal.clusters <- readRDS('/data0/RDS_Files/Med_Submission/scRNA_Martin-Ileal.rds')

#multiome
multiome <- readRDS('/data0/RDS_Files/Med_Submission/multiome_annotated.RDS')

##########################
#### Figures 1, S1 #######
##########################

#Figure 1a
vsd.CD <- vst(dds.CD, blind = FALSE)
pca.CD.data <- plotPCA(vsd.CD, intgroup = 'BiopsyRegion', returnData = TRUE)
percentVar.CD <- round(100 * attr(pca.CD.data, "percentVar"))
ggplot(pca.CD.data, aes(PC1, PC2, color = BiopsyRegion)) +
  geom_point(size = 2) +
  xlab(paste0("PC1: ",percentVar.CD[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar.CD[2],"% variance")) +
  theme_light() +
  theme(legend.position="none", legend.title = element_blank()) +
  theme(axis.text = element_text(colour = "black", family = "Arial", size = 14),
        axis.title = element_text(size = 14)) +
  scale_color_manual(values = c('slategray3', 'lightgreen', 'khaki', 'burlywood4', 'coral', 'deeppink', 'navy'),
                     breaks = c('Ileum', 'Cecum', 'Right_colon', 'Transverse', 'Left_colon', 'Sigmoid', 'Rectum')) +
  coord_fixed()

#Figure 1B
ileal.wilcox <- compare_means(PC1 ~ BiopsyRegion,  data = pca.CD.data, ref.group = "Ileum",
              method = "wilcox")
ggplot(data = pca.CD.data, aes(x = BiopsyRegion, y = PC1, color = BiopsyRegion)) +
  geom_boxplot() +
  scale_color_manual(values = c('slategray3', 'lightgreen', 'khaki', 'burlywood4', 'coral', 'deeppink', 'navy'),
                     breaks = c('Ileum', 'Cecum', 'Right_colon', 'Transverse', 'Left_colon', 'Sigmoid', 'Rectum')) +
  scale_x_discrete(limits=c('Ileum', 'Cecum', 'Right_colon', 'Transverse', 'Left_colon', 'Sigmoid', 'Rectum')) +
  theme_classic() +
  theme(axis.text = element_text(colour = "black", family = "Arial", size = 12),
        axis.title = element_text(size = 13)) +
  theme(legend.position = 'none') +
  xlab(label = NULL) +
  stat_pvalue_manual(
    ileal.wilcox, 
    y.position = 35, step.increase = 0.1,
    label = "p.adj"
  ) 

#Figure 1C
pca.CD.data2 <- plotPCA(vsd.CD, intgroup = 'BiopsyType', returnData = TRUE)
percentVar.CD2 <- round(100 * attr(pca.CD.data2, "percentVar"))
ggplot(pca.CD.data2, aes(PC1, PC2, color = BiopsyType)) +
  geom_point(size = 2) +
  xlab(paste0("PC1: ",percentVar.CD2[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar.CD2[2],"% variance")) +
  theme_light() +
  theme(legend.position="none", legend.title = element_blank()) +
  theme(axis.text = element_text(colour = "black", family = "Arial", size = 14),
        axis.title = element_text(size = 14)) +
  scale_color_manual(values = c('red', 'blue')) +
  coord_fixed()

#Figure 1D
compare_means(BiopsyType ~ PC2, data = pca.CD.data2, method = 'wilcox.test')
ggplot(data = pca.CD.data2, aes(x = BiopsyType, y = PC2, color = BiopsyType)) +
  geom_boxplot() +
  scale_color_manual(values = c('red', 'blue')) +
  theme_classic() +
  theme(axis.text = element_text(colour = "black", family = "Arial", size = 12),
        axis.title = element_text(size = 13)) +
  theme(legend.position="none") +
  xlab(label = NULL) +
  stat_compare_means(method = 't.test', label.x = 1.4)

#Figure S1B
res.ileum <- results(dds.ileum, contrast = c('BiopsyType', 'Inflamed', 'Non_inflamed'), alpha = 0.05)
#summary(res.ileum) #4071 up, 13% ; 3593 down, 11%

custom.cols.ileum <- ifelse(
  res.ileum$log2FoldChange < 0, 'blue', 'red')
names(custom.cols.ileum)[custom.cols.ileum == 'blue'] <- 'Non_inflamed'
names(custom.cols.ileum)[custom.cols.ileum == 'red'] <- 'Inflamed'

EnhancedVolcano(res.ileum,
                lab = rownames(res.ileum),
                legendPosition = 'none',
                x = 'log2FoldChange',
                y = 'padj',
                xlim = c(-6,6),
                ylim = c(0, 40),
                pCutoff = 5e-10,
                pointSize = 1,
                #col = c('black', 'orange', 'skyblue', 'lightcoral'),
                colCustom = custom.cols.ileum,
                labSize = 5,
                ylab = bquote(~-Log[10]~italic(Padj)),
                title = 'Inflamed Vs. Non-inflamed',
                subtitle = 'Ileum Biopsies',
                selectLab = c('ANXA10', 'GALNTL6', 'MUC5AC', 'S100A8', 'TCN1'),
                boxedLabels = TRUE,
                drawConnectors = TRUE,
                colConnectors = 'black',
                widthConnectors = 0.5
)+ ggplot2::labs(title = NULL, subtitle = NULL, caption = NULL)

res.rectum <- results(dds.rectum, contrast = c('BiopsyType', 'Inflamed', 'Non_inflamed'), alpha = 0.05)
#summary(res.rectum) #6222 up, 20% ; 4466 down, 14%

custom.cols.rectum <- ifelse(
  res.rectum$log2FoldChange < 0, 'blue', 'red')
names(custom.cols.rectum)[custom.cols.rectum == 'blue'] <- 'Non_inflamed'
names(custom.cols.rectum)[custom.cols.rectum == 'red'] <- 'Inflamed'

EnhancedVolcano(res.rectum,
                lab = rownames(res.rectum),
                legendPosition = 'none',
                x = 'log2FoldChange',
                y = 'padj',
                xlim = c(-6,6),
                ylim = c(0, 80),
                pCutoff = 5e-10,
                pointSize = 1,
                #col = c('black', 'orange', 'skyblue', 'lightcoral'),
                colCustom = custom.cols.rectum,
                labSize = 5,
                ylab = bquote(~-Log[10]~italic(Padj)),
                title = 'Inflamed Vs. Non-inflamed',
                subtitle = 'Rectum Biopsies',
                selectLab = c('IDO1','TNIP3','CHI3L1','CXCL11','S100A9'),
                boxedLabels = TRUE,
                drawConnectors = TRUE,
                colConnectors = 'black',
                widthConnectors = 0.5
)+ ggplot2::labs(title = NULL, subtitle = NULL, caption = NULL)

#Figure S1G
library(ComplexHeatmap)
library(circlize)

region <- c('ileum', 'rectum')
anno = HeatmapAnnotation(region = region, col=list(region = c('ileum' = 'slategray3', 'rectum' = 'navy')),which = 'row', show_annotation_name = FALSE)

means <- read.csv('/data0/RDS_Files/Med_Submission/matrix_means.csv', row.names = 1)
means <- as.matrix(means)
means <- t(means)
means_pal <- colorRamp2(c(0,10,100,5000,6000), c('white', 'lightcyan', 'dodgerblue','cornflowerblue','darkblue'))
means.hm <- Heatmap(means, name = 'Mean Exp.', col = means_pal, right_annotation = anno,
        cluster_rows = FALSE, cluster_columns = TRUE, border = TRUE, 
        width = ncol(means)*unit(3, 'mm'),
        height = nrow(means)*unit(5, "mm"),
        row_names_gp = grid::gpar(fontsize = 10),
        column_names_gp = grid::gpar(fontsize = 8),
        heatmap_legend_param = list(direction = "horizontal"))
means.hm

fc <- read.csv('/data0/RDS_Files/Med_Submission/matrix_fc.csv', row.names = 1)
fc <- as.matrix(fc)
fc <- t(fc)
fc_pal <- colorRamp2(c(-0.1,0.5,1,3,5), c('white', 'cornsilk', 'wheat2','wheat3','tan4'))
fc.hm <- Heatmap(fc, name = 'Log2 Fold Change', col = fc_pal, right_annotation = anno,
        cluster_rows = FALSE, cluster_columns = TRUE, border = TRUE,
        width = ncol(means)*unit(3, 'mm'),
        height = nrow(means)*unit(5, "mm"),
        row_names_gp = grid::gpar(fontsize = 10),
        column_names_gp = grid::gpar(fontsize = 8),
        heatmap_legend_param = list(direction = "horizontal"))
fc.hm

padj <- read.csv('/data0/RDS_Files/Med_Submission/matrix_padj.csv', row.names = 1)
padj <- as.matrix(padj)
#log transform p-values
padj <- -log10(padj)
summary(padj)
#new range: 0-53
padj <- t(padj)
#padj_pal <- colorRamp2(c(0,0.00001,1), c('red', 'white', 'blue'))
padj_pal <- colorRamp2(c(0, 5,10, 15, 25,55), c("white", "pink",  "pink1", "pink2", "pink3", "pink4"))
padj.hm <- Heatmap(padj, name = '-Log10(P-adj.)', col = padj_pal, right_annotation = anno,
        cluster_rows = FALSE, cluster_columns = TRUE, border = TRUE,
        width = ncol(means)*unit(3, 'mm'),
        height = nrow(means)*unit(5, "mm"),
        row_names_gp = grid::gpar(fontsize = 10),
        column_names_gp = grid::gpar(fontsize = 8),
        heatmap_legend_param = list(direction = "horizontal"))
padj.hm

hm_list = fc.hm %v% padj.hm %v% means.hm
draw(hm_list, heatmap_legend_side='top', annotation_legend_side='top')

#Figure S1D
DimPlot(combined.cd, label=TRUE, label.box=TRUE, repel=TRUE)
#Figure S1E
DimPlot(combined.cd, group.by = 'orig.ident') + ggplot2::theme(legend.position = 'bottom') + ggtitle(label=NULL)
#Figure S1F
library(tidyseurat)
tibble <- structure(combined.cd)
tibble <- as_tibble(tibble)
#split into groups
tibble.martin <- subset(tibble, subset = orig.ident=='Martin')
tibble.martin <- tibble.martin[order(tibble.martin$ClusterAnnotation, decreasing=TRUE),]
table(tibble.martin$PubID_Obj)

tibble.perianal <- subset(tibble, subset = orig.ident=='PerianalCD')
tibble.perianal <- tibble.perianal[order(tibble.perianal$ClusterAnnotation, decreasing=TRUE),]
table(tibble.perianal$PubID_Obj)

#color pal
bar.pal <- c('#F8766D', '#CD9600',  '#E68613', '#FF68A1',  '#ABA300',  '#00A9FF',  '#00C19A', '#0CB702','#FF61CC' ,  '#7CAE00' ,'#00BE67', '#8494FF',  '#C77CFF', '#00BFC4' , '#00B8E7', '#ED68ED')

#orders to set x-axis
patient.list.order1 <- c(
    'Patient6_Inv',
    'Patient6_Non',
    'Patient7_Inv',
    'Patient7_Non',
    'Patient8_Inv',
    'Patient8_Non',
    'Patient10_Inv',
    'Patient10_Non',
    'Patient14_Inv',
    'Patient14_Non',
    'Patient15_Inv',
    'Patient15_Non',
    'Patient16_Inv',
    'Patient16_Non'
    )

patient.list.order2 <- c(
  "Perianal1_Transverse", "Perianal1_Sig",
  "Perianal2_Sig", "Perianal2_Rect",
  "Perianal3_Sig", "Perianal3_Rect",
  "Perianal4_TI", "Perianal4_Rect",
  "Perianal5_Sig", "Perianal5_Rect",
  "Perianal6_Sig", "Perianal6_Rect",
  "Perianal7_Sig", "Perianal7_Rect",
  "Perianal8_Sig", "Perianal8_Rect",
  "Perianal9_Sig", "Perianal9_Rect",
  "Perianal10_Sig", "Perianal10_Rect",
  "Perianal11_Fistula", "Perianal11_Rect",
  "Perianal12_Sig", "Perianal12_Rect",
  "Perianal13_Sig", "Perianal13_Rect",
  "Perianal14_Fistula", "Perianal14_Rect"
)

proportion.martin <- ggplot(tibble.martin, aes(x = PubID_Obj, fill = ClusterAnnotation)) +
  geom_bar(position = "fill", width = 0.98) + 
  theme_classic() +
  theme(legend.position = 'right', legend.title = element_blank()) +
  theme(axis.text.x = element_text(colour = "black", family = "Arial", size = 12, angle=90),
        axis.text.y = element_text(colour = "black", family = "Arial", size = 12)) +
  xlab(label = NULL) +
  ylab(label = NULL) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) + 
  scale_fill_manual(values=bar.pal) +
  scale_x_discrete(limits = patient.list.order1)
proportion.martin

proportion.perianal <- ggplot(tibble.perianal, 
  aes(x = PubID_Obj, fill = ClusterAnnotation)) +
  geom_bar(position = "fill", width = 0.98) + 
  theme_classic() +
  theme(legend.position = 'right', legend.title = element_blank()) +
  theme(axis.text.x = element_text(colour = "black", family = "Arial", size = 12, angle=90),
        axis.text.y = element_text(colour = "black", family = "Arial", size = 12)) +
  xlab(label = NULL) +
  ylab(label = NULL) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) + 
  scale_fill_manual(values=bar.pal) +
  scale_x_discrete(limits = patient.list.order2)
proportion.perianal

proportion.martin + proportion.perianal + plot_layout(guides='collect')

#Figure 1F
DimPlot(perianal.cd, label=FALSE) + ggplot2::theme(legend.position = 'bottom')

#Figure 1G
umap.pal <- c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69",
  "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#1B9E77", "#D95F02" )
DimPlot(perianal.cd, group.by='PubID', cols=umap.pal) + ggplot2::theme(legend.position = 'bottom') + ggtitle(label=NULL)

#Figure 1H
vln.iga <- VlnPlot(combined.cd, features = c('IGHA1'), split.by = 'orig.ident', split.plot = TRUE, idents = c('IgA Plasma'), pt.size=0)
iga.res <- prop.test(x=c(5047,39180), n=c(49331,130961)) #proportions z-test

#Figure 1I
combined.IRS <- subset(combined.cd, subset = Region == c('Rectum', 'Sigmoid', 'Ileum'))
do_DotPlot(integrated.IRS, features = c('CXCL3', 'CCL2', 'CCL7', 'CCR2', 'CCR4', 'IL1RN', 'IL10', 'IL11', 'IL6', 'IL6ST'), 
   split.by='Region', cluster.idents = TRUE, flip = F, 
   use_viridis=TRUE, viridis_color_map='viridis', viridis_direction=1) + 
   ggplot2::theme(legend.position='none') +
   ggplot2::theme(axis.text.y=element_text(size=10))


##########################
#### Figures 2, S3 #######
##########################

#Figure S3A
DimPlot(fistula.clusters, label = TRUE, repel = TRUE, label.box = TRUE) +
  ggplot2::labs(title = NULL) + ggplot2::theme(legend.position = 'none', axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) + 
  xlab(NULL) + ylab(NULL)
#Figure S3C
do_DotPlot(clusters, features = c('TGFB1', 'SNAI1', 'SNAI2', 'ETS1', 'DKK1', 'IL13', 'ITGB6', 'MMP3', 'MMP9', 'MMP13'), 
           colors.use = c('blue', 'red'), group.by = 'Region', legend.length = 8,
           cluster.idents = T)

#generate myeloid-stromal subclusters
#load objects from n = 12 colorectal biopsy patients
myeloid <- readRDS('/data0/RDS_Files/Med_Submission/scRNA_myeloid_subset.rds')
stromal <- readRDS('/data0/RDS_Files/Med_Submission/scRNA_stromal_subset.rds')

#subset fistula clusters for stromal, myeloid 
stromal.subset <- c('Endothelial, Lymphatic, Blood Vessel', 'Pericytes, Contractile', 'Endothelial ', 'Fibroblasts', 'Myofibroblasts')
myeloid.subset <- c('moMac, moDC', 'Inflammatory Macrophages', 'Plasmacytoid DC', 'Mast Cells', 'Dendritic Cells')
fistula.sub <- subset(fistula.clusters, idents = c(stromal.subset, myeloid.subset))
table(fistula.sub@active.ident) #n = 5699

#re-cluster all stromal/myeoloid cells together
integration.list <- list(stromal, myeloid, fistula.sub)
anchors <- FindIntegrationAnchors(object.list = integration.list, dims = 1:15)  
integrated <- IntegrateData(anchorset = anchors, dims = 1:15)  
DefaultAssay(integrated) <- "integrated"
integrated <- ScaleData(integrated)
integrated <- RunPCA(integrated, npcs = 20, verbose = FALSE)
ElbowPlot(integrated)
integrated <- RunUMAP(integrated, reduction = "pca", dims = 1:18)
integrated <- FindNeighbors(integrated, reduction = "pca", dims = 1:18)
integrated <- FindClusters(integrated, resolution = 0.8)
DefaultAssay(integrated) <- "RNA"
DimPlot(integrated, label = TRUE, label.size = 8) + ggplot2::theme(legend.position = "none")
write_clip(table(integrated@active.ident))

#assign fistula metadata
DimPlot(object = integrated, group.by = 'Region')
integrated$FistulaAnnotation <- ifelse(integrated$Region == 'Fistula', 'Fistula', 'Tissue')

#Figure S3D
all.genes <- rownames(integrated)
integrated <- ScaleData(integrated, features = all.genes)
levels.hm <- c('0',	'1',	'6',	'25',	'2',	'7',	'17',	'23',	'12',	'14',	'21',	'10',	'11',	'16',	'22',	'18',	'5',	'3',	'20',	'15',	'26',	'4',	'8',	'9',	'24',	'13',	'19')
integrated@active.ident <- factor(integrated@active.ident, levels = levels.hm)
gene.list.hm <- c('ADAMDEC1',	'CCL8',	'PTN',	'CCL2',	'STMN2',			
	'ABCA8',	'ADH1B',	'FN1',	'DCN',			
'CFD',	'APOE',	'PTGDS',	'CCL11',	'CXCL14',			
'COL6A2',		'PLAT',	'SOX6',	'MT1A',	'PDGFRA',		
	' POSTN',	' SOX6',	' F3',	' COL6A2',			
'POSTN',		'HSD17B2',					
'MTRNR2L8',	'CXCL13',	'FDCSP',	'CHI3L1',	'CHI3L2',			
'CSF3',	'MMP1',	'ADM',	'TUBB4B',	'MT1E',			
'TAGLN',	'HHIP',	'SOSTDC1',	'ACTG2',	'MYH11',	'ACTA2',		
'TPM2',			'TPM1',				
'RERGL',	'MUSTN1',	'PLN',		'NET1',			
'RGS5',	'NOTCH3',	'HIGD1B',	'STEAP4',	'MGP',			
'PLVAP',	'FLT1',	'PECAM1',	'PODXL',	'AQP1',	'CD36',		
	'CLDN5',	'CD320',	'FABP5',				
'ACKR1',	'CCL14',	'VWF',		'CPE',			
'NRXN1',	'CRYAB',	'PLP1',	'GPM6B',	'CDH19',			
'IL1B',	'CCL3',	'S100A8',	'EREG',	'CCL3L1',	'S100A9',	'SOD2',	
'CCL18',	'C1QB',	'C1QA',	'C1QC',	'LYZ',	'CD14',		
	'CYBA',			'TYROBP',	'CST3',	'MERTK',	'CD206',
'HLA-DQB1',	'HLA-DQA1',	'HLA-DPB1',	'HLA-DRA',	'HLA-DPA1',	'CLEC9A',	'CLEC10A',	
'CCR7',	'CCL17',	'IL7R',	'CSF2RA',	'DAPP1',	'CFP',	'CD1E',	
'MS4A2',	'CTSG',	'HPGDS',	'CPA3',	'CD69',			
'TPSAB1',	'TPSB2',	'LTC4S',		'JUND',			
	'ADCYAP1',	'HDC',		'IL1RL1',			
	'BATF',		'HSP90AA1',	'CACYBP',			
'MZB1',	'CD79A',	'MS4A1',	'CXCR4',	'PLCG2',			
	'KLRB1',	'CD3D',	'CCL5',	'TRBC2')
unanno.hm <- DoHeatmap(integrated, features = gene.list.hm, ) + NoLegend()

#Figure 2D
DimPlot(ms.int, label = T, repel = T, label.size = 5) +
  ggplot2::labs(title = NULL) + ggplot2::theme(legend.position = 'none') 
#Figure 2E
ms.int$Anno <- ms.int@active.ident

tibble2 <- structure(ms.int)
tibble2 <- as_tibble(tibble2)

ggplot(tibble2, aes(x = Anno, fill = FistulaAnnotation)) +
  geom_bar(position = "fill", width = 0.98) + 
  theme_classic() +
  theme(legend.position = 'none') +
  #guides(fill=guide_legend(nrow=5, byrow=TRUE)) + #https://stackoverflow.com/questions/27130610/legend-on-bottom-two-rows-wrapped-in-ggplot2-in-r
  theme(axis.text.x = element_text(colour = "black", family = "Arial", size = 12),
        axis.text.y = element_text(colour = "black", family = "Arial", size = 12)) +
  xlab(label = NULL) +
  ylab(label = NULL) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  scale_fill_manual(values=c('maroon3', 'navy')) +
  coord_flip()
#Figure S3B
ms.int <- SetIdent(ms.int, value = 'Anno')
VlnPlot(ms.int, features = c('ACTA2', 'MYH11'), split.by = 'FistulaAnnotation', split.plot = TRUE,
        idents = c('Myofibroblasts'), cols = c('maroon3', 'navy'))

#Figure 2F
do_DotPlot(ms.int, features=c('CXCL14', 'CXCR4', 'CSF3', 'CSF3R','CCL2', 'NR3C1', 'THY1', 'ITGAM', 'ITGAX',
                              'PTGES2', 'PTGER4', 'TNF', 'TNFRSF1A', 'CSF1', 'CSF1R', 'FTL', 'FTH1', 'SCARA5'),
           cluster.idents = TRUE, flip = F, 
           use_viridis=TRUE, viridis_color_map='viridis', viridis_direction=1
           ) + ggplot2::theme(legend.position = 'bottom')
#Figure 2G
stromal.subset <- c('CHI3L1hi Fibroblasts', 'ADAMDEC1hi Fibroblasts', 'PDGFRAhi Fibroblasts', 'CCL11hi Fibroblasts', 'Pericytes',
                    'CD36hi Endothelial', 'Myofibroblasts', 'Enteric Neurons', 'ACKR1hi Endothelial, Lymphatic')
stromal.int <- subset(ms.int, idents = stromal.subset)
DEG.stromal.n14 <- FindMarkers(stromal.int, ident.1 = 'Fistula', group.by = 'FistulaAnnotation', logfc.threshold = 0)

custom.volcano.cols3 <- ifelse(
  DEG.stromal.n14$avg_log2FC < 0, 'navy', 'maroon3')
names(custom.volcano.cols3)[custom.volcano.cols3 == 'navy'] <- 'Rectum'
names(custom.volcano.cols3)[custom.volcano.cols3 == 'maroon3'] <- 'Fistula'

DEG.stromal.n14 %>%
  top_n(n = 20, wt = avg_log2FC) -> top20.stromal.up
DEG.stromal.n14 %>%
  top_n(n = -20, wt = avg_log2FC) -> top20.stromal.down
stromal.list1 <- rownames(top20.stromal.up)
stromal.list2 <- rownames(top20.stromal.down)
stromal.list.volcano <- c(stromal.list1, stromal.list2)

stromal.volcano.plot <- EnhancedVolcano(DEG.stromal.n14,
                                        lab = rownames(DEG.stromal.n14),
                                        x = 'avg_log2FC',
                                        y = 'p_val_adj',
                                        xlim = c(-6,5),
                                        ylim = c(0, 300),
                                        pCutoff = 5e-5,
                                        pointSize = 1,
                                        colCustom = custom.volcano.cols3,
                                        axisLabSize = 11,
                                        labSize = 8,
                                        ylab = NULL,
                                        xlab = NULL,
                                        title = 'Fistula Tract vs. Rectal Mucosa Origin',
                                        subtitle = 'stromal',
                                        selectLab = stromal.list.volcano,
                                        drawConnectors = TRUE,
                                        colConnectors = 'black',
                                        widthConnectors = 0.5,
                                        legendPosition = 'none') + ggplot2::labs(title = NULL, subtitle = NULL, caption = NULL)
#Figure 2H
myeloid.subset <- c('CD14hi moMac, moDC', 'Mast Cells', 'Inflammatory Macrophages', 'Dendritic Cells', 'Activated Dendritic Cells')
myeloid.int <- subset(ms.int, idents = myeloid.subset)
DEG.myeloid.n14 <- FindMarkers(myeloid.int, ident.1 = 'Fistula', group.by = 'FistulaAnnotation', logfc.threshold = 0)

custom.volcano.cols2 <- ifelse(
  DEG.myeloid.n14$avg_log2FC < 0, 'navy', 'maroon3')
names(custom.volcano.cols2)[custom.volcano.cols2 == 'navy'] <- 'Rectum'
names(custom.volcano.cols2)[custom.volcano.cols2 == 'maroon3'] <- 'Fistula'

DEG.myeloid.n14 %>%
  top_n(n = 20, wt = avg_log2FC) -> top20.myeloid.up
DEG.myeloid.n14 %>%
  top_n(n = -20, wt = avg_log2FC) -> top20.myeloid.down
myeloid.list1 <- rownames(top20.myeloid.up)
myeloid.list2 <- rownames(top20.myeloid.down)
myeloid.list.volcano <- c(myeloid.list1, myeloid.list2)

myeloid.volcano.plot <- EnhancedVolcano(DEG.myeloid.n14,
                                        lab = rownames(DEG.myeloid.n14),
                                        x = 'avg_log2FC',
                                        y = 'p_val_adj',
                                        xlim = c(-4.5,4.5),
                                        ylim = c(0, 350),
                                        pCutoff = 5e-50,
                                        pointSize = 1,
                                        colCustom = custom.volcano.cols2,
                                        axisLabSize = 11,
                                        labSize = 8,
                                        ylab = NULL,
                                        xlab = NULL,
                                        title = 'Fistula Tract vs. Rectal Mucosa Origin',
                                        subtitle = 'Myeloid',
                                        selectLab = myeloid.list.volcano,
                                        boxedLabels = FALSE,
                                        drawConnectors = TRUE,
                                        colConnectors = 'black',
                                        widthConnectors = 0.5,
                                        legendPosition = 'none') + ggplot2::labs(title = NULL, subtitle = NULL, caption = NULL)
#Figure 2I
feat.CHI3L1 <- FeaturePlot(ms.int, features = c('CHI3L1'), order=TRUE)
feat.CD14 <- FeaturePlot(ms.int, features=c('CD14'), order=TRUE)
feat.CD14 / feat.CHI3L1

#Figure S2E
DEG.myeloid.n14 %>%
  top_n(n = 30, wt = avg_log2FC) -> top30.myeloid.up
myeloid.list <- rownames(top30.myeloid.up)
DEG.stromal.n14 %>%
  top_n(n = 30, wt = avg_log2FC) -> top30.stromal.up
stromal.list <- rownames(top30.stromal.up)

myeloid.int <- AddModuleScore(myeloid.int, features = stromal.list, name = 'StromalUp')
FeaturePlot(myeloid.int, features = 'StromalUp1', order = T)
myeloid.module <- do_FeaturePlot(myeloid.int, features = 'StromalUp1', order = T, pt.size = 0.5, enforce_symmetry = T, 
                                 split.by = 'FistulaAnnotation',legend.title = 'Stromal Module Score', font.size = 12, legend.length = 6, legend.position = 'right', legend.width = 0.5)
myeloid.module

stromal.int<- AddModuleScore(stromal.int, features = myeloid.list, name = 'MyeloidUp')
FeaturePlot(stromal.int, features = 'MyeloidUp1')
stromal.module <- do_FeaturePlot(stromal.int, features = 'MyeloidUp1', enforce_symmetry = T, order = T, 
                                 pt.size = 0.5, split.by = 'FistulaAnnotation', legend.title = 'Myeloid Module Score', font.size = 12, legend.length = 6, legend.position = 'right', legend.width = 0.5)
stromal.module

myeloid.module / stromal.module

##########################
#### Figures 3, S4 #######
##########################

#Figure 3A
fibrotic.genes <- c('COL1A1', 'ACTA2', 'CCL5', 'EDN1', 'VCAM1', 'CXCL8', 'ENPP2', 
                    'CTGF', 'COL3A1', 'THBS1', 'CALU', 'ITGAV', 'RECK', 'COL1A2')
destructive.genes <- c('FAP', 'PDPN', 'THY1',
                       'MMP3', 'MMP9', 'MMP13', 'IL18', 'CCL9', 'TNFSF11', 
                       'CXCL1', 'PRG4', 'CLIC5', 'TSPAN15', 'COL22A1')

enrichment.genes <- list("Fibrotic"=fibrotic.genes,
                         "Destructive"=destructive.genes)
do_EnrichmentHeatmap(ms.int, input_gene_list = enrichment.genes, 
                     cluster_rows = T, symmetrical_scale = T, use_viridis = F)

#ileal & fistula fibroblasts recluster
fibroblasts.ileal <- subset(ileal.clusters, idents = c('Activated fibroblasts', 'Fibroblasts'))
fibroblasts.involved <- subset(fibroblasts.ileal, subset = group == 'inflamed')

fibroblasts.rectal <- subset(ms.int, idents = c('CHI3L1hi Fibroblasts', 'ADAMDEC1hi Fibroblasts', 'PDGFRAhi Fibroblasts', 'CCL11hi Fibroblasts'))
fibroblasts.fistula <- subset(fibroblasts.rectal, subset = FistulaAnnotation == 'Fistula')

fibroblasts.involved$Original <- 'Ileal Involved'
fibroblasts.fistula$Original <- 'Fistula'

integration.list2 <- list(fibroblasts.stricture, fibroblasts.fistula)
anchors2 <- FindIntegrationAnchors(object.list = integration.list2, dims = 1:15)  
integrated2 <- IntegrateData(anchorset = anchors2, dims = 1:15)  
DefaultAssay(integrated2) <- "integrated"
integrated2 <- ScaleData(integrated2)
integrated2 <- RunPCA(integrated2, npcs = 20, verbose = FALSE)
ElbowPlot(integrated2)
integrated2 <- RunUMAP(integrated2, reduction = "pca", dims = 1:18)
integrated2 <- FindNeighbors(integrated2, reduction = "pca", dims = 1:18)
integrated2 <- FindClusters(integrated2, resolution = 0.8)
DefaultAssay(integrated2) <- "RNA"

#Figure S4A
DimPlot(integrated2, label = FALSE, group.by = 'Original', cols = c('maroon3', 'slategray3')) +
  ggplot2::labs(title = NULL) + ggplot2::theme(legend.position = 'bottom', axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) + 
  xlab(NULL) + ylab(NULL)

#Figure S4B
integrated2.df = data.frame(CHI3L1 = integrated2[['RNA']]@data['CHI3L1',],
                      OriginalRegion = integrated2$Original)
mean.CHI3L1 <- aggregate(CHI3L1 ~ OriginalRegion, integrated2.df, mean)
CHI3L1.vln <- VlnPlot(integrated2, features = 'CHI3L1', group.by = 'Original', cols = c('maroon3', 'slategray3'))
CHI3L1.vln + stat_compare_means(method = 'wilcox.test', label = 'p.format') + 
  ggplot2::labs(title = NULL) +
  stat_summary(fun=mean, geom="point", shape=16, size=3, color='black') +
  ylab(label = 'CHI3L1 Expression') +
  theme(axis.text.x  = element_text(colour = "black", family = "Arial", size = 12, angle = 0, hjust = 0.5),
        axis.text.y  = element_text(colour = "black", family = "Arial", size = 12, angle = 0),
        axis.title.x = element_blank()) +
  theme(legend.position = 'none')

#Figure 3B
CHI3L1.pos <- subset(integrated2, subset = CHI3L1 > 0)
table(CHI3L1.pos$Original) #305 fistula, 280 ileal involved
DEG.CHI3L1pos <- FindMarkers(CHI3L1.pos, ident.1 = 'Fistula', group.by = 'Original', logfc.threshold = 0)

DEG.CHI3L1pos %>%
  top_n(n = 10, wt = avg_log2FC) -> top10.up
top10.up <- top10.up[order(-top10.up$avg_log2FC),]

DEG.CHI3L1pos %>%
  top_n(n = -10, wt = avg_log2FC) -> top10.down
top10.down <- top10.down[order(-top10.down$avg_log2FC),]

top10.up <- rownames(top10.up)
top10.down <- rownames(top10.down)
DEG.CHI3L1.list <- c(top10.up, top10.down)

all.genes <- rownames(CHI3L1.pos)
CHI3L1.pos <- ScaleData(CHI3L1.pos, features = all.genes)
DoHeatmap(CHI3L1.pos, features = DEG.CHI3L1.list, group.by = 'Original', 
          group.colors = c('maroon3', 'slategray3'), label = FALSE, ) + RotatedAxis() +
  scale_fill_gradientn(colors = c("blue", "white", "red")) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=12)) +
  theme(legend.position = 'none')

#Figure S4C
m2.markers <- c('CD163', 'CD200R1', 'CLEC10A', 'CXCR1', 'CXCR2', 'MRC1', #surface markers, MRC1 = CD206, CD163 = mac/mono marker
                'PPARG', 'STAT6',  #transcription factors
                'CCL1', 'CCL2', 'TGFB1', 'IL10', 'IL1RN', 'CCL14', 'CCL17', 'CCL18', 'CCL22', 'CCL23' #secreted
                )

ms.int <- AddModuleScore(ms.int, features = c(m2.markers), name = 'M2_Module')
m2.vln.1 <- VlnPlot(ms.int, features = 'M2_Module1', split.by = 'FistulaAnnotation', 
                    idents = c('Inflammatory Macrophages'), cols = c('maroon3', 'navy'), pt.size = 1)
m2.vln.2 <- VlnPlot(ms.int, features = 'M2_Module1', split.by = 'FistulaAnnotation', 
                    idents = c('CD14hi moMac, moDC'), cols = c('maroon3', 'navy'), pt.size = 1)
vlndf.m2 = data.frame(M2_Module = ms.int@meta.data$M2_Module1,
                      FistulaAnnotation = ms.int$FistulaAnnotation,
                      FistulaCurrent = ms.int$FistulaCurrent)
mean.m2.fistula <- aggregate(M2_Module ~ FistulaAnnotation, vlndf.m2, mean)
m2.vln.1 <- m2.vln.1 + stat_compare_means(method = 'wilcox.test', label = 'p.format', bracket.size = 1) + 
  ggplot2::labs(title = NULL) +
  theme(axis.text.x  = element_text(colour = "black", family = "Arial", size = 12, angle = 0, hjust = 0.5),
        axis.text.y  = element_text(colour = "black", family = "Arial", size = 12, angle = 0),
        axis.title = element_blank()) +
  theme(legend.position = 'none') +
  geom_violin(alpha=0.1)
m2.vln.2 <- m2.vln.2 + stat_compare_means(method = 'wilcox.test', label = 'p.format') + 
  ggplot2::labs(title = NULL) +
  theme(axis.text.x  = element_text(colour = "black", family = "Arial", size = 12, angle = 0, hjust = 0.5),
        axis.text.y  = element_text(colour = "black", family = "Arial", size = 12, angle = 0),
        axis.title = element_blank()) +
  theme(legend.position = 'none')
(m2.vln.1 | m2.vln.2)

#Figure 3C
#to convert mouse to human genes from Sanin 2022
mouse_human_genes = read.csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")
convert_mouse_to_human <- function(gene_list){
  
  output = c()
  
  for(gene in gene_list){
    class_key = (mouse_human_genes %>% filter(Symbol == gene & Common.Organism.Name=="mouse, laboratory"))[['DB.Class.Key']]
    if(!identical(class_key, integer(0)) ){
      human_genes = (mouse_human_genes %>% filter(DB.Class.Key == class_key & Common.Organism.Name=="human"))[,"Symbol"]
      for(human_gene in human_genes){
        output = append(output,human_gene)
      }
    }
  }
  
  return (output)
}

#p1
phagocytic.mouse <- c('Pf41',	'Ltc4s',	'C1qa',	'Ctsb',	'C1qc',	'Apoe',	'C1qb',	'Csf1r',	'Dab2',	'Gas6',	'Ninj1',	'Fcgrt',	'Timp2',	'Lgmn',	'Pltp',	'Itm2b',	'Wfdc17',	'Blvrb',	'Serinc3',	'Rnase4',	'Grn',	'Ctsd1',	'Tgfbr2',	'C3ar1',	'Nrp1',	'Trf',	'Cfh',	'Cd63',	'C5ar1',	'Ang',	'Dhrs3',	'Ctsl',	'Glul',	'Pla2g15',	'Emp1',	'Cd33',	'Cltc1',	'Rhob',	'Egr1',	'Lamp1',	'Ptpn18',	'Ccl6',	'Aplp2',	'Atp13a2',	'Ier3',	'Ms4a7',	'Ptger4',	'Aldh2',	'Mafb',	'Hexa',	'Ftl1',	'Adgre1',	'Zfp36l1','St6galnac4',	'Tpp1',	'Plk2',	'Tcn2',	'AI467606',	'Fcgr2b',	'Hmox1',	'Man2b1',	'Camk1',	'Marcks',	'Pigk',	'Clta',	'Tmem37',	'Adam15',	'Tcf3',	'Prkcb',	'Qk',	'Tbc1d17',	'Cfp1',	'Cd81',	'Vps33a',	'Ccl2',	'Gramd1a',	'Snx2',	'Galk2',	'Canx',	'Tbxas1',	'Ap2a2',	'Dusp6',	'Snx6',	'Golph3',	'Calm2',	'Plod1',	'Trib1',	'Fam3c',	'Aldh9a1',	'Gnl3',	'Cd68',	'Ptprj',	'Rab31',	'Herpud21',	'Fam214b',	'Dok3',	'Sec14l1',	'Man1a',	'Dnajb9',	'Stard5',	'Cst3',	'Ggnbp2',	'Qrich1',	'Cln8',	'Ap2m1',	'Arl11',	'Fkbp1a',	'Mfsd1',	'Rit1',	'Pea15a',	'Ccl9',	'Ctnnb1',	'Ptp4a3',	'Atp6v1a',	'Pepd',	'Naglu',	'Snx5',	'Rgs18',	'Cd37',	'Dennd5a',	'Ehd4',	'Tmem86a1',	'Trem2',	'Gas7',	'Lipa',	'Gpx4',	'Cysltr1',	'Lgals1',	'Atp6ap1',	'Rin2',	'Dpm3',	'Dgkz',	'Hpgds1',	'Bin1',	'Nucb1',	'Mat2a',	'Tspan4',	'As3mt1',	'Ubc',	'Pink1',	'Nfic',	'Tm9sf2',	'Mef2c1',	'Sertad1',	'Tm2d3',	'Unc93b1',	'Vat1',	'Setd3',	'Ctla2b',	'Il6ra',	'Jmjd1c1',	'Atp6v0b',	'Snx3',	'Pdcd5',	'Smap1',	'Trappc5',	'Tom1',	'Rtn41',	'Ubn1',	'Bmp2k',	'Hmgn1',	'Get4',	'Wsb11',	'Mrfap11',	'Slc43a2',	'Agfg1',	'Akirin1',	'Ppil4',	'Ppp5c',	'Rgs10',	'Acat1',	'Ctsa1',	'Ankrd13a',	'Mcfd2',	'Stom',	'Brd2',	'Cat',	'Helz',	'Atraid',	'Bri3',	'Txnip',	'Mpp11',	'Abcf1',	'Tgfbi1',	'Idh1',	'Fli1',	'Sh3bp5',	'Slc29a1',	'Snx4',	'Nenf',	'Nmd3',	'Smim11',	'Arhgap17',	'Rragc',	'Egln2',	'Ddx3x',	'Arrdc1',	'Foxp1',	'Nisch',	'Anxa5',	'Frmd4b',	'Trim47',	'Ddx5',	'Leprot',	'Ssh2',	'Hspa9',	'Gns1',	'Gnpat',	'Asah1',	'Dpp7',	'Ncaph2',	'Cenpb',	'Clec4a1',	'Spag7',	'Nagpa',	'Inpp5d',	'Jund',	'Tcf4',	'P2ry6',	'Sgpl11',	'Cebpg',	'Clk3',	'Luc7l2',	'Pkig',	'Atf4',	'Slc16a6',	'Rab3il1',	'Klf2',	'Pld3',	'Tgoln1',	'App',	'Gpr107',	'Sypl',	'Gabarapl1',	'Ost4',	'Scp2',	'Ivns1abp',	'Abhd121',	'Hacd4',	'Ifi27',	'Txndc5',	'Idh2',	'Dstn',	'Eps8',	'U2af2',	'Ring1',	'Stau1',	'Hsbp1',	'Top1',	'Fbxw4',	'Gatm',	'Tmem106a',	'Atp6ap2',	'Nsmce1',	'Rfk',	'Nsun2',	'Dusp22',	'Plbd2',	'Scarb2',	'Naa50',	'Anp32a',	'Tsc22d3',	'Otulin',	'Ndufc1',	'Nop56',	'Sptlc2',	'Slc11a1',	'Hnrnph1',	'Mafg',	'Dnajc1',	'Mef2a',	'Galnt1',	'Rabac1',	'Sh3bgrl',	'Ubl4a',	'Nr3c1',	'Elk3',	'Tm9sf3',	'Ndufa2',	'Rab11a',	'Sgpp1',	'Plin3',	'Kctd121',	'Uchl3',	'Slc25a5',	'Stat3',	'Cpne3',	'Tmx1',	'Map7d1',	'Ech1',	'Swap70',	'Cmtm6',	'Zmiz1',	'Lamtor1',	'Caml',	'Dnase2a1',	'Ykt6',	'Mrpl34',	'Igsf8',	'Smim15',	'Comt',	'Tmem50a',	'Mapk3',	'Lamtor3',	'Rab24',	'Fos1',	'Stx4a',	'Plxnb2',	'Ncbp2',	'Zfand5',	'Cndp2',	'Gna12',	'Lypla2',	'Pnrc1',	'Eif4g2',	'Mrpl27',	'Tspan3',	'Fam234a',	'Gaa',	'Epn1',	'Txndc12',	'Dynlt3',	'Syap1',	'Hist1h2bc',	'Gnpda1',	'Il10rb',	'Eif5b1',	'Pon2',	'Tnfsf12',	'Clptm1',	'Hnrnpu',	'Kras',	'Psmd1')
phagocytic <- convert_mouse_to_human(phagocytic.mouse) 
rm.P1 <- c('FCGR2C', 'SMIM11', 'ANP32C')
phagocytic <- phagocytic[!phagocytic %in% rm.P1] 
ms.int <- AddModuleScore(ms.int, features = c(phagocytic), name = 'Phagocytic_Module')

p1.feat <- do_FeaturePlot(ms.int, features = 'Phagocytic_Module1', order = T, 
                          pt.size = 0.5, split.by = 'FistulaAnnotation', ncol = 1,
                          legend.title = '', font.size = 10, legend.length = 6, legend.position = 'left')
vlndf.p1 = data.frame(Phagocytic_Module = ms.int@meta.data$Phagocytic_Module1,
                      FistulaAnnotation = ms.int$FistulaAnnotation,
                      FistulaCurrent = ms.int$FistulaCurrent)
mean.p1.fistula <- aggregate(Phagocytic_Module ~ FistulaAnnotation, vlndf.p1, mean)

p1.wilcox <- compare_means(Phagocytic_Module ~ FistulaAnnotation,  data = vlndf.p1, method = "wilcox.test")
p1.bar <- ggplot(data=vlndf.p1, aes(x=FistulaAnnotation, y=Phagocytic_Module, color=FistulaAnnotation)) +
  geom_boxplot() +
  scale_color_manual(values = c('maroon3', 'navy')) +
  theme_classic() +
  theme(axis.text = element_text(colour = "black", family = "Arial", size = 12),
        axis.title = element_text(size = 13)) +
  theme(legend.position="none") +
  xlab(label = NULL) +
  ylab(label = 'Module Score') +
  stat_pvalue_manual(
    p1.wilcox, 
    y.position = 4.5, step.increase = 0.1,
    label = "p.adj"
  )+
  geom_text(data = mean.p1.fistula, aes(label = round(Phagocytic_Module, digits = 4), y = 3.5), nudge_x = -0.22)

#p2
oxidative.mouse <- c('Fn12',	'Emilin2',	'Thbs1',	'Acly',	'Itgam',	'Laptm5',	'Gngt2',	'Ltc4s',	'Man2b1',	'Ctsd2',	'Wfdc17',	'Ssr4',	'Ccl22',	'Tcn2',	'Cd9',	'Gpx4',	'Itgb2',	'Tyrobp',	'Adgre2')
oxidative <- convert_mouse_to_human(oxidative.mouse)
ms.int <- AddModuleScore(ms.int, features = c(oxidative), name = 'Oxidative_Module')
p2.feat <- do_FeaturePlot(ms.int, features = 'Oxidative_Module1', order = T, 
                          pt.size = 0.5, split.by = 'FistulaAnnotation',ncol = 1,
                          legend.title = '', font.size = 10, legend.length = 6, legend.position = 'left')
vlndf.p2 = data.frame(Oxidative_Module = ms.int@meta.data$Oxidative_Module1,
                      FistulaAnnotation = ms.int$FistulaAnnotation,
                      FistulaCurrent = ms.int$FistulaCurrent)
mean.p2.fistula <- aggregate(Oxidative_Module ~ FistulaAnnotation, vlndf.p2, mean)
p2.wilcox <- compare_means(Oxidative_Module ~ FistulaAnnotation,  data = vlndf.p2, method = "wilcox.test")

p2.bar <- ggplot(data=vlndf.p2, aes(x=FistulaAnnotation, y=Oxidative_Module, color=FistulaAnnotation)) +
  geom_boxplot() +
  scale_color_manual(values = c('maroon3', 'navy')) +
  theme_classic() +
  theme(axis.text = element_text(colour = "black", family = "Arial", size = 12),
        axis.title = element_text(size = 13)) +
  theme(legend.position="none") +
  xlab(label = NULL) +
  ylab(label = 'Module Score') +
  stat_pvalue_manual(
    p2.wilcox, 
    y.position = 4.5, step.increase = 0.1,
    label = "p.adj"
  ) +
  geom_text(data = mean.p2.fistula, aes(label = round(Oxidative_Module, digits = 4), y = 3.65))

#p3
inflammatory.mouse <- c('Msrb1',	'Plac8',	'Ifitm3',	'Prdx5',	'Lyz2',	'Samhd1',	'Gsr',	'Tyrobp',	'Ifitm6',	'Smpdl3a',	'Lst1',	'Gngt2',	'S100a4',	'Emilin2',	'Tmsb10',	'Coro1a',	'Tpd5',	'Ly6e',	'Thbs1',	'Fxyd5',	'Napsa',	'Cybb',	'Irgm1',	'Taldo1',	'Clec4e',	'Ptpn11',	'Arpc3',	'Wtap',	'Cib1',	'Snx1',	'Fyb2',	'Snx18',	'Fcer1g',	'Klf13',	'Pkm',	'Alox5ap',	'Mrpl30',	'Pira2',	'Ifitm2',	'Oas1a',	'Msr1',	'Sf3b6',	'Ostf',	'Gyg',	'Stat1',	'Fcgr1',	'Prr13',	'Dok3',	'Sec61g',	'Btg1',	'Tnfrsf1a',	'Rbms1',	'B4galnt1',	'Tln1',	'Aprt',	'Cytip',	'Fam49b',	'Srgn',	'Ncf2',	'Ppp1r21',	'Lyn',	'Gda',	'Sirpb1c',	'Dck',	'Pla2g7',	'Cdc42se1',	'Itgam',	'Hck',	'Arhgdib',	'Sat1',	'Ndufb8',	'Ddx24',	'Actg1',	'Rtp4',	'Ube2d3',	'Mrpl4',	'Rnf149',	'Zcrb1',	'Ptp4a1',	'Echs1',	'Fgr',	'Rbm7',	'Cyba',	'Rasgrp2',	'Ppp1ca',	'Mpeg1',	'Ifi47',	'Emb',	'Gch1',	'Mpp1',	'Ptpn6',	'H3f3a',	'Tifab',	'G6pdx',	'Plaur',	'Vdac3',	'Serp1',	'Rap1b',	'Rassf3',	'Samsn1',	'Blvra',	'Smim7',	'Ube2b',	'Slfn1',	'Psma7',	'Ube2a',	'Ptbp3',	'Stk24',	'Spi1',	'Gpx1',	'Pycard',	'Eno1',	'Mxd1',	'Gtf2b',	'Clec4a1',	'Uqcrh',	'Lamtor5',	'Tor3a',	'Igsf6',	'Cdkn1a',	'Clec4a3',	'Rab8a',	'Epsti1',	'Gpr35',	'Eif3k',	'Zeb2',	'Eif3e',	'Capzb',	'Plgrkt',	'Tspan13',	'Xdh',	'Tsc22d',	'Nadk',	'Flna',	'Tpgs1',	'Sppl2a',	'Add3',	'Ccl4',	'Glud1',	'Naca',	'Smox',	'S100a11',	'Tnfrsf1b',	'Itgb2',	'Plxnb2',	'Scand1',	'Metrnl',	'Shisa5',	'Rras',	'Milr1',	'Osbpl9',	'Degs1',	'Ugp2',	'Itga4',	'Cdk2ap2',	'Xpnpep1',	'Prdx6',	'Mcemp1',	'Plin2',	'Dcaf13',	'Dusp3',	'Klf3',	'Myd88',	'Ifngr1',	'Dr1',	'Fndc3a',	'Dctn3',	'Stard3',	'Csnk1e',	'Tnfaip2',	'Adssl1',	'Cd14',	'Zfp622',	'Sra1',	'Tmem51',	'Il17ra',	'Pnrc2',	'Rap1a',	'Gsto1',	'Mafb',	'Abi3',	'Bak1',	'Stk17b')
inflammatory <- convert_mouse_to_human(inflammatory.mouse) #n =184 to n = 60, then to 194 after fixing end numbers
rm.p3 <- c('LILRA3', 'FCGR1BP', 'H3-3A', 'H3-3B', 'RAP1BL', 'GPX1', 'CCL4L1', 'ADSS1')
inflammatory <- inflammatory[!inflammatory %in% rm.p3]
ms.int <- AddModuleScore(ms.int, features = c(inflammatory), name = 'Inflammatory_Module')
p3.feat <- do_FeaturePlot(ms.int, features = 'Inflammatory_Module1', order = T, 
                          pt.size = 0.5, split.by = 'FistulaAnnotation', ncol = 1,
                          legend.title = '', font.size = 10, legend.length = 6, legend.position = 'left')
vlndf.p3 = data.frame(Inflammatory_Module = ms.int@meta.data$Inflammatory_Module1,
                      FistulaAnnotation = ms.int$FistulaAnnotation,
                      FistulaCurrent = ms.int$FistulaCurrent)
mean.p3.fistula <- aggregate(Inflammatory_Module ~ FistulaAnnotation, vlndf.p3, mean)
p3.wilcox <- compare_means(Inflammatory_Module ~ FistulaAnnotation,  data = vlndf.p3, method = "wilcox.test")

p3.bar <- ggplot(data=vlndf.p3, aes(x=FistulaAnnotation, y=Inflammatory_Module, color=FistulaAnnotation)) +
  geom_boxplot() +
  scale_color_manual(values = c('maroon3', 'navy')) +
  theme_classic() +
  theme(axis.text = element_text(colour = "black", family = "Arial", size = 12),
        axis.title = element_text(size = 13)) +
  theme(legend.position="none") +
  xlab(label = NULL) +
  ylab(label = 'Module Score') +
  stat_pvalue_manual(
    p3.wilcox, 
    y.position = 4.5, step.increase = 0.1,
    label = "p.adj"
  ) +
  geom_text(data = mean.p3.fistula, aes(label = round(Inflammatory_Module, digits = 4), y = 3.5))

phago <- p1.feat + p1.bar + plot_layout(ncol=2) 
oxi <- p2.feat + p2.bar + plot_layout(ncol=2) 
infl <- p3.feat + p3.bar + plot_layout(ncol=2) 

phago / oxi / infl

##########################
### Figures 4, S5, S6 ####
##########################

#Figure 4B
DefaultAssay(multiome) <- 'peaks'
DimPlot(multiome, reduction = 'wnn.umap', label=TRUE, repel = TRUE, label.size=8) + 
  ggplot2::theme(legend.position = 'none')

#Figure 4C
DefaultAssay(multiome) <- 'SCT'
DotPlot(multiome, features=c('KRT8', 'PIGR', 'LAMA3','COL3A1', 'FN1', 'CXCL14'), cols='RdBu', cluster.idents=TRUE) + ggplot2::theme(legend.position = 'bottom')

#Figure S5A
tibble3 <- structure(multiome)
tibble3 <- as_tibble(tibble3)
bar.pal <- do_ColorPalette(colors.use='steelblue', n=19)
set.seed(872335)
bar.pal <- sample(bar.pal)
multiome.order <- c('Perianal4', 'Perianal5', 'Perianal9', 'Perianal10',  'Perianal12', 'Perianal13')

ggplot(tibble3, aes(x = Publication_ID, fill = Anno)) +
  geom_bar(position = "fill", width = 0.98) + 
  theme_classic() +
  theme(legend.position = 'top', legend.title = element_blank(), legend.text=element_text(size=20)) +
  guides(fill=guide_legend(nrow=5, byrow=TRUE)) + #https://stackoverflow.com/questions/27130610/legend-on-bottom-two-rows-wrapped-in-ggplot2-in-r
  theme(axis.text.x = element_text(colour = "black", family = "Arial", size = 12),
        axis.text.y = element_text(colour = "black", family = "Arial", size = 12)) +
  xlab(label = NULL) +
  ylab(label = NULL) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  scale_fill_manual(values=bar.pal) +
  scale_x_discrete(limits = multiome.order) +
  coord_flip()

#Figure S5C
DefaultAssay(multiome) <- 'peaks'
multiome <- SetIdent(multiome, value='Bin')
CoveragePlot(multiome, region = c('CD8A'), group.by = 'Bin', extend.downstream=1000) #Cytotoxic T
CoveragePlot(multiome, region = c('XCL1'), group.by = 'Bin', extend.upstream=1000)  #NK/ILC
CoveragePlot(multiome, region = c('CD4'), group.by = 'Bin', extend.upstream=1000) #T 
CoveragePlot(multiome, region = c('CTLA4'), group.by = 'Bin', extend.upstream=1000) #TReg
CoveragePlot(multiome, region = c('MS4A1'), group.by = 'Bin', extend.upstream=1000) #B Cell
CoveragePlot(multiome, region = c('CD80'), group.by = 'Bin', extend.downstream=1000) #DC
CoveragePlot(multiome, region = c('CD68'), group.by = 'Bin', extend.upstream=1000, extend.downstream=500) #Myeloid
CoveragePlot(multiome, region = c('JCHAIN'), group.by = 'Bin', extend.downstream=1000) #plasma
CoveragePlot(multiome, region = c('EPCAM'), group.by = 'Bin', extend.upstream=1000) #epithelial
CoveragePlot(multiome, region = c('LUM'), group.by = 'Bin', extend.downstream=1000) #fibroblasts
CoveragePlot(multiome, region = c('ACKR1'), group.by = 'Bin') #stromal; no clear promoter either direction
CoveragePlot(multiome, region = c('KDR'), group.by = 'Bin', extend.downstream=1000) #stromal

#Figure 4F, 4G: bar plots
DefaultAssay(perianal.cd) <- 'SCT'
ap1.df <- data.frame(
  JUN = perianal.cd[['SCT']]@data["JUN",],
  JUNB = perianal.cd[['SCT']]@data["JUNB",],
  JUND = perianal.cd[['SCT']]@data["JUND",],
  FOS = perianal.cd[['SCT']]@data["FOS",],
  FOSL1 = perianal.cd[['SCT']]@data["FOSL1",],
  FOSL2 = perianal.cd[['SCT']]@data["FOSL2",],
  cluster = perianal.cd$ClusterAnnotation,
  Bin = perianal.cd$Bin,
  ancestry = perianal.cd$Ancestry
)
box.JUN <- ggplot(data = ap1.df, aes(x = Bin, y = JUN, color = Bin)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.title.y = element_text(size = 16)) +
  theme(legend.position = 'bottom') +
  xlab(label = NULL) +
  ylab(label = 'Normalized JUN Expression')
box.FOS <- ggplot(data = ap1.df, aes(x = Bin, y = FOS, color = Bin)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.title.y = element_text(size = 16)) +
  theme(legend.position = 'bottom') +
  xlab(label = NULL) +
  ylab(label = 'Normalized FOS Expression')


box.JUN + box.FOS + plot_layout(guides='collect') & theme(legend.position = 'bottom')

#Figure S6A
infl.chi3l1 <- CoveragePlot(multiome, region = c('CHI3L1'), group.by = 'Inflammation', 
extend.upstream=1000, extend.downstream=1000) & scale_fill_manual(values = c('red', 'blue'))

infl.osm <- CoveragePlot(multiome, region = c('OSM'), group.by = 'Inflammation',
extend.upstream=1000, extend.downstream=1000) & scale_fill_manual(values = c('red', 'blue'))

infl.chi3l1 / infl.osm

FindMarkers(
    object = multiome,
    features=c('CHI3L1', 'OSM'),
    ident.1 = 'Inflamed',
    group.by = 'Inflammation',
    test.use = 'LR',
    min.pct=0,
    latent.vars = 'nCount_peaks'
)

#Figures S6C, D
#MSCCR
DEGs.msccr <- read.csv('/data0/RDS_Files/Med_Submission/DEGs_MSCCR_Ancestry_AA-EA.csv', row.names=1)
custom.cols <- ifelse(
  DEGs.msccr$log2FoldChange < 0, 'darkgoldenrod1', 'lightskyblue')
names(custom.cols)[custom.cols == 'darkgoldenrod1'] <- 'EA'
names(custom.cols)[custom.cols == 'lightskyblue'] <- 'AA'
volcano.msccr <- EnhancedVolcano(DEGs.msccr,
                    lab = rownames(DEGs.msccr),
                    x = 'log2FoldChange',
                    y = 'padj',
                    xlim = c(-2, 2),
                    ylim=c(0,5),
                    pCutoff = 0.05,
                    pointSize = 1,
                    #col = c('black', 'orange', 'navy', 'magenta'),
                    colCustom = custom.cols,
                    axisLabSize = 11,
                    labSize = 4.5,
                    ylab = bquote(~-Log[10]~italic(Padj)),
                    #ylab = NULL,
                    #xlab = NULL,
                    title = 'AA vs. EA',
                    subtitle = 'Bulk RNAseq',
                    selectLab = c('JUN', 'JUNB', 'JUND', 'FOS'),
                    boxedLabels = TRUE,
                    drawConnectors = TRUE,
                    colConnectors = 'black',
                    widthConnectors = 0.5,
                    legendPosition = 'none') + ggplot2::labs(title = NULL, caption = NULL)

#scRNA
#all bx
DEGs.AAvsEA <- read.csv(file = '/data0/RDS_Files/Med_Submission/DEGs_scRNA_AAvsEA_All.csv', row.names=1)
custom.cols.all <- ifelse(
  DEGs.AAvsEA$avg_log2FC < 0, 'darkgoldenrod1', 'lightskyblue')
names(custom.cols.all)[custom.cols.all == 'darkgoldenrod1'] <- 'EA'
names(custom.cols.all)[custom.cols.all == 'lightskyblue'] <- 'AA'

#infl
DEGs.AAvsEA.infl <- read.csv(file = '/data0/RDS_Files/Med_Submission/DEGs_scRNA_AAvsEA_Infl.csv', row.names=1)
custom.cols.infl <- ifelse(
  DEGs.AAvsEA.infl$avg_log2FC < 0, 'darkgoldenrod1', 'lightskyblue')
names(custom.cols.infl)[custom.cols.infl == 'darkgoldenrod1'] <- 'EA'
names(custom.cols.infl)[custom.cols.infl == 'lightskyblue'] <- 'AA'

#non
DEGs.AAvsEA.non <- read.csv(file = '/data0/RDS_Files/Med_Submission/DEGs_scRNA_AAvsEA_Non.csv', row.names=1)
custom.cols.non <- ifelse(
  DEGs.AAvsEA.non$avg_log2FC < 0, 'darkgoldenrod1', 'lightskyblue')
names(custom.cols.non)[custom.cols.non == 'darkgoldenrod1'] <- 'EA'
names(custom.cols.non)[custom.cols.non == 'lightskyblue'] <- 'AA'

volcano.all <- EnhancedVolcano(DEGs.AAvsEA,
                                    lab = rownames(DEGs.AAvsEA),
                                    x = 'avg_log2FC',
                                    y = 'p_val_adj',
                                    xlim = c(-2,2),
                                    ylim = c(0, 350),
                                    pCutoff = 5e-50,
                                    pointSize = 1,
                                    #col = c('black', 'orange', 'navy', 'magenta'),
                                    colCustom = custom.cols.all,
                                    axisLabSize = 11,
                                    labSize = 4.5,
                                    ylab = bquote(~-Log[10]~italic(Padj)),
                                    #ylab = NULL,
                                    #xlab = NULL,
                                    title = 'AA vs. EA',
                                    subtitle = 'All scRNAseq',
                                    selectLab = c('JUN', 'JUNB', 'JUND', 'FOS'),
                                    boxedLabels = TRUE,
                                    drawConnectors = TRUE,
                                    colConnectors = 'black',
                                    widthConnectors = 0.5,
                                    legendPosition = 'none') + ggplot2::labs(title = NULL, caption = NULL)

volcano.infl <- EnhancedVolcano(DEGs.AAvsEA.infl,
                                    lab = rownames(DEGs.AAvsEA.infl),
                                    x = 'avg_log2FC',
                                    y = 'p_val_adj',
                                    xlim = c(-2,2),
                                    ylim = c(0, 350),
                                    pCutoff = 5e-50,
                                    pointSize = 1,
                                    #col = c('black', 'orange', 'navy', 'magenta'),
                                    colCustom = custom.cols.infl,
                                    axisLabSize = 11,
                                    labSize = 4.5,
                                    ylab = bquote(~-Log[10]~italic(Padj)),
                                    #ylab = NULL,
                                    #xlab = NULL,
                                    title = 'AA vs. EA',
                                    subtitle = 'Inflamed scRNAseq',
                                    selectLab = c('JUN', 'JUNB', 'JUND', 'FOS'),
                                    boxedLabels = TRUE,
                                    drawConnectors = TRUE,
                                    colConnectors = 'black',
                                    widthConnectors = 0.5,
                                    legendPosition = 'none') + ggplot2::labs(title = NULL, caption = NULL)

volcano.non <- EnhancedVolcano(DEGs.AAvsEA.non,
                                    lab = rownames(DEGs.AAvsEA.non),
                                    x = 'avg_log2FC',
                                    y = 'p_val_adj',
                                    xlim = c(-2,2),
                                    ylim = c(0, 350),
                                    pCutoff = 5e-50,
                                    pointSize = 1,
                                    #col = c('black', 'orange', 'navy', 'magenta'),
                                    colCustom = custom.cols.non,
                                    axisLabSize = 11,
                                    labSize = 4.5,
                                    ylab = bquote(~-Log[10]~italic(Padj)),
                                    #ylab = NULL,
                                    #xlab = NULL,
                                    title = 'AA vs. EA',
                                    subtitle = 'Non inflamed - scRNAseq',
                                    selectLab = c('JUN', 'JUNB', 'JUND', 'FOS'),
                                    boxedLabels = TRUE,
                                    drawConnectors = TRUE,
                                    colConnectors = 'black',
                                    widthConnectors = 0.5,
                                    legendPosition = 'none') + ggplot2::labs(title = NULL, caption = NULL)

(volcano.all + volcano.msccr) / (volcano.infl + volcano.non)

##########################
### Figures 5, 6, S7 #####
##########################

#Figure 5D
library(ComplexUpset)
input <- read.csv('/data0/RDS_Files/Med_Submission/UpSet_input.csv')
input$IBD_Locus <- as.factor(input$IBD_Locus)

upset(
   input, intersect = c('B_Cells', 'Plasma_Cells', 'T_Cells', 'CytoT_Innate', 'Myeloid', 'Fibroblasts', 'Stromal', 'Epithelial'), 
   name = 'Nucleus Cluster Bin', min_size = 4, 
   base_annotations = list(
      'Intersection size'=intersection_size(
         mapping = aes(fill=IBD_Locus),
         text_colors = c(on_background='black', on_bar='white'))
         + scale_fill_manual(values=c('0' = 'darkgrey', '1' = 'turquoise2'))
         + ylab('Number of TF Footprints per Intersection')
      ),
   set_sizes = (
      upset_set_size(
         geom = geom_bar(
            aes(fill=IBD_Locus)), position = 'left'
      ) 
      + scale_fill_manual(values=c('0' = 'darkgrey', '1' = 'turquoise2'))
      + ylab('Number of TF Footprints per Set')
   ), guides='over'
)

#Figures 6B, 6E
FeaturePlot(perianal.CD, features=c('PTGER4', 'LACC1'), order=TRUE, cols=c('darkgrey', 'red'))

#Figure S7C
CoveragePlot(multiome, region = c('chr5-40409000-40411600')) # PTGER4 SNP locus = chr5:40410482

#Figure S7D
CoveragePlot(multiome, region = 'LACC1') # LACC1 SNP locus = chr13:43883789