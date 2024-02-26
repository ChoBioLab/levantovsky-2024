#code to generate figures in Levantovsky 2023
#revised: 09/06/2023; 11/28/23

##########################
####### Libraries ########
##########################
library(future)
plan(
  multicore,
  workers = 8
) # parallelization
options(future.globals.maxSize = 12 * 1024^2 * 1000)
library(dplyr)
library(Seurat)
library(harmony)
library(ggplot2)
library(pals)
library(patchwork)
library(tidyseurat)
library(Scillus)
library(ggrepel)
options(ggrepel.max.overlaps = Inf) #https://ggrepel.slowkow.com/articles/examples.html
library(EnhancedVolcano)
library(ggpubr)
library(DESeq2)
library(SCpubr)
library(viridis)
library(Nebulosa)
library(RColorBrewer)
library(Signac)
library(BSgenome.Hsapiens.UCSC.hg38)

##########################
####### RDS Files ########
##########################
#scRNAseq
CombinedCD <- readRDS('/data0/RDS_Files/Med_Resubmission/CombinedCD.RDS')
fistula.clusters <- readRDS('/data0/RDS_Files/Med_Resubmission/Fistula_n-3_Anno.rds')
PerianalCD <- readRDS('/data0/RDS_Files/Med_Resubmission/PerianalCD_n-15_Anno.rds')
ms.int <- readRDS('/data0/RDS_Files/Med_Resubmission/MS-int_Anno.rds')
ms.int <- readRDS('/data0/RDS_Files/2023-11-26_MS-int_Anno.rds')

ms.int_abbr <- RenameIdents(ms.int, 
  "T Cells" = "T Cells",
  "B Cells" = "B Cells", 
  "Macrophages" = "Macrophages", 
  "cDC2" = "cDC2", 
  "mregDCs" = "mregDCs", 
  "Monocytes" = "Monocytes", 
  "Monocytes/Macrophages" = "MonoMacs", 
  "Pericytes" = "Pericytes", 
  "Enteric Neurons"= "Enteric Neurons", 
  "Myofibroblasts" = "Myofibroblasts", 
  "Endothelial" = "Endothelial",
  "Fibroblastic Reticular Cells" = "FRCs", 
  "CHI3L1hi Fibroblasts" = "CHI3L1hi Fibro.", 
  "PDGFRAhi Fibroblasts" = "PDGFRAhi Fibro.", 
  "ADAMDEC1hi Fibroblasts" = "ADAMDEC1hi Fibro.", 
  "CCL11hi Fibroblasts" = "CCL11hi Fibro.")
fibroblasts <- readRDS('/data0/RDS_Files/Med_Resubmission/FibroblastSubset_UnAnno.rds')

#bulk RNA seq
msccr.rectum <- readRDS('/data0/RDS_Files/Med_Resubmission/MSCCR_dds_Rectum_Inf-Non.rds')
msccr.ancestry <- readRDS('/data0/RDS_Files/Med_Resubmission/MSCCR_dds_Ancestry.rds') 

#multiome
multiome <- readRDS('/data0/RDS_Files/Med_Resubmission/multiome_annotated.RDS')

#upsetPlot
input <- readRDS('/data0/RDS_Files/Med_Resubmission/UpSet-TOBIAS-Input.csv')

##########################
####### Gene Lists #######
##########################
saul_sen_mayo <- c('ACVR1B',	'ANG',	'ANGPT1',	'ANGPTL4',	'AREG',	'AXL',	'BEX3',	'BMP2',	'BMP6',	'C3',	'CCL1',	'CCL13',	'CCL16',	'CCL2',	'CCL20',	'CCL24',	'CCL26',	'CCL3',	'CCL3L1',	'CCL4',	'CCL5',	'CCL7',	'CCL8',	'CD55',	'CD9',	'CSF1',	'CSF2',	'CSF2RB',	'CST4',	'CTNNB1',	'CTSB',	'CXCL1',	'CXCL10',	'CXCL12',	'CXCL16',	'CXCL2',	'CXCL3',	'CXCL8',	'CXCR2',	'DKK1',	'EDN1',	'EGF',	'EGFR',	'EREG',	'ESM1',	'ETS2',	'FAS',	'FGF1',	'FGF2',	'FGF7',	'GDF15',	'GEM',	'GMFG',	'HGF',	'HMGB1',	'ICAM1',	'ICAM3',	'IGF1',	'IGFBP1',	'IGFBP2',	'IGFBP3',	'IGFBP4',	'IGFBP5',	'IGFBP6',	'IGFBP7',	'IL10',	'IL13',	'IL15',	'IL18',	'IL1A',	'IL1B',	'IL2',	'IL32',	'IL6',	'IL6ST',	'IL7',	'INHA',	'IQGAP2',	'ITGA2',	'ITPKA',	'JUN',	'KITLG',	'LCP1',	'MIF',	'MMP1',	'MMP10',	'MMP12',	'MMP13',	'MMP14',	'MMP2',	'MMP3',	'MMP9',	'NAP1L4',	'NRG1',	'PAPPA',	'PECAM1',	'PGF',	'PIGF',	'PLAT',	'PLAU',	'PLAUR',	'PTBP1',	'PTGER2',	'PTGES',	'RPS6KA5',	'SCAMP4',	'SELPLG',	'SEMA3F',	'SERPINB4',	'SERPINE1',	'SERPINE2',	'SPP1',	'SPX',	'TIMP2',	'TNF',	'TNFRSF10C',	'TNFRSF11B',	'TNFRSF1A',	'TNFRSF1B',	'TUBGCP2',	'VEGFA',	'VEGFC',	'VGF',	'WNT16',	'WNT2')
LR_pairs <- c('FTL', 'FTH1', 'SCARA5', 'CCL2', 'CXCL8', 'NR3C1', 'ACKR1','CSF1', 'CSF3', 'CSF1R', 'TNF', 'TNFRSF1A', 'CXCL14', 'CXCR4', 'PLAU', 'PLAUR') #related to figure 2D

#fibroblast modules
fibrotic.genes <- c('COL1A1', 'ACTA2', 'CCL5', 'EDN1', 'VCAM1', 'CXCL8', 'ENPP2', 'CTGF', 'COL3A1', 'THBS1', 'CALU', 'ITGAV', 'RECK', 'COL1A2')
destructive.genes <- c('FAP', 'PDPN', 'THY1', 'MMP3', 'MMP9', 'MMP13', 'IL18', 'CCL9', 'TNFSF11', 'CXCL1', 'PRG4', 'CLIC5', 'TSPAN15', 'COL22A1')
enrichment.genes <- list("Fibrotic"=fibrotic.genes, "Destructive"=destructive.genes)

#macrophage activation paths
phagocytic <- c('LTC4S',
  'C1QA',
  'CTSB',
  'C1QC',
  'APOE',
  'C1QB',
  'CSF1R',
  'DAB2',
  'GAS6',
  'NINJ1',
  'FCGRT',
  'TIMP2',
  'LGMN',
  'PLTP',
  'ITM2B',
  'BLVRB',
  'SERINC3',
  'RNASE4',
  'GRN',
  'TGFBR2',
  'C3AR1',
  'NRP1',
  'TF',
  'CFH',
  'CD63',
  'C5AR1',
  'ANG',
  'DHRS3',
  'CTSL',
  'CTSV',
  'GLUL',
  'PLA2G15',
  'EMP1',
  'CD33',
  'SIGLEC6',
  'RHOB',
  'EGR1',
  'LAMP1',
  'PTPN18',
  'CCL15',
  'CCL23',
  'APLP2',
  'ATP13A2',
  'IER3',
  'MS4A7',
  'PTGER4',
  'ALDH2',
  'MAFB',
  'HEXA',
  'FTL',
  'ADGRE1',
  'ZFP36L1',
  'ST6GALNAC4',
  'TPP1',
  'PLK2',
  'TCN2',
  'C16orf54',
  'FCGR2A',
  'FCGR2B',
  'HMOX1',
  'MAN2B1',
  'CAMK1',
  'MARCKS',
  'PIGK',
  'CLTA',
  'TMEM37',
  'ADAM15',
  'TCF3',
  'PRKCB',
  'TBC1D17',
  'CD81',
  'VPS33A',
  'CCL13',
  'CCL2',
  'GRAMD1A',
  'SNX2',
  'GALK2',
  'CANX',
  'TBXAS1',
  'AP2A2',
  'DUSP6',
  'SNX6',
  'GOLPH3',
  'CALM2',
  'PLOD1',
  'TRIB1',
  'FAM3C',
  'ALDH9A1',
  'GNL3',
  'CD68',
  'PTPRJ',
  'RAB31',
  'DOK3',
  'SEC14L1',
  'MAN1A1',
  'DNAJB9',
  'STARD5',
  'CST3',
  'GGNBP2',
  'QRICH1',
  'CLN8',
  'AP2M1',
  'ARL11',
  'FKBP1A',
  'FKBP1C',
  'MFSD1',
  'RIT1',
  'PEA15',
  'CCL15',
  'CCL23',
  'CTNNB1',
  'PTP4A3',
  'ATP6V1A',
  'PEPD',
  'NAGLU',
  'SNX5',
  'RGS18',
  'CD37',
  'DENND5A',
  'EHD4',
  'TREM2',
  'GAS7',
  'LIPA',
  'GPX4',
  'CYSLTR1',
  'LGALS1',
  'ATP6AP1',
  'RIN2',
  'DPM3',
  'DGKZ',
  'BIN1',
  'NUCB1',
  'MAT2A',
  'TSPAN4',
  'UBC',
  'PINK1',
  'NFIC',
  'TM9SF2',
  'SERTAD1',
  'TM2D3',
  'UNC93B1',
  'VAT1',
  'SETD3',
  'IL6R',
  'ATP6V0B',
  'SNX3',
  'PDCD5',
  'SMAP1',
  'TRAPPC5',
  'TOM1',
  'UBN1',
  'BMP2K',
  'HMGN1',
  'GET4',
  'SLC43A2',
  'AGFG1',
  'AKIRIN1',
  'PPIL4',
  'PPP5C',
  'RGS10',
  'ACAT1',
  'ANKRD13A',
  'MCFD2',
  'STOM',
  'BRD2',
  'CAT',
  'HELZ',
  'ATRAID',
  'BRI3',
  'TXNIP',
  'ABCF1',
  'IDH1',
  'FLI1',
  'SH3BP5',
  'SLC29A1',
  'SNX4',
  'NENF',
  'NMD3',
  'ARHGAP17',
  'RRAGC',
  'EGLN2',
  'DDX3X',
  'ARRDC1',
  'FOXP1',
  'NISCH',
  'ANXA5',
  'FRMD4B',
  'TRIM47',
  'DDX5',
  'LEPROT',
  'SSH2',
  'HSPA9',
  'GNPAT',
  'ASAH1',
  'DPP7',
  'NCAPH2',
  'CENPB',
  'CLEC4A',
  'ZNF705A',
  'SPAG7',
  'NAGPA',
  'INPP5D',
  'JUND',
  'TCF4',
  'P2RY6',
  'CEBPG',
  'CLK3',
  'LUC7L2',
  'PKIG',
  'ATF4',
  'SLC16A6',
  'RAB3IL1',
  'KLF2',
  'PLD3',
  'TGOLN2',
  'APP',
  'GPR107',
  'SYPL1',
  'GABARAPL1',
  'OST4',
  'SCP2',
  'IVNS1ABP',
  'HACD4',
  'IFI27',
  'IFI27L1',
  'IFI27L2',
  'TXNDC5',
  'IDH2',
  'DSTN',
  'EPS8',
  'U2AF2',
  'RING1',
  'STAU1',
  'HSBP1',
  'TOP1',
  'FBXW4',
  'GATM',
  'TMEM106A',
  'ATP6AP2',
  'NSMCE1',
  'RFK',
  'NSUN2',
  'DUSP22',
  'PLBD2',
  'SCARB2',
  'NAA50',
  'ANP32A',
  'TSC22D3',
  'OTULIN',
  'NDUFC1',
  'NOP56',
  'SPTLC2',
  'SLC11A1',
  'HNRNPH1',
  'MAFG',
  'DNAJC1',
  'MEF2A',
  'GALNT1',
  'RABAC1',
  'SH3BGRL',
  'UBL4A',
  'NR3C1',
  'ELK3',
  'TM9SF3',
  'NDUFA2',
  'RAB11A',
  'SGPP1',
  'PLIN3',
  'UCHL3',
  'SLC25A5',
  'STAT3',
  'CPNE3',
  'TMX1',
  'MAP7D1',
  'ECH1',
  'SWAP70',
  'CMTM6',
  'ZMIZ1',
  'LAMTOR1',
  'CAMLG',
  'YKT6',
  'MRPL34',
  'IGSF8',
  'SMIM15',
  'COMT',
  'TMEM50A',
  'MAPK3',
  'LAMTOR3',
  'RAB24',
  'STX4',
  'PLXNB2',
  'NCBP2',
  'NCBP2L',
  'ZFAND5',
  'CNDP2',
  'GNA12',
  'LYPLA2',
  'PNRC1',
  'EIF4G2',
  'MRPL27',
  'TSPAN3',
  'FAM234A',
  'GAA',
  'EPN1',
  'TXNDC12',
  'DYNLT3',
  'SYAP1',
  'GNPDA1',
  'IL10RB',
  'PON2',
  'TNFSF12',
  'TNFSF12-TNFSF13',
  'CLPTM1',
  'HNRNPU',
  'KRAS',
  'PSMD1')
oxidative <- c('EMILIN2',
  'THBS1',
  'ACLY',
  'ITGAM',
  'LAPTM5',
  'GNGT2',
  'LTC4S',
  'MAN2B1',
  'SSR4',
  'CCL22',
  'TCN2',
  'CD9',
  'GPX4',
  'ITGB2',
  'TYROBP')
inflammatory <- c('MSRB1',
  'PLAC8',
  'IFITM3',
  'PRDX5',
  'LYZ',
  'SAMHD1',
  'GSR',
  'TYROBP',
  'SMPDL3A',
  'LST1',
  'GNGT2',
  'S100A4',
  'EMILIN2',
  'TMSB10',
  'CORO1A',
  'LY6E',
  'THBS1',
  'FXYD5',
  'NAPSA',
  'CYBB',
  'IRGM',
  'TALDO1',
  'CLEC4E',
  'PTPN11',
  'ARPC3',
  'WTAP',
  'CIB1',
  'SNX1',
  'FYB2',
  'SNX18',
  'FCER1G',
  'KLF13',
  'PKM',
  'ALOX5AP',
  'MRPL30',
  'LILRA1',
  'LILRA2',
  'LILRA4',
  'LILRA6',
  'LILRB3',
  'IFITM2',
  'IFITM3',
  'OAS1',
  'MSR1',
  'SF3B6',
  'GYG1',
  'STAT1',
  'FCGR1A',
  'PRR13',
  'DOK3',
  'SEC61G',
  'BTG1',
  'TNFRSF1A',
  'RBMS1',
  'B4GALNT1',
  'TLN1',
  'APRT',
  'CYTIP',
  'SRGN',
  'NCF2',
  'PPP1R21',
  'LYN',
  'GDA',
  'SIRPB1',
  'SIRPG',
  'DCK',
  'PLA2G7',
  'CDC42SE1',
  'ITGAM',
  'HCK',
  'ARHGDIB',
  'SAT1',
  'NDUFB8',
  'DDX24',
  'ACTG1',
  'RTP4',
  'UBE2D3',
  'MRPL4',
  'RNF149',
  'ZCRB1',
  'PTP4A1',
  'ECHS1',
  'FGR',
  'RBM7',
  'CYBA',
  'RASGRP2',
  'PPP1CA',
  'MPEG1',
  'EMB',
  'GCH1',
  'MPP1',
  'PTPN6',
  'TIFAB',
  'G6PD',
  'PLAUR',
  'VDAC3',
  'SERP1',
  'RAP1B',
  'RASSF3',
  'SAMSN1',
  'BLVRA',
  'SMIM7',
  'UBE2B',
  'SLFN12',
  'SLFN12L',
  'PSMA7',
  'UBE2A',
  'PTBP3',
  'STK24',
  'SPI1',
  'PYCARD',
  'ENO1',
  'MXD1',
  'GTF2B',
  'CLEC4A',
  'ZNF705A',
  'UQCRH',
  'UQCRHL',
  'LAMTOR5',
  'TOR3A',
  'IGSF6',
  'CDKN1A',
  'CLEC4A',
  'RAB8A',
  'EPSTI1',
  'GPR35',
  'EIF3K',
  'ZEB2',
  'EIF3E',
  'CAPZB',
  'PLGRKT',
  'TSPAN13',
  'XDH',
  'NADK',
  'FLNA',
  'TPGS1',
  'SPPL2A',
  'ADD3',
  'CCL4',
  'CCL4L2',
  'GLUD1',
  'GLUD2',
  'NACA',
  'SMOX',
  'S100A11',
  'TNFRSF1B',
  'ITGB2',
  'PLXNB2',
  'SCAND1',
  'METRNL',
  'SHISA5',
  'RRAS',
  'MILR1',
  'OSBPL9',
  'DEGS1',
  'UGP2',
  'ITGA4',
  'CDK2AP2',
  'XPNPEP1',
  'PRDX6',
  'MCEMP1',
  'PLIN2',
  'DCAF13',
  'DUSP3',
  'KLF3',
  'MYD88',
  'IFNGR1',
  'DR1',
  'FNDC3A',
  'DCTN3',
  'STARD3',
  'CSNK1E',
  'TPTEP2-CSNK1E',
  'TNFAIP2',
  'CD14',
  'ZNF622',
  'SRA1',
  'TMEM51',
  'IL17RA',
  'PNRC2',
  'RAP1A',
  'GSTO1',
  'MAFB',
  'ABI3',
  'BAK1',
  'STK17B')
#M2
m2.markers <- c('CD163', 'CD200R1', 'CLEC10A', 'CXCR2', 'MRC1', #surface markers, MRC1 = CD206, CD163 = mac/mono marker
                'PPARG', 'STAT6',  #transcription factors
                'CCL2', 'TGFB1', 'IL10', 'IL1RN', 'CCL14', 'CCL17', 'CCL18', 'CCL22', 'CCL23') #secreted

##########################
### Figures 1, S1, S2 ####
##########################
#Figure 1C
DimPlot(fistula.clusters, label = FALSE, cols='alphabet') +
  ggplot2::labs(title = NULL) + ggplot2::theme(legend.position = 'left', axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) + 
  xlab(NULL) + ylab(NULL)

#Figure S1B
DotPlot(fistula.clusters, features = c('TGFB1', 'SNAI1', 'SNAI2', 'ETS1', 'MMP3', 'MMP9', 'MMP13', 'TNF', 'TNFRSF1A', 'IL13', 'IL13RA1'), 
  idents = c('Fibroblasts', 'Stromal Epithelial Mix', 'Endothelial', 'Monocytes, Macrophages'), cols = c('blue', 'red'), split.by = 'Region') + xlab(NULL) + ylab(NULL) + theme(legend.position='top')

#figure 1D
table(fistula.clusters$Region)
table(fistula.clusters$FistulaClustersAnno, fistula.clusters$Region)
#dendritic cells
DC.prop <- prop.test(x=c(251, 62), n=c(15580,20865)) 
#monocytes/macrophages
MM.prop <- prop.test(x=c(1166, 349), n=c(15580,20865)) 

#Re: Figure 1 DEG of fistula for IPA
#DEG.FvR <- FindMarkers(fistula.clusters, ident.1 = 'Fistula', group.by = 'Region', logfc.threshold = 0)
#write.csv(DEG.FvR, file='/data0/2023-09-14_FistulaClusters_DEGs.csv')

#Figure S1C
#feat.FTL <- FeaturePlot(fistula.clusters, features = c('FTL'), order=TRUE)
feat.FTL <- plot_density(fistula.clusters, features = c('FTL'))
#feat.SCARA5 <- FeaturePlot(fistula.clusters, features=c('SCARA5'), order=TRUE)
feat.SCARA5 <- plot_density(fistula.clusters, features=c('SCARA5'))
feat.FTL + feat.SCARA5 + plot_layout(ncol=2)

#Figure S1D 
DotPlot(fistula.clusters, features = senescence, split.by = 'Region',cols='RdBu')
scale.fistula <- rownames(fistula.clusters)
fistula.clusters <- ScaleData(fistula.clusters, features = scale.fistula)
fistula.clusters.pal <- names(pals::alphabet(n=17)) #colors used below w substitutions

plot_heatmap(
  dataset = fistula.clusters,
  markers = saul_sen_mayo,
  sort_var = c('FistulaClustersAnno', 'Region'),
  anno_var = c('FistulaClustersAnno', 'Region'),
  anno_colors = list(c('plum1', 'blue', 'tan4', 'purple4', 'gray0', 'seagreen4', 'springgreen3', 'bisque', 
                      'gray', 'darkseagreen1', 'goldenrod4', 'olivedrab1', 'magenta', 'dodgerblue4', 'orange', 'pink', '#426600'),
                      c('maroon3', 'navy')
                      )
)

#Figure 1H
DimPlot(PerianalCD, label = FALSE, cols='alphabet2') +
  ggplot2::labs(title = NULL) + ggplot2::theme(legend.position = 'left', axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) + 
  xlab(NULL) + ylab(NULL)

#Figure S2C
PerianalCD <- SetIdent(PerianalCD, value = 'PubID_Obj')
patient.list.1 <- c(
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
  "Perianal14_Fistula", "Perianal14_Rect",
  "Perianal16_Fistula", 'Perianal16_Rect'
)
PerianalCD$PubID_Obj <- factor(PerianalCD$PubID_Obj, levels = patient.list.1)

qc.count <- VlnPlot(PerianalCD, features = c("nCount_RNA"), pt.size = 0) + 
  ggplot2::theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + xlab(NULL) +
  ggplot2::theme(legend.position = 'bottom', legend.justification = 'center', legend.key.size = unit(3, 'mm')) 
qc.feature <- VlnPlot(PerianalCD, features = c("nFeature_RNA"), pt.size = 0) +
  ggplot2::theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + xlab(NULL) +
  ggplot2::theme(legend.position = 'bottom', legend.justification = 'center', legend.key.size = unit(3, 'mm')) 
qc.sct <- VlnPlot(PerianalCD, features = c("nCount_SCT"), pt.size = 0) +
  ggplot2::theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + xlab(NULL) +
  ggplot2::theme(legend.position = 'bottom', legend.justification = 'center', legend.key.size = unit(3, 'mm')) 
qc.mt <- VlnPlot(PerianalCD, features = c("percent.mt"), pt.size = 0) +
  ggplot2::theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + xlab(NULL) +
  ggplot2::theme(legend.position = 'bottom', legend.justification = 'center', legend.key.size = unit(3, 'mm')) 

qc.count + qc.feature + qc.sct + qc.mt + plot_layout(ncol=4, guides = "collect")  & theme(legend.position = 'bottom')

#Figure S2D
tibble.PerianalCD <- structure(PerianalCD)
tibble.PerianalCD <- as_tibble(tibble.PerianalCD)

proportion.perianal <- ggplot(tibble.PerianalCD, aes(x = PubID_Obj, fill = PerianalCDAnno)) +
  geom_bar(position = "fill", width = 0.98) + 
  theme_classic() +
  theme(legend.position = 'right', legend.title = element_blank()) +
  theme(axis.text.x = element_text(colour = "black", family = "Arial", size = 12, angle=90),
        axis.text.y = element_text(colour = "black", family = "Arial", size = 12)) +
  xlab(label = NULL) +
  ylab(label = NULL) +
  scale_fill_manual(values=as.vector(alphabet2(14))) + #https://github.com/kwstat/pals/issues/3
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0))  +
  scale_x_discrete(limits = patient.list.1)
proportion.perianal

#Figure S2E
DimPlot(CombinedCD, label = FALSE, cols='alphabet') +
  ggplot2::labs(title = NULL) + ggplot2::theme(legend.position = 'right', axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) + 
  xlab(NULL) + ylab(NULL)

#Figure S2F
DimPlot(CombinedCD, group.by = 'orig.ident', cols=c('firebrick', 'deepskyblue')) + ggplot2::theme(legend.position = 'top', legend.justification = 'center', axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) + 
  ggtitle(label=NULL) + xlab(NULL) + ylab(NULL)

#Figure S2G
table(CombinedCD$orig.ident, CombinedCD$ClusterAnnotation)
#proportion z-test: c(n cell type Martin, n cell type PerianalCD), c(n total Martin, n total PerianalCD)

#iga
iga.prop <- prop.test(x=c(10442, 42415), n=c(86126,142708)) #proportions z-test
iga.vln <- VlnPlot(CombinedCD, features = c('IGHA1'), split.by = 'orig.ident', split.plot = TRUE,
        idents = c('IgA Plasma'), pt.size=0, cols = c('firebrick', 'deepskyblue')) +
        ggplot2::theme(legend.position = 'none', axis.text.x = element_blank(), axis.ticks.x = element_blank()) + xlab('IgA Plasma')
#igg
igg.prop <- prop.test(x=c(2642, 8643), n=c(86126,142708)) #proportions z-test
igg.vln <- VlnPlot(CombinedCD, features = c('IGHG1'), split.by = 'orig.ident', split.plot = TRUE,
        idents = c('IgG Plasma'), pt.size=0, cols = c('firebrick', 'deepskyblue')) +
        ggplot2::theme(legend.position = 'none', axis.text.x = element_blank(), axis.ticks.x = element_blank()) + xlab('IgG Plasma')

#fibroblast
fb.prop <- prop.test(x=c(2089, 8409), n=c(86126,142708)) #proportions z-test
fibro.vln <- VlnPlot(CombinedCD, features = c('COL1A1'), split.by = 'orig.ident', split.plot = TRUE,
        idents = c('Fibroblasts'), pt.size=0, cols = c('firebrick', 'deepskyblue')) +
        ggplot2::theme(legend.position = 'none', axis.text.x = element_blank(), axis.ticks.x = element_blank()) + xlab('Fibroblasts')

#cd4
cd4.prop <- prop.test(x=c(29147, 27277), n=c(86126,142708)) #proportions z-test
cd4.vln <- VlnPlot(CombinedCD, features = c('IL7R'), split.by = 'orig.ident', split.plot = TRUE,
        idents = c('CD4 T Cells'), pt.size=0, cols = c('firebrick', 'deepskyblue')) +
        ggplot2::theme(legend.position = 'none', axis.text.x = element_blank(), axis.ticks.x = element_blank()) + xlab('CD4 T')

#cd8
cd8.prop <- prop.test(x=c(9705, 8985), n=c(86126,142708)) #proportions z-test
cd8.vln <- VlnPlot(CombinedCD, features = c('CD8A'), split.by = 'orig.ident', split.plot = TRUE,
        idents = c('CD8 T Cells'), pt.size=0, cols = c('firebrick', 'deepskyblue')) +
        ggplot2::theme(legend.position = 'none', axis.text.x = element_blank(), axis.ticks.x = element_blank()) + xlab('CD8 T')

#cytotoxic cd8/nk
cyto.prop <- prop.test(x=c(4979, 4336), n=c(86126,142708)) #proportions z-test
cytotoxic.vln <- VlnPlot(CombinedCD, features = c('GZMA'), split.by = 'orig.ident', split.plot = TRUE,
        idents = c('Cytotoxic CD8 T, NK'), pt.size=0, cols = c('firebrick', 'deepskyblue')) +
        ggplot2::theme(legend.position = 'none', axis.text.x = element_blank(), axis.ticks.x = element_blank()) + xlab('Cytotoxic CD8/NK')

iga.vln + fibro.vln + igg.vln + cytotoxic.vln + cd8.vln + cd4.vln + plot_layout(ncol=6)

##########################
#### Figures 2, S3 #######
##########################
#Figure 2A
DimPlot(ms.int, label = FALSE, cols='polychrome') +
  ggplot2::labs(title = NULL) + ggplot2::theme(legend.position = 'right', axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) + 
  xlab(NULL) + ylab(NULL)

#Figure 2B
DimPlot(ms.int, label = FALSE, group.by = 'FistulaAnnotation', cols = c('maroon3', 'navy')) +
  ggplot2::labs(title = NULL) + ggplot2::theme(legend.position = 'bottom', legend.justification = "center", axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) + 
  xlab(NULL) + ylab(NULL)

#Figure 2D
DotPlot(ms.int_abbr, features = LR_pairs, split.by = 'FistulaAnnotation', cols='RdBu') + 
  ggplot2::theme(axis.text.x = element_text(angle=45, hjust=1 ,vjust=1)) +
  xlab(NULL) + ylab(NULL)

#Figure 2E
stromal.subset <- c('ADAMDEC1hi Fibroblasts', 'PDGFRAhi Fibroblasts', 'Endothelial',
  'Myofibroblasts', 'Pericytes', 'Fibroblastic Reticular Cells', 'CHI3L1hi Fibroblasts',
  'Enteric Neurons', 'CCL11hi Fibroblasts')
stromal.int <- subset(ms.int, idents = stromal.subset)
DEG.stromal <- FindMarkers(stromal.int, ident.1 = 'Fistula', group.by = 'FistulaAnnotation', logfc.threshold = 0)
write.csv(DEG.stromal, file = '/data0/2023-11-14_DEGstromal_n-15.csv')

volcano.cols.s <- ifelse(
  DEG.stromal$avg_log2FC < 0, 'navy', 'maroon3')
names(volcano.cols.s)[volcano.cols.s == 'navy'] <- 'Rectum'
names(volcano.cols.s)[volcano.cols.s == 'maroon3'] <- 'Fistula'

DEG.stromal %>%
  top_n(n = 10, wt = avg_log2FC) -> top10.stromal.up
DEG.stromal %>%
  top_n(n = -10, wt = avg_log2FC) -> top10.stromal.down
stromal.list1 <- rownames(top10.stromal.up)
stromal.list2 <- rownames(top10.stromal.down)
stromal.list.volcano <- c(stromal.list1, stromal.list2)

EnhancedVolcano(DEG.stromal,
 lab = rownames(DEG.stromal),
  x = 'avg_log2FC',
  y = 'p_val_adj',
  xlim = c(-1.5,1.5),
  ylim = c(0, 320),
  pCutoff = 5e-50,
  pointSize = 1,
  colCustom = volcano.cols.s,
  axisLabSize = 12,
  labSize = 6,
  ylab = 'Adj. P-val',
  xlab = 'Log2FC',
  title = 'Fistula Tract vs. Rectal Mucosa Origin',
  subtitle = 'Stromal',
  selectLab = c(stromal.list.volcano, 'MMP1', 'MMP3', 'MMP13'),
  boxedLabels = FALSE,
  drawConnectors = TRUE,
  colConnectors = 'black',
  widthConnectors = 0.5,
  legendPosition = 'none') +  ggplot2::labs(title = NULL,caption = NULL)

#Figure S3A
res.rectum <- results(msccr.rectum, contrast = c('BiopsyType', 'Inflamed', 'Non_inflamed'), alpha = 0.05)
summary(res.rectum) #6222 up, 20% ; 4466 down, 14%

custom.cols.rectum <- ifelse(
  res.rectum$log2FoldChange < 0, 'blue', 'red')
names(custom.cols.rectum)[custom.cols.rectum == 'blue'] <- 'Non_inflamed'
names(custom.cols.rectum)[custom.cols.rectum == 'red'] <- 'Inflamed'

EnhancedVolcano(res.rectum,
                lab = rownames(res.rectum),
                legendPosition = 'none',
                x = 'log2FoldChange',
                y = 'padj',
                xlim = c(-3,6),
                ylim = c(0, 80),
                pCutoff = 5e-10,
                pointSize = 1,
                #col = c('black', 'orange', 'skyblue', 'lightcoral'),
                colCustom = custom.cols.rectum,
                labSize = 5,
                ylab = bquote(~-Log[10]~italic(Padj)),
                title = 'Inflamed Vs. Non-inflamed',
                subtitle = 'Bulk RNA: Rectum Biopsies',
                selectLab = c('IDO1','TNIP3','CHI3L1','CXCL11','S100A9'),
                boxedLabels = TRUE,
                drawConnectors = TRUE,
                colConnectors = 'black',
                widthConnectors = 0.5
)+ ggplot2::labs(caption = NULL)

#Figure S3B
CHI3L1_Meta = data.frame(CHI3L1 = ms.int[['SCT']]@data['CHI3L1',],
  ClusterAnnotation = ms.int$MyeloidStromalAnno,
  Patient = ms.int$PubID,
  FistulaOrigin = ms.int$FistulaAnnotation,
  Ancestry = ms.int$Ancestry,
  Race = ms.int$Race,
  TNF_Status = ms.int$TNFAtProcedure,
  Inflammation = ms.int$Inflammation,
  FistulaCurrent = ms.int$FistulaCurrent
  )
#write.csv(CHI3L1_Meta, file='/data0/2023-09-18_MS-int_CHI3L1exp-Meta.csv')
ggplot(data=CHI3L1_Meta, aes(x=Patient, y=CHI3L1, color=Inflammation)) +
  geom_boxplot() 
ggplot(data=CHI3L1_Meta, aes(x=ClusterAnnotation, y=CHI3L1, color=Patient)) +
  geom_boxplot() 
 
#endoscopy subset
CHI3L1_Meta_endo <- CHI3L1_Meta %>%
  filter(!Patient %in%  c('Perianal11', 'Perianal14', 'Perianal16'))

#inflammation
endo.infl <- ggplot(CHI3L1_Meta_endo, aes(x = Inflammation, y = CHI3L1, color=Inflammation)) +
  geom_violin(trim = TRUE) + 
  scale_color_manual(values=c('red', 'blue'))+
  stat_compare_means(method='wilcox.test', label='p.format', label.x = 1.5)+
  theme_minimal() +
  theme(legend.position = 'none', axis.title = element_text(size=16), axis.text = element_text(size=14))
aggregate(CHI3L1 ~ Inflammation, CHI3L1_Meta_endo, mean)
#antiTNF current
endo.tnf <- ggplot(CHI3L1_Meta_endo, aes(x = TNF_Status, y = CHI3L1, color=TNF_Status)) +
  geom_violin(trim = TRUE) + 
  scale_color_manual(values=c('darkgreen', 'coral'))+
  stat_compare_means(method='wilcox.test', label='p.format', label.x = 1.5)+
  theme_minimal() +
  theme(legend.position = 'none', axis.title = element_text(size=16), axis.text = element_text(size=14))
aggregate(CHI3L1 ~ TNF_Status, CHI3L1_Meta_endo, mean)
#endo.ancestry+endo.race+endo.infl+endo.tnf+ plot_layout(ncol=4)
endo.infl+endo.tnf+ plot_layout(ncol=1)

#Figure 2F
do_EnrichmentHeatmap(ms.int, input_gene_list = enrichment.genes, 
                     cluster_rows = T, symmetrical_scale = T, use_viridis = F, heatmap.legend.length = 50)


#Figure S3C
DimPlot(fibroblasts, label = FALSE, group.by = 'Region', cols = c('maroon3', 'slategray3')) +
  ggplot2::labs(title = NULL) + ggplot2::theme(legend.position = 'none',axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) + 
  xlab(NULL) + ylab(NULL)

#Figure S3D
fibroblasts.df = data.frame(CHI3L1 = fibroblasts[['RNA']]@data['CHI3L1',],
                      OriginalRegion = fibroblasts$Region)
mean.CHI3L1 <- aggregate(CHI3L1 ~ OriginalRegion, fibroblasts.df, mean)

ggplot(data=fibroblasts.df, aes(x=OriginalRegion, y=CHI3L1, color=OriginalRegion)) +
  geom_boxplot() +
  scale_color_manual(values = c('maroon3', 'slategray3')) +
  theme_classic() +
  theme(axis.text = element_text(colour = "black", family = "Arial", size = 12),
        axis.title = element_text(size = 13)) +
  theme(legend.position="none") +
  xlab(label = NULL) +
  ylab(label = 'CHI3L1 Expression') +
  stat_compare_means(label.x = 1.6, label.y = 100)

#Figure 2G
CHI3L1.cells <- WhichCells(fibroblasts, expression = CHI3L1 > 0)
fibroblasts$CHI3L1_score <- ifelse(colnames(fibroblasts) %in% CHI3L1.cells, 'Pos', 'Neg') #https://www.biostars.org/p/408922/#443884

CHI3L1.pos <- subset(fibroblasts, subset = CHI3L1 > 0)
table(CHI3L1.pos$orig.ident) #281 perianal, 108 martin
DEG.CHI3L1pos <- FindMarkers(CHI3L1.pos, ident.1 = 'PerianalCD', group.by = 'orig.ident', logfc.threshold = 0)
DEG.CHI3L1pos %>%
  top_n(n = 25, wt = avg_log2FC) -> top25.up
top25.up <- top25.up[order(-top25.up$avg_log2FC),]
DEG.CHI3L1pos %>%
  top_n(n = -25, wt = avg_log2FC) -> top25.down
top25.down <- top25.down[order(-top25.down$avg_log2FC),]
top25.up <- rownames(top25.up)
top25.down <- rownames(top25.down)
DEG.CHI3L1.list <- c(top25.up, top25.down)

scale_fibroblasts <- rownames(fibroblasts)
fibroblasts <- ScaleData(fibroblasts, features = scale_fibroblasts)
plot_heatmap(
  dataset = fibroblasts,
  markers = DEG.CHI3L1.list,
  sort_var = c('Region', 'CHI3L1_score'),
  anno_var = c('Region', 'CHI3L1_score'),
  anno_colors = list(c('maroon3', 'slategray3'),
                      c('darkgreen', 'coral'))
) 

#Figure S3E - DEGs input into IPA
CHI3L1.DEGs <- FindMarkers(fibroblasts, ident.1 = 'Pos', group.by = 'CHI3L1_score')
#write.csv(CHI3L1.DEGs, file = '/data0/2023-09-11_Fibroblasts_CHI3L1pos_DEGs.csv')

#Figure 2H
myeloid.subset <- c('Macrophages', 'Monocytes', 'Monocytes/Macrophages', 'mregDCs','cDC2')
myeloid.int <- subset(ms.int, idents = myeloid.subset)
DEG.myeloid <- FindMarkers(myeloid.int, ident.1 = 'Fistula', group.by = 'FistulaAnnotation', logfc.threshold = 0)
#write.csv(DEG.myeloid, file = '/data0/2023-11-14_DEGmyeloid_n-15.csv')

volcano.cols.m <- ifelse(
  DEG.myeloid$avg_log2FC < 0, 'navy', 'maroon3')
names(volcano.cols.m)[volcano.cols.m == 'navy'] <- 'Rectum'
names(volcano.cols.m)[volcano.cols.m == 'maroon3'] <- 'Fistula'

DEG.myeloid %>%
  top_n(n = 10, wt = avg_log2FC) -> top10.myeloid.up
DEG.myeloid %>%
  top_n(n = -10, wt = avg_log2FC) -> top10.myeloid.down
myeloid.list1 <- rownames(top10.myeloid.up)
myeloid.list2 <- rownames(top10.myeloid.down)
myeloid.list.volcano <- c(myeloid.list1, myeloid.list2)

EnhancedVolcano(DEG.myeloid,
 lab = rownames(DEG.myeloid),
  x = 'avg_log2FC',
  y = 'p_val_adj',
  xlim = c(-1.5,1.5),
  ylim = c(0, 320),
  pCutoff = 5e-50,
  pointSize = 1,
  colCustom = volcano.cols.m,
  axisLabSize = 12,
  labSize = 6,
  ylab = 'Adj. P-val',
  xlab = 'Log2FC',
  title = 'Fistula Tract vs. Rectal Mucosa Origin',
  subtitle = 'Myeloid',
  selectLab = c(myeloid.list.volcano, 'SOD2'),
  boxedLabels = FALSE,
  drawConnectors = TRUE,
  colConnectors = 'black',
  widthConnectors = 0.5,
  legendPosition = 'none') + 
  ggplot2::labs(title = NULL,caption = NULL)

myeloid.int <- SetIdent(myeloid.int, value = 'FistulaAnnotation')
CD14.cells <- WhichCells(myeloid.int, expression = CD14 > 0)
myeloid.int$CD14_score <- ifelse(colnames(myeloid.int) %in% CD14.cells, 'Pos', 'Neg') 
table(myeloid.int$CD14_score) #neg 3121, pos 1637

#Figure 2I
#phagocytic module
rm.phago <- c('FKBP1C', 'ZNF705A', 'ANP32CP', 'NCBP2L', 'TNFSF12-TNFSF13') #genes in original list not in SCT assay
phagocytic <- phagocytic[!phagocytic %in% rm.phago] 
myeloid.int <- AddModuleScore(myeloid.int, features = c(phagocytic), name = 'Phagocytic_Module')
phago.feat <- plot_density(myeloid.int,'Phagocytic_Module1')  +
  ggplot2::labs(title = 'Phagocytic Module') + ggplot2::theme(legend.position = 'bottom', legend.justification = 'center', legend.text = element_text(angle=45, vjust = 0.7))
vlndf.phago = data.frame(Phagocytic_Module = myeloid.int$Phagocytic_Module1,
                      FistulaAnnotation = myeloid.int$FistulaAnnotation,
                      CD14_score = myeloid.int$CD14_score)
mean.phago.fistula <- aggregate(Phagocytic_Module ~ FistulaAnnotation, vlndf.phago, mean)
phago.vln <- ggplot(data=vlndf.phago, aes(x=FistulaAnnotation, y=Phagocytic_Module, color=FistulaAnnotation)) +
  geom_violin() +
  scale_color_manual(values = c('maroon3', 'navy')) +
  theme_classic() +
  theme(axis.text = element_text(colour = "black", family = "Arial", size = 12),
        axis.title = element_text(size = 13)) +
  theme(legend.position="none") +
  xlab(label = 'Phagocytic Module') +
  ylab(label = 'Module Score') +
  stat_compare_means(method = 'wilcox.test', label = 'p.format', label.x=1.5, label.y=1.5) +
  geom_text(data = mean.phago.fistula, aes(label = round(Phagocytic_Module, digits = 4), y = 0.9), nudge_x = -0.22)
compare_means(Phagocytic_Module ~ CD14_score, data = vlndf.phago, method = 'wilcox.test') #p = 0.091

#oxidative module
myeloid.int <- AddModuleScore(myeloid.int, features = c(oxidative), name = 'Oxidative_Module')
oxi.feat <- plot_density(myeloid.int,'Oxidative_Module1')  +
  ggplot2::labs(title = 'Oxidative Module') + ggplot2::theme(legend.position = 'bottom', legend.justification = 'center', legend.text = element_text(angle=45, vjust = 0.7))

vlndf.oxi = data.frame(Oxidative_Module = myeloid.int$Oxidative_Module1,
                      FistulaAnnotation = myeloid.int$FistulaAnnotation, 
                      CD14_score = myeloid.int$CD14_score)
mean.oxi.fistula <- aggregate(Oxidative_Module ~ FistulaAnnotation, vlndf.oxi, mean)
oxi.vln <- ggplot(data=vlndf.oxi, aes(x=FistulaAnnotation, y=Oxidative_Module, color=FistulaAnnotation)) +
  geom_violin() +
  scale_color_manual(values = c('maroon3', 'navy')) +
  theme_classic() +
  theme(axis.text = element_text(colour = "black", family = "Arial", size = 12),
        axis.title = element_text(size = 13)) +
  theme(legend.position="none") +
  xlab(label = 'Oxidative Module') +
  ylab(label = 'Module Score') +
  stat_compare_means(method = 'wilcox.test', label = 'p.format', label.y=3, label.x=1.5) +
  geom_text(data = mean.oxi.fistula, aes(label = round(Oxidative_Module, digits = 4), y = 1.9))
compare_means(Oxidative_Module ~ CD14_score, data = vlndf.oxi, method = 'wilcox.test') #p = 0.18

#inflammatory module
rm.infl <- c('NAPSA', 'FYB2', 'ZNF705A')
inflammatory <- inflammatory[!inflammatory %in% rm.infl]
myeloid.int <- AddModuleScore(myeloid.int, features = c(inflammatory), name = 'Inflammatory_Module')

infl.feat <- plot_density(myeloid.int,'Inflammatory_Module1')  +
  ggplot2::labs(title = 'Inflammatory Module') + ggplot2::theme(legend.position = 'bottom', legend.justification = 'center', legend.text = element_text(angle=45, vjust = 0.7))
vlndf.infl = data.frame(Inflammatory_Module = myeloid.int$Inflammatory_Module1,
                      FistulaAnnotation = myeloid.int$FistulaAnnotation,
                      CD14_score = myeloid.int$CD14_score)
mean.infl.fistula <- aggregate(Inflammatory_Module ~ FistulaAnnotation, vlndf.infl, mean)
infl.vln <- ggplot(data=vlndf.infl, aes(x=FistulaAnnotation, y=Inflammatory_Module, color=FistulaAnnotation)) +
  geom_violin() +
  scale_color_manual(values = c('maroon3', 'navy')) +
  theme_classic() +
  theme(axis.text = element_text(colour = "black", family = "Arial", size = 12),
        axis.title = element_text(size = 13)) +
  theme(legend.position="none") +
  xlab(label = 'Inflammatory Module') +
  ylab(label = 'Module Score') +
  stat_compare_means(method = 'wilcox.test', label = 'p.format', label.y=2, label.x=1.5) +
  geom_text(data = mean.infl.fistula, aes(label = round(Inflammatory_Module, digits = 4), y = 1.5))
compare_means(Inflammatory_Module ~ CD14_score, data = vlndf.infl, method = 'wilcox.test') #p = 2.4E-46
aggregate(Inflammatory_Module ~ CD14_score, vlndf.infl, mean)

phago.feat + phago.vln + oxi.feat + oxi.vln + infl.feat + infl.vln + plot_layout(ncol=6)

#Figure S3F 
DEG.myeloid %>%
  top_n(n = 30, wt = avg_log2FC) -> top30.myeloid.up
myeloid.list <- rownames(top30.myeloid.up)
DEG.stromal %>%
  top_n(n = 30, wt = avg_log2FC) -> top30.stromal.up
stromal.list <- rownames(top30.stromal.up)

myeloid.int <- AddModuleScore(myeloid.int, features = stromal.list, name = 'StromalUp')
myeloid.module <- do_FeaturePlot(myeloid.int, features = 'StromalUp1', order = T, pt.size = 0.5, enforce_symmetry = TRUE, plot_cell_borders = FALSE, 
                                 split.by = 'FistulaAnnotation',legend.title = 'Stromal Module', font.size = 12, legend.length = 6, legend.position = 'bottom', legend.width = 0.5)

stromal.int<- AddModuleScore(stromal.int, features = myeloid.list, name = 'MyeloidUp')
stromal.module <- do_FeaturePlot(stromal.int, features = 'MyeloidUp1', enforce_symmetry = T, order = T, plot_cell_borders = FALSE, 
                                 pt.size = 0.5, split.by = 'FistulaAnnotation', legend.title = 'Myeloid Module', font.size = 12, legend.length = 6, legend.position = 'bottom', legend.width = 0.5)

myeloid.module / stromal.module

#Figure S3G
ms.int <- AddModuleScore(ms.int, features = c(m2.markers), name = 'M2_Module')
m2.vln.1 <- VlnPlot(ms.int, features = 'M2_Module1', split.by = 'FistulaAnnotation', 
    idents = c('Macrophages'), cols = c('maroon3', 'navy'), pt.size = 1)
m2.vln.2 <- VlnPlot(ms.int, features = 'M2_Module1', split.by = 'FistulaAnnotation', 
    idents = c('Monocytes'), cols = c('maroon3', 'navy'), pt.size = 1)
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
m2.vln.1 + m2.vln.2 + plot_layout(ncol = 2) 


#Figure 2J
ms.int <- SetIdent(ms.int, value = 'FistulaAnnotation')
plot_density(ms.int, c("CD14", "CHI3L1"), joint = TRUE)

#Figure S3H
#counts in myeloid
table(myeloid.int$FistulaAnnotation) #F: 1354, T: 3404
myeloid.F <- subset(myeloid.int, subset = FistulaAnnotation == 'Fistula')
myeloid.T <- subset(myeloid.int, subset = FistulaAnnotation == 'Tissue')
summary(WhichCells(myeloid.T, expression = CHI3L1 > 0)) #21
summary(WhichCells(myeloid.F, expression = CHI3L1 > 0)) #137
summary(WhichCells(myeloid.T, expression = CD14 > 0)) #1114
summary(WhichCells(myeloid.F, expression = CD14 > 0)) #523
summary(WhichCells(myeloid.T, expression = CD14 > 0 & CHI3L1 > 0)) #3
summary(WhichCells(myeloid.F, expression = CD14 > 0 & CHI3L1 > 0)) #70
#counts in stromal
table(stromal.int$FistulaAnnotation) #F: 904, T: 12083
stromal.F <- subset(stromal.int, subset = FistulaAnnotation == 'Fistula')
stromal.T <- subset(stromal.int, subset = FistulaAnnotation == 'Tissue')
summary(WhichCells(stromal.T, expression = CHI3L1 > 0)) #254
summary(WhichCells(stromal.F, expression = CHI3L1 > 0)) #372
summary(WhichCells(stromal.T, expression = CD14 > 0)) #298
summary(WhichCells(stromal.F, expression = CD14 > 0)) #26
summary(WhichCells(stromal.T, expression = CD14 > 0 & CHI3L1 > 0)) #5
summary(WhichCells(stromal.F, expression = CD14 > 0 & CHI3L1 > 0)) #12
#figure in Prism9

##########################
###### Figure S4 #########
##########################

#Figure S4B
VlnPlot(ms.int, features = c('ACTA2'), split.by = 'FistulaAnnotation', split.plot = TRUE,
        idents = c('Myofibroblasts'), cols = c('maroon3', 'navy')) +
        ggplot2::theme(legend.position = 'bottom', legend.justification = 'center', axis.text.x= element_text(angle=0, hjust=0.5, size=14),
          axis.text.y= element_text(size=14)) + 
        xlab(NULL)

#Figure S4E
DotPlot(ms.int_abbr, features=c('CD14', 'CCL18', 'MRC1'), cols='RdBu', split.by = 'FistulaAnnotation') + 
  xlab(NULL) + ylab(NULL) +
  ggplot2::theme(legend.position = 'right') +
  ggplot2::theme(axis.text.x = element_text(colour = "black", family = "Arial"))

##########################
#### Figures 5, S6 #######
##########################
#Figure 5B
DefaultAssay(multiome) <- 'peaks'
DimPlot(multiome, reduction = 'wnn.umap',label = FALSE, cols='alphabet') + 
  ggplot2::labs(title = NULL) + ggplot2::theme(legend.position = 'right', legend.justification = 'center',axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) + 
  xlab(NULL) + ylab(NULL)

#Figure S6A
tibble.multiome <- structure(multiome)
tibble.multiome <- as_tibble(tibble.multiome)
multiome.order <- c('Perianal4', 'Perianal5', 'Perianal9', 'Perianal10',  'Perianal12', 'Perianal13')
names(pals::alphabet(n=19)) #colors used below w substitutions
multiome.pal <- c('plum1', 'blue', 'tan4', 'purple4', 'gray0', 'seagreen4', 'springgreen3', 'bisque', 'gray', 'darkseagreen1', 
  'goldenrod4', 'olivedrab1', 'magenta', 'dodgerblue4', 'orange', 'pink', '#426600', 'red', 'turquoise1')

ggplot(tibble.multiome, aes(x = Publication_ID, fill = Anno)) +
  geom_bar(position = "fill", width = 0.98) + 
  theme_classic() +
  theme(legend.position = 'right', legend.title = element_blank(), legend.key.size = unit(3, 'mm')) +
  theme(axis.text.x = element_text(colour = "black", family = "Arial", size = 12),
        axis.text.y = element_text(colour = "black", family = "Arial", size = 12)) +
  xlab(label = NULL) +
  ylab(label = NULL) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0))  +
  scale_fill_manual(values=multiome.pal) +
  scale_x_discrete(limits = multiome.order) +
  coord_flip()

#Re: Supp. Table
multiome.celltypist <- read.csv('/data0/2023-11-06_Multiome_CellTypist.csv', row.names = 1)
str(multiome.celltypist)
multiome <- AddMetaData(multiome, metadata = multiome.celltypist)
MajorityVotingxAnno <- table(multiome$majority_voting, multiome$ClusterAnnotation)
MajorityVotingxBin <- table(multiome$majority_voting, multiome$Bin)
#write.csv(MajorityVotingxAnno, file = '/data0/2023-11-06_Multiome_MajorityxAnno.csv')
#write.csv(MajorityVotingxBin, file = '/data0/2023-11-06_Multiome_MajorityxBin.csv')
PredictedLabelsxAnno <- table(multiome$predicted_labels, multiome$ClusterAnnotation)
MajorityVotingxPredictedLabels <- table(multiome$majority_voting, multiome$predicted_labels)
write.csv(PredictedLabelsxAnno, file = '/data0/2024-02-19_Multiome_PredictedLabelsxAnno.csv')
write.csv(MajorityVotingxPredictedLabels, file = '/data0/2024-02-19_Multiome_MajorityVotingxPredictedLabels.csv')

multiome <- SetIdent(multiome, value = 'majority_voting')
DefaultAssay(multiome) <- 'peaks'
multiome.cols <- c(rev(pals::alphabet(13)), pals::glasbey(30))
DimPlot(multiome, reduction = 'wnn.umap',label = TRUE, label.size=5, repel=TRUE) + 
  ggplot2::labs(title = NULL) + ggplot2::theme(legend.position = 'none', axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) + 
  xlab(NULL) + ylab(NULL)

#Figure S6C
DefaultAssay(multiome) <- 'peaks'
multiome <- SetIdent(multiome, value='Bin')
bin.cols = c("#F8766D", "#CD9600", "#7CAE00", "#00BE67", "#00BFC4", "#00A9FF", 
"#C77CFF", "#FF61CC")

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

#Figure S6D
qc.count_RNA <- VlnPlot(multiome, features = c("nCount_RNA"), pt.size = 0) + 
  ggplot2::theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + xlab(NULL) +
  ggplot2::theme(legend.position = 'bottom', legend.justification = 'center', legend.key.size = unit(3, 'mm')) 
qc.feature_RNA <- VlnPlot(multiome, features = c("nFeature_RNA"), pt.size = 0) +
  ggplot2::theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + xlab(NULL) +
  ggplot2::theme(legend.position = 'bottom', legend.justification = 'center', legend.key.size = unit(3, 'mm')) 
qc.mt <- VlnPlot(multiome, features = c("percent.mt"), pt.size = 0) +
  ggplot2::theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + xlab(NULL) +
  ggplot2::theme(legend.position = 'bottom', legend.justification = 'center', legend.key.size = unit(3, 'mm')) 
qc.atac_frag <- VlnPlot(multiome, features = c("atac_fragments"), pt.size = 0) +
  ggplot2::theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + xlab(NULL) +
  ggplot2::theme(legend.position = 'bottom', legend.justification = 'center', legend.key.size = unit(3, 'mm')) 
qc.atac_TSS <- VlnPlot(multiome, features = c("atac_TSS_fragments"), pt.size = 0) +
  ggplot2::theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + xlab(NULL) +
  ggplot2::theme(legend.position = 'bottom', legend.justification = 'center', legend.key.size = unit(3, 'mm')) 
qc.atac_region <- VlnPlot(multiome, features = c("atac_peak_region_fragments"), pt.size = 0) +
  ggplot2::theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + xlab(NULL) +
  ggplot2::theme(legend.position = 'bottom', legend.justification = 'center', legend.key.size = unit(3, 'mm')) 
qc.count_atac <- VlnPlot(multiome, features = c("nCount_ATAC"), pt.size = 0) +
  ggplot2::theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + xlab(NULL) +
  ggplot2::theme(legend.position = 'bottom', legend.justification = 'center', legend.key.size = unit(3, 'mm')) 
qc.feature_atac <- VlnPlot(multiome, features = c("nFeature_ATAC"), pt.size = 0) +
  ggplot2::theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + xlab(NULL) +
  ggplot2::theme(legend.position = 'bottom', legend.justification = 'center', legend.key.size = unit(3, 'mm')) 
qc.count_peaks <- VlnPlot(multiome, features = c("nCount_peaks"), pt.size = 0) +
  ggplot2::theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + xlab(NULL) +
  ggplot2::theme(legend.position = 'bottom', legend.justification = 'center', legend.key.size = unit(3, 'mm'))   
qc.feature_peaks <- VlnPlot(multiome, features = c("nFeature_peaks"), pt.size = 0) +
  ggplot2::theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + xlab(NULL) +
  ggplot2::theme(legend.position = 'bottom', legend.justification = 'center', legend.key.size = unit(3, 'mm')) 
qc.reads_peaks <- VlnPlot(multiome, features = c("pct_reads_in_peaks"), pt.size = 0) +
  ggplot2::theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + xlab(NULL) +
  ggplot2::theme(legend.position = 'bottom', legend.justification = 'center', legend.key.size = unit(3, 'mm'))   
qc.nuc_sig <- VlnPlot(multiome, features = c("nucleosome_signal"), pt.size = 0) +
  ggplot2::theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + xlab(NULL) +
  ggplot2::theme(legend.position = 'bottom', legend.justification = 'center', legend.key.size = unit(3, 'mm'))  
qc.nuc_pct <- VlnPlot(multiome, features = c("nucleosome_percentile"), pt.size = 0) +
  ggplot2::theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + xlab(NULL) +
  ggplot2::theme(legend.position = 'bottom', legend.justification = 'center', legend.key.size = unit(3, 'mm'))  
qc.tss_enrich <- VlnPlot(multiome, features = c("TSS.enrichment"), pt.size = 0) +
  ggplot2::theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + xlab(NULL) +
  ggplot2::theme(legend.position = 'bottom', legend.justification = 'center', legend.key.size = unit(3, 'mm'))  
qc.tss_pct <- VlnPlot(multiome, features = c("TSS.percentile"), pt.size = 0) +
  ggplot2::theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + xlab(NULL) +
  ggplot2::theme(legend.position = 'bottom', legend.justification = 'center', legend.key.size = unit(3, 'mm'))  

qc.count_RNA + qc.feature_RNA + qc.mt + qc.atac_frag +  qc.atac_TSS + qc.atac_region + qc.count_atac + qc.feature_atac + qc.count_peaks + qc.feature_peaks + qc.reads_peaks +
  qc.nuc_sig + qc.nuc_pct + qc.tss_enrich + qc.tss_pct + guide_area() +
  plot_layout(ncol=4, guides = "collect")  & theme(legend.position = 'right')
qc.count_RNA + qc.feature_RNA + qc.mt + qc.atac_frag +  qc.atac_TSS + qc.atac_region + qc.count_atac + qc.feature_atac + qc.count_peaks + qc.feature_peaks + qc.reads_peaks +
  qc.nuc_sig + qc.nuc_pct + qc.tss_enrich + qc.tss_pct + guide_area() +
  plot_layout(ncol=5, guides = "collect")  & theme(legend.position = 'bottom', legend.justification=c('left','center'))

#Figure S6E
infl.chi3l1 <- CoveragePlot(multiome, region = c('CHI3L1'), group.by = 'Inflammation', 
extend.upstream=1000, extend.downstream=1000) & scale_fill_manual(values = c('red', 'blue'))

infl.osm <- CoveragePlot(multiome, region = c('OSM'), group.by = 'Inflammation',
extend.upstream=1000, extend.downstream=1000) & scale_fill_manual(values = c('red', 'blue'))

infl.plau <- CoveragePlot(multiome, region = c('PLAU'), group.by = 'Inflammation',
extend.upstream=1000, extend.downstream=1000) & scale_fill_manual(values = c('red', 'blue'))

infl.plaur <- CoveragePlot(multiome, region = c('PLAUR'), group.by = 'Inflammation',
extend.upstream=1000, extend.downstream=1000) & scale_fill_manual(values = c('red', 'blue'))

infl.chi3l1 + infl.osm + infl.plau + infl.plaur + plot_layout(ncol=4)

#Figure 5G
DefaultAssay(PerianalCD) <- 'SCT'
ap1.df <- data.frame(
  JUN = PerianalCD[['SCT']]@data["JUN",],
  JUNB = PerianalCD[['SCT']]@data["JUNB",],
  JUND = PerianalCD[['SCT']]@data["JUND",],
  FOS = PerianalCD[['SCT']]@data["FOS",],
  FOSL1 = PerianalCD[['SCT']]@data["FOSL1",],
  FOSL2 = PerianalCD[['SCT']]@data["FOSL2",],
  cluster = PerianalCD$ClusterAnnotation,
  Bin = PerianalCD$Bin,
  ancestry = PerianalCD$Ancestry,
  inflammation = PerianalCD$Inflammation
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

box.JUN + box.FOS + plot_layout(guides='collect') & theme(legend.position = 'right')

#Figure S6G
#scRNA - ancestry
DEGs.AAvsEA <- FindMarkers(PerianalCD, group.by='Ancestry', ident.1='AA', ident.2='EA', logfc.threshold=0)
#write.csv(DEGs.AAvsEA, file = '/data0/2023-11-12_PerianalCD_AAvsEA-DEGs.csv')

custom.cols.all <- ifelse(
  DEGs.AAvsEA$avg_log2FC < 0, 'darkgoldenrod1', 'lightskyblue')
names(custom.cols.all)[custom.cols.all == 'darkgoldenrod1'] <- 'EA'
names(custom.cols.all)[custom.cols.all == 'lightskyblue'] <- 'AA'

volcano.scrna <- EnhancedVolcano(DEGs.AAvsEA,
                                    lab = rownames(DEGs.AAvsEA),
                                    x = 'avg_log2FC',
                                    y = 'p_val_adj',
                                    xlim = c(-1,1),
                                    ylim = c(0, 100),
                                    pCutoff = 0.05,
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
volcano.scrna                                    

mean.jun <- aggregate(JUN ~ inflammation, ap1.df, mean)
ggplot(data = ap1.df, aes(x = inflammation, y = JUN, color = inflammation)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(size = 14),
        axis.title.y = element_text(size = 16)) +
  theme(legend.position = 'none') +
  scale_color_manual(values=c('red', 'blue')) +
  xlab(label = NULL) +
  ylab(label = 'Normalized JUN Expression') +
  stat_compare_means(method='wilcox.test', label.y = 6, label= 'p.format')

#bulk RNA - ancestry
DEGs.msccr <- results(msccr.ancestry, contrast = c('ancestry_genetic', 'Black', 'White'), alpha = 0.05)
summary(DEGs.msccr) #1617 up, 5% ; 2281 down, 7%

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
volcano.msccr

volcano.scrna + volcano.msccr

##########################
#### Figures 6, 7, S7 ####
##########################

#Figure 6
library(ComplexUpset)
input <- read.csv('Codes/UpSet/2023-02-14_UpSet_TOBIAS-AllCluster_Input.csv')
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

#Figure 7D
plot_density(PerianalCD, features = c('PTGER4')) + 
  ggplot2::theme(legend.position = 'bottom', legend.justification = 'center', legend.text = element_text(angle=45, vjust = 0.7))

#Figure 7G
plot_density(PerianalCD, features = c('LACC1')) +
  ggplot2::theme(legend.position = 'bottom', legend.justification = 'center', legend.text = element_text(angle=45, vjust = 0.7))


#Figure S7F
CoveragePlot(multiome, region = c('chr5-40409000-40411600')) # PTGER4 SNP locus = chr5:40410482

#Figure S7G
CoveragePlot(multiome, region = 'LACC1') # LACC1 SNP locus = chr13:43883789