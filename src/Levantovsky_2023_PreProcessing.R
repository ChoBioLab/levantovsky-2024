#code to generate .rds object and other pre-processing for Levantovsky 2023
library(future)
plan(
  multicore,
  workers = 8
) # parallelization
options(future.globals.maxSize = 12 * 1024^2 * 1000)
library(dplyr)
library(Seurat)
library(harmony)

#starting object: Combined PerianalCD and Martin 2019 cohorts
CombinedCD <- readRDS('/data0/RDS_Files/2023-09-06_CombinedCD-v3.RDS')

#fistula.clusters
CombinedCD <- SetIdent(integrated, value = 'PubID')
fistula.clusters <- subset(integrated, idents = c('Perianal11', 'Perianal14', 'Perianal16'))

# re normalize, harmonize, reduce, and cluster the subset object
# make sure to select appropriate values for PCs, dimensions, and resolution
fistula.clusters <- SCTransform(
  object = fistula.clusters,
  vst.flavor = "v2"
) %>%
  RunPCA(
    npcs = 20
  )
fistula.clusters <- RunHarmony(
  object = fistula.clusters,
  group.by.vars = "object",
  reduction = "pca",
  assay.use = "SCT",
  reduction.save = "harmony"
)
fistula.clusters <- RunUMAP(
  object = fistula.clusters,
  reduction = "harmony",
  dims = 1:20
) %>%
  FindNeighbors(
    reduction = "harmony",
    dims = 1:20
  ) %>%
  FindClusters(
    resolution = 1.0 #0.8: fibroblasts not separated out
  )

fistula.clusters <- RenameIdents(fistula.clusters, '0' = 'B Cells',
'1' = 'CD4 T Cells',
'2' = 'B Cells',
'3' = 'CD8 T, NK Cells',
'4' = 'Fibroblasts',
'5' = 'CD4 T Cells, Tregs',
'6' = 'Monocytes, Macrophages',
'7' = 'CD4 T Cells',
'8' = 'CD8 T, NK Cells',
'9' = 'B Cells',
'10' = 'IgA Plasma',
'11' = 'NK Cells, ILC1',
'12' = 'Mast Cells',
'13' = 'Stromal',
'14' = 'B Cells',
'15' = 'CD4 T Cells',
'16' = 'IgG Plasma',
'17' = 'IgA Plasma',
'18' = 'IgA Plasma',
'19' = 'IgG Plasma',
'20' = 'B Cells',
'21' = 'B Cells',
'22' = 'Endothelial',
'23' = 'IgA Plasma',
'24' = 'IgG Plasma',
'25' = 'IgA Plasma',
'26' = 'CD4 T Cells',
'27' = 'IgG Plasma',
'28' = 'IgG Plasma',
'29' = 'B Cells',
'30' = 'Stromal Epithelial Mix',
'31' = 'Dendritic Cells',
'32' = 'IgG Plasma',
'33' = 'IgG Plasma',
'34' = 'Epithelial, TA',
'35' = 'B Cells',
'36' = 'CD4 Tfh Like',
'37' = 'IgA Plasma',
'38' = 'Enteric Neurons',
'39' = 'IgG Plasma')
fistula.clusters$FistulaClustersAnno <- fistula.clusters@active.ident

#PerianalCD cohort subset
CombinedCD <- SetIdent(CombinedCD, value = 'orig.ident')
PerianalCD <- subset(CombinedCD, idents = c('PerianalCD'))

# re normalize, harmonize, reduce, and cluster the subset object
# make sure to select appropriate values for PCs, dimensions, and resolution
PerianalCD <- SCTransform(
  object = PerianalCD,
  vst.flavor = "v2"
) %>%
  RunPCA(
    npcs = 20
  )
PerianalCD <- RunHarmony(
  object = PerianalCD,
  group.by.vars = "object",
  reduction = "pca",
  assay.use = "SCT",
  reduction.save = "harmony"
)
PerianalCD <- RunUMAP(
  object = PerianalCD,
  reduction = "harmony",
  dims = 1:20
) %>%
  FindNeighbors(
    reduction = "harmony",
    dims = 1:20
  ) %>%
  FindClusters(
    resolution = 0.8 
  )

PerianalCD <- RenameIdents(PerianalCD, '0' = 'CD4 T Cells',
'1' = 'B Cells',
'2' = 'CD4 T Cells',
'3' = 'B Cells',
'4' = 'CD8 T Cells',
'5' = 'IgA Plasma',
'6' = 'Fibroblasts',
'7' = 'Myeloid',
'8' = 'CD8 T, NK, ILC1',
'9' = 'IgA Plasma',
'10' = 'IgA Plasma',
'11' = 'Mast Cells',
'12' = 'IgA Plasma',
'13' = 'IgA Plasma',
'14' = 'Epithelial',
'15' = 'IgA Plasma',
'16' = 'IgG Plasma',
'17' = 'Endothelial, Enteric Neurons',
'18' = 'NK, ILC1',
'19' = 'Myofibroblasts, Pericytes',
'20' = 'IgA Plasma',
'21' = 'IgA Plasma',
'22' = 'IgA Plasma',
'23' = 'IgA Plasma',
'24' = 'IgG Plasma',
'25' = 'IgA Plasma',
'26' = 'CD4 T, T Regs',
'27' = 'IgG Plasma',
'28' = 'B Cells',
'29' = 'IgG Plasma',
'30' = 'IgA Plasma',
'31' = 'B Cells',
'32' = 'IgG Plasma',
'33' = 'IgG Plasma',
'34' = 'IgA Plasma',
'35' = 'IgA Plasma',
'36' = 'IgA Plasma',
'37' = 'IgA Plasma',
'38' = 'IgG Plasma',
'39' = 'B Cells',
'40' = 'IgA Plasma',
'41' = 'IgA Plasma',
'42' = 'CD4 T Cells',
'43' = 'IgA Plasma',
'44' = 'CD4 T Cells',
'45' = 'IgA Plasma',
'46' = 'Epithelial',
'47' = 'IgA Plasma',
'48' = 'IgA Plasma',
'49' = 'IgA Plasma',
'50' = 'B Cells',
'51' = 'Epithelial',
'52' = 'CD4 T Cells')
PerianalCD$PerianalCDAnno <- PerianalCD@active.ident

#myeloid-stromal sub-clustering and integration
PerianalCD <- readRDS('/data0/RDS_Files/2023-09-08_PerianalCD_n-15_Anno.rds')
myeloid.stromal.idents <- c('Myeloid', 'Endothelial, Enteric Neurons', 'Myofibroblasts, Pericytes', 'Fibroblasts')
ms.int <- subset(PerianalCD, idents = myeloid.stromal.idents)

# re normalize, harmonize, reduce, and cluster the subset object
# make sure to select appropriate values for PCs, dimensions, and resolution
ms.int <- SCTransform(
  object = ms.int,
  vst.flavor = "v2"
) %>%
  RunPCA(
    npcs = 20
  )
ms.int <- RunHarmony(
  object = ms.int,
  group.by.vars = "object",
  reduction = "pca",
  assay.use = "SCT",
  reduction.save = "harmony"
)
ms.int <- RunUMAP(
  object = ms.int,
  reduction = "harmony",
  dims = 1:20
) %>%
  FindNeighbors(
    reduction = "harmony",
    dims = 1:20
  ) %>%
  FindClusters(
    resolution = 0.8 
  )
ms.int$FistulaAnnotation <- ifelse(ms.int$Region == 'Fistula', 'Fistula', 'Tissue')
#assign ms-int annotations
ms.int <- RenameIdents(ms.int, '0' = 'ADAMDEC1hi Fibroblasts',
'1' = 'PDGFRAhi Fibroblasts',
'2' = 'Endothelial',
'3' = 'Macrophages',
'4' = 'Myofibroblasts',
'5' = 'cDC2',
'6' = 'Monocytes',
'7' = 'Pericytes',
'8' = 'Fibroblastic Reticular Cells',
'9' = 'ADAMDEC1hi Fibroblasts',
'10' = 'CHI3L1hi Fibroblasts',
'11' = 'Enteric Neurons',
'12' = 'Monocytes/Macrophages',
'13' = 'B Cells',
'14' = 'ADAMDEC1hi Fibroblasts',
'15' = 'mregDCs',
'16' = 'T Cells',
'17' = 'CCL11hi Fibroblasts',
'18' = 'CCL11hi Fibroblasts'
)
ms.int$MyeloidStromalAnno <- ms.int@active.ident
#remove B & T cells
ms.int <- subset(ms.int, idents = c('B Cells', 'T Cells'), invert = TRUE)

#subset ileal and fistula fibroblasts
fibroblasts <- subset(CombinedCD, idents=c('Fibroblasts'))
fibroblasts <- SetIdent(fibroblasts, value = c('Region'))
fibroblasts2 <- subset(fibroblasts, idents = c('Fistula', 'Ileum'))
fibroblasts2 <- subset(fibroblasts2, subset = Inflammation == 'Inflamed')

# re normalize, harmonize, reduce, and cluster the subset object
# make sure to select appropriate values for PCs, dimensions, and resolution
fibroblasts2 <- SCTransform(
  object = fibroblasts2,
  vst.flavor = "v2"
) %>%
  RunPCA(
    npcs = 20
  )
fibroblasts2 <- RunHarmony(
  object = fibroblasts2,
  group.by.vars = "object",
  reduction = "pca",
  assay.use = "SCT",
  reduction.save = "harmony"
)
fibroblasts2 <- RunUMAP(
  object = fibroblasts2,
  reduction = "harmony",
  dims = 1:20
) %>%
  FindNeighbors(
    reduction = "harmony",
    dims = 1:20
  ) %>%
  FindClusters(
    resolution = 0.8 
  )

#gene lists for macrophage activation paths
# to convert mouse to human genes from Sanin 2022
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

phagocytic.mouse <- c('Pf41',	'Ltc4s',	'C1qa',	'Ctsb',	'C1qc',	'Apoe',	'C1qb',	'Csf1r',	'Dab2',	'Gas6',	'Ninj1',	'Fcgrt',	'Timp2',	'Lgmn',	'Pltp',	'Itm2b',	'Wfdc17',	'Blvrb',	'Serinc3',	'Rnase4',	'Grn',	'Ctsd1',	'Tgfbr2',	'C3ar1',	'Nrp1',	'Trf',	'Cfh',	'Cd63',	'C5ar1',	'Ang',	'Dhrs3',	'Ctsl',	'Glul',	'Pla2g15',	'Emp1',	'Cd33',	'Cltc1',	'Rhob',	'Egr1',	'Lamp1',	'Ptpn18',	'Ccl6',	'Aplp2',	'Atp13a2',	'Ier3',	'Ms4a7',	'Ptger4',	'Aldh2',	'Mafb',	'Hexa',	'Ftl1',	'Adgre1',	'Zfp36l1','St6galnac4',	'Tpp1',	'Plk2',	'Tcn2',	'AI467606',	'Fcgr2b',	'Hmox1',	'Man2b1',	'Camk1',	'Marcks',	'Pigk',	'Clta',	'Tmem37',	'Adam15',	'Tcf3',	'Prkcb',	'Qk',	'Tbc1d17',	'Cfp1',	'Cd81',	'Vps33a',	'Ccl2',	'Gramd1a',	'Snx2',	'Galk2',	'Canx',	'Tbxas1',	'Ap2a2',	'Dusp6',	'Snx6',	'Golph3',	'Calm2',	'Plod1',	'Trib1',	'Fam3c',	'Aldh9a1',	'Gnl3',	'Cd68',	'Ptprj',	'Rab31',	'Herpud21',	'Fam214b',	'Dok3',	'Sec14l1',	'Man1a',	'Dnajb9',	'Stard5',	'Cst3',	'Ggnbp2',	'Qrich1',	'Cln8',	'Ap2m1',	'Arl11',	'Fkbp1a',	'Mfsd1',	'Rit1',	'Pea15a',	'Ccl9',	'Ctnnb1',	'Ptp4a3',	'Atp6v1a',	'Pepd',	'Naglu',	'Snx5',	'Rgs18',	'Cd37',	'Dennd5a',	'Ehd4',	'Tmem86a1',	'Trem2',	'Gas7',	'Lipa',	'Gpx4',	'Cysltr1',	'Lgals1',	'Atp6ap1',	'Rin2',	'Dpm3',	'Dgkz',	'Hpgds1',	'Bin1',	'Nucb1',	'Mat2a',	'Tspan4',	'As3mt1',	'Ubc',	'Pink1',	'Nfic',	'Tm9sf2',	'Mef2c1',	'Sertad1',	'Tm2d3',	'Unc93b1',	'Vat1',	'Setd3',	'Ctla2b',	'Il6ra',	'Jmjd1c1',	'Atp6v0b',	'Snx3',	'Pdcd5',	'Smap1',	'Trappc5',	'Tom1',	'Rtn41',	'Ubn1',	'Bmp2k',	'Hmgn1',	'Get4',	'Wsb11',	'Mrfap11',	'Slc43a2',	'Agfg1',	'Akirin1',	'Ppil4',	'Ppp5c',	'Rgs10',	'Acat1',	'Ctsa1',	'Ankrd13a',	'Mcfd2',	'Stom',	'Brd2',	'Cat',	'Helz',	'Atraid',	'Bri3',	'Txnip',	'Mpp11',	'Abcf1',	'Tgfbi1',	'Idh1',	'Fli1',	'Sh3bp5',	'Slc29a1',	'Snx4',	'Nenf',	'Nmd3',	'Smim11',	'Arhgap17',	'Rragc',	'Egln2',	'Ddx3x',	'Arrdc1',	'Foxp1',	'Nisch',	'Anxa5',	'Frmd4b',	'Trim47',	'Ddx5',	'Leprot',	'Ssh2',	'Hspa9',	'Gns1',	'Gnpat',	'Asah1',	'Dpp7',	'Ncaph2',	'Cenpb',	'Clec4a1',	'Spag7',	'Nagpa',	'Inpp5d',	'Jund',	'Tcf4',	'P2ry6',	'Sgpl11',	'Cebpg',	'Clk3',	'Luc7l2',	'Pkig',	'Atf4',	'Slc16a6',	'Rab3il1',	'Klf2',	'Pld3',	'Tgoln1',	'App',	'Gpr107',	'Sypl',	'Gabarapl1',	'Ost4',	'Scp2',	'Ivns1abp',	'Abhd121',	'Hacd4',	'Ifi27',	'Txndc5',	'Idh2',	'Dstn',	'Eps8',	'U2af2',	'Ring1',	'Stau1',	'Hsbp1',	'Top1',	'Fbxw4',	'Gatm',	'Tmem106a',	'Atp6ap2',	'Nsmce1',	'Rfk',	'Nsun2',	'Dusp22',	'Plbd2',	'Scarb2',	'Naa50',	'Anp32a',	'Tsc22d3',	'Otulin',	'Ndufc1',	'Nop56',	'Sptlc2',	'Slc11a1',	'Hnrnph1',	'Mafg',	'Dnajc1',	'Mef2a',	'Galnt1',	'Rabac1',	'Sh3bgrl',	'Ubl4a',	'Nr3c1',	'Elk3',	'Tm9sf3',	'Ndufa2',	'Rab11a',	'Sgpp1',	'Plin3',	'Kctd121',	'Uchl3',	'Slc25a5',	'Stat3',	'Cpne3',	'Tmx1',	'Map7d1',	'Ech1',	'Swap70',	'Cmtm6',	'Zmiz1',	'Lamtor1',	'Caml',	'Dnase2a1',	'Ykt6',	'Mrpl34',	'Igsf8',	'Smim15',	'Comt',	'Tmem50a',	'Mapk3',	'Lamtor3',	'Rab24',	'Fos1',	'Stx4a',	'Plxnb2',	'Ncbp2',	'Zfand5',	'Cndp2',	'Gna12',	'Lypla2',	'Pnrc1',	'Eif4g2',	'Mrpl27',	'Tspan3',	'Fam234a',	'Gaa',	'Epn1',	'Txndc12',	'Dynlt3',	'Syap1',	'Hist1h2bc',	'Gnpda1',	'Il10rb',	'Eif5b1',	'Pon2',	'Tnfsf12',	'Clptm1',	'Hnrnpu',	'Kras',	'Psmd1')
phagocytic <- convert_mouse_to_human(phagocytic.mouse) 
oxidative.mouse <- c('Fn12',	'Emilin2',	'Thbs1',	'Acly',	'Itgam',	'Laptm5',	'Gngt2',	'Ltc4s',	'Man2b1',	'Ctsd2',	'Wfdc17',	'Ssr4',	'Ccl22',	'Tcn2',	'Cd9',	'Gpx4',	'Itgb2',	'Tyrobp',	'Adgre2')
oxidative <- convert_mouse_to_human(oxidative.mouse)
inflammatory.mouse <- c('Msrb1',	'Plac8',	'Ifitm3',	'Prdx5',	'Lyz2',	'Samhd1',	'Gsr',	'Tyrobp',	'Ifitm6',	'Smpdl3a',	'Lst1',	'Gngt2',	'S100a4',	'Emilin2',	'Tmsb10',	'Coro1a',	'Tpd5',	'Ly6e',	'Thbs1',	'Fxyd5',	'Napsa',	'Cybb',	'Irgm1',	'Taldo1',	'Clec4e',	'Ptpn11',	'Arpc3',	'Wtap',	'Cib1',	'Snx1',	'Fyb2',	'Snx18',	'Fcer1g',	'Klf13',	'Pkm',	'Alox5ap',	'Mrpl30',	'Pira2',	'Ifitm2',	'Oas1a',	'Msr1',	'Sf3b6',	'Ostf',	'Gyg',	'Stat1',	'Fcgr1',	'Prr13',	'Dok3',	'Sec61g',	'Btg1',	'Tnfrsf1a',	'Rbms1',	'B4galnt1',	'Tln1',	'Aprt',	'Cytip',	'Fam49b',	'Srgn',	'Ncf2',	'Ppp1r21',	'Lyn',	'Gda',	'Sirpb1c',	'Dck',	'Pla2g7',	'Cdc42se1',	'Itgam',	'Hck',	'Arhgdib',	'Sat1',	'Ndufb8',	'Ddx24',	'Actg1',	'Rtp4',	'Ube2d3',	'Mrpl4',	'Rnf149',	'Zcrb1',	'Ptp4a1',	'Echs1',	'Fgr',	'Rbm7',	'Cyba',	'Rasgrp2',	'Ppp1ca',	'Mpeg1',	'Ifi47',	'Emb',	'Gch1',	'Mpp1',	'Ptpn6',	'H3f3a',	'Tifab',	'G6pdx',	'Plaur',	'Vdac3',	'Serp1',	'Rap1b',	'Rassf3',	'Samsn1',	'Blvra',	'Smim7',	'Ube2b',	'Slfn1',	'Psma7',	'Ube2a',	'Ptbp3',	'Stk24',	'Spi1',	'Gpx1',	'Pycard',	'Eno1',	'Mxd1',	'Gtf2b',	'Clec4a1',	'Uqcrh',	'Lamtor5',	'Tor3a',	'Igsf6',	'Cdkn1a',	'Clec4a3',	'Rab8a',	'Epsti1',	'Gpr35',	'Eif3k',	'Zeb2',	'Eif3e',	'Capzb',	'Plgrkt',	'Tspan13',	'Xdh',	'Tsc22d',	'Nadk',	'Flna',	'Tpgs1',	'Sppl2a',	'Add3',	'Ccl4',	'Glud1',	'Naca',	'Smox',	'S100a11',	'Tnfrsf1b',	'Itgb2',	'Plxnb2',	'Scand1',	'Metrnl',	'Shisa5',	'Rras',	'Milr1',	'Osbpl9',	'Degs1',	'Ugp2',	'Itga4',	'Cdk2ap2',	'Xpnpep1',	'Prdx6',	'Mcemp1',	'Plin2',	'Dcaf13',	'Dusp3',	'Klf3',	'Myd88',	'Ifngr1',	'Dr1',	'Fndc3a',	'Dctn3',	'Stard3',	'Csnk1e',	'Tnfaip2',	'Adssl1',	'Cd14',	'Zfp622',	'Sra1',	'Tmem51',	'Il17ra',	'Pnrc2',	'Rap1a',	'Gsto1',	'Mafb',	'Abi3',	'Bak1',	'Stk17b')
inflammatory <- convert_mouse_to_human(inflammatory.mouse)

