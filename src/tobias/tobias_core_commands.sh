# bias correct for 8 bins 

TOBIAS ATACorrect --bam /mnt/cho_lab/disk1/jiayuzh/tmp/RL_Hu_Epithelial.bam --genome /mnt/cho_lab/disk2/cache/ref-genomes/10x/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa --peaks ../seurat-manual/Epithelial_peaks.bed --blacklist resources/hg38-blacklist.v2.bed --outdir results/nonbcell_bcell/ --cores 20
TOBIAS ATACorrect --bam /mnt/cho_lab/disk1/jiayuzh/tmp/RL_Hu_CytotoxicT_InnateLymphoid.bam --genome /mnt/cho_lab/disk2/cache/ref-genomes/10x/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa --peaks ../seurat-manual/CytotoxicT_InnateLymphoid_peaks.bed --blacklist resources/hg38-blacklist.v2.bed --outdir results/noncyto_cyto/ --cores 20
TOBIAS ATACorrect --bam /mnt/cho_lab/disk1/jiayuzh/tmp/RL_Hu_Fibroblasts.bam --genome /mnt/cho_lab/disk2/cache/ref-genomes/10x/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa --peaks ../seurat-manual/Fibroblasts_peaks.bed --blacklist resources/hg38-blacklist.v2.bed --outdir results/nonepi_epi/ --cores 20
TOBIAS ATACorrect --bam /mnt/cho_lab/disk1/jiayuzh/tmp/RL_Hu_Fibroblasts.bam --genome /mnt/cho_lab/disk2/cache/ref-genomes/10x/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa --peaks ../seurat-manual/Fibroblasts_peaks.bed --blacklist resources/hg38-blacklist.v2.bed --outdir results/nonfibro_fibro/ --cores 20
TOBIAS ATACorrect --bam /mnt/cho_lab/disk1/jiayuzh/tmp/RL_Hu_Myeloid.bam --genome /mnt/cho_lab/disk2/cache/ref-genomes/10x/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa --peaks ../seurat-manual/Myeloid_peaks.bed --blacklist resources/hg38-blacklist.v2.bed --outdir results/nonmyeloid_myeloid/ --cores 20
TOBIAS ATACorrect --bam /mnt/cho_lab/disk1/jiayuzh/tmp/RL_Hu_Plasma_Cells.bam --genome /mnt/cho_lab/disk2/cache/ref-genomes/10x/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa --peaks ../seurat-manual/Plasma_Cells_peaks.bed --blacklist resources/hg38-blacklist.v2.bed --outdir results/nonplasma_plasma/ --cores 20
TOBIAS ATACorrect --bam /mnt/cho_lab/disk1/jiayuzh/tmp/RL_Hu_Stromal.bam --genome /mnt/cho_lab/disk2/cache/ref-genomes/10x/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa --peaks ../seurat-manual/Stromal_peaks.bed --blacklist resources/hg38-blacklist.v2.bed --outdir results/nonstromal_stromal/ --cores 20
TOBIAS ATACorrect --bam /mnt/cho_lab/disk1/jiayuzh/tmp/RL_Hu_T_Lymphoid.bam --genome /mnt/cho_lab/disk2/cache/ref-genomes/10x/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa --peaks ../seurat-manual/T_Lymphoid_peaks.bed --blacklist resources/hg38-blacklist.v2.bed --outdir results/nontcell_tcell/ --cores 20

# footprint score calculation for 8 bins

TOBIAS FootprintScores --signal results/nonbcell_bcell/RL_Hu_Epithelial_corrected.bw --regions ../seurat-manual/Epithelial_peaks.bed --output results/nonbcell_bcell/Epithelial_footprints.bw --cores 20
TOBIAS FootprintScores --signal results/noncyto_cyto/RL_Hu_CytotoxicT_InnateLymphoid_corrected.bw --regions ../seurat-manual/CytotoxicT_InnateLymphoid_peaks.bed --output results/noncyto_cyto/CytotoxicT_InnateLymphoid_footprints.bw --cores 20
TOBIAS FootprintScores --signal results/nonepi_epi/RL_Hu_Epithelial_corrected.bw --regions ../seurat-manual/Epithelial_peaks.bed --output results/nonepi_epi/Epithelial_footprints.bw --cores 20
TOBIAS FootprintScores --signal results/nonepi_epi/RL_Hu_Fibroblasts_corrected.bw --regions ../seurat-manual/Fibroblasts_peaks.bed --output results/nonfibro_fibro/Fibroblasts_footprints.bw --cores 20
TOBIAS FootprintScores --signal results/nonepi_epi/RL_Hu_Myeloid_corrected.bw --regions ../seurat-manual/Myeloid_peaks.bed --output results/nonmyeloid_myeloid/Myeloid_footprints.bw --cores 20
TOBIAS FootprintScores --signal results/nonepi_epi/RL_Hu_Plasma_Cells_corrected.bw --regions ../seurat-manual/Plasma_Cells_peaks.bed --output results/nonplasma_plasma/Plasma_Cells_footprints.bw --cores 20
TOBIAS FootprintScores --signal results/nonepi_epi/RL_Hu_Stromal_corrected.bw --regions ../seurat-manual/Stromal_peaks.bed --output results/nonstromal_stromal/Stromal_footprints.bw --cores 20
TOBIAS FootprintScores --signal results/nonepi_epi/RL_Hu_T_Lymphoid_corrected.bw --regions ../seurat-manual/T_Lymphoid_peaks.bed --output results/nontcell_tcell/T_Lymphoid_footprints.bw --cores 20

# annotate peak files (bed)

uropa --bed Epithelial_peaks.bed --gtf genes.gtf
uropa --bed CytotoxicT_InnateLymphoid_peaks.bed --gtf genes.gtf
uropa --bed Epithelial_peaks.bed --gtf genes.gtf
uropa --bed Fibroblasts_peaks.bed --gtf genes.gtf
uropa --bed Myeloid_peaks.bed --gtf genes.gtf
uropa --bed Plasma_Cells_peaks.bed --gtf genes.gtf
uropa --bed Stromal_peaks.bed --gtf genes.gtf
uropa --bed T_Lymphoid_peaks.bed --gtf genes.gtf

# generate peak header file

head -n 1 ../seurat-manual/Epithelial_peaks_allhits.txt > resources/peak_header.txt

# bind detect for 8 bins

TOBIAS BINDetect --motifs ../../data/JASPAR2022_CORE_vertebrates_non-redundant_pfms_jaspar/*.jaspar --signals results/nonbcell_bcell/B_Cells_footprints.bw --genome /mnt/cho_lab/disk2/cache/ref-genomes/10x/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa --peaks ../seurat-manual/B_Cells_peaks_allhits.bed --peak_header resources/peak_header.txt --outdir results/nonbcell_bcell/ --cores 20 --bound-pvalue 0.05
TOBIAS BINDetect --motifs ../../data/JASPAR2022_CORE_vertebrates_non-redundant_pfms_jaspar/*.jaspar --signals results/nonbcell_bcell/CytotoxicT_InnateLymphoid_footprints.bw --genome /mnt/cho_lab/disk2/cache/ref-genomes/10x/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa --peaks ../seurat-manual/CytotoxicT_InnateLymphoid_peaks_allhits.bed --peak_header resources/peak_header.txt --outdir results/noncyto_cyto/ --cores 20 --bound-pvalue 0.05
TOBIAS BINDetect --motifs ../../data/JASPAR2022_CORE_vertebrates_non-redundant_pfms_jaspar/*.jaspar --signals results/nonbcell_bcell/Epithelial_footprints.bw --genome /mnt/cho_lab/disk2/cache/ref-genomes/10x/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa --peaks ../seurat-manual/Epithelial_peaks_allhits.bed --peak_header resources/peak_header.txt --outdir results/nonepi_epi/ --cores 20 --bound-pvalue 0.05
TOBIAS BINDetect --motifs ../../data/JASPAR2022_CORE_vertebrates_non-redundant_pfms_jaspar/*.jaspar --signals results/nonbcell_bcell/Fibroblasts_footprints.bw --genome /mnt/cho_lab/disk2/cache/ref-genomes/10x/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa --peaks ../seurat-manual/Fibroblasts_peaks_allhits.bed --peak_header resources/peak_header.txt --outdir results/nonfibro_fibro/ --cores 20 --bound-pvalue 0.05
TOBIAS BINDetect --motifs ../../data/JASPAR2022_CORE_vertebrates_non-redundant_pfms_jaspar/*.jaspar --signals results/nonbcell_bcell/Plasma_Cells_footprints.bw --genome /mnt/cho_lab/disk2/cache/ref-genomes/10x/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa --peaks ../seurat-manual/Plasma_Cells_peaks_allhits.bed --peak_header resources/peak_header.txt --outdir results/nonplasma_plasma/ --cores 20 --bound-pvalue 0.05
TOBIAS BINDetect --motifs ../../data/JASPAR2022_CORE_vertebrates_non-redundant_pfms_jaspar/*.jaspar --signals results/nonbcell_bcell/Stromal_footprints.bw --genome /mnt/cho_lab/disk2/cache/ref-genomes/10x/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa --peaks ../seurat-manual/Stromal_peaks_allhits.bed --peak_header resources/peak_header.txt --outdir results/nonstromal_stromal/ --cores 20 --bound-pvalue 0.05
TOBIAS BINDetect --motifs ../../data/JASPAR2022_CORE_vertebrates_non-redundant_pfms_jaspar/*.jaspar --signals results/nonbcell_bcell/T_Lymphoid_footprints.bw --genome /mnt/cho_lab/disk2/cache/ref-genomes/10x/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa --peaks ../seurat-manual/T_Lymphoid_peaks_allhits.bed --peak_header resources/peak_header.txt --outdir results/nontcell_tcell/ --cores 20 --bound-pvalue 0.05
