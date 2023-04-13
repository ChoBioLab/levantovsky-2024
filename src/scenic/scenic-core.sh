pyscenic grn \--num_workers 20 \--output PerianalCD.adjacencies.tsv \--method grnboost2 \PerianalCD_filtered.loom \allTFs_hg38.txt

pyscenic ctx \PerianalCD.adjacencies.tsv \hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather\--annotations_fname motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl \--expression_mtx_fname PerianalCD_filtered.loom \--mode "dask_multiprocessing" \--output PerianalCD.motifs.csv \--num_workers 20 \--mask_dropouts

pyscenic aucell \PerianalCD_filtered.loom \PerianalCD.motifs.csv \--output PerianalCD_scenic.loom \--num_workers 20