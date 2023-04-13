# Activate the environment
source activate cpdb

# Running with statistical methods
cellphonedb method statistical_analysis meta_table.txt count_martix.txt --counts-data gene_name

# Plotting with statistical results
cellphonedb plot dot_plot --rows rows.txt --columns cols.txt
cellphonedb plot heatmap_plot meta_table.txt
