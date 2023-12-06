projdir = '/home/ubuntu/perianal-cd/'
outdir = projdir + 'results/'
tmpdir = projdir + 'tmp/'
import os
import pandas as pd

path_to_blacklist = '/home/ubuntu/perianal-cd/data/hg38-blacklist.v2.bed'
# Create cisTopic object
from pycisTopic.cistopic_class import *
count_matrix=pd.read_csv(projdir+'data/PerianalCD_MultiomeReRun_ATAC.tsv', sep='\t')
cistopic_obj = create_cistopic_object(fragment_matrix=count_matrix, path_to_blacklist=path_to_blacklist)
cell_data = pd.read_csv(projdir+'data/PerianalCD_MultiomeReRun_metadata.csv', index_col=0)

cistopic_obj.add_cell_data(cell_data)

from pycisTopic.cistopic_class import *
models = run_cgs_models(cistopic_obj, n_topics=[2,5,10,15,30,40], n_cpu=6, n_iter=500, _temp_dir = tmpdir)
model=evaluate_models(models, select_model = None, return_model = True, metrics = 
                      ['Arun_2010','Cao_Juan_2009', 'Minmo_2011', 'loglikelihood'], plot_metrics = False, 
                      save = outdir + 'model_selection.pdf')
cistopic_obj.add_LDA_model(model)

from pycisTopic.clust_vis import *
find_clusters(cistopic_obj, target = 'cells', k = 17, res = 0.6, prefix = 'pycisTopic_', scale = True, split_pattern = '-')
run_umap(cistopic_obj, target  = 'cell', scale=True)
os.mkdir(outdir+'visualization')
plot_metadata(cistopic_obj, reduction_name='UMAP', variables = ['ClusterAnnotation', 'Bin'], target='cell', num_columns=2, 
              text_size=10, dot_size=5, figsize=(10,5), save = outdir + 'visualization/UMAP_label.pdf')
plot_topic(cistopic_obj, reduction_name = 'UMAP', target = 'cell', num_columns=5, 
           save= outdir + 'visualization/UMAP_topic_contribution.pdf')

from pycisTopic.topic_binarization import *
region_bin_topics_otsu = binarize_topics(cistopic_obj, method='otsu')
from pycisTopic.diff_features import *
imputed_acc_obj = impute_accessibility(cistopic_obj, selected_cells=None, selected_regions=None, scale_factor=10**6)
normalized_imputed_acc_obj = normalize_scores(imputed_acc_obj, scale_factor=10**4)
variable_regions = find_highly_variable_features(normalized_imputed_acc_obj, plot = True, 
                                                 save = outdir + 'visualization/DAR_HVR_plot.pdf')
markers_dict_bin = find_diff_features(cistopic_obj, imputed_acc_obj, variable='Bin', var_features=variable_regions, 
                                      split_pattern = '-')
markers_dict_majvote = find_diff_features(cistopic_obj, imputed_acc_obj, variable='majority_voting', 
                                          var_features=variable_regions, split_pattern = '-')
os.makedirs(outdir + 'candidate_enhancers')
pickle.dump(region_bin_topics_otsu, open(outdir+ 'candidate_enhancers/region_bin_topics_otsu.pkl', 'wb'))
pickle.dump(markers_dict_bin, open(outdir+ 'candidate_enhancers/markers_dict_bin.pkl', 'wb'))
pickle.dump(markers_dict_majvote, open(outdir+ 'candidate_enhancers/markers_dict_majvote.pkl', 'wb'))

import pyranges as pr
from pycistarget.utils import region_names_to_coordinates
region_sets = {}
region_sets['topics_otsu'] = {}
region_sets['DARs_bin'] = {}
region_sets['DARs_majvote'] = {}
for topic in region_bin_topics_otsu.keys():
	regions = region_bin_topics_otsu[topic].index[region_bin_topics_otsu[topic].index.str.startswith('chr')] #only keep regions on known chromosomes
	region_sets['topics_otsu'][topic] = pr.PyRanges(region_names_to_coordinates(regions))
for DAR in markers_dict_bin.keys():
	regions = markers_dict_bin[DAR].index[markers_dict_bin[DAR].index.str.startswith('chr')]
	region_sets['DARs_bin'][DAR] = pr.PyRanges(region_names_to_coordinates(regions))
for DAR in markers_dict_majvote.keys():
	regions = markers_dict_majvote[DAR].index[markers_dict_majvote[DAR].index.str.startswith('chr')]
	region_sets['DARs_majvote'][DAR] = pr.PyRanges(region_names_to_coordinates(regions))
for key in region_sets.keys():
	print(f'{key}: {region_sets[key].keys()}')
 
rankings_db = projdir + 'data/hg38_screen_v10_clust.regions_vs_motifs.rankings.feather'
scores_db = projdir + 'data/hg38_screen_v10_clust.regions_vs_motifs.scores.feather'
motif_annotation = projdir + 'data/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl'
from scenicplus.wrappers.run_pycistarget import run_pycistarget
run_pycistarget(region_sets = region_sets,species = 'homo_sapiens', 
                save_path = os.path.join(outdir, 'motifs'), ctx_db_path = rankings_db, 
                dem_db_path = scores_db, path_to_motif_annotations = motif_annotation, run_without_promoters = True,
                n_cpu = 30, _temp_dir = os.path.join(tmpdir, 'ray_spill'), annotation_version = 'v10nr_clust')

import scanpy as sc
adata = sc.read_h5ad(projdir + 'data/PerianalCD_Multi_RNA.h5ad')
menr = dill.load(open(outdir + 'motifs/menr.pkl', 'rb'))
from scenicplus.scenicplus_class import create_SCENICPLUS_object
import numpy as np
scplus_obj = create_SCENICPLUS_object(GEX_anndata = adata.raw.to_adata(), cisTopic_obj = cistopic_obj, menr = menr, 
                                      bc_transform_func = lambda x: x + '___cisTopic')
scplus_obj.X_EXP = np.array(scplus_obj.X_EXP.todense())
highly_variable_gene_names = adata.var_names[adata.var['highly_variable']]
highly_variable_gene_names
scplus_obj.subset(genes = highly_variable_gene_names)

from scenicplus.cistromes import *
merge_cistromes(scplus_obj)

from scenicplus.enhancer_to_gene import get_search_space, calculate_regions_to_genes_relationships, GBM_KWARGS
tf_file = projdir + 'data/TF_names_v_1.01.txt'
get_search_space(scplus_obj, biomart_host="http://sep2019.archive.ensembl.org/", species = 'hsapiens',
                 assembly = 'hg38',upstream = [1000, 150000],downstream = [1000, 150000])
calculate_regions_to_genes_relationships(scplus_obj, ray_n_cpu = 50, _temp_dir = tmpdir + 'ray_spill', 
                                         importance_scoring_method = 'GBM', importance_scoring_kwargs = GBM_KWARGS)
from scenicplus.TF_to_gene import *
calculate_TFs_to_genes_relationships(scplus_obj,tf_file = tf_file,ray_n_cpu =60,method = 'GBM',_temp_dir =tmpdir + 'ray_spill',
                                     key= 'TF2G_adj')

from scenicplus.grn_builder.gsea_approach import build_grn
build_grn(scplus_obj,min_target_genes = 10,adj_pval_thr = 1,min_regions_per_gene = 0,quantiles = (0.85, 0.90, 0.95),
          top_n_regionTogenes_per_gene = (5, 10, 15),top_n_regionTogenes_per_region = (),binarize_using_basc = True,
          rho_dichotomize_tf2g = True,rho_dichotomize_r2g = True,rho_dichotomize_eregulon = True,rho_threshold = 0.05,
          keep_extended_motif_annot = True,merge_eRegulons = True,order_regions_to_genes_by = 'importance',
          order_TFs_to_genes_by = 'importance',key_added = 'eRegulons_importance',cistromes_key = 'Unfiltered',
          disable_tqdm = False,ray_n_cpu = 20, _temp_dir =tmpdir + 'ray_spill')

from scenicplus.utils import format_egrns
format_egrns(scplus_obj, eregulons_key = 'eRegulons_importance', TF2G_key = 'TF2G_adj', key_added = 'eRegulon_metadata')
get_eRegulons_as_signatures(scplus_obj, eRegulon_metadata_key='eRegulon_metadata', key_added='eRegulon_signatures')

from scenicplus.cistromes import *
region_ranking = make_rankings(scplus_obj, target='region')
score_eRegulons(scplus_obj,
                ranking = region_ranking,
                eRegulon_signatures_key = 'eRegulon_signatures',
                key_added = 'eRegulon_AUC',
                enrichment_type= 'region',
                auc_threshold = 0.05,
                normalize = False,
                n_cpu = 12)

gene_ranking = make_rankings(scplus_obj, target='gene')
# Score gene regulons
score_eRegulons(scplus_obj,
                gene_ranking,
                eRegulon_signatures_key = 'eRegulon_signatures',
                key_added = 'eRegulon_AUC',
                enrichment_type = 'gene',
                auc_threshold = 0.05,
                normalize= False,
                n_cpu = 12)

generate_pseudobulks(scplus_obj,
                         variable = 'ACC_majority_voting',
                         auc_key = 'eRegulon_AUC',
                         signature_key = 'Gene_based',
                         nr_cells = 1,
                         nr_pseudobulks = 20,
                         seed=555)
generate_pseudobulks(scplus_obj,
                         variable = 'ACC_majority_voting',
                         auc_key = 'eRegulon_AUC',
                         signature_key = 'Region_based',
                         nr_cells = 1,
                         nr_pseudobulks = 20,
                         seed=555)
TF_cistrome_correlation(scplus_obj,
                        variable = 'ACC_majority_voting',
                        auc_key = 'eRegulon_AUC',
                        signature_key = 'Gene_based',
                        out_key = 'ACC_majority_voting_eGRN_gene_based')
TF_cistrome_correlation(scplus_obj,
                        variable = 'ACC_majority_voting',
                        auc_key = 'eRegulon_AUC',
                        signature_key = 'Region_based',
                        out_key = 'ACC_majority_voting_eGRN_region_based')

df1 = scplus_obj.uns['eRegulon_AUC']['Gene_based'].copy()
df2 = scplus_obj.uns['eRegulon_AUC']['Region_based'].copy()
df1.columns = [x.split('_(')[0] for x in df1.columns]
df2.columns = [x.split('_(')[0] for x in df2.columns]
correlations = df1.corrwith(df2, axis = 0)
correlations
correlations = correlations[abs(correlations) > 0.5]
correlations
keep = [x for x in correlations.index if '+_+' in x] + [x for x in correlations.index if '-_+' in x]
# Keep extended if not direct
extended = [x for x in keep if 'extended' in x]
direct = [x for x in keep if not 'extended' in x]
keep_extended = [x for x in extended if not x.replace('extended_', '') in direct]
keep = direct + keep_extended
keep
len(keep)
keep_gene = [x for x in scplus_obj.uns['eRegulon_AUC']['Gene_based'].columns if x.split('_(')[0] in keep]
keep_gene = [x for x in keep_gene if (int(x.split('_(')[1].replace('g)', '')) > 10)]
keep_all = [x.split('_(')[0] for x in keep_gene]
keep_all
keep_region = [x for x in scplus_obj.uns['eRegulon_AUC']['Region_based'].columns if x.split('_(')[0] in keep]
scplus_obj.uns['selected_eRegulons'] = {}
scplus_obj.uns['selected_eRegulons']['Gene_based'] = keep_gene
scplus_obj.uns['selected_eRegulons']['Region_based'] = keep_region



from scenicplus.dimensionality_reduction import *
run_eRegulons_umap(scplus_obj,
                   scale=True, signature_keys=['Gene_based', 'Region_based'], selected_regulons=scplus_obj.uns['selected_eRegulons']['Gene_based'])
run_eRegulons_umap(scplus_obj,
                   scale=True, signature_keys=['Gene_based'],
                   reduction_name='eRegulons_UMAP_gb', selected_regulons=scplus_obj.uns['selected_eRegulons']['Gene_based'])
run_eRegulons_umap(scplus_obj,
                   scale=True, signature_keys=['Region_based'],
                   reduction_name='eRegulons_UMAP_rb', selected_regulons=scplus_obj.uns['selected_eRegulons']['Region_based'])

from scenicplus.RSS import *
regulon_specificity_scores(scplus_obj, 'ACC_Bin', signature_keys=['Gene_based'],
                           selected_regulons=scplus_obj.uns['selected_eRegulons']['Gene_based'],out_key_suffix='_gene_based',
                           scale=False)
regulon_specificity_scores(scplus_obj, 'ACC_Bin', signature_keys=['Region_based'], 
                           selected_regulons=scplus_obj.uns['selected_eRegulons']['Region_based'],out_key_suffix='_region_based',
                           scale=False)

gb_rss = scplus_obj.uns['RSS']['Gene_based']
rb_rss = scplus_obj.uns['RSS']['Region_based']

def _plot_rss_internal(rss, cell_type, top_n=5, max_n=None, ax=None):
    if ax is None:
        _, ax = plt.subplots(1, 1, figsize=(4, 4))
    if max_n is None:
        max_n = rss.shape[1]
    data = rss.T[cell_type].sort_values(ascending=False)[0:max_n]
    ax.plot(np.arange(len(data)), data, '.')
    ax.set_ylim([np.floor(data.min() * 100.0) / 100.0,np.ceil(data.max() * 100.0) / 100.0])
    ax.set_title(cell_type)
    ax.set_xticklabels([])
    font = {'color': 'black', 'weight': 'normal'}
    x_offsets = np.linspace(0, (max_n / 25) * (top_n - 1), top_n)
    for idx, (regulon_name, rss_val) in enumerate(zip(data[0:top_n].index, data[0:top_n].values)):
        ax.plot([idx, idx], [rss_val, rss_val], 'r.')
        ax.text(idx + (max_n / 25),rss_val, regulon_name, fontdict=font, horizontalalignment='left', verticalalignment='center')

# region based
data_mat = rb_rss
cats = ['Myeloid', 'Fibroblasts', 'Stromal', 'Epithelial', 'CytotoxicT_InnateLymphoid', 'T_Lymphoid', 'Plasma_Cells', 'B_Cells']
i = 1
fig = plt.figure(figsize=(16,8))
pdf = None
for c in cats:
    x = data_mat.T[c]
    ax = fig.add_subplot(2, 4, i)
    i = i+1
    _plot_rss_internal(data_mat, c, top_n=5,max_n=None, ax=ax)
    ax.set_ylim(x.min()-(x.max()-x.min())*0.05,x.max()+(x.max()-x.min())*0.05)
    for t in ax.texts:
            t.set_fontsize(13)
    ax.set_ylabel('')
    ax.set_xlabel('')
    adjust_text(ax.texts, autoalign='xy', ha='right', va='bottom', arrowprops=dict(arrowstyle='-', color='lightgrey'), precision=0.001)
fig.text(0.5, 0.0, 'eRegulon rank', ha='center',va='center', size='x-large')
fig.text(0.00, 0.5, 'Region Based eRegulon specificity score (eRSS)',ha='center', va='center', rotation='vertical', size='x-large')
plt.tight_layout()
plt.rcParams.update({
    'figure.autolayout': True,
    'figure.titlesize': 'large',
    'axes.labelsize': 'medium',
    'axes.titlesize': 'large',
    'xtick.labelsize': 'medium',
    'ytick.labelsize': 'medium'
})
fig.savefig(projdir + 'results/1204_region_based_rss.pdf', bbox_inches='tight')
plt.show()

# gene based
data_mat = gb_rss
cats = ['Myeloid', 'Fibroblasts', 'Stromal', 'Epithelial', 'CytotoxicT_InnateLymphoid', 'T_Lymphoid', 'Plasma_Cells', 'B_Cells']
i = 1
fig = plt.figure(figsize=(16,8))
pdf = None
for c in cats:
    x = data_mat.T[c]
    ax = fig.add_subplot(2, 4, i)
    i = i+1
    _plot_rss_internal(data_mat, c, top_n=5,max_n=None, ax=ax)
    ax.set_ylim(x.min()-(x.max()-x.min())*0.05,x.max()+(x.max()-x.min())*0.05)
    for t in ax.texts:
            t.set_fontsize(13)
    ax.set_ylabel('')
    ax.set_xlabel('')
    adjust_text(ax.texts, autoalign='xy', ha='right', va='bottom', arrowprops=dict(arrowstyle='-', color='lightgrey'), precision=0.001)
fig.text(0.5, 0.0, 'eRegulon rank', ha='center',va='center', size='x-large')
fig.text(0.00, 0.5, 'Gene Based eRegulon specificity score (eRSS)',ha='center', va='center', rotation='vertical', size='x-large')
plt.tight_layout()
plt.rcParams.update({
    'figure.autolayout': True,
    'figure.titlesize': 'large',
    'axes.labelsize': 'medium',
    'axes.titlesize': 'large',
    'xtick.labelsize': 'medium',
    'ytick.labelsize': 'medium'
})
fig.savefig(projdir + 'results/1204_gene_based_rss.pdf', bbox_inches='tight')
plt.show()

