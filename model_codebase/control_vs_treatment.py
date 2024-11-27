# -*- coding: utf-8 -*-
import os
import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("--control_data_path", type=str, required=True, help='Path of h5ad file for control')
parser.add_argument("--condition1_data_path", type=str, required=True, help='Path of h5ad file for condition 1')
parser.add_argument("--condition2_data_path", type=str, default='', help='Path of h5ad file for condition 2')
parser.add_argument("--output_folder", type=str, default='./results', help='output folder')

args = parser.parse_args()
os.makedirs(args.output_folder, exist_ok=True)

# load control data in h5ad format
adata = sc.read_h5ad(args.control_data_path)
print(f'cells: {adata.X.shape[0]} genes: {adata.X.shape[1]}')
# check data quality to ensure the data is in correct format/range
if adata.X.max() > 20.0:
    # it is raw count
    print(f'adata.X.max:{adata.X.max()}, data is raw count')
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
if adata.X.min() < 0.0:
    # it is scaled log1p data
    print(f'min:{adata.X.min()}, data is already scaled')
    if adata.raw is not None:
        print(f'restore from adata.raw')
        # if we want to keep original adata, assign to a different object instead of adata
        adata = adata.raw.to_adata()
        assert adata.X.min() >= 0.0, f"min:{adata.X.min()} seems also scaled, please use unscaled data"
        if adata.X.max() < 20.0:
            print(f"max is {adata.X.max()}, this is not raw count, seems log-normalized")
        else:
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)
    else:
        assert False, 'adata.raw does not exist, please fix your data'

# load condition1 data in h5ad format
adata1 = sc.read_h5ad(args.condition1_data_path)
print(f'cells: {adata1.X.shape[0]} genes: {adata1.X.shape[1]}')
# check data quality to ensure the data is in correct format/range
if adata1.X.max() > 20.0:
    # it is raw count
    print(f'adata1.X.max:{adata1.X.max()}, data is raw count')
    sc.pp.normalize_total(adata1, target_sum=1e4)
    sc.pp.log1p(adata1)
if adata1.X.min() < 0.0:
    # it is scaled log1p data
    print(f'min:{adata1.X.min()}, data is already scaled')
    if adata1.raw is not None:
        print(f'restore from adata1.raw')
        # if we want to keep original adata, assign to a different object instead of adata
        adata1 = adata1.raw.to_adata()
        assert adata1.X.min() >= 0.0, f"min:{adata1.X.min()} seems also scaled, please use unscaled data"
        if adata1.X.max() < 20.0:
            print(f"max is {adata1.X.max()}, this is not raw count, seems log-normalized")
        else:
            sc.pp.normalize_total(adata1, target_sum=1e4)
            sc.pp.log1p(adata1)
    else:
        assert False, 'adata1.raw does not exist, please fix your data'

# load condition2 data in h5ad format if applicable
adata2 = None
if args.condition2_data_path != '':
    adata2 = sc.read_h5ad(args.condition2_data_path)
    print(f'cells: {adata2.X.shape[0]} genes: {adata2.X.shape[1]}')
    # check data quality to ensure the data is in correct format/range
    if adata2.X.max() > 20.0:
        # it is raw count
        print(f'adata2.X.max:{adata2.X.max()}, data is raw count')
        sc.pp.normalize_total(adata2, target_sum=1e4)
        sc.pp.log1p(adata2)
    if adata2.X.min() < 0.0:
        # it is scaled log1p data
        print(f'min:{adata2.X.min()}, data is already scaled')
        if adata2.raw is not None:
            print(f'restore from adata2.raw')
            # if we want to keep original adata, assign to a different object instead of adata
            adata2 = adata2.raw.to_adata()
            assert adata2.X.min() >= 0.0, f"min:{adata2.X.min()} seems also scaled, please use unscaled data"
            if adata2.X.max() < 20.0:
                print(f"max is {adata2.X.max()}, this is not raw count, seems log-normalized")
            else:
                sc.pp.normalize_total(adata2, target_sum=1e4)
                sc.pp.log1p(adata2)
        else:
            assert False, 'adata2.raw does not exist, please fix your data'

adata.obs['tag'] = 'control'
adata1.obs['tag'] = 'condition1'

#adata_combined = adata.concatenate(adata1, batch_key='tag', batch_categories=["control", "condition1"])
adata_combined = ad.concat([adata, adata1], label='tag', keys=["control", "condition1"])

# find marker genes
sc.tl.rank_genes_groups(adata_combined, groupby='tag', method='wilcoxon')

result = adata_combined.uns["rank_genes_groups"]
groups = result["names"].dtype.names
# write out full marker gene information
dat = pd.DataFrame({group + '_' + key[:1]: result[key][group] for group in groups for key in ['names', 'logfoldchanges','scores','pvals']})
dat.to_csv(f"{args.output_folder}/control_vs_condition1_marker_genes.csv")
print(f'full maker genes list is saved to {args.output_folder}/control_vs_condition1_marker_genes.csv')

# find significant markers for control group and condition 1 group, respectively
control1_sig_genes = dat.loc[(dat['control_l'] > 1) & (dat['control_p'] < 0.05), 'control_n'].tolist()
condition1_sig_genes = dat.loc[(dat['condition1_l'] > 1) & (dat['condition1_p'] < 0.05), 'condition1_n'].tolist()

# write out top 25 markers genes for control vs. condition 1
with open(f'{args.output_folder}/control_vs_condition1_top25_markers.txt', 'w') as f:
    top_features = {}
    n_top_genes = 25
    for group in groups:
        top_features[group] = result["names"][group][:n_top_genes]
    # print the top features for each cluster
    for group, features in top_features.items():
        f.write(f"Cluster {group} top features:\n")
        for feature in features:
            f.write(feature + '\n')
        f.write('\n')
print(f'top25 maker genes for each cell type group is saved to {args.output_folder}/control_vs_condition1_top25_markers.txt')

# dotplot top 10 genes for each cluster
top_genes = pd.DataFrame(result['names']).iloc[:10].melt()['value'].unique()
sc.pl.dotplot(adata_combined, var_names=top_genes, groupby='tag', show=False, figsize=(10, 6))
plt.savefig(f'{args.output_folder}/control_vs_condition1_top10_genes_dotplot.pdf', bbox_inches='tight')

if adata2 is not None:
    adata2.obs['tag'] = 'condition2'
#    adata_combined = adata.concatenate(adata2, batch_key='tag', batch_categories=["control", "condition2"])
    adata_combined = ad.concat([adata, adata2], label='tag', keys=["control", "condition2"])
    # find marker genes
    sc.tl.rank_genes_groups(adata_combined, groupby='tag', method='wilcoxon')

    result = adata_combined.uns["rank_genes_groups"]
    groups = result["names"].dtype.names
    # write out full marker gene information
    dat = pd.DataFrame({group + '_' + key[:1]: result[key][group] for group in groups for key in ['names', 'logfoldchanges','scores','pvals']})
    dat.to_csv(f"{args.output_folder}/control_vs_condition2_marker_genes.csv")
    print(f'full maker genes list is saved to {args.output_folder}/control_vs_condition2_marker_genes.csv')

    control2_sig_genes = dat.loc[(dat['control_l'] > 1) & (dat['control_p'] < 0.05), 'control_n'].tolist()
    condition2_sig_genes = dat.loc[(dat['condition2_l'] > 1) & (dat['condition2_p'] < 0.05), 'condition2_n'].tolist()

    # write out top 25 markers genes for each predicted cell type group
    with open(f'{args.output_folder}/control_vs_condition2_top25_markers.txt', 'w') as f:
        top_features = {}
        n_top_genes = 25
        for group in groups:
            top_features[group] = result["names"][group][:n_top_genes]
        # print the top features for each cluster
        for group, features in top_features.items():
            f.write(f"Cluster {group} top features:\n")
            for feature in features:
                f.write(feature + '\n')
            f.write('\n')
    print(f'top25 maker genes for each cell type group is saved to {args.output_folder}/control_vs_condition2_top25_markers.txt')

    # dotplot top 10 genes for each cluster
    top_genes = pd.DataFrame(result['names']).iloc[:10].melt()['value'].unique()
    sc.pl.dotplot(adata_combined, var_names=top_genes, groupby='tag', show=False, figsize=(10, 6))
    plt.savefig(f'{args.output_folder}/control_vs_condition2_top10_genes_dotplot.pdf', bbox_inches='tight')

    # find common significant marker genes between control_vs_condition1 and control_vs_condition2
    common_sig_markers1 = list(set(control1_sig_genes) & set(control2_sig_genes))
    print(f'common down-regulated significant genes are saved to {args.output_folder}/control_vs_conditions_common_sig_markers.txt')
    np.savetxt(f'{args.output_folder}/control_vs_conditions_common_sig_markers.txt', common_sig_markers1, fmt='%s')
    common_sig_markers2 = list(set(condition1_sig_genes) & set(condition2_sig_genes))
    print(f'common up-regulated significant genes are saved to {args.output_folder}/conditions_vs_control_common_sig_markers.txt')
    np.savetxt(f'{args.output_folder}/conditions_vs_control_common_sig_markers.txt', common_sig_markers2, fmt='%s')
    
