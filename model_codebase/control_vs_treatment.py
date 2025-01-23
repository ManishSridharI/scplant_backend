# -*- coding: utf-8 -*-
import os
import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt

def load_data(data_path, data_type):
    if data_type == 'h5ad':
        if os.path.isfile(data_path):
            adata = sc.read_h5ad(data_path)
        else:
            assert False, f'{data_path} does not exist or is not a file'
    elif data_type == '10x':
        if os.path.isdir(data_path):
            adata = sc.read_10x_mtx(data_path, var_names="gene_symbols")
        else:
            assert False, f'{data_path} does not exist or is not a directory'
    else:
        assert False, f"unsupported data type: {data_type}"

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

    return adata

parser = argparse.ArgumentParser()
parser.add_argument("--control_data_path", type=str, required=True, help='Path of h5ad file for control')
parser.add_argument("--condition1_data_path", type=str, required=True, help='Path of h5ad file for condition 1')
parser.add_argument("--condition2_data_path", type=str, default='', help='Path of h5ad file for condition 2')
parser.add_argument("--data_type", type=str, default='h5ad', help='type of input data (h5ad|10x)')
parser.add_argument("--output_folder", type=str, default='./results', help='output folder')

args = parser.parse_args()
os.makedirs(args.output_folder, exist_ok=True)

# load control data
adata = load_data(args.control_data_path, args.data_type)
print(f'cells: {adata.X.shape[0]} genes: {adata.X.shape[1]}')
adata.obs['tag'] = 'control'

# load condition1 data
adata1 = load_data(args.condition1_data_path, args.data_type)
print(f'cells: {adata1.X.shape[0]} genes: {adata1.X.shape[1]}')
adata1.obs['tag'] = 'condition1'

# load condition2 data if provided
adata2 = None
if args.condition2_data_path != '':
    adata2 = load_data(args.condition2_data_path, args.data_type)
    print(f'cells: {adata2.X.shape[0]} genes: {adata2.X.shape[1]}')

#adata_combined = adata.concatenate(adata1, batch_key='tag', batch_categories=["control", "condition1"])
adata_combined = ad.concat([adata, adata1], label='tag', keys=["control", "condition1"])

# find marker genes
sc.tl.rank_genes_groups(adata_combined, groupby='tag', method='wilcoxon')
rank_genes_groups = adata_combined.uns["rank_genes_groups"]
groups = rank_genes_groups["names"].dtype.names

# write out full marker gene information
dat = pd.DataFrame({group + '_' + key: rank_genes_groups[key][group] for group in groups for key in ['names', 'logfoldchanges','scores','pvals']})
dat.to_csv(f"{args.output_folder}/control_vs_condition1_marker_genes.csv")
print(f'full maker genes list is saved to {args.output_folder}/control_vs_condition1_marker_genes.csv')

# write out all and top marker genes for control/condition1
output_marker_folder=f'{args.output_folder}/control_vs_condition1'
os.makedirs(output_marker_folder, exist_ok=True)
# dump all marker genes for control/condition1
for group in groups:
    group_df = pd.DataFrame({
        "gene": rank_genes_groups["names"][group],
        "logfoldchange": rank_genes_groups["logfoldchanges"][group],
        "score": rank_genes_groups["scores"][group],
        "p_values": rank_genes_groups["pvals"][group]
    })
    file_name = f"{output_marker_folder}/{group}_all_marker_genes.csv"
    group_df.to_csv(file_name, index=False)

for n_top_genes in (25, 10, 5):
    for group in groups:
        group_df = pd.DataFrame({
            "gene": rank_genes_groups["names"][group][:n_top_genes],
            "logfoldchange": rank_genes_groups["logfoldchanges"][group][:n_top_genes],
            "score": rank_genes_groups["scores"][group][:n_top_genes],
            "p_values": rank_genes_groups["pvals"][group][:n_top_genes]
        })
        file_name = f"{output_marker_folder}/{group}_top{n_top_genes}_marker_genes.csv"
        group_df.to_csv(file_name, index=False)

# find significant markers for control group and condition 1 group, respectively
control1_sig_genes = dat.loc[(dat['control_logfoldchanges'] > 1) & (dat['control_pvals'] < 0.05), 'control_names'].tolist()
condition1_sig_genes = dat.loc[(dat['condition1_logfoldchanges'] > 1) & (dat['condition1_pvals'] < 0.05), 'condition1_names'].tolist()

# dotplot top 10 genes for each condition
top_genes = pd.DataFrame(rank_genes_groups['names']).iloc[:10].melt()['value'].unique()
sc.pl.dotplot(adata_combined, var_names=top_genes, groupby='tag', show=False, figsize=(10, 6))
plt.savefig(f'{args.output_folder}/control_vs_condition1_top10_genes_dotplot.pdf', bbox_inches='tight')
#plt.savefig(f'{args.output_folder}/control_vs_condition1_top10_genes_dotplot.png', bbox_inches='tight', dpi=300)

if adata2 is not None:
    adata2.obs['tag'] = 'condition2'
#    adata_combined = adata.concatenate(adata2, batch_key='tag', batch_categories=["control", "condition2"])
    adata_combined = ad.concat([adata, adata2], label='tag', keys=["control", "condition2"])
    # find marker genes
    sc.tl.rank_genes_groups(adata_combined, groupby='tag', method='wilcoxon')

    rank_genes_groups = adata_combined.uns["rank_genes_groups"]
    groups = rank_genes_groups["names"].dtype.names
    # write out full marker gene information
    dat = pd.DataFrame({group + '_' + key: rank_genes_groups[key][group] for group in groups for key in ['names', 'logfoldchanges','scores','pvals']})
    dat.to_csv(f"{args.output_folder}/control_vs_condition2_marker_genes.csv")
    print(f'full maker genes list is saved to {args.output_folder}/control_vs_condition2_marker_genes.csv')

    control2_sig_genes = dat.loc[(dat['control_logfoldchanges'] > 1) & (dat['control_pvals'] < 0.05), 'control_names'].tolist()
    condition2_sig_genes = dat.loc[(dat['condition2_logfoldchanges'] > 1) & (dat['condition2_pvals'] < 0.05), 'condition2_names'].tolist()

    # write out top marker genes for control/condition2
    output_marker_folder=f'{args.output_folder}/control_vs_condition2'
    os.makedirs(output_marker_folder, exist_ok=True)
    for n_top_genes in (25, 10, 5):
        for group in groups:
            group_df = pd.DataFrame({
                "gene": rank_genes_groups["names"][group][:n_top_genes],
                "logfoldchange": rank_genes_groups["logfoldchanges"][group][:n_top_genes],
                "score": rank_genes_groups["scores"][group][:n_top_genes],
                "p_values": rank_genes_groups["pvals"][group][:n_top_genes]
            })
            group = group.replace(" ", "_")
            group = group.replace("/", "_")
            file_name = f"{output_marker_folder}/{group}_top{n_top_genes}_marker_genes.csv"
            group_df.to_csv(file_name, index=False)

    # dotplot top 10 genes for each condition
    top_genes = pd.DataFrame(rank_genes_groups['names']).iloc[:10].melt()['value'].unique()
    sc.pl.dotplot(adata_combined, var_names=top_genes, groupby='tag', show=False, figsize=(10, 6))
    plt.savefig(f'{args.output_folder}/control_vs_condition2_top10_genes_dotplot.pdf', bbox_inches='tight')
    #plt.savefig(f'{args.output_folder}/control_vs_condition2_top10_genes_dotplot.png', bbox_inches='tight', dpi=300)

    # find common significant marker genes between control_vs_condition1 and control_vs_condition2
    common_sig_markers1 = list(set(control1_sig_genes) & set(control2_sig_genes))
    print(f'common down-regulated significant genes are saved to {args.output_folder}/control_vs_conditions_common_sig_markers.txt')
    np.savetxt(f'{args.output_folder}/control_vs_conditions_common_sig_markers.txt', common_sig_markers1, fmt='%s')
    common_sig_markers2 = list(set(condition1_sig_genes) & set(condition2_sig_genes))
    print(f'common up-regulated significant genes are saved to {args.output_folder}/conditions_vs_control_common_sig_markers.txt')
    np.savetxt(f'{args.output_folder}/conditions_vs_control_common_sig_markers.txt', common_sig_markers2, fmt='%s')
