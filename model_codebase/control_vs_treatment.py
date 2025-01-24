# -*- coding: utf-8 -*-
import os
import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt

celltype_column="scPlantAnnotate_celltype"
cellid_column="Cell_Name"

def find_and_dump_DEGs(adata1, adata2, common_cell_types, df1, df2, IDs, n_top_genes, output_folder):

    group_df = {}
    sig_genes = {}
    for group in IDs:
        group_df[group] = {}
        sig_genes[group] = {}

    for ctype in common_cell_types:
        # map 
        sheet_name = ctype.replace(" ", "_")
        sheet_name = sheet_name.replace("/", "_")

        # Select cells of the given cell type from data1
        matching_cells = df1[df1[celltype_column] == ctype][cellid_column].tolist()
        adata1_subset = adata1[adata1.obs_names.isin(matching_cells)].copy()
        if adata1_subset.shape[0] < 3:
            print(f'control has only {adata1_subset.shape[0]} cells in cell type {ctype}, not sufficient for DEG calculation, skip it' )
            continue

        # Select cells of the given cell type from data2
        matching_cells = df2[df2[celltype_column] == ctype][cellid_column].tolist()
        adata2_subset = adata2[adata2.obs_names.isin(matching_cells)].copy()
        if adata2_subset.shape[0] < 3:
            print(f'treatment has only {adata2_subset.shape[0]} cells in cell type {ctype}, not sufficient for DEG calculation, skip it' )
            continue

        # combine them into one object
        adata_combined = ad.concat([adata1_subset, adata2_subset], label='tag', keys=IDs)  # ["control", "condition1"]

        sc.tl.rank_genes_groups(adata_combined, groupby='tag', method='wilcoxon')
        rank_genes_groups = adata_combined.uns["rank_genes_groups"]
        groups = rank_genes_groups["names"].dtype.names

        dat = pd.DataFrame({group + '_' + key: rank_genes_groups[key][group] for group in groups for key in ['names', 'logfoldchanges','scores','pvals']})

        # write out top DEGs for each group: control/condition1
        for g_idx, group in enumerate(groups):  # group = control | condition1
            assert group == IDs[g_idx], f"g_idx={g_idx} group={group} does not match IDs[{g_idx}]={IDs[g_idx]}"
            if n_top_genes == 0:
                group_df[group][sheet_name] = pd.DataFrame({
                    "gene": rank_genes_groups["names"][group],
                    "logfoldchange": rank_genes_groups["logfoldchanges"][group],
                    "score": rank_genes_groups["scores"][group],
                    "p_values": rank_genes_groups["pvals"][group]
                })
            else:
                group_df[group][sheet_name] = pd.DataFrame({
                    "gene": rank_genes_groups["names"][group][:n_top_genes],
                    "logfoldchange": rank_genes_groups["logfoldchanges"][group][:n_top_genes],
                    "score": rank_genes_groups["scores"][group][:n_top_genes],
                    "p_values": rank_genes_groups["pvals"][group][:n_top_genes]
                })

            sig_genes[group][ctype] = dat.loc[(dat[f'{group}_logfoldchanges'] > 1) & (dat[f'{group}_pvals'] < 0.05), f'{group}_names'].tolist()

    # write out DEGs in xlsx for each group (specified in IDs)
    output_path = f'{output_folder}/{IDs[0]}_vs_{IDs[1]}'
    os.makedirs(output_path, exist_ok=True)
    for group in groups:
        if n_top_genes == 0:
            excel_file = f'{output_path}/{group}_all_DEGs.xlsx'
        else:
            excel_file = f'{output_path}/{group}_top{n_top_genes}_DEGs.xlsx'
        # Write DataFrames to Excel, each in its own sheet
        with pd.ExcelWriter(excel_file, engine="xlsxwriter") as writer:
            for sheet_name, df in group_df[group].items():
                df.to_excel(writer, sheet_name=sheet_name, index=False)
        print(f"saved DEGs to {excel_file}")

    return sig_genes

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
parser.add_argument("--data_type", type=str, default='h5ad', help='type of input data (h5ad|10x)')
parser.add_argument("--control_data_path", type=str, required=True, help='Path of h5ad file for control')
parser.add_argument("--control_pred_file", type=str, required=True, help='Path of cell type prediction of control data')
parser.add_argument("--condition1_data_path", type=str, required=True, help='Path of h5ad file for condition 1')
parser.add_argument("--condition1_pred_file", type=str, required=True, help='Path of cell type prediction of condition 1 data')
parser.add_argument("--condition2_data_path", type=str, default='', help='Path of h5ad file for condition 2')
parser.add_argument("--condition2_pred_file", type=str, default='', help='Path of cell type prediction of condition 2 data')
parser.add_argument("--output_folder", type=str, default='./results', help='output folder')

args = parser.parse_args()
os.makedirs(args.output_folder, exist_ok=True)

# load control data and its prediction
adata = load_data(args.control_data_path, args.data_type)
print(f'cells: {adata.X.shape[0]} genes: {adata.X.shape[1]}')
adata.obs['tag'] = 'control'
df_ctr = pd.read_csv(args.control_pred_file, header=0)
vc_ctr = df_ctr[celltype_column].value_counts(normalize=True)

# load condition1 data and its prediction
adata1 = load_data(args.condition1_data_path, args.data_type)
print(f'cells: {adata1.X.shape[0]} genes: {adata1.X.shape[1]}')
adata1.obs['tag'] = 'condition1'
df_cond1 = pd.read_csv(args.condition1_pred_file, header=0)
vc_cond1 = df_cond1[celltype_column].value_counts(normalize=True)

# load condition2 data and its prediction if provided
adata2 = None
if args.condition2_data_path != '':
    adata2 = load_data(args.condition2_data_path, args.data_type)
    print(f'cells: {adata2.X.shape[0]} genes: {adata2.X.shape[1]}')
    df_cond2 = pd.read_csv(args.condition2_pred_file, header=0)
    vc_cond2 = df_cond2[celltype_column].value_counts(normalize=True)

# Combine value counts into a single DataFrame
if df_cond2 is None:
    combined_counts = pd.concat([vc_ctr, vc_cond1], axis=1, keys=['control', 'condition1']).fillna(0)
else:
    combined_counts = pd.concat([vc_ctr, vc_cond1, vc_cond2], axis=1, keys=['control', 'condition1', 'condition2']).fillna(0)
combined_counts = combined_counts.astype(float)
print(combined_counts)
# Plot bar chart
combined_counts.plot(kind='bar', figsize=(10, 6))
plt.title('Comparison of cell type ditributions across condtions')
plt.xlabel('cell type')
plt.ylabel('cell type distribution within each condition (%)')
#plt.show()
plt.savefig(f'{args.output_folder}/compare_celltype_distributions.pdf', bbox_inches='tight')

# find common cell types in control and condition_1
common_cell_types = list(set(df_ctr[celltype_column]).intersection(set(df_cond1[celltype_column])))
if len(common_cell_types) > 0:
    # find DEGs for each common cell type
    group_df_list = []
    n_top_genes_list = [0, 25, 10, 5]   # 0 means 'all'
    for n_top_genes in (n_top_genes_list):
        sig_genes = find_and_dump_DEGs(adata, adata1, common_cell_types, df_ctr, df_cond1, ["control", "condition1"], n_top_genes, args.output_folder)
else:
    print("there is no common cell types in control and condition_1")

if adata2 is not None:
    adata2.obs['tag'] = 'condition2'
    # find common cell types in control and condition_2
    common_cell_types_2 = list(set(df_ctr[celltype_column]).intersection(set(df_cond2[celltype_column])))
    if len(common_cell_types_2) > 0:
        # find DEGs for each common cell type
        group_df_list = []
        n_top_genes_list = [0, 25, 10, 5]   # 0 means 'all'
        for n_top_genes in (n_top_genes_list):
            sig_genes2 = find_and_dump_DEGs(adata, adata2, common_cell_types_2, df_ctr, df_cond2, ["control", "condition2"], n_top_genes, args.output_folder)
    else:
        print("there is no common cell types in control and condition_2")

    # find common cell types in all groups (control/condition1/condition2)
    common_all = list(set(common_cell_types) & set(common_cell_types_2))

    # find and dump common DOWN-regulated signficant genes (for each cell type)
    group_df = {}
    for ctype in common_all:
        if ctype not in sig_genes['control']:
            print(f'did not find {ctype} in sig_genes[control], skip')
            continue
        if ctype not in sig_genes2['control']:
            print(f'did not find {ctype} in sig_genes2[control], skip')
            continue

        # find common significant marker genes between control_vs_condition1 and control_vs_condition2
        common_sig_markers = list(set(sig_genes["control"][ctype]) & set(sig_genes2["control"][ctype]))

        sheet_name = ctype.replace(" ", "_")
        sheet_name = sheet_name.replace("/", "_")
        group_df[sheet_name] = pd.DataFrame({
            "gene": common_sig_markers,
        })
    excel_file = f'{args.output_folder}/control_vs_conditions_common_sig_markers.xlsx'
    with pd.ExcelWriter(excel_file, engine="xlsxwriter") as writer:
        for sheet_name, df in group_df.items():
            df.to_excel(writer, sheet_name=sheet_name, index=False, header=False)
    print(f"common down-regulated significant genes are saved to {excel_file}")

    # find and dump common UP-regulated signficant genes (for each cell type)
    group_df = {}
    for ctype in common_all:
        if ctype not in sig_genes['condition1']:
            print(f'did not find {ctype} in sig_genes[condition1], skip')
            continue
        if ctype not in sig_genes2['condition2']:
            print(f'did not find {ctype} in sig_genes2[condition2], skip')
            continue

        # find common significant marker genes between condition1_vs_condition and condition2_vs_control
        common_sig_markers = list(set(sig_genes["condition1"][ctype]) & set(sig_genes2["condition2"][ctype]))

        sheet_name = ctype.replace(" ", "_")
        sheet_name = sheet_name.replace("/", "_")
        group_df[sheet_name] = pd.DataFrame({
            "gene": common_sig_markers,
        })
    excel_file = f'{args.output_folder}/conditions_vs_control_common_sig_markers.xlsx'
    with pd.ExcelWriter(excel_file, engine="xlsxwriter") as writer:
        for sheet_name, df in group_df.items():
            df.to_excel(writer, sheet_name=sheet_name, index=False, header=False)
    print(f"common up-regulated significant genes are saved to {excel_file}")
