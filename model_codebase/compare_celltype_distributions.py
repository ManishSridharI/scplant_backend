# -*- coding: utf-8 -*-
import os
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("--control_pred_file", type=str, required=True, help='Path of cell type prediction of control data')
parser.add_argument("--condition1_pred_file", type=str, required=True, help='Path of cell type prediction of condition 1 data')
parser.add_argument("--condition2_pred_file", type=str, default='', help='Path of cell type prediction of condition 2 data')
parser.add_argument("--output_folder", type=str, default='./results', help='output folder')

args = parser.parse_args()
os.makedirs(args.output_folder, exist_ok=True)

df_ctr = pd.read_csv(args.control_pred_file, header=0)
vc_ctr = df_ctr['scPlantAnnotate_celltype'].value_counts(normalize=True)

df_cond1 = pd.read_csv(args.condition1_pred_file, header=0)
vc_cond1 = df_cond1['scPlantAnnotate_celltype'].value_counts(normalize=True)

df_cond2 = None
if args.condition2_pred_file != '':
    df_cond2 = pd.read_csv(args.condition2_pred_file, header=0)
    vc_cond2 = df_cond2['scPlantAnnotate_celltype'].value_counts(normalize=True)

# Combine value counts into a single DataFrame
if df_cond2 is None:
    combined_counts = pd.concat([vc_ctr, vc_cond1], axis=1, keys=['control', 'condition1']).fillna(0)
else:
    combined_counts = pd.concat([vc_ctr, vc_cond1, vc_cond2], axis=1, keys=['control', 'condition1', 'condition2']).fillna(0)

# Ensure counts are floats
combined_counts = combined_counts.astype(float)
print(combined_counts)

# Plot bar chart
combined_counts.plot(kind='bar', figsize=(10, 6))
plt.title('Comparison of cell type ditributions across condtions')
plt.xlabel('cell type')
plt.ylabel('cell type distribution within each condition (%)')
#plt.show()
plt.savefig(f'{args.output_folder}/compare_celltype_distributions.pdf', bbox_inches='tight')