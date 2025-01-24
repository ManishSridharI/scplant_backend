Install packages:

pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118
pip install pandas scanpy anndata scipy performer_pytorch scikit-learn
pip install xlsxwriter

(1) To run inference and plot results in tSNE/UMAP and marker genes dot plots.

Note that input model are changed to use a new version that include genes list so that the gene number is inferred from the model. No need to specify gene_num on the command line.

python annotate_and_plot.py --model_path /data/SHARE/models/arabidopsis_hvg20k_add_genes.ckpt --data_path /data/SHARE/data/arabidopsis/SRP171_hvg20k.h5ad --data_type h5ad --output_folder logs/171
python annotate_and_plot.py --model_path /data/SHARE/models/arabidopsis_hvg20k_add_genes.ckpt --data_path /data/SHARE/data/arabidopsis/SRP235_hvg20k.h5ad --data_type h5ad --output_folder logs/235
python annotate_and_plot.py --model_path /data/SHARE/models/arabidopsis_hvg20k_add_genes.ckpt --data_path /data/SHARE/data/arabidopsis/SRP330_hvg20k.h5ad --data_type h5ad --output_folder logs/330

If input format is 10x (cellranger output), use this:
python annotate_and_plot.py --model_path /data/SHARE/models/arabidopsis_hvg20k_add_genes.ckpt --data_path /path/to/cellranger/output/folder --data_type 10x --output_folder logs/171

note that --data_path must be a folder which include 3 files (barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz) which are default filenames produced by cellranger

Output files are (use logs/171 as example):
  prediction.csv:                           predicted cell types for each cell
  stats.csv and stats.pdf:                  counts for predicted cell types and plot
  marker_genes.csv:                         full list of marker genes for all predicted cell type group
  annotate_tsne.pdf and annotate_umap.pdf:  t-SNE and UMAP plot of all cells grouped by predicted cell types
  top3_genes_dotplot.pdf:                   dotplot for top-3 marker genes for each cell group
  all_marker_genes_by_group:                folder containing full marker gene list for each cell type group
  top25_marker_genes_by_group:              folder containing top 25 marker genes for each cell type group
  top10_marker_genes_by_group:              folder containing top 10 marker genes for each cell type group
  top5_marker_genes_by_group:               folder containing top 5 marker genes for each cell type group

(2) To compare data across different conditions: control and condition 1 (treatment 1), and (optionally) condition 2 (treatment 2)

python control_vs_treatment.py --control_data_path /data/SHARE/data/arabidopsis/SRP171_hvg20k.h5ad --condition1_data_path /data/SHARE/data/arabidopsis/SRP235_hvg20k.h5ad --condition2_data_path /data/SHARE/data/arabidopsis/SRP330_hvg20k.h5ad --output_folder results

if input data format is 10x (cellranger output), use this:

python control_vs_treatment.py --data_type 10x --control_data_path /path/to/cellranger/control_data_folder --condition1_data_path /path/to/cellranger/condition1_data_folder --condition2_data_path /path/to/cellranger/condition2_data_folder --output_folder results

Outut files are:
  control_vs_condition1_marker_genes.csv:           full list of marker genes comparing control data vs. condition 1 data
  control_vs_condition1_top10_genes_dotplot.pdf:    dot plot figure of top 10 marker genes in control data and condition 1 data, respectively
  control_vs_condition1:                            subfolder containing all, top 25, 10, 5 marker genes in control data and condition 1 data, respectively

  control_vs_condition2_marker_genes.csv:           full list of marker genes comparing control data vs. condition 2 data
  control_vs_condition2_top10_genes_dotplot.pdf:    dot plot figure of top 10 marker genes in control data and condition 2 data, respectively
  control_vs_condition2:                            subfolder containing all, top 25, 10, 5 marker genes in control data and condition 2 data, respectively

  control_vs_conditions_common_sig_markers.txt:     the shared genes in significant marker genes in control_vs_condition1 comparison and control_vs_condition2 comparison, which means the list of DOWN-regulated genes in both conditions 
  conditions_vs_control_common_sig_markers.txt:     the shared genes in significant marker genes in condition1_vs_control comparison and condition2_vs_control comparison, which means the list of UP-regulated genes in both conditions 

(3) To compare cell type distributions in different data: control and condition 1 (treatment 1), and (optionally) condition 2 (treatment 2)

python compare_celltype_distributions.py --control_pred_file logs/171/prediction.csv --condition1_pred_file logs/235/prediction.csv --condition2_pred_file logs/330/prediction.csv
