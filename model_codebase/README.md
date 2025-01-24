# Install packages:

`pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118`  
`pip install pandas scanpy anndata scipy performer_pytorch scikit-learn`  
`pip install xlsxwriter`  

# 1. Predict cell types and find marker genes for each predicted cell type

Note that input model are changed to use a new version that include genes list so that the gene number is inferred from the model. No need to specify gene_num on the command line.

`python annotate_and_plot.py --model_path /data/SHARE/models/arabidopsis_hvg20k_add_genes.ckpt --data_path /data/SHARE/data/arabidopsis/SRP171_hvg20k.h5ad --data_type h5ad --output_folder logs/171`  
`python annotate_and_plot.py --model_path /data/SHARE/models/arabidopsis_hvg20k_add_genes.ckpt --data_path /data/SHARE/data/arabidopsis/SRP235_hvg20k.h5ad --data_type h5ad --output_folder logs/235`  
`python annotate_and_plot.py --model_path /data/SHARE/models/arabidopsis_hvg20k_add_genes.ckpt --data_path /data/SHARE/data/arabidopsis/SRP330_hvg20k.h5ad --data_type h5ad --output_folder logs/330`  

If input format is 10x (cellranger output), use this:  
`python annotate_and_plot.py --model_path /data/SHARE/models/arabidopsis_hvg20k_add_genes.ckpt --data_path /path/to/cellranger/output/folder --data_type 10x --output_folder logs/171  `

note that --data_path must be a folder which include 3 files (barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz) which are default filenames produced by cellranger

Output files are (use logs/171 as example):  
| output file(s)  | Description     |
|--------------|---------------------|
| prediction.csv | predicted cell types for each cell |
| stats.csv and stats.pdf | counts for predicted cell types and plot |
| marker_genes.csv | full list of marker genes for all predicted cell type group |
| annotate_tsne.pdf and annotate_umap.pdf | t-SNE and UMAP plot of all cells grouped by predicted cell types |
| top3_genes_dotplot.pdf | dotplot for top-3 marker genes for each cell group | 
| all_marker_genes.xlsx |                   excel file containing full marker genes for each cell type, each in a separate sheet |
| top25_marker_genes.xlsx |                  excel file containing top 25 marker genes for each cell type, each in a separate sheet |
| top10_marker_genes.xlsx |                 excel file containing top 10 marker genes for each cell type, each in a separate sheet |
| top5_marker_genes.xlsx |                 excel file containing top  5 marker genes for each cell type, each in a separate sheet |

# 2. Compare data across different conditions: control and condition 1 (treatment 1), and (optionally) condition 2 (treatment 2)

NOTE: this tool requires prediction results from annotate_and_plot.py tool. So make sure run that tool beforehand.

`python control_vs_treatment.py --control_data_path /data/SHARE/data/arabidopsis/SRP171_hvg20k.h5ad --condition1_data_path /data/SHARE/data/arabidopsis/SRP235_hvg20k.h5ad --condition2_data_path /data/SHARE/data/arabidopsis/SRP330_hvg20k.h5ad --control_pred_file logs/171.new/prediction.csv --condition1_pred_file logs/235/prediction.csv --condition2_pred_file logs/330/prediction.csv --output_folder results`  

if input data format is 10x (cellranger output), use this:

`python control_vs_treatment.py --data_type 10x --control_data_path /path/to/cellranger/control_data_folder --condition1_data_path /path/to/cellranger/condition1_data_folder --condition2_data_path /path/to/cellranger/condition2_data_folder --control_pred_file path_to_control_data_pred_csv_file --condition1_pred_file path_to_condition1_data_pred_csv_file --condition2_pred_file path_to_control2_data_pred_csv_file --output_folder results`  

Outut files are:
| output file(s)  | Description     |
|--------------|---------------------|
| compare_celltype_distributions.pdf |             plot of cell type distributions in all conditions (control, condition1, and condition2 if provided) |
|  control_vs_condition1 |                          subfolder containing all, top 25, 10, 5 DEGs (differentially expressed genes) in control data and condition 1 data, respectively. Each excel file contains multiple sheets, one for each cell type. |
  
  If condition2 data and predictions are provided, a few more output files:  

  | output file(s)  | Description     |
|--------------|---------------------|
  |  control_vs_condition2 |                          subfolder containing all, top 25, 10, 5 DEGs (differentially expressed genes) in control data and condition 2 data, respectively. Each excel file contains multiple sheets, one for each cell type. |
  | control_vs_conditions_common_sig_markers.xlsx | shared genes in significant DEGs in control_vs_condition1 comparison and control_vs_condition2 comparison, which means the list of DOWN-regulated genes in both conditions, grouped by cell type, each one on a separate sheet. |
 | conditions_vs_control_common_sig_markers.txt: | shared genes in significant DEGs in condition1_vs_control comparison and condition2_vs_control comparison, which means the list of UP-regulated genes in both conditions, grouped by cell type, each one on a separate sheet. |
