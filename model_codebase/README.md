Install packages:

pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118

pip install pandas scanpy anndata scipy performer_pytorch scikit-learn

(1) To run inference:

python inference.py --gene_num 20000 --model_path /data/share/models/arabidopsis_hvg20k.ckpt --data_path /data/share/data/arabidopsis/SRP171_hvg20k.h5ad --data_type h5ad --log_file log_arabidopsis_SRP171.txt --prediction_file arabidopsis_SRP171_pred.csv --stats_file arabidopsis_SRP171_stats.csv

python inference.py --gene_num 10000 --model_path /data/share/models/corn_hvg10k.ckpt --data_path /data/share/data/corn_SRP335_hvg10k.h5ad --data_type h5ad --log_file log_corn_SRP335.txt --prediction_file corn_SRP335_pred.csv --stats_file corn_SRP335_stats.csv

python inference.py --gene_num 20000 --model_path /data/share/models/rice_hvg20k.ckpt --data_path /data/share/data/rice_SRP286_hvg20k.h5ad --data_type h5ad --log_file log_rice_SRP286.txt --prediction_file rice_SRP286_pred.csv --stats_file rice_SRP286_stats.csv

python inference.py --gene_num 20000 --model_path /data/share/models/soybean_hvg20k.ckpt --data_path /data/share/data/soybean_flowerbud_hvg20k.h5ad --data_type h5ad --log_file log_soybean_flowerbud.txt --prediction_file soybean_flowerbud_pred.csv --stats_file soybean_flowerbud_stats.csv

If input format is 10x (cellranger output), use this: 
python inference.py --gene_num 20000 --model_path /data/share/models/arabidopsis_hvg20k.ckpt --data_path /path/to/cellranger/output/folder --data_type 10x --log_file log_arabidopsis_SRP171.txt --prediction_file arabidopsis_SRP171_pred.csv --stats_file arabidopsis_SRP171_stats.csv

(2) To run inference and plot results in tSNE/UMAP and marker genes dot plots.

python annotate_and_plot.py --gene_num 20000 --model_path /data/share/models/arabidopsis_hvg20k.ckpt --data_path /data/share/data/arabidopsis/SRP171_hvg20k.h5ad --data_type h5ad --output_folder logs/171
python annotate_and_plot.py --gene_num 20000 --model_path /data/share/models/arabidopsis_hvg20k.ckpt --data_path /data/share/data/arabidopsis/SRP235_hvg20k.h5ad --data_type h5ad --output_folder logs/235
python annotate_and_plot.py --gene_num 20000 --model_path /data/share/models/arabidopsis_hvg20k.ckpt --data_path /data/share/data/arabidopsis/SRP330_hvg20k.h5ad --data_type h5ad --output_folder logs/330

if input format is 10x (cellranger output), use this:
python annotate_and_plot.py --gene_num 20000 --model_path /data/share/models/arabidopsis_hvg20k.ckpt --data_path /path/to/cellranger/output/folder --data_type 10x --output_folder logs/171

(3) To compare data across different conditions: control and condition 1 (treatment 1), and (optionally) condition 2 (treatment 2)

python control_vs_treatment.py --control_data_path /data/share/data/arabidopsis/SRP171_hvg20k.h5ad --condition1_data_path /data/share/data/arabidopsis/SRP235_hvg20k.h5ad --condition2_data_path /data/share/data/arabidopsis/SRP330_hvg20k.h5ad --output_folder results


(4) To compare cell type distributions in different data: control and condition 1 (treatment 1), and (optionally) condition 2 (treatment 2)

python compare_celltype_distributions.py --control_pred_file logs/171/prediction.csv --condition1_pred_file logs/235/prediction.csv --condition2_pred_file logs/330/prediction.csv
