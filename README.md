Install packages:

pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118

pip install pandas scanpy anndata scipy performer_pytorch scikit-learn

To run inference:

python inference.py --gene_num 20000 --model_path /data/share/models/arabidopsis_hvg20k.ckpt --data_path /data/share/data/arabidopsis_SRP171_hvg20k.h5ad --log_file log_arabidopsis_SRP171.txt --prediction_file arabidopsis_SRP171_pred.csv

python inference.py --gene_num 10000 --model_path /data/share/models/corn_hvg10k.ckpt --data_path /data/share/data/corn_SRP335_hvg10k.h5ad --log_file log_corn_SRP335.txt --prediction_file corn_SRP335_pred.csv

python inference.py --gene_num 20000 --model_path /data/share/models/rice_hvg20k.ckpt --data_path /data/share/data/rice_SRP286_hvg20k.h5ad --log_file log_rice_SRP286.txt --prediction_file rice_SRP286_pred.csv

python inference.py --gene_num 20000 --model_path /data/share/models/soybean_hvg20k.ckpt --data_path /data/share/data/soybean_flowerbud_hvg20k.h5ad --log_file log_soybean_flowerbud.txt --prediction_file soybean_flowerbud_pred.csv

