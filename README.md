Install packages:
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118
pip install pandas scanpy anndata scipy performer_pytorch scikit-learn

To run inference:
python inference.py --model_path model-checkpoint-path --data_path test-file-in-h5ad-format --log_file log-file --prediction_file output-prediction-file
