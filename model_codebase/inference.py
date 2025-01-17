# -*- coding: utf-8 -*-
import os

os.environ['LD_LIBRARY_PATH'] = '/home/clcdp/data/miniconda3_space/envs/jupyterlab/lib'

import gc
import argparse
import json
import random
import math
import random
from functools import reduce
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.model_selection import train_test_split, ShuffleSplit, StratifiedShuffleSplit, StratifiedKFold
from sklearn.metrics import accuracy_score, f1_score, confusion_matrix, precision_recall_fscore_support, classification_report
import torch
from torch import nn
from torch.optim import Adam, SGD, AdamW
from torch.nn import functional as F
from torch.optim.lr_scheduler import StepLR, CosineAnnealingWarmRestarts, CyclicLR
from torch.utils.data import DataLoader, Dataset
from torch.utils.data.distributed import DistributedSampler
from torch.nn.parallel import DistributedDataParallel as DDP
#import torch.distributed as dist

from performer_pytorch import PerformerLM
import scanpy as sc
import anndata as ad
from utils import *
import pickle as pkl

class FlushFileHandler(logging.FileHandler):
    def emit(self, record):
        super().emit(record)
        self.flush()

def configure_logging(log_file):
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    # Remove any existing handlers
    for handler in logger.handlers[:]:
        logger.removeHandler(handler)

    # Add custom file handler that flushes immediately
    file_handler = FlushFileHandler(log_file)
    file_handler.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

parser = argparse.ArgumentParser()
parser.add_argument("--local_rank", type=int, default=0, help='Local process rank.')
parser.add_argument("--bin_num", type=int, default=5, help='Number of bins.')
parser.add_argument("--gene_num", type=int, default=20000, help='Number of genes.')
parser.add_argument("--embed_dim", type=int, default=150, help='embedding dimension.')
parser.add_argument("--num_layer", type=int, default=5, help='Number of transformer layers.')
parser.add_argument("--batch_size", type=int, default=1, help='Number of batch size.')
parser.add_argument("--pos_embed", action='store_true', help='Using Gene2vec encoding or not.')
parser.add_argument("--mask_class", action='store_true', help='use only cell types in the test data.')
parser.add_argument("--data_path", type=str, default='./data/test/SRP171040_hvg20k_test.h5ad', help='Path of data for evaluation.')
parser.add_argument("--data_type", type=str, default='h5ad', help='type of input data (h5ad|10x)')
parser.add_argument("--celltype_column", type=str, default='', help='celltype column in h5ad if available')
parser.add_argument("--model_path", type=str, default='./ckpts_arabidopsis_ft/pt_on_all_ft_on_all_noval_best.pth', help='Path of best finetuned model.')
parser.add_argument("--log_file", type=str, default='./logs/log.txt', help='log file path')
parser.add_argument("--prediction_file", type=str, default='./logs/pred.csv', help='cell type prediction file path')
parser.add_argument("--stats_file", type=str, default='./logs/stats.csv', help='prediction stats file path')

args = parser.parse_args()
configure_logging(args.log_file)
logging.info("Starting the evaluation process")
logging.info(args)

#rank = int(os.environ["RANK"])
local_rank = args.local_rank
BATCH_SIZE = args.batch_size
SEQ_LEN = args.gene_num + 1
UNASSIGN_THRES = 0.0
CLASS = args.bin_num + 2
POS_EMBED_USING = args.pos_embed

torch.cuda.set_device(local_rank)
device = torch.device("cuda", local_rank)

class SCDataset(Dataset):
    def __init__(self, data, label, names):
        super().__init__()
        self.data = data
        self.label = label
        self.names = names

    def __getitem__(self, index):
        index = index % self.data.shape[0]
        full_seq = self.data[index].toarray()[0]
        full_seq[full_seq > (CLASS - 2)] = CLASS - 2
        full_seq = torch.from_numpy(full_seq).long()
        full_seq = torch.cat((full_seq, torch.tensor([0]))).to(device)
        seq_label = self.label[index]
        name = self.names[index]
        return full_seq, seq_label, name

    def __len__(self):
        return self.data.shape[0]

class Identity(torch.nn.Module):
    def __init__(self, seq_len, embed_dim, dropout = 0., h_dim = 100, out_dim = 10):
        super(Identity, self).__init__()
        self.conv1 = nn.Conv2d(1, 1, (1, embed_dim))
        self.act = nn.ReLU()
        self.fc1 = nn.Linear(in_features=seq_len, out_features=512, bias=True)
        self.act1 = nn.ReLU()
        self.dropout1 = nn.Dropout(dropout)
        self.fc2 = nn.Linear(in_features=512, out_features=h_dim, bias=True)
        self.act2 = nn.ReLU()
        self.dropout2 = nn.Dropout(dropout)
        self.fc3 = nn.Linear(in_features=h_dim, out_features=out_dim, bias=True)

    def forward(self, x):
        x = x[:,None,:,:]
        x = self.conv1(x)
        x = self.act(x)
        x = x.view(x.shape[0],-1)
        x = self.fc1(x)
        x = self.act1(x)
        x = self.dropout1(x)
        x = self.fc2(x)
        x = self.act2(x)
        x = self.dropout2(x)
        x = self.fc3(x)
        return x

if args.data_type == 'h5ad':
    if os.path.isfile(args.data_path):
        adata = sc.read_h5ad(args.data_path)
    else:
        assert False, f'{args.data_path} is not a file'
elif args.data_type == '10x':
    if os.path.isdir(args.data_path):
        adata = sc.read_10x_mtx(args.data_path, var_names="gene_symbols")
    else:
        assert False, f'{args.data_path} is not a directory'
else:
    assert False, f"unsupported data type: {args.data_type}"

num_cells = adata.X.shape[0]
#print(adata)
path = args.model_path
ckpt = torch.load(path, map_location='cpu')

if 'output_node_names' in ckpt:
    label_dict_list = ckpt['output_node_names']
else:
    logging.error(f"did not find output_node_names in {args.model_path}")
    assert False, f"did not find output_node_names in {args.model_path}"

reverse_dict = {ctype: idx for idx, ctype in enumerate(label_dict_list)}
has_label = False
if args.celltype_column != '':
    print(f'reference cell type is in {args.celltype_column} column')
    if args.celltype_column in adata.obs:
        adata.obs['CelltypeID'] = adata.obs[args.celltype_column].map(reverse_dict)
        if adata.obs['CelltypeID'].isnull().any():
            adata.obs['CelltypeID'] = adata.obs['CelltypeID'].fillna(-1).astype(int)
        else:
            adata.obs['CelltypeID'] = adata.obs['CelltypeID'].astype(int)
        label_list = list(range(len(label_dict_list)))
        label = torch.from_numpy(adata.obs['CelltypeID'].to_numpy())
        has_label = True
    else:
        print(f"did not find {args.celltype_column} column in {args.data_path}, skip it")
        label = torch.zeros(num_cells)
else:
    # no celltype available, fake some
    label = torch.zeros(num_cells)

cell_names = adata.obs.index.tolist()
#print(cell_names)

data = adata.X

acc = []
f1 = []
f1w = []
val_dataset = SCDataset(data, label, cell_names)
val_loader = DataLoader(val_dataset, batch_size=BATCH_SIZE)
print(f"Number of samples in val_dataset: {len(val_dataset)}")
logging.info(f"Number of samples in val_dataset: {len(val_dataset)}")

model = PerformerLM(
    num_tokens = CLASS,
    dim = args.embed_dim,
    depth = args.num_layer,
    max_seq_len = SEQ_LEN,
    heads = 10,
    local_attn_heads = 0,
    g2v_position_emb = POS_EMBED_USING
)
nclass = len(label_dict_list)
model.to_out = Identity(args.gene_num + 1, args.embed_dim, dropout=0., h_dim=128, out_dim=nclass)
model.load_state_dict(ckpt['model_state_dict'])
model = model.to(device)
loss_fn = nn.CrossEntropyLoss(weight=None).to(device)
softmax = nn.Softmax(dim=-1)

if args.mask_class:
    val_clist = adata.obs['Celltype'].unique().tolist()
    mask = [1 if i in val_clist else 0 for i in label_dict_list]
    mask = torch.tensor(mask).to(device)
    nclass_masked = len(val_clist)
    print(f'mask prediction classes from {nclass} to {nclass_masked}')

model.eval()
running_loss = 0.0
predictions = []
truths = []
nskip = 0
cells = []
pred_ctypes = []
with torch.no_grad():
    for index, (data_v, labels_v, names_v) in enumerate(val_loader):
        index += 1
        data_v, labels_v = data_v.to(device), labels_v.to(device)
#        if labels_v.item() < 0 or labels_v.item() >= nclass:
#            print(f'index={index} label={labels_v.item()} is out of bound, skip this sample')
#            nskip += 1
#            continue
        logits = model(data_v)
        final_prob = softmax(logits)
        if args.mask_class:
            final_prob = final_prob * mask
        final = final_prob.argmax(dim=-1)
        final[np.amax(np.array(final_prob.cpu()), axis=-1) < UNASSIGN_THRES] = -1
        predictions.append(final)
        truths.append(labels_v)
        cells.extend(names_v)
        pred_ctypes.extend([label_dict_list[i] for i in final])
        if index % 1000 == 0:
            print(f'index={index}')
    del data_v, labels_v, names_v, logits, final_prob, final
    # gather
    no_drop = predictions != -1
    if no_drop:
        predictions = np.array([tensor.cpu().detach().numpy() for tensor in predictions])
        truths = np.array([tensor.cpu().detach().numpy() for tensor in truths])
    else:
        predictions = np.array((predictions[no_drop]).cpu())
        truths = np.array((truths[no_drop]).cpu())
    if has_label:
        cur_acc = accuracy_score(truths, predictions)
        f1 = f1_score(truths, predictions, average='macro')
        print(f'number of skipped samples: {nskip}')
        print(f'    ==  Accuracy: {cur_acc:.6f} | F1 Score: {f1:.6f}  ==')
        print(confusion_matrix(truths, predictions))
        print(classification_report(truths, predictions, labels=label_list, target_names=label_dict_list, digits=4))

        logging.info(f'number of skipped samples: {nskip}')
        logging.info(f'    ==  Accuracy: {cur_acc:.6f} | F1 Score: {f1:.6f}  ==')
        logging.info(confusion_matrix(truths, predictions))
        logging.info(classification_report(truths, predictions, labels=label_list, target_names=label_dict_list, digits=4))

logging.info("Finished Evaluation")

# Create a DataFrame with cell names and predictions
df = pd.DataFrame({
    'Cell_Name': cells,
    'Prediction': pred_ctypes
})

# Save to CSV
df.to_csv(args.prediction_file, index=False, header=False)

# write counts of each cell type
counts = df['Prediction'].value_counts()
counts.to_csv(args.stats_file)
