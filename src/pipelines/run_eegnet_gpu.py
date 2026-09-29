"""EEGNet exploratória: blocos de aquisição e exclusão de 120 s entre bordas."""
import sys
from pathlib import Path
if __package__ in (None, ''):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import argparse
import numpy as np
import pandas as pd
import torch
from torch.utils.data import DataLoader, TensorDataset
from src.models.eegnet import EEGNet
from src.decoding import load_sessions, normalized_tensor, metrics, fold_arrays
from src.validation import temporal_folds
from src.study import RESULTS, provenance


def train_and_predict(xtr, ytr, xte, seed=42, epochs=45, device='cpu'):
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
    model = EEGNet(samples=xtr.shape[-1]).to(device)
    counts = np.bincount(ytr, minlength=2)
    if np.any(counts == 0):
        raise ValueError('Treino requer duas classes')
    weights = torch.tensor(len(ytr)/(2*counts), dtype=torch.float32, device=device)
    loss_fn = torch.nn.CrossEntropyLoss(weight=weights)
    optimizer = torch.optim.AdamW(model.parameters(), lr=.001, weight_decay=.01)
    data = TensorDataset(torch.from_numpy(xtr).unsqueeze(1), torch.from_numpy(ytr))
    loader = DataLoader(data, batch_size=32, shuffle=True, generator=torch.Generator().manual_seed(seed))
    model.train()
    for _ in range(epochs):
        for x, y in loader:
            optimizer.zero_grad()
            loss = loss_fn(model(x.to(device)), y.to(device))
            loss.backward()
            optimizer.step()
    model.eval()
    with torch.no_grad():
        return np.concatenate([model(torch.from_numpy(batch).unsqueeze(1).to(device)).softmax(1)[:,1].cpu().numpy()
                               for batch in np.array_split(xte, max(1, (len(xte)+63)//64))])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--epochs', type=int, default=45)
    parser.add_argument('--device', choices=['auto','cpu','cuda'], default='auto')
    args = parser.parse_args()
    if args.epochs < 1:
        parser.error('--epochs deve ser positivo')
    device = ('cuda' if torch.cuda.is_available() else 'cpu') if args.device == 'auto' else args.device
    torch.set_num_threads(4)
    sessions = load_sessions()[:2]
    arrays = [normalized_tensor(ep) for _, ep, _ in sessions]
    ids = [i for _, _, i in sessions]
    rows = []
    for fold, (train, test, windows) in enumerate(temporal_folds(ids), 1):
        xtr, ytr, xte, yte = fold_arrays(arrays, train, test)
        seed = 2026+fold*47
        pred = train_and_predict(xtr,ytr,xte,seed,args.epochs,device)
        # Controle diagnóstico único, não uma distribuição nula inferencial.
        permuted = np.random.default_rng(seed).permutation(ytr)
        null = train_and_predict(xtr,permuted,xte,seed,args.epochs,device)
        row = dict(fold=fold, n_train=len(ytr), n_test=len(yte),
                   base_test_start_s=windows[0][0]*2, bb_test_start_s=windows[1][0]*2,
                   gap_s=120, epochs=args.epochs, **metrics(yte,pred),
                   **{'null_'+k:v for k,v in metrics(yte,null).items()},
                   majority_accuracy=max(np.mean(yte==0),np.mean(yte==1)))
        rows.append(row)
        print(row, flush=True)
    if not rows:
        raise ValueError('Nenhum fold válido com gap especificado')
    RESULTS.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(RESULTS/'eegnet_temporal.csv', index=False)
    provenance('eegnet', dict(epochs=args.epochs, device=device, seed_rule='2026+fold*47',
               gap_s=120, block_s=60, threshold_uv=500, class_weights='inverse frequency',
               null='one training-label shuffle per fold; no inference',
               limitation='session confounded; whole-session zero-phase filtering'))
    print(pd.DataFrame(rows)[['accuracy','balanced_accuracy','null_balanced_accuracy']].mean())


if __name__ == '__main__':
    main()
