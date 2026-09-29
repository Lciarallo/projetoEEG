"""Dados e métricas comuns aos classificadores exploratórios."""
import numpy as np
from sklearn.metrics import accuracy_score, balanced_accuracy_score, roc_auc_score
from src.preprocessing import load_raw_openbci, apply_filters, segment_epochs
from src.study import SESSIONS, session_path


def load_sessions():
    sessions = []
    for label, stamp in SESSIONS:
        clean = apply_filters(load_raw_openbci(session_path(stamp)))
        epochs, ids = segment_epochs(clean, threshold_frontal_p2p=500)
        sessions.append((label, epochs, ids))
    return sessions


def normalized_tensor(epochs):
    scaled = (epochs-epochs.mean(axis=1,keepdims=True))/(epochs.std(axis=1,keepdims=True)+1e-8)
    return scaled.transpose(0,2,1).astype(np.float32)


def metrics(y, probs):
    return dict(accuracy=accuracy_score(y, probs >= .5),
                balanced_accuracy=balanced_accuracy_score(y, probs >= .5),
                auc=roc_auc_score(y, probs))


def fold_arrays(arrays, train, test):
    xtr = np.concatenate([a[idx] for a, idx in zip(arrays, train)])
    xte = np.concatenate([a[idx] for a, idx in zip(arrays, test)])
    ytr = np.concatenate([np.full(len(idx), label, dtype=np.int64) for label, idx in enumerate(train)])
    yte = np.concatenate([np.full(len(idx), label, dtype=np.int64) for label, idx in enumerate(test)])
    return xtr, ytr, xte, yte
