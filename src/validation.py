"""Partições por tempo de aquisição, preservando lacunas de épocas rejeitadas."""
import numpy as np


def temporal_folds(indices_by_session, block_epochs=30, gap_epochs=60):
    """Pareia blocos cobrindo toda a gravação; gap medido entre bordas das épocas.

    Folds compartilham treino e não constituem réplicas independentes.
    """
    if block_epochs < 1 or gap_epochs < 0 or len(indices_by_session) != 2:
        raise ValueError('Esperadas duas sessões e durações válidas')
    arrays = [np.asarray(ids, dtype=int) for ids in indices_by_session]
    if any(len(ids) == 0 or np.any(np.diff(ids) <= 0) for ids in arrays):
        raise ValueError('Índices devem ser não vazios e crescentes')
    blocks = [np.unique(ids // block_epochs) for ids in arrays]
    count = min(map(len, blocks))
    chosen = [b[np.linspace(0, len(b)-1, count).astype(int)] for b in blocks]
    for fold in range(count):
        train, test, windows = [], [], []
        for ids, selected in zip(arrays, chosen):
            start = int(selected[fold]*block_epochs)
            stop = start+block_epochs
            test.append(np.flatnonzero((ids >= start) & (ids < stop)))
            train.append(np.flatnonzero((ids+1 <= start-gap_epochs) | (ids >= stop+gap_epochs)))
            windows.append((start, stop))
        if all(len(x) >= 2 for x in train+test):
            yield train, test, windows


def holm_adjust(pvalues):
    pvalues = np.asarray(pvalues, dtype=float)
    if pvalues.ndim != 1 or not np.isfinite(pvalues).all() or np.any((pvalues < 0) | (pvalues > 1)):
        raise ValueError('p-valores inválidos')
    order = np.argsort(pvalues)
    adjusted = np.empty_like(pvalues)
    adjusted[order] = np.minimum(1, np.maximum.accumulate(pvalues[order]*np.arange(len(order), 0, -1)))
    return adjusted
