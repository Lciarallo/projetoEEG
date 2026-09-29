"""Random Forest com gap real e transferência descritiva para sessões posteriores."""
import sys
from pathlib import Path
if __package__ in (None, ''):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from src.features import extract_neuro_features
from src.decoding import load_sessions, fold_arrays, metrics
from src.validation import temporal_folds
from src.study import RESULTS, provenance


def classifier(seed=42):
    return RandomForestClassifier(n_estimators=150, max_depth=6, class_weight='balanced', random_state=seed, n_jobs=4)


def main():
    sessions = load_sessions()
    features = [extract_neuro_features(ep)[0] for _,ep,_ in sessions]
    names = extract_neuro_features(sessions[0][1])[1]
    rows = []
    for fold, (train,test,windows) in enumerate(temporal_folds([s[2] for s in sessions[:2]]),1):
        xtr,ytr,xte,yte = fold_arrays(features[:2],train,test)
        model = classifier(42+fold).fit(xtr,ytr)
        row = dict(fold=fold, gap_s=120, n_train=len(ytr), n_test=len(yte),
                   base_test_start_s=windows[0][0]*2, bb_test_start_s=windows[1][0]*2,
                   **metrics(yte,model.predict_proba(xte)[:,1]),
                   majority_accuracy=max(np.mean(yte==0),np.mean(yte==1)))
        null = classifier(42+fold).fit(xtr,np.random.default_rng(42+fold).permutation(ytr))
        row.update({'null_'+k:v for k,v in metrics(yte,null.predict_proba(xte)[:,1]).items()})
        rows.append(row)
    if not rows:
        raise ValueError('Nenhum fold válido')
    RESULTS.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(RESULTS/'random_forest_temporal.csv',index=False)
    x = np.concatenate(features[:2])
    y = np.concatenate([np.full(len(a),i) for i,a in enumerate(features[:2])])
    model = classifier().fit(x,y)
    transfer = []
    for (label,_,ids),feats in zip(sessions,features):
        probs = model.predict_proba(feats)[:,1]
        transfer.append(dict(sessao=label,epocas=len(feats),prob_classe_bb=float(probs.mean()),
                             fracao_classe_bb=float(np.mean(probs>=.5)),
                             avaliacao='in_sample' if label in ('Base','BB') else 'out_of_session_sem_rotulo_validado'))
    pd.DataFrame(transfer).to_csv(RESULTS/'transferencia_rf.csv',index=False)
    pd.DataFrame(dict(feature=names,importance=model.feature_importances_)).sort_values('importance',ascending=False).to_csv(RESULTS/'importancias_rf.csv',index=False)
    provenance('random_forest',dict(gap_s=120,block_s=60,threshold_uv=500,features=24,
               seeds='42+fold; full fit 42',class_weight='balanced',n_estimators=150,max_depth=6,
               limitation='session confounded; no causal or population inference; A/A historical only'))
    print(pd.DataFrame(rows).to_string(index=False))
    print(pd.DataFrame(transfer).to_string(index=False))


if __name__ == '__main__':
    main()
