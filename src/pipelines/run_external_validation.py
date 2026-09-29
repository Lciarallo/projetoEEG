"""Reanálise local orientada a relaxamento e duas bases externas independentes."""
import argparse
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.metrics import balanced_accuracy_score, roc_auc_score
from sklearn.model_selection import LeaveOneGroupOut
from src.relaxation import METRICS, SEED, frontal_features, paired_summary, local_block_difference
from src.preprocessing import load_raw_openbci
from src.study import SESSIONS, session_path
from src.validation import holm_adjust

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'results'/'validation_2026_09_28'
DATA=ROOT/'data'/'external'


def paired_contrasts(frame,condition_col,positive,negative,metrics):
    """Nunca usa épocas como réplicas de participantes."""
    rows=[]
    for metric in metrics:
        wide=frame.pivot(index='subject',columns=condition_col,values=metric)
        pairs=wide[[positive,negative]].dropna()
        rows.append(dict(metric=metric,**paired_summary(pairs[positive]-pairs[negative])))
    result=pd.DataFrame(rows)
    result['p_holm']=holm_adjust(result.p)
    return result


def eegmat():
    import mne
    rows,quality=[],[]
    for subject in range(36):
        for code,condition in [(1,'rest'),(2,'arithmetic')]:
            path=DATA/'eegmat'/f'Subject{subject:02d}_{code}.edf'
            recording=mne.io.read_raw_edf(path,verbose=False)
            fs=recording.info['sfreq']
            samples=int(60*fs)
            if recording.n_times < samples:
                raise ValueError(f'Registro curto: {path}')
            start=recording.n_times-samples if code==1 else 0
            raw=recording.get_data(picks=['EEG F3','EEG F4'],start=start,stop=start+samples).T*1e6
            frame,q=frontal_features(raw,fs)
            row=dict(subject=subject,condition=condition,duration_s=recording.n_times/fs,
                     window_start_s=start/fs,fs=fs,**q)
            quality.append(row)
            if q['accepted']>=14:
                rows.append(dict(subject=subject,condition=condition,**frame[METRICS].median().to_dict()))
        print(f'EEGMAT {subject+1}/36',flush=True)
    pd.DataFrame(quality).to_csv(OUT/'eegmat_quality.csv',index=False)
    frame=pd.DataFrame(rows)
    counts=frame.groupby('subject').size()
    frame=frame[frame.subject.isin(counts[counts==2].index)].sort_values(['subject','condition'])
    frame.to_csv(OUT/'eegmat_subject_features.csv',index=False)
    # Sinal positivo = repouso menos cálculo, não "relaxamento confirmado".
    paired_contrasts(frame,'condition','rest','arithmetic',METRICS).to_csv(OUT/'eegmat_paired.csv',index=False)
    return frame


def local_study():
    all_frames=[]
    for label,stamp in SESSIONS:
        raw=load_raw_openbci(session_path(stamp))[:,:2]
        windows=[(f'BB {i+1}',i*60000,(i+1)*60000) for i in range(4)] if label=='BB' else [(label,0,60000)]
        for condition,start,stop in windows:
            for threshold in (150.,500.):
                for notch in (False,True):
                    frame,q=frontal_features(raw[start:stop],200,threshold,notch)
                    frame['condition']=condition
                    frame['threshold_uv']=threshold
                    frame['notch20']=notch
                    frame['candidate']=q['candidate']
                    all_frames.append(frame)
    epochs=pd.concat(all_frames,ignore_index=True)
    epochs.to_csv(OUT/'local_epochs.csv',index=False)
    summary=epochs.groupby(['condition','threshold_uv','notch20'],sort=False)[METRICS].median().reset_index()
    counts=epochs.groupby(['condition','threshold_uv','notch20'],sort=False).size().to_numpy()
    summary['accepted']=counts
    summary.to_csv(OUT/'local_features.csv',index=False)
    contrasts=[]
    for (threshold,notch),group in epochs.groupby(['threshold_uv','notch20']):
        base=group[group.condition=='Base']
        for condition in ['BB 1','BB 2','BB 3','BB 4','Pós 1','Pós 2']:
            for metric in METRICS:
                for block_s in (30,60):
                    contrasts.append(dict(condition=condition,metric=metric,threshold_uv=threshold,notch20=notch,
                        block_s=block_s,**local_block_difference(base,group[group.condition==condition],metric,block_s)))
    pd.DataFrame(contrasts).to_csv(OUT/'local_block_contrasts.csv',index=False)
    return summary


def fit_model():
    return make_pipeline(StandardScaler(),LogisticRegression(C=1.,solver='liblinear',random_state=SEED))


def grouped_predictions(x,y,groups):
    predictions=np.empty(len(y))
    for train,test in LeaveOneGroupOut().split(x,y,groups):
        if set(groups[train]) & set(groups[test]):
            raise AssertionError('Vazamento de participantes')
        model=fit_model().fit(x[train],y[train])
        predictions[test]=model.predict_proba(x[test])[:,1]
    return predictions


def decoding(external,local,n_permutations):
    x=external[METRICS].to_numpy()
    y=(external.condition=='rest').to_numpy(int)
    groups=external.subject.to_numpy()
    pred=grouped_predictions(x,y,groups)
    output=external[['subject','condition']].copy()
    output['prob_rest']=pred
    output.to_csv(OUT/'eegmat_out_of_subject_predictions.csv',index=False)
    score=balanced_accuracy_score(y,pred>=.5)
    auc=roc_auc_score(y,pred)
    rng=np.random.default_rng(SEED)
    null=[]
    for iteration in range(n_permutations):
        shuffled=y.copy()
        for subject in np.unique(groups):
            if rng.integers(2):
                mask=groups==subject
                shuffled[mask]=1-shuffled[mask]
        p=grouped_predictions(x,shuffled,groups)
        null.append(balanced_accuracy_score(shuffled,p>=.5))
        if (iteration+1)%25==0:
            print(f'Permutações por participante: {iteration+1}/{n_permutations}',flush=True)
    participant_acc=np.array([np.mean((pred[groups==s]>=.5)==y[groups==s]) for s in np.unique(groups)])
    boot=rng.choice(participant_acc,(5000,len(participant_acc))).mean(axis=1)
    result=dict(n_subjects=len(np.unique(groups)),balanced_accuracy=score,auc=auc,
        ci_low=float(np.quantile(boot,.025)),ci_high=float(np.quantile(boot,.975)),
        n_permutations=n_permutations,permutation_p=float((1+np.sum(np.array(null)>=score))/(1+n_permutations)),
        null_mean=float(np.mean(null)),features=METRICS)
    (OUT/'eegmat_model.json').write_text(json.dumps(result,indent=2)+'\n')
    pd.DataFrame(dict(balanced_accuracy=null)).to_csv(OUT/'eegmat_permutation_scores.csv',index=False)
    model=fit_model().fit(x,y)
    transfer=local.copy()
    transfer['prob_external_rest']=model.predict_proba(transfer[METRICS].to_numpy())[:,1]
    lower,upper=np.quantile(x,[.025,.975],axis=0)
    transfer['features_outside_external95']=np.sum((transfer[METRICS].to_numpy()<lower)|(transfer[METRICS].to_numpy()>upper),axis=1)
    transfer.to_csv(OUT/'local_external_model_transfer.csv',index=False)
    print('Modelo externo:',result,flush=True)


def sudre():
    questionnaire,eeg_results,correlations,coverage,paired_rows=[],[],[],[],[]
    for group in ('G6','G40'):
        relax=pd.read_csv(DATA/'sudre2024'/f'relax_{group}.csv',sep=';')
        relax.columns=relax.columns.str.strip()
        if relax.duplicated(['subject','Condition']).any():
            raise ValueError('Questionário duplicado')
        for metric in ('PA','CT','PT'):
            wide=relax.pivot(index='subject',columns='Condition',values=metric)
            pairs=wide[['C0','C1']].dropna()
            d=pairs.C1-pairs.C0
            questionnaire.append(dict(group=group,metric=metric,**paired_summary(d)))
            for subject,val in d.items():
                paired_rows.append(dict(group=group,subject=subject,metric=metric,difference=val))
        eeg=pd.read_csv(DATA/'sudre2024'/f'eeg_{group}.csv',sep=';')
        duplicates=int(eeg.duplicated().sum())
        eeg=eeg.drop_duplicates()
        if eeg.duplicated(['subject','conditions','epoch','electrode']).any():
            raise ValueError('Duplicatas conflitantes de EEG')
        roi=eeg[eeg.electrode.isin(['F3','Fz','F4']) & eeg.conditions.isin(['C0','C1'])].copy()
        roi['phase']=roi.epoch.str.extract(r'^(Bruit_1|Stim1|Stim2|Stim3)(?:_|$)',expand=False)
        roi=roi.dropna(subset=['phase'])
        # Peso igual por canal e por bloco, não pelo número de épocas retidas.
        channels=roi.groupby(['subject','conditions','phase','electrode']).relativeAlpha.mean().reset_index()
        phases=channels.groupby(['subject','conditions','phase']).relativeAlpha.agg(['mean','count']).reset_index()
        phases=phases[phases['count']==3]
        wide=phases.pivot(index=['subject','conditions'],columns='phase',values='mean')
        complete=wide.reindex(columns=['Bruit_1','Stim1','Stim2','Stim3']).dropna()
        delta=complete[['Stim1','Stim2','Stim3']].mean(axis=1)-complete.Bruit_1
        paired=delta.unstack('conditions').reindex(columns=['C0','C1']).dropna()
        difference=paired.C1-paired.C0
        eeg_results.append(dict(group=group,metric='frontal_alpha_delta_C1_minus_C0',**paired_summary(difference)))
        for subject,val in difference.items():
            paired_rows.append(dict(group=group,subject=subject,metric='frontal_alpha_delta',difference=val))
        ct=relax.pivot(index='subject',columns='Condition',values='CT')
        joint=pd.concat([difference.rename('eeg'),(ct.C1-ct.C0).rename('ct')],axis=1).dropna()
        rho,p=stats.spearmanr(joint.eeg,joint.ct)
        correlations.append(dict(group=group,n=len(joint),rho=rho,p=p))
        coverage.append(dict(group=group,questionnaire_subjects=int(relax.subject.nunique()),
            eeg_subjects=int(eeg.subject.nunique()),frontal_complete_pairs=len(paired),
            eeg_exact_duplicates_removed=duplicates,subject_ids_complete=','.join(map(str,paired.index)),
            eeg_ids_not_in_questionnaire=','.join(map(str,set(eeg.subject)-set(relax.subject))),
            questionnaire_ids_no_complete_roi=','.join(map(str,sorted(set(relax.subject)-set(paired.index))))))
    q=pd.DataFrame(questionnaire);q['p_holm']=holm_adjust(q.p)
    e=pd.DataFrame(eeg_results);e['p_holm']=holm_adjust(e.p)
    c=pd.DataFrame(correlations);c['p_holm']=holm_adjust(c.p)
    for name,frame in [('sudre_questionnaire',q),('sudre_frontal_eeg',e),('sudre_eeg_questionnaire_correlation',c),
                       ('sudre_coverage',pd.DataFrame(coverage)),('sudre_paired_differences',pd.DataFrame(paired_rows))]:
        frame.to_csv(OUT/f'{name}.csv',index=False)
    print(q.to_string(index=False),flush=True)
    print(e.to_string(index=False),flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--permutations',type=int,default=199)
    args=parser.parse_args()
    if args.permutations<1:
        parser.error('Número de permutações deve ser positivo')
    OUT.mkdir(parents=True,exist_ok=True)
    sudre()
    external=eegmat()
    local=local_study()
    decoding(external,local,args.permutations)
    paths=list((ROOT/'src').rglob('*.py'))+list((ROOT/'data'/'raw').glob('*.txt'))
    manifest=dict(generated_utc=datetime.now(timezone.utc).isoformat(),seed=SEED,
        versions={name:importlib.metadata.version(name) for name in ['numpy','scipy','pandas','mne','scikit-learn']},
        code_and_local_data_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        external_manifest_sha256=hashlib.sha256((OUT/'download_manifest.json').read_bytes()).hexdigest(),
        parameters=dict(metrics=METRICS,fs=200,band=[1,30],epoch_s=2,edge_s=2,primary_threshold_uv=150,
                        primary_notch20=False,minimum_external_epochs=14,permutations=args.permutations,
                        eegmat_window_s=60,local_window_s=300,local_bootstrap_block_s=[30,60]))
    (OUT/'analysis_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')


if __name__=='__main__':
    main()
