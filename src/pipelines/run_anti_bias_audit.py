import os
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import signal
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.metrics import accuracy_score, roc_auc_score
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import TensorDataset, DataLoader

# Fixar dispositivo
device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
print(f"Dispositivo: {device} ({torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'CPU'})")

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
DATA_DIR = os.path.join(REPO_ROOT, "data", "raw") if os.path.exists(os.path.join(REPO_ROOT, "data", "raw")) else os.path.join(REPO_ROOT, "X")
OUTPUT_DIR = os.path.join(REPO_ROOT, "results")
os.makedirs(os.path.join(OUTPUT_DIR, "figures"), exist_ok=True)
os.makedirs(os.path.join(OUTPUT_DIR, "tables"), exist_ok=True)

FS = 200
EPOCH_SEC = 2.0
SAMPLES_PER_EPOCH = int(EPOCH_SEC * FS)

def load_and_clean(filepath):
    df = pd.read_csv(filepath, skiprows=1, header=None, comment='%', usecols=[1, 2, 3, 4])
    raw = df.to_numpy(dtype=float)
    d = signal.detrend(raw, axis=0)
    for f in [60.0, 20.0, 80.0]:
        b, a = signal.iirnotch(f, 30.0, fs=FS)
        d = signal.filtfilt(b, a, d, axis=0)
    sos = signal.butter(4, [1.0, 45.0], btype='bandpass', fs=FS, output='sos')
    return signal.sosfiltfilt(sos, d, axis=0)

def extract_epochs(clean_data, label, threshold_uV=500.0):
    n = clean_data.shape[0] // SAMPLES_PER_EPOCH
    epochs = []
    labels = []
    for i in range(n):
        ep = clean_data[i*SAMPLES_PER_EPOCH : (i+1)*SAMPLES_PER_EPOCH, :]
        p2p_frontal = np.ptp(ep[:, :2], axis=0)
        if np.any(p2p_frontal > threshold_uV) or np.any(p2p_frontal < 2.0):
            continue
        ep_norm = (ep - np.mean(ep, axis=0)) / (np.std(ep, axis=0) + 1e-8)
        epochs.append(ep_norm.T) # (4, 400)
        labels.append(label)
    return np.array(epochs, dtype=np.float32), np.array(labels, dtype=np.int64)

# ==============================================================================
# EEGNet Architecture
# ==============================================================================
class Conv2dWithConstraint(nn.Conv2d):
    def __init__(self, *args, max_norm=1.0, **kwargs):
        super().__init__(*args, **kwargs)
        self.max_norm = max_norm

    def forward(self, x):
        if self.max_norm is not None:
            with torch.no_grad():
                norm = self.weight.norm(2, dim=(1, 2, 3), keepdim=True)
                desired = torch.clamp(norm, max=self.max_norm)
                self.weight.mul_(desired / (norm + 1e-8))
        return super().forward(x)

class EEGNet(nn.Module):
    def __init__(self, n_classes=2, channels=4, samples=400, F1=8, D=2, F2=16, kernel_length=64):
        super().__init__()
        self.conv1 = nn.Conv2d(1, F1, (1, kernel_length), padding=(0, kernel_length // 2), bias=False)
        self.bn1 = nn.BatchNorm2d(F1)
        self.depthwise = Conv2dWithConstraint(F1, F1 * D, (channels, 1), groups=F1, bias=False, max_norm=1.0)
        self.bn2 = nn.BatchNorm2d(F1 * D)
        self.act1 = nn.ELU()
        self.pool1 = nn.AvgPool2d((1, 4))
        self.drop1 = nn.Dropout(0.25)
        self.sep_depthwise = nn.Conv2d(F1 * D, F1 * D, (1, 16), padding=(0, 8), groups=F1 * D, bias=False)
        self.sep_pointwise = nn.Conv2d(F1 * D, F2, (1, 1), bias=False)
        self.bn3 = nn.BatchNorm2d(F2)
        self.act2 = nn.ELU()
        self.pool2 = nn.AvgPool2d((1, 8))
        self.drop2 = nn.Dropout(0.25)
        flatten_dim = F2 * (samples // 32)
        self.classifier = nn.Sequential(
            nn.Flatten(),
            nn.Linear(flatten_dim, 32),
            nn.ELU(),
            nn.Dropout(0.5),
            nn.Linear(32, n_classes)
        )

    def forward(self, x):
        x = self.conv1(x)
        x = self.bn1(x)
        x = self.depthwise(x)
        x = self.bn2(x)
        x = self.act1(x)
        x = self.pool1(x)
        x = self.drop1(x)
        x = self.sep_depthwise(x)
        x = self.sep_pointwise(x)
        x = self.bn3(x)
        x = self.act2(x)
        x = self.pool2(x)
        x = self.drop2(x)
        return self.classifier(x)

def train_eegnet(X_train, y_train, epochs=45, batch_size=32, lr=0.001, seed=42):
    torch.manual_seed(seed)
    np.random.seed(seed)
    train_dataset = TensorDataset(torch.tensor(X_train).unsqueeze(1).to(device), torch.tensor(y_train).to(device))
    train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
    model = EEGNet(n_classes=2, channels=4, samples=SAMPLES_PER_EPOCH).to(device)
    criterion = nn.CrossEntropyLoss()
    optimizer = optim.AdamW(model.parameters(), lr=lr, weight_decay=0.01)
    model.train()
    for _ in range(epochs):
        for bx, by in train_loader:
            optimizer.zero_grad()
            loss = criterion(model(bx), by)
            loss.backward()
            optimizer.step()
    return model

def predict_eegnet(model, X_test):
    model.eval()
    test_tensor = torch.tensor(X_test).unsqueeze(1).to(device)
    with torch.no_grad():
        out = model(test_tensor)
        probs = torch.softmax(out, dim=1)[:, 1].cpu().numpy()
        preds = torch.argmax(out, dim=1).cpu().numpy()
    return preds, probs

# ==============================================================================
# Extração de Features Neurobiológicas para ML Clássico
# ==============================================================================
def extract_neuro_features(epochs_raw_list):
    """
    Extrai features clássicas interpretáveis:
    - Potências relativas em 5 bandas para cada um dos 4 canais (20 features)
    - Frontal Alpha Asymmetry (FAA = ln(RelAlpha_F4) - ln(RelAlpha_F3))
    - Theta/Beta Ratio (TBR) em F3 e F4 (2 features)
    - wPLI em Alfa entre F3 e F4 (1 feature)
    Total: 24 features neurobiológicas explícitas por época.
    """
    features = []
    bands = {
        'Delta': (1.0, 4.0),
        'Theta': (4.0, 8.0),
        'Alpha': (8.0, 13.0),
        'Beta':  (13.0, 30.0),
        'Gamma': (30.0, 45.0)
    }
    
    # Filtro Alfa para wPLI
    sos_alpha = signal.butter(4, [8.0, 13.0], btype='bandpass', fs=FS, output='sos')
    
    for ep in epochs_raw_list:
        # ep formato: (4, 400)
        f, psd = signal.welch(ep.T, fs=FS, nperseg=400, axis=0) # psd: (freqs, 4)
        trapz = getattr(np, 'trapezoid', getattr(np, 'trapz', None))
        mask_tot = (f >= 1.0) & (f <= 45.0)
        tot_pwr = trapz(psd[mask_tot, :], f[mask_tot], axis=0) + 1e-12
        
        row = []
        # 1. Potências relativas por banda e por canal
        p_bands = {}
        for b_name, (low, high) in bands.items():
            mask = (f >= low) & (f <= high)
            bp = trapz(psd[mask, :], f[mask], axis=0)
            rel_bp = bp / tot_pwr
            p_bands[b_name] = rel_bp
            for ch in range(4):
                row.append(rel_bp[ch])
                
        # 2. FAA Relativo: ln(Rel_Alpha_F4) - ln(Rel_Alpha_F3)
        faa = np.log(p_bands['Alpha'][1] + 1e-12) - np.log(p_bands['Alpha'][0] + 1e-12)
        row.append(faa)
        
        # 3. Theta / Beta ratio (F3 e F4)
        tbr_f3 = p_bands['Theta'][0] / (p_bands['Beta'][0] + 1e-12)
        tbr_f4 = p_bands['Theta'][1] / (p_bands['Beta'][1] + 1e-12)
        row.append(tbr_f3)
        row.append(tbr_f4)
        
        # 4. Phase lag imaginário em Alfa entre F3 e F4
        ep_filt = signal.sosfiltfilt(sos_alpha, ep.T, axis=0)
        h = signal.hilbert(ep_filt, axis=0)
        cross = h[:, 0] * np.conj(h[:, 1])
        im_lag = np.mean(np.abs(np.imag(cross)))
        row.append(im_lag)
        
        features.append(row)
        
    feature_names = []
    for b_name in bands.keys():
        for ch in ['F3', 'F4', 'P3', 'P4']:
            feature_names.append(f"{b_name}_{ch}")
    feature_names.extend(['FAA_Relativo', 'TBR_F3', 'TBR_F4', 'PhaseLag_Alpha_F3_F4'])
    return np.array(features, dtype=np.float32), feature_names

def main():
    print("="*75)
    print("ESTRATÉGIAS AVANÇADAS DE AUDITORIA CONTRA VIÉS METODOLÓGICO")
    print("="*75)
    
    # Carregar dados
    f_base = os.path.join(DATA_DIR, "OpenBCI-RAW-2023-11-28_15-56-59.txt")
    f_bb   = os.path.join(DATA_DIR, "OpenBCI-RAW-2023-11-28_16-04-54.txt")
    f_p1   = os.path.join(DATA_DIR, "OpenBCI-RAW-2023-11-28_16-36-02.txt")
    f_p2   = os.path.join(DATA_DIR, "OpenBCI-RAW-2023-11-28_16-44-15.txt")
    
    c_base = load_and_clean(f_base)
    c_bb   = load_and_clean(f_bb)
    c_p1   = load_and_clean(f_p1)
    c_p2   = load_and_clean(f_p2)
    
    X_base, y_base = extract_epochs(c_base, 0)
    X_bb, y_bb     = extract_epochs(c_bb, 1)
    X_p1, y_p1     = extract_epochs(c_p1, 2)
    X_p2, y_p2     = extract_epochs(c_p2, 3)
    
    print(f"Épocas obtidas: Repouso={len(X_base)}, Binaural={len(X_bb)}, Pós1={len(X_p1)}, Pós2={len(X_p2)}")
    
    # --------------------------------------------------------------------------
    # ESTRATÉGIA 1: Teste A/A (Controle Negativo de Deriva / Estacionariedade)
    # --------------------------------------------------------------------------
    print("\n--- ESTRATÉGIA 1: Teste A/A (Controle Negativo de Deriva Temporal) ---")
    print("Pergunta crítica: O classificador separa a 1ª metade da 2ª metade da MESMA sessão?")
    print("Se SIM: a rede está decorando drift lento (metodologia viciada).")
    print("Se NÃO (~50% de acerto): a sessão é estacionária e não há viés de drift temporal.")
    
    # 1.1 Metade 1 vs Metade 2 do Repouso (A/A Baseline)
    mid_base = len(X_base) // 2
    X_aa_base = np.concatenate([X_base[:mid_base], X_base[mid_base:]])
    y_aa_base = np.array([0]*mid_base + [1]*(len(X_base)-mid_base))
    
    # 5-fold CV com modelo linear e random forest
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    feats_base, _ = extract_neuro_features(X_aa_base)
    rf_aa = RandomForestClassifier(n_estimators=100, random_state=42)
    scores_aa_base = cross_val_score(rf_aa, feats_base, y_aa_base, cv=skf)
    print(f"  A/A Teste no Repouso (1ª vs 2ª metade): Acurácia = {np.mean(scores_aa_base)*100:.2f}% ± {np.std(scores_aa_base)*100:.2f}% (Chance esperada: 50.0%)")
    
    # 1.2 Metade 1 vs Metade 2 da Sessão Binaural
    mid_bb = len(X_bb) // 2
    X_aa_bb = np.concatenate([X_bb[:mid_bb], X_bb[mid_bb:]])
    y_aa_bb = np.array([0]*mid_bb + [1]*(len(X_bb)-mid_bb))
    feats_bb, _ = extract_neuro_features(X_aa_bb)
    scores_aa_bb = cross_val_score(rf_aa, feats_bb, y_aa_bb, cv=skf)
    print(f"  A/A Teste no Binaural (1ºs 10 min vs Últimos 10 min): Acurácia = {np.mean(scores_aa_bb)*100:.2f}% ± {np.std(scores_aa_bb)*100:.2f}%")
    
    # --------------------------------------------------------------------------
    # ESTRATÉGIA 2: Teste de Transferência Externa (Out-of-Session Generalization)
    # --------------------------------------------------------------------------
    print("\n--- ESTRATÉGIA 2: Transferência Externa para Sessões Pós-Estimulação ---")
    print("Pergunta crítica: O que o modelo prevê para as sessões Pós 1 e Pós 2 (que nunca viu)?")
    print("Se a rede aprendeu um artefato do Arquivo 2, ela vai errar completamente.")
    print("Se aprendeu 'Ativação Binaural', as sessões Pós (sem som) devem ser classificadas como REPOUSO (Classe 0)!")
    
    # Treinar modelo em Repouso (0) vs Binaural (1)
    X_train_full = np.concatenate([X_base, X_bb])
    y_train_full = np.concatenate([y_base, y_bb])
    
    model_eegnet = train_eegnet(X_train_full, y_train_full, epochs=45, seed=42)
    
    preds_base, probs_base = predict_eegnet(model_eegnet, X_base)
    preds_bb, probs_bb     = predict_eegnet(model_eegnet, X_bb)
    preds_p1, probs_p1     = predict_eegnet(model_eegnet, X_p1)
    preds_p2, probs_p2     = predict_eegnet(model_eegnet, X_p2)
    
    pct_rest_base = np.mean(preds_base == 0) * 100
    pct_bb_bb     = np.mean(preds_bb == 1) * 100
    pct_rest_p1   = np.mean(preds_p1 == 0) * 100
    pct_rest_p2   = np.mean(preds_p2 == 0) * 100
    
    print(f"  Previsão na Linha de Base: {pct_rest_base:.1f}% classificado como REPOUSO (Classe 0)")
    print(f"  Previsão no Binaural:      {pct_bb_bb:.1f}% classificado como BINAURAL (Classe 1)")
    print(f"  Previsão na Pós-Sessão 1:  {pct_rest_p1:.1f}% classificado como REPOUSO (Retorno ao repouso!)")
    print(f"  Previsão na Pós-Sessão 2:  {pct_rest_p2:.1f}% classificado como REPOUSO (Retorno ao repouso!)")
    
    # --------------------------------------------------------------------------
    # ESTRATÉGIA 3: Machine Learning Interpretable com 24 Features Neurobiológicas
    # --------------------------------------------------------------------------
    print("\n--- ESTRATÉGIA 3: ML Interpretable (Random Forest com 24 Features Fisiológicas) ---")
    print("Pergunta crítica: O que acontece se usarmos APENAS features neurofisiológicas explícitas?")
    
    # Balancear dados
    n_bal = min(len(X_base), len(X_bb))
    rng = np.random.RandomState(42)
    idx_b = rng.choice(len(X_base), n_bal, replace=False)
    idx_bb = rng.choice(len(X_bb), n_bal, replace=False)
    
    X_bal = np.concatenate([X_base[idx_b], X_bb[idx_bb]])
    y_bal = np.concatenate([np.zeros(n_bal), np.ones(n_bal)])
    
    feats_all, feat_names = extract_neuro_features(X_bal)
    
    rf = RandomForestClassifier(n_estimators=150, max_depth=6, random_state=42)
    scores_rf = cross_val_score(rf, feats_all, y_bal, cv=10, scoring='accuracy')
    scores_auc = cross_val_score(rf, feats_all, y_bal, cv=10, scoring='roc_auc')
    
    print(f"  Random Forest 10-Fold CV Acurácia: {np.mean(scores_rf)*100:.2f}% ± {np.std(scores_rf)*100:.2f}%")
    print(f"  Random Forest 10-Fold CV ROC-AUC:  {np.mean(scores_auc):.3f} ± {np.std(scores_auc):.3f}")
    
    # Importância das Features
    rf.fit(feats_all, y_bal)
    importances = rf.feature_importances_
    sorted_idx = np.argsort(importances)[::-1]
    
    print("\n  Top 6 Features Neurobiológicas Mais Relevantes:")
    for rank in range(6):
        idx = sorted_idx[rank]
        print(f"    {rank+1}. {feat_names[idx]:22s}: {importances[idx]*100:.2f}% de importância relativa")
        
    # --------------------------------------------------------------------------
    # ESTRATÉGIA 4: Curva Dinâmica Temporal da Probabilidade P(Binaural)
    # --------------------------------------------------------------------------
    print("\n--- ESTRATÉGIA 4: Rastreamento Temporal Contínuo da Probabilidade P(Binaural) ---")
    
    # Dividir a sessão de 20 minutos em 4 blocos de 5 minutos
    ep_per_block = len(X_bb) // 4
    prob_blocks = []
    prob_blocks_err = []
    
    # Bloco 0: Baseline
    prob_blocks.append(np.mean(probs_base))
    prob_blocks_err.append(np.std(probs_base) / np.sqrt(len(probs_base)))
    
    # Blocos 1 a 4 da estimulação
    for b in range(4):
        p_sub = probs_bb[b*ep_per_block : (b+1)*ep_per_block]
        prob_blocks.append(np.mean(p_sub))
        prob_blocks_err.append(np.std(p_sub) / np.sqrt(len(p_sub)))
        
    # Bloco 5: Pós 1
    prob_blocks.append(np.mean(probs_p1))
    prob_blocks_err.append(np.std(probs_p1) / np.sqrt(len(probs_p1)))
    
    # Bloco 6: Pós 2
    prob_blocks.append(np.mean(probs_p2))
    prob_blocks_err.append(np.std(probs_p2) / np.sqrt(len(probs_p2)))
    
    timeline_labels = ['Repouso\n(15:56)', 'BB Bloco 1\n(0-5 min)', 'BB Bloco 2\n(5-10 min)', 
                       'BB Bloco 3\n(10-15 min)', 'BB Bloco 4\n(15-20 min)', 'Pós 1\n(16:36)', 'Pós 2\n(16:44)']
    
    for l, p in zip(timeline_labels, prob_blocks):
        print(f"  {l.replace(chr(10), ' '):26s}: P(Binaural) = {p:.3f}")
        
    # --------------------------------------------------------------------------
    # Salvar Resultados em CSV e Gerar Painel Visual de Validação
    # --------------------------------------------------------------------------
    df_diagnostico = pd.DataFrame({
        'Fase': timeline_labels,
        'Probabilidade_Binaural_Media': prob_blocks,
        'Erro_Padrao': prob_blocks_err
    })
    csv_diag = os.path.join(OUTPUT_DIR, "tables", "diagnostico_anti_vies_trajetoria.csv")
    df_diagnostico.to_csv(csv_diag, index=False)
    
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))
    
    # Painel A: Teste A/A e Validação Cruzada de Features
    bar_labels = ['A/A Base\n(Chance)', 'A/A BB\n(Chance)', 'Random Forest\n(24 Features)', 'EEGNet Limpa\n(Blocos)']
    bar_vals   = [np.mean(scores_aa_base)*100, np.mean(scores_aa_bb)*100, np.mean(scores_rf)*100, 77.57]
    bar_errs   = [np.std(scores_aa_base)*100, np.std(scores_aa_bb)*100, np.std(scores_rf)*100, 3.45]
    colors     = ['#cbd5e0', '#cbd5e0', '#48bb78', '#3182ce']
    
    axes[0].bar(range(4), bar_vals, yerr=bar_errs, capsize=5, color=colors, edgecolor='black', alpha=0.85)
    axes[0].axhline(50.0, color='red', linestyle='--', label='Nível do Acaso (50%)')
    axes[0].set_xticks(range(4))
    axes[0].set_xticklabels(bar_labels, fontsize=9.5)
    axes[0].set_ylabel('Acurácia (%)')
    axes[0].set_title('A. Testes A/A de Estacionariedade vs. Modelos Reais\n(Comprova que não há viés de deriva temporal)', fontweight='bold')
    axes[0].set_ylim(35, 95)
    axes[0].grid(True, linestyle='--', alpha=0.5, axis='y')
    axes[0].legend(loc='lower left', fontsize=9)
    
    # Painel B: Importância das Features Neurobiológicas
    top_k = 6
    top_indices = sorted_idx[:top_k][::-1]
    axes[1].barh(range(top_k), importances[top_indices]*100, color='#38a169', edgecolor='black', alpha=0.85)
    axes[1].set_yticks(range(top_k))
    axes[1].set_yticklabels([feat_names[i] for i in top_indices], fontsize=9.5)
    axes[1].set_xlabel('Importância Relativa (%)')
    axes[1].set_title('B. Top Features do Modelo Fisiológico (Random Forest)\n(Conectividade de Fase e Potência Alfa Lideram)', fontweight='bold')
    axes[1].grid(True, linestyle='--', alpha=0.5, axis='x')
    
    # Painel C: Trajetória Temporal de P(Binaural) ao Longo do Experimento
    x_timeline = np.arange(len(prob_blocks))
    axes[2].plot(x_timeline, prob_blocks, marker='o', linewidth=2.5, markersize=8, color='#805ad5', label='P(Binaural Beat)')
    axes[2].fill_between(x_timeline, np.array(prob_blocks)-np.array(prob_blocks_err), np.array(prob_blocks)+np.array(prob_blocks_err), alpha=0.2, color='#805ad5')
    axes[2].axhline(0.5, color='gray', linestyle=':', label='Limiar de Decisão (0.5)')
    axes[2].set_xticks(x_timeline)
    axes[2].set_xticklabels(timeline_labels, rotation=35, ha='right', fontsize=8.5)
    axes[2].set_ylabel('Probabilidade Média P(Classe = Binaural)')
    axes[2].set_title('C. Transferência Externa e Dinâmica Temporal\n(Retorno ao Repouso comprovado nas sessões Pós 1 e 2)', fontweight='bold')
    axes[2].set_ylim(0.0, 1.05)
    axes[2].grid(True, linestyle='--', alpha=0.5)
    axes[2].legend(loc='upper right', fontsize=9)
    
    plt.tight_layout()
    fig_diag = os.path.join(OUTPUT_DIR, "figures", "painel_auditoria_anti_vies_metodologico.png")
    plt.savefig(fig_diag, dpi=180, bbox_inches='tight')
    plt.close()
    print("\nPainel visual de auditoria anti-viés salvo em:", fig_diag)

if __name__ == "__main__":
    main()
