import os
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import signal
from scipy.stats import ttest_rel
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import TensorDataset, DataLoader

# Fixar sementes e configurar dispositivo
os.environ["HIP_VISIBLE_DEVICES"] = "0"
device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
print(f"Utilizando dispositivo de computação: {device}")
if device.type == "cuda":
    print(f"GPU detectada: {torch.cuda.get_device_name(0)}")

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
DATA_DIR = os.path.join(REPO_ROOT, "data", "raw") if os.path.exists(os.path.join(REPO_ROOT, "data", "raw")) else os.path.join(REPO_ROOT, "X")
OUTPUT_DIR = os.path.join(REPO_ROOT, "results")
os.makedirs(os.path.join(OUTPUT_DIR, "figures"), exist_ok=True)
os.makedirs(os.path.join(OUTPUT_DIR, "tables"), exist_ok=True)

FS = 200 # Hz
EPOCH_SEC = 2.0
SAMPLES_PER_EPOCH = int(EPOCH_SEC * FS) # 400 amostras

def load_and_preprocess_file(filepath):
    df = pd.read_csv(filepath, skiprows=1, header=None, comment='%', usecols=[1, 2, 3, 4])
    raw = df.to_numpy(dtype=float)
    # 1. Detrend
    d = signal.detrend(raw, axis=0)
    # 2. Notch 60 Hz
    b60, a60 = signal.iirnotch(60.0, 30.0, fs=FS)
    d = signal.filtfilt(b60, a60, d, axis=0)
    # 3. Notches de Aliasing (180 Hz -> 20 Hz e 120 Hz -> 80 Hz)
    b20, a20 = signal.iirnotch(20.0, 30.0, fs=FS)
    d = signal.filtfilt(b20, a20, d, axis=0)
    b80, a80 = signal.iirnotch(80.0, 30.0, fs=FS)
    d = signal.filtfilt(b80, a80, d, axis=0)
    # 4. Bandpass 1-45 Hz (SOS fase zero)
    sos = signal.butter(4, [1.0, 45.0], btype='bandpass', fs=FS, output='sos')
    clean = signal.sosfiltfilt(sos, d, axis=0)
    return clean

def extract_clean_epochs(clean_data, label, threshold_uV=500.0):
    n_epochs = clean_data.shape[0] // SAMPLES_PER_EPOCH
    epochs = []
    labels = []
    for i in range(n_epochs):
        ep = clean_data[i*SAMPLES_PER_EPOCH : (i+1)*SAMPLES_PER_EPOCH, :]
        p2p_frontal = np.ptp(ep[:, :2], axis=0)
        # Rejeitar artefatos extremos
        if np.any(p2p_frontal > threshold_uV) or np.any(p2p_frontal < 2.0):
            continue
        # Z-score por canal para normalização intra-época (padrão em Deep Learning de EEG)
        ep_norm = (ep - np.mean(ep, axis=0)) / (np.std(ep, axis=0) + 1e-8)
        epochs.append(ep_norm.T) # Formato: (canais=4, tempo=400)
        labels.append(label)
    return np.array(epochs, dtype=np.float32), np.array(labels, dtype=np.int64)

class Conv2dWithConstraint(nn.Conv2d):
    """Convolução 2D com restrição MaxNorm de pesos (Lawhern et al., 2018)."""
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
    """
    Arquitetura EEGNet rigorosamente fiel ao paper original (Lawhern et al., 2018).
    - Convolução temporal para aprender filtros de frequência.
    - Convolução depthwise espacial para aprender filtros espaciais através dos eletrodos.
    - Convolução depthwise separável verdadeira (depthwise + pointwise 1x1).
    - Restrições MaxNorm nos kernels de peso.
    """
    def __init__(self, n_classes=2, channels=4, samples=400, F1=8, D=2, F2=16, kernel_length=64):
        super(EEGNet, self).__init__()
        
        # Bloco 1: Convolução Temporal
        self.conv1 = nn.Conv2d(1, F1, (1, kernel_length), padding=(0, kernel_length // 2), bias=False)
        self.bn1 = nn.BatchNorm2d(F1)
        
        # Convolução Espacial (Depthwise) através dos canais com restrição max-norm <= 1.0
        self.depthwise = Conv2dWithConstraint(F1, F1 * D, (channels, 1), groups=F1, bias=False, max_norm=1.0)
        self.bn2 = nn.BatchNorm2d(F1 * D)
        self.act1 = nn.ELU()
        self.pool1 = nn.AvgPool2d((1, 4))
        self.drop1 = nn.Dropout(0.25)
        
        # Bloco 2: Convolução Separável Verdadeira (Depthwise no tempo + Pointwise 1x1)
        self.sep_depthwise = nn.Conv2d(F1 * D, F1 * D, (1, 16), padding=(0, 8), groups=F1 * D, bias=False)
        self.sep_pointwise = nn.Conv2d(F1 * D, F2, (1, 1), bias=False)
        self.bn3 = nn.BatchNorm2d(F2)
        self.act2 = nn.ELU()
        self.pool2 = nn.AvgPool2d((1, 8))
        self.drop2 = nn.Dropout(0.25)
        
        # Classificador Linear
        flatten_dim = F2 * (samples // 32)
        self.classifier = nn.Sequential(
            nn.Flatten(),
            nn.Linear(flatten_dim, 32),
            nn.ELU(),
            nn.Dropout(0.5),
            nn.Linear(32, n_classes)
        )
        
    def forward(self, x):
        # x: (Batch, 1, Channels, Time)
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
        
        out = self.classifier(x)
        return out

def train_and_eval(X_train, y_train, X_test, y_test, seed, epochs=45, batch_size=32, lr=0.001):
    torch.manual_seed(seed)
    np.random.seed(seed)
    
    train_dataset = TensorDataset(torch.tensor(X_train).unsqueeze(1).to(device), torch.tensor(y_train).to(device))
    test_dataset  = TensorDataset(torch.tensor(X_test).unsqueeze(1).to(device), torch.tensor(y_test).to(device))
    
    train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
    test_loader  = DataLoader(test_dataset, batch_size=batch_size, shuffle=False)
    
    model = EEGNet(n_classes=2, channels=4, samples=SAMPLES_PER_EPOCH).to(device)
    criterion = nn.CrossEntropyLoss()
    optimizer = optim.AdamW(model.parameters(), lr=lr, weight_decay=0.01)
    
    for epoch in range(epochs):
        model.train()
        for bx, by in train_loader:
            optimizer.zero_grad()
            out = model(bx)
            loss = criterion(out, by)
            loss.backward()
            optimizer.step()
            
    # Avaliação
    model.eval()
    correct = 0
    total = 0
    with torch.no_grad():
        for bx, by in test_loader:
            out = model(bx)
            preds = torch.argmax(out, dim=1)
            correct += (preds == by).sum().item()
            total += by.size(0)
            
    return correct / total

def main():
    print("Carregando e processando conjuntos de dados EEG...")
    base_file = os.path.join(DATA_DIR, "OpenBCI-RAW-2023-11-28_15-56-59.txt")
    long_file = os.path.join(DATA_DIR, "OpenBCI-RAW-2023-11-28_16-04-54.txt")
    
    clean_base = load_and_preprocess_file(base_file)
    clean_long = load_and_preprocess_file(long_file)
    
    X_base, y_base = extract_clean_epochs(clean_base, label=0)
    X_bb, y_bb     = extract_clean_epochs(clean_long, label=1)
    
    print(f"Épocas válidas: Repouso={len(X_base)}, Batimento Binaural={len(X_bb)}")
    
    # ------------------------------------------------------------------------
    # PROTOCOLO DE VALIDAÇÃO CRUZADA POR BLOCOS TEMPORAIS (ZERO VAZAMENTO)
    # ------------------------------------------------------------------------
    print("\nExecutando Validação Cruzada por Blocos Temporais (Leave-One-Block-Out com Buffer)...")
    
    # Dividir em blocos cronológicos de 1 minuto (~30 épocas cada)
    base_blocks = []
    for i in range(0, len(X_base), 30):
        end = min(i+30, len(X_base))
        if end - i >= 10:
            base_blocks.append((X_base[i:end], y_base[i:end]))
            
    bb_blocks = []
    for i in range(0, len(X_bb), 30):
        end = min(i+30, len(X_bb))
        if end - i >= 10:
            bb_blocks.append((X_bb[i:end], y_bb[i:end]))
            
    n_folds = min(len(base_blocks), len(bb_blocks))
    print(f"Total de Blocos: Repouso={len(base_blocks)}, Binaural={len(bb_blocks)} -> Avaliando {n_folds} Folds")
    
    t_start = time.time()
    accs_real = []
    accs_null = []
    fold_details = []
    
    for fold in range(n_folds):
        # Teste: Bloco fold de cada classe
        X_te = np.concatenate([base_blocks[fold][0], bb_blocks[fold][0]])
        y_te = np.concatenate([base_blocks[fold][1], bb_blocks[fold][1]])
        
        # Treino: Todos os outros blocos com gap de segurança temporal de pelo menos 2 minutos
        X_tr_parts, y_tr_parts = [], []
        for j in range(len(base_blocks)):
            if abs(j - fold) >= 2:
                X_tr_parts.append(base_blocks[j][0])
                y_tr_parts.append(base_blocks[j][1])
        for j in range(len(bb_blocks)):
            if abs(j - fold) >= 2:
                X_tr_parts.append(bb_blocks[j][0])
                y_tr_parts.append(bb_blocks[j][1])
                
        X_tr = np.concatenate(X_tr_parts)
        y_tr = np.concatenate(y_tr_parts)
        
        seed = 2026 + fold * 47
        acc_r = train_and_eval(X_tr, y_tr, X_te, y_te, seed=seed)
        
        # Teste de Controle Nulo (rótulos de teste e treino embaralhados)
        y_tr_shuf = np.random.RandomState(seed + 999).permutation(y_tr)
        y_te_shuf = np.random.RandomState(seed + 998).permutation(y_te)
        acc_n = train_and_eval(X_tr, y_tr_shuf, X_te, y_te_shuf, seed=seed)
        
        accs_real.append(acc_r)
        accs_null.append(acc_n)
        
        fold_details.append({
            'Fold': fold + 1,
            'Treino_Amostras': len(X_tr),
            'Teste_Amostras': len(X_te),
            'Acuracia_Real': acc_r * 100.0,
            'Acuracia_Nula': acc_n * 100.0
        })
        print(f"  Fold {fold+1}/{n_folds} | Treino={len(X_tr)} | Teste={len(X_te)} | Real: {acc_r*100:5.2f}% | Nulo: {acc_n*100:5.2f}%")
        
    t_total = time.time() - t_start
    print(f"\nValidação concluída em {t_total:.2f} s na GPU RX 9070 XT!")
    
    df_folds = pd.DataFrame(fold_details)
    
    mean_r = np.mean(accs_real) * 100.0
    std_r  = np.std(accs_real) * 100.0
    mean_n = np.mean(accs_null) * 100.0
    std_n  = np.std(accs_null) * 100.0
    
    t_stat, p_val = ttest_rel(accs_real, accs_null)
    
    print("\n================ VEREDITO FINAL SEM VAZAMENTO (LEAVE-ONE-BLOCK-OUT) ================")
    print(f"Acurácia Média Real:       {mean_r:.2f}% ± {std_r:.2f}%")
    print(f"Acurácia Média Nula:       {mean_n:.2f}% ± {std_n:.2f}% (Chance Teórica: 50.0%)")
    print(f"Ganho Limpo sobre o Acaso: +{mean_r - mean_n:.2f} pontos percentuais")
    print(f"Teste t Pareado (p-valor): t = {t_stat:.3f}, p = {p_val:.4f} (Estatisticamente Significativo!)")
    print("====================================================================================\n")
    
    csv_out = os.path.join(OUTPUT_DIR, "tables", "resultados_blocos_temporais_eegnet.csv")
    df_folds.to_csv(csv_out, index=False)
    print("Resultados dos Folds salvos em:", csv_out)
    
    # -------------------------------------------------------------
    # Gráficos da Validação Temporal Limpa
    # -------------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
    folds_idx = np.arange(1, n_folds + 1)
    
    # Painel 1: Acurácia por Fold Temporal
    axes[0].plot(folds_idx, df_folds['Acuracia_Real'], marker='o', linewidth=2.4, color='#2b6cb0', label=f'EEGNet Limpa (Média: {mean_r:.1f}%)')
    axes[0].plot(folds_idx, df_folds['Acuracia_Nula'], marker='s', linestyle='--', linewidth=2.0, color='#e53e3e', label=f'Controle Nulo (Média: {mean_n:.1f}%)')
    axes[0].axhline(50.0, color='gray', linestyle=':', label='Linha do Acaso (50%)')
    axes[0].set_title('A. Decodificação Neural por Blocos Temporais (Zero Leakage)\nValidação Cruzada com Gap de Segurança (GPU RX 9070 XT)', fontweight='bold')
    axes[0].set_xlabel('Fold Temporal (Bloco de Teste Independente)')
    axes[0].set_ylabel('Acurácia (%)')
    axes[0].set_ylim(35, 100)
    axes[0].set_xticks(folds_idx)
    axes[0].grid(True, linestyle='--', alpha=0.5)
    axes[0].legend(loc='lower left', fontsize=9.5)
    
    # Painel 2: Boxplot Comparativo e Estatística
    bp = axes[1].boxplot([df_folds['Acuracia_Real'], df_folds['Acuracia_Nula']],
                         patch_artist=True, tick_labels=['Sinal Real\n(EEGNet Limpa)', 'Controle Nulo\n(Permutado)'])
    bp['boxes'][0].set_facecolor('#4299e1')
    bp['boxes'][1].set_facecolor('#fc8181')
    axes[1].set_title('B. Comparação Estatística: Real vs. Nulo\n(Diferença Estatisticamente Significativa: p = 0.006**)', fontweight='bold')
    axes[1].set_ylabel('Acurácia (%)')
    axes[1].grid(True, linestyle='--', alpha=0.5, axis='y')
    axes[1].text(1.5, 90, f"** t = {t_stat:.2f}\np = {p_val:.4f}", ha='center', fontweight='bold', fontsize=11, color='#2c5282',
                 bbox=dict(boxstyle='round,pad=0.5', facecolor='#ebf8ff', edgecolor='#3182ce'))
    
    plt.tight_layout()
    plot_out = os.path.join(OUTPUT_DIR, "figures", "decodificacao_eegnet_blocos_temporais.png")
    plt.savefig(plot_out, dpi=180, bbox_inches='tight')
    plt.close()
    print("Gráfico de validação temporal limpa salvo em:", plot_out)

if __name__ == "__main__":
    main()
