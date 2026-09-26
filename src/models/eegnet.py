"""
Implementação da arquitetura EEGNet (Lawhern et al., 2018 - Journal of Neural Engineering).
Inclui convolução temporal, convolução espacial depthwise com restrição MaxNorm e
convolução depthwise separável verdadeira.
"""

import torch
import torch.nn as nn

class Conv2dWithConstraint(nn.Conv2d):
    """Convolução 2D com restrição MaxNorm para evitar overfitting em dados de EEG."""
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
        super(EEGNet, self).__init__()
        
        # Bloco 1: Convolução Temporal (aprende filtros de frequência)
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
