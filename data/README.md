# Estrutura e Protocolo dos Dados de EEG

## Hardware de Aquisição
* **Dispositivo:** OpenBCI Ganglion (4 canais analógicos diferenciais)
* **Taxa de Amostragem ($f_s$):** $200\text{ Hz}$ (amostras transmitidas em pares via Bluetooth Low Energy)
* **Headset:** OpenBCI Ultracortex Mark IV (impresso em 3D) com eletrodos secos de pinos (*dry comb electrodes*)
* **Referência e Terra:** Clipes auriculares nos lóbulos A1 e A2 (referência neutra compartilhada)

## Mapeamento de Canais (Sistema Internacional 10–20)
| Canal OpenBCI | Eletrodo 10–20 | Região Cortical |
| :---: | :---: | :--- |
| **EXG Channel 0** | **F3** | Córtex Pré-frontal Dorsolateral Esquerdo |
| **EXG Channel 1** | **F4** | Córtex Pré-frontal Dorsolateral Direito |
| **EXG Channel 2** | **P3 / O1** | Região Centro-Parietal / Occipital Esquerda |
| **EXG Channel 3** | **P4 / O2** | Região Centro-Parietal / Occipital Direita |

## Arquivos dos Registros (`data/raw/`)
1. `OpenBCI-RAW-2023-11-28_15-56-59.txt` (5.8 min, 69.786 linhas):
   - **Condição:** Linha de Base em repouso pré-estímulo.
2. `OpenBCI-RAW-2023-11-28_16-04-54.txt` (20.4 min, 245.118 linhas):
   - **Condição:** Estimulação Contínua com Batimento Binaural (dividida analiticamente em 4 blocos de 5 minutos).
3. `OpenBCI-RAW-2023-11-28_16-36-02.txt` (5.2 min, 62.216 linhas):
   - **Condição:** Pós-Estimulação 1 (recuperação imediata após desligamento do áudio).
4. `OpenBCI-RAW-2023-11-28_16-44-15.txt` (5.3 min, 63.776 linhas):
   - **Condição:** Pós-Estimulação 2 (consolidação tardia).
