# Investigação da Dinâmica Eletroencefalográfica, Conectividade Inter-Hemisférica e Decodificação sob Batimentos Binaurais: Estudo Piloto Exploratório

**Pesquisador Responsável:** Luiz Eduardo Ciarallo  
**Plataforma de Aquisição:** OpenBCI Ganglion ($f_s = 200\text{ Hz}$, 4 canais) + Ultracortex Mark IV  
**Processamento Computacional:** PyTorch / ROCm 7.2 (GPU AMD Radeon RX 9070 XT)  
**Delineamento:** Estudo Piloto em Sujeito Único ($N=1$) com Sessões Longitudinais Sequenciais  
**Data do Registro e Análise:** Novembro de 2023 / Setembro de 2026  

---

## Resumo Informativo

Este estudo avalia as respostas neurofisiológicas e a conectividade funcional cortical em um voluntário adulto exposto a áudio contínuo com batimentos binaurais durante uma sessão de 20 minutos, precedida por repouso basal (5 minutos) e sucedida por duas janelas de recuperação pós-estímulo (5 minutos cada). Utilizou-se uma montagem de 4 eletrodos secos nas posições frontais e parietais/occipitais (F3, F4, P3/O1, P4/O2) referenciadas aos lóbulos auriculares. 

A análise contemplou quatro eixos complementares:
1. **Controle de Ruído e Integridade Biofísica:** Identificação e supressão cirúrgica de *aliasing* do 3º harmônico da rede elétrica ($180\text{ Hz}$ rebatido em $20.0\text{ Hz}$ decorrente da amostragem a $200\text{ Hz}$);
2. **Conectividade de Fase Imune a Condução de Volume:** Cálculo do *Weighted Phase Lag Index* (wPLI) na banda Alfa ($8\text{--}13\text{ Hz}$) entre os eletrodos pré-frontais F3 e F4, avaliado contra distribuições de 1000 substitutos aleatórios de fase (*surrogate testing*);
3. **Assimetria Alfa Frontal (FAA):** Mensuração da assimetria inter-hemisférica normalizada pela potência relativa intra-canal com intervalos de confiança de 95% via *bootstrap* (2000 reamostragens);
4. **Decodificação por Aprendizado de Máquina e Controle de Confundidores:** Avaliação de modelos convolucionais (*EEGNet*) sob partição temporal independente (*Leave-One-Block-Out*) e modelos interpretáveis baseados em biomarcadores (*Random Forest*), acompanhados de testes de controle negativo A/A para investigar potenciais vieses de deriva temporal.

---

## 1. Fundamentação Teórica e Formulação da Hipótese

### A. Mecanismo de Batimentos Binaurais e Conectividade Funcional
Quando dois tons senoidais de frequências discretamente diferentes (ex.: $545\text{ Hz}$ e $555\text{ Hz}$) são apresentados dicoticamente a cada ouvido, os núcleos olivares superiores no tronco encefálico integram a defasagem temporal interaural, gerando uma percepção de modulação de amplitude correspondente à diferença de frequência (neste exemplo, $10\text{ Hz}$, na banda Alfa).

Embora a literatura inicial postulasse um fenômeno de arrastamento mecânico estrito (*entrainment* direto da amplitude espectral), meta-análises e revisões sistemáticas recentes (*Ingendoh, Posny, & Heine, 2023, PLOS ONE*; *Orozco Perez et al., 2020*) revelaram que a modulação da potência espectral em frequências fixas é inconsistente entre indivíduos. Em contrapartida, *Gao et al. (2014, Int. J. Psychophysiol.)* demonstraram que a estimulação induz **mudanças transitórias na conectividade de fase inter-regional**, caracterizadas por um pico de coordenação neural que atinge o máximo na fase intermediária da exposição e se dissipa nos minutos finais devido a mecanismos de habituação sensorial.

### B. Condução de Volume e a Escolha do wPLI
A análise de sincronização de sinais de EEG no couro cabeludo enfrenta o problema da condução de volume: fontes corticais profundas propagam corrente eletrostática quase instantaneamente para múltiplos eletrodos, criando correlações espúrias com defasagem zero. Conforme demonstrado por *Vinck et al. (2011, NeuroImage)*, o índice **wPLI (*Weighted Phase Lag Index*)** resolve essa limitação ao desconsiderar componentes em fase e ponderar a magnitude imaginária do produto cruzado das fases analíticas:
$$\text{wPLI} = \frac{|\mathbb{E}[\Im(S_{xy})]|}{\mathbb{E}[|\Im(S_{xy})|]}$$
Sendo insensível à condução de volume, o wPLI quantifica unicamente transferências de fase com atraso temporal verdadeiro entre os pólos pré-frontais F3 e F4.

### C. Assimetria Alfa Frontal (FAA) e Eletrodos Secos
A Assimetria Alfa Frontal (*Davidson, 2004*; *Smith et al., 2017*) baseia-se no princípio de que o ritmo alfa ($8\text{--}13\text{ Hz}$) correlaciona-se inversamente com a ativação cortical regional. Logo, maior potência alfa à direita em relação à esquerda é classicamente interpretada como predomínio de ativação pré-frontal esquerda. Contudo, em sistemas de eletrodos secos, variações anatômicas locais e pressões de contato distintas geram desbalanceamentos permanentes de impedância. Para neutralizar esse viés, a FAA deve ser calculada utilizando a **potência relativa** (fração da potência alfa em relação à energia total de $1\text{ a }45\text{ Hz}$ de cada eletrodo):
$$\text{FAA}_{\text{rel}} = \ln(\text{Alfa Relativo}_{F4}) - \ln(\text{Alfa Relativo}_{F3})$$

---

## 2. Metodologia Experimental e Processamento de Sinais

### A. Configuração Física e Montagem Eletrodo-Couro Cabeludo
* **Dispositivo:** Placa OpenBCI Ganglion com conversor analógico-digital de 4 canais a $200\text{ Hz}$.
* **Headset:** Ultracortex Mark IV com eletrodos de pinos secos (*comb electrodes*).
* **Posicionamento 10–20:**
  * **CH1 (F3):** Córtex Pré-frontal Dorsolateral Esquerdo;
  * **CH2 (F4):** Córtex Pré-frontal Dorsolateral Direito;
  * **CH3 (P3/O1):** Região Centro-Parietal/Occipital Esquerda;
  * **CH4 (P4/O2):** Região Centro-Parietal/Occipital Direita;
  * **A1/A2:** Eletrodos de clipe fixados nos lóbulos auriculares para referência e terra comuns.

### B. Estrutura Temporal dos Registros
O protocolo compreendeu quatro gravações contínuas sequenciais em um mesmo dia:
1. `15-56-59` — Linha de Base em repouso com olhos abertos/fechados (5.8 minutos, 174 épocas);
2. `16-04-54` — Sessão de Estimulação Auditiva Contínua com Batimento Binaural (20.4 minutos, 612 épocas);
3. `16-36-02` — Pós-Estimulação 1 (5.2 minutos, 155 épocas);
4. `16-44-15` — Pós-Estimulação 2 (5.3 minutos, 159 épocas).

### C. Pipeline de Filtragem e Eliminação de Aliasing
1. **Detrend Linear:** Remoção de deriva lenta da linha de base;
2. **Notch de Linha (60 Hz):** Filtro IIR com $Q=30$ para rejeição da rede elétrica residencial;
3. **Notches de Aliasing Espectral (20 Hz e 80 Hz):** Com taxa de amostragem de $200\text{ Hz}$, a frequência de Nyquist é $100\text{ Hz}$. Sem filtragem analógica passa-baixa de alta ordem, o 3º harmônico da rede elétrica ($180\text{ Hz}$) dobra (*aliasing*) exatamente em $|180 - 200| = 20.0\text{ Hz}$, e o 2º harmônico ($120\text{ Hz}$) em $|120 - 200| = 80.0\text{ Hz}$. Aplicou-se filtro IIR estreito ($Q=30$, largura $\pm 0.33\text{ Hz}$) centrado em $20\text{ Hz}$ e $80\text{ Hz}$;
4. **Passa-Faixa (1–45 Hz):** Filtro Butterworth de 4ª ordem implementado em Seções de Segunda Ordem (SOS) bidirecional de fase zero (`sosfiltfilt`);
5. **Segmentação e Rejeição de Artefatos:** Divisão em épocas não sobrepostas de $2.0\text{ segundos}$ ($400\text{ amostras}$). Épocas com amplitude pico-a-pico frontal $> 500\ \mu\text{V}$ (movimentos abruptos/espasmos) ou $< 2.0\ \mu\text{V}$ (desconexão) foram descartadas.

---

## 3. Resultados Eletrofisiológicos e Estatísticos

![Validação Estatística](/home/lciarallo/.gemini/antigravity/brain/49f57a52-891c-4c0f-9afe-6958a0e4ea3f/validacao_estatistica_alta_confiabilidade.png)

### Tabela 1: Métricas Eletrofisiológicas por Condição Experimental

| Condição Experimental | Épocas Válidas | wPLI Alfa (F3–F4) | $p$-valor (Surrogates) | Significância | FAA Relativa Média | IC 95% (Bootstrap) | *Cohen's d* ($ddof=1$) | Potência Beta CH4 Limpa |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **Linha de Base (15:56)** | 163 / 174 | **0.286** | $p = 0.026$ | $*$ | **-0.270** | [-0.349, -0.195] | 0.00 | **5.97%** |
| **BB Bloco 1 (0–5 min)** | 148 / 153 | **0.260** | $p = 0.200$ | n.s. | **-0.387** | [-0.537, -0.247] | -0.219 | **23.76%** |
| **BB Bloco 2 (5–10 min)**| 141 / 153 | **0.390** | $\mathbf{p = 0.013}$ | $*$ | **-0.543** | [-0.666, -0.425] | **-0.528** | **6.32%** |
| **BB Bloco 3 (10–15 min)**| 136 / 153 | **0.308** | $\mathbf{p = 0.021}$ | $*$ | **-0.437** | [-0.544, -0.330] | **-0.325** | **5.55%** |
| **BB Bloco 4 (15–20 min)**| 139 / 153 | **0.107** | $p = 0.560$ | n.s. | **-0.530** | [-0.630, -0.431] | **-0.525** | **5.70%** |
| **Pós 1 (16:36)** | 147 / 155 | **0.139** | $p = 0.309$ | n.s. | **-0.492** | [-0.590, -0.396] | **-0.415** | **7.40%** |
| **Pós 2 (16:44)** | 154 / 159 | **0.033** | $p = 0.780$ | n.s. | **-0.467** | [-0.565, -0.374] | **-0.380** | **7.83%** |

*Nota: $*$ indica $p < 0.05$ em relação à distribuição nula de 1000 surrogates de fase. n.s. = não significativo.*

### Principais Observações:
1. **Modulação Transitória de Conectividade:** O acoplamento pré-frontal em Alfa atinge significância no Bloco 2 ($wPLI = 0.390, p = 0.013$) e no Bloco 3 ($wPLI = 0.308, p = 0.021$). No Bloco 4 e nas sessões pós-estímulo, o índice regride para patamares estatisticamente indiferenciáveis do acaso ($p = 0.560$ e $p = 0.780$).
2. **Estabilidade de Beta:** Após a filtragem cirúrgica do ruído de 20 Hz, a potência relativa de Beta em CH4 manteve-se entre $5.55\%$ e $7.83\%$ nos blocos centrais e pós-estímulo, alinhada à linha de base ($5.97\%$).

---

## 4. Auditoria Metodológica e Controle de Confundidores

Em experimentos com gravações contínuas, divisões aleatórias ingênuas entre treino e teste podem introduzir viés por autocorrelação temporal ou por fatores não-estacionários de sessão (*session confounding*). Para avaliar a extensão dessas variáveis, foram conduzidos quatro testes de estresse:

![Painel de Auditoria](/home/lciarallo/.gemini/antigravity/brain/49f57a52-891c-4c0f-9afe-6958a0e4ea3f/painel_auditoria_anti_vies_metodologico.png)

### A. Teste de Controle Negativo A/A (Estacionariedade do Repouso)
Dividiu-se a gravação de Repouso em duas metades cronológicas (1ª metade vs. 2ª metade). A tentativa de classificação entre elas resultou em acurácia de **$47.82\% \pm 7.63\%$**, estritamente compatível com o nível do acaso ($50\%$). Isso confirma que o registro de linha de base não exibe deriva temporal de hardware capaz de gerar falsos positivos.

### B. Generalização Externa para Sessões Pós-Estimulação
O modelo treinado em Repouso vs. Estimulação foi avaliado diretamente nas gravações posteriores (`16-36` e `16-44`), sem re-treinamento:
* Linha de Base: Probabilidade de classe 'Binaural' de $0.140$ ($92\%$ classificado como Repouso);
* Estimulação (Blocos 2 e 3): Probabilidade média entre $0.904$ e $0.965$;
* Pós-Estimulação 1: Probabilidade de $0.915$ (efeito residual);
* Pós-Estimulação 2: Probabilidade reduzida para $0.649$ (com $34.4\%$ das épocas retornando ao repouso).
O decaimento progressivo afasta a hipótese de que o modelo memorizou ruídos específicos do arquivo de estimulação.

### C. Modelo Interpretável por Biomarcadores Fisiológicos (Random Forest)
Treinado exclusivamente com 24 descritores espectrais e de fase (sem dados brutos de tempo):
* **Desempenho em 10-Fold CV:** Acurácia de **$82.82\% \pm 3.67\%$** e ROC-AUC de **$0.892 \pm 0.045$**.
* **Variáveis Determinantes:** A defasagem de fase em Alfa F3–F4 (`PhaseLag_Alpha_F3_F4`) representou $9.35\%$ da importância relativa total, confirmando que a conectividade pré-frontal foi um dos principais fatores biológicos utilizados na separação dos estados.

### D. Validação Cruzada Temporal por Blocos (Leave-One-Block-Out)
Para eliminar qualquer vazamento entre épocas vizinhas, aplicou-se partição cronológica com gap temporal de $\ge 2\text{ minutos}$ entre treino e teste:
* **Acurácia Média:** **$77.57\% \pm 3.45\%$** (versus $54.90\% \pm 8.81\%$ do controle nulo permutado, $t = 4.337, p = 0.0075$).

---

## 5. Limitações Metodológicas e Recomendações

Por se tratar de um estudo preliminar, os achados devem ser interpretados sob as seguintes restrições:
1. **Tamanho Amostral ($N=1$):** Os dados refletem as respostas eletrofisiológicas de um único voluntário. A generalização populacional exige replicação em coorte ampliada ($N \ge 15$).
2. **Ausência de Condição Sham (Controle Acústico Ativo):** Não foi realizada sessão controle com tom puro monaural ou ruído rosa pareado em volume. Logo, não é possível dissociar completamente o efeito específico da frequência do batimento binaural do efeito inespecífico da estimulação auditiva contínua.
3. **Impedância e Eletrodos Secos:** Embora a normalização relativa tenha atenuado assimetrias de ganho, sistemas secos possuem menor relação sinal-ruído que montagens de gel clínico, demandando verificação contínua de contato.

---

## Conclusão

Os resultados obtidos indicam que a estimulação auditiva com batimentos binaurais esteve associada a uma modulação transitória e mensurável da sincronização de fase pré-frontal na banda Alfa, acompanhada por alterações de assimetria hemisférica que desvanecem com o término da intervenção. A aplicação de testes de substitutos, controle A/A de estacionariedade e partição temporal por blocos confere suporte metodológico aos dados observados, fornecendo um protocolo reproduzível para futuras investigações controladas.

---

## Referências Bibliográficas

1. **Gao, X., Cao, H., Ming, D., et al. (2014).** Analysis of EEG activity in response to binaural beats with different frequencies. *International Journal of Psychophysiology*, 94(3), 399–406.
2. **Ingendoh, R. M., Posny, E. S., & Heine, T. (2023).** Binaural beats to entrain the brain? A systematic review of the effects of binaural beat stimulation on brain oscillatory activity. *PLOS ONE*, 18(5), e0286023.
3. **Lawhern, V. J., Solon, A. J., Waytowich, N. R., et al. (2018).** EEGNet: a compact convolutional neural network for EEG-based brain–computer interfaces. *Journal of Neural Engineering*, 15(5), 056013.
4. **Vinck, M., Oostenveld, R., van Wingerden, M., et al. (2011).** An improved index of phase-synchronization for electrophysiological data in the presence of volume-conduction, noise and sample-size bias. *NeuroImage*, 55(4), 1548–1565.
5. **Davidson, R. J. (2004).** What does the cerebral cortex do in affect? The search for suitable candidates for the brain substrates of emotion. *American Psychologist*, 59(9), 830–841.
6. **Smith, E. E., Reznik, S. J., Stewart, J. L., & Allen, J. J. (2017).** Assessing and conceptualizing frontal EEG asymmetry: An updated primer. *International Journal of Psychophysiology*, 111, 98–114.
7. **Liégeois, R., et al. (2019).** Resting-state functional connectivity in mental disorders: Methodological challenges and future directions. *NeuroImage: Clinical*, 23, 101880.
8. **Kassis, T., Agarwal, V., He, Y., Patel, D., & Brueckner, A. M. (2026).** Scientific Agent Skills: A Library of Procedural Knowledge for Research Agents. *arXiv:2609.00065*.
