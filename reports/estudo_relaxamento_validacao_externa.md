# Áudio binaural e relaxamento: estudo do EEG com comparação externa

**Relatório de 28/09/2026.** Pergunta: os registros locais apresentam mudanças de EEG que sustentem aumento de relaxamento durante o áudio? Este estudo acrescenta duas bases independentes à auditoria de 26/09/2026. É uma reanálise exploratória, sem novas coletas e sem escala subjetiva de relaxamento disponível no piloto.

## Resultado principal

Os dados locais **não confirmaram aumento de relaxamento pelos indicadores avaliados**. O alfa frontal relativo foi menor em todos os quatro blocos do áudio que na base, nas quatro configurações de notch/limiar examinadas. Na EEGMAT, o mesmo descritor foi, em média, maior no repouso que no cálculo mental. Essa discordância não prova piora de relaxamento: repouso não é uma medida subjetiva de relaxamento e as gravações têm diferenças de aquisição e de seleção de épocas.

A base externa com áudio e questionário permitiu testar a especificidade binaural: nenhum dos seis contrastes de questionário C1-C0 apresentou p de Holm menor que 0,05. Os dados de EEG frontal também não estabeleceram diferença após ajuste. Esses resultados não são prova de equivalência ou ausência de efeito. A cobertura de canais, especialmente para EEG, foi pequena e diferente da amostra descrita no artigo.

## 1. Fontes e alcance da validação

- **Piloto local:** quatro registros de um participante, rotulados Base, BB, Pós 1 e Pós 2. O usuário esclareceu que o objetivo era relaxamento. Os canais posteriores têm posição ambígua; por isso a comparação com sinais externos usa apenas F3/F4.
- **EEGMAT 1.0.0:** 36 participantes, repouso e cálculo mental, 72 arquivos EDF. É uma referência de contraste entre estados experimentais, não uma base rotulada por relaxamento subjetivo. Fonte: [PhysioNet EEGMAT](https://physionet.org/content/eegmat/1.0.0/).
- **Sudre et al. (2024):** arquivos derivados de EEG e questionários disponibilizados no [OSF stzfd](https://osf.io/stzfd/), com condições binaural C1 e monaural C0, grupos G6/G40. O estudo publicado fornece o protocolo e a interpretação das escalas: [artigo e dicionário de dados](https://doi.org/10.1371/journal.pone.0306427).

Validar externamente significa aqui testar descritores em uma amostra independente e reanalisar contrastes de uma base de áudio. Não significa criar um grupo controle para o participante local ou demonstrar que o efeito de outra população ocorreu nele. Não foram combinados todos os participantes como se tivessem recebido o mesmo protocolo.

## 2. Plano e integridade dos dados

O plano local `reports/plano_validacao_externa.md` foi escrito antes do cálculo dos resultados externos nesta execução. Não é pré-registro público. Incluiu todos os 36 participantes EEGMAT, quatro descritores definidos previamente e contrastes binaural-monaural separados por frequência. URLs, tamanhos e hashes dos 81 arquivos baixados estão em `download_manifest.json`; os EDF conferem com os SHA-256 oficiais. Os sinais locais não foram enviados a serviços externos.

As durações reais EEGMAT foram verificadas pelos cabeçalhos: repouso de 80 a 188 s e tarefa de 62 s. Não se assumiu que todos os registros tinham a duração simplificada apresentada na página. Foi utilizada uma janela de 60 s finais do repouso e 60 s iniciais da tarefa; após bordas, todos os participantes retiveram 28 épocas por condição.

| Grupo | IDs de questionário | IDs de EEG | Pares frontais completos |
| --- | --- | --- | --- |
| G6 | 22 | 20 | 7 |
| G40 | 19 | 19 | 5 |

Nos CSVs de Sudre, G6 tem 22 IDs de questionário e 20 de EEG, enquanto o artigo descreve 19 por grupo. Não havia, nos arquivos inspecionados, lista suficiente para reconstruir todas as exclusões do artigo; esta reanálise usa casos completos documentados, sem selecionar IDs para chegar a 19. Em G40 foram identificadas e retiradas em memória 1.092 duplicatas exatas, todas na condição C2, fora do contraste principal. A exigência de F3/Fz/F4, base e três blocos de estimulação em C0 e C1 deixou apenas 7 e 5 pares de EEG. Os IDs incluídos/excluídos estão em `sudre_coverage.csv`. Portanto, esta não é uma reprodução exata da análise publicada.

## 3. Método comum do piloto e EEGMAT

Foram usados F3 e F4 em µV, detrend e Butterworth de quarta ordem 1-30 Hz por janela, reamostragem para 200 Hz, remoção de 2 s de cada borda e épocas de 2 s sem sobreposição. Rejeição principal: amplitude pico-a-pico frontal fora de 2-150 µV. Os quatro descritores foram:

1. Alfa relativo: potência 8-13 Hz dividida pela potência 1-30 Hz, média de F3/F4.
2. Theta relativo: potência 4-8 Hz dividida pela potência 1-30 Hz, média de F3/F4.
3. Razão logarítmica alfa/beta baixa: média de ln(potência 8-13 Hz / potência 13-18 Hz). Evita incluir diretamente o pico de 20 Hz na banda beta.
4. FAA: ln(alfa relativo F4) - ln(alfa relativo F3).

A potência usa Welch, janela de 400 amostras. Na EEGMAT cada pessoa/condição foi resumida pela mediana das épocas. As diferenças pareadas foram testadas com Wilcoxon bilateral e Holm sobre quatro indicadores; ICs percentis de 95% da diferença média usam 5.000 reamostragens de participantes. Semente: 20260928. Nenhum desses descritores é uma escala validada de relaxamento por si só.

No piloto, cada janela de cinco minutos gerou 148 épocas candidatas após remoção das bordas. As contagens diferem da auditoria anterior porque a filtragem e as bordas agora são definidas por janela. A sensibilidade inclui 500 µV e notch de 20 Hz ligado/desligado. O processamento comum não remove diferenças de referência, eletrodos, olhos ou limpeza prévia pelos autores da EEGMAT.

## 4. O que aparece no EEG local

| Condição | Épocas / 148 | Alfa (%) | Theta (%) | ln(alfa/beta baixa) | FAA |
| --- | --- | --- | --- | --- | --- |
| Base | 142/148 | 11.19 | 15.78 | 0.703 | -0.247 |
| BB 1 | 74/148 | 9.57 | 14.79 | 0.383 | -0.389 |
| BB 2 | 80/148 | 10.53 | 15.94 | 0.624 | -0.514 |
| BB 3 | 106/148 | 10.19 | 16.38 | 0.461 | -0.412 |
| BB 4 | 103/148 | 7.51 | 15.15 | 0.549 | -0.614 |
| Pós 1 | 143/148 | 9.71 | 12.38 | 0.556 | -0.537 |
| Pós 2 | 124/148 | 12.72 | 15.01 | 0.600 | -0.376 |

A mediana do alfa relativo passa de 11.19% na base para 7.51% no BB 4, recuperando-se para 12.72% na Pós 2. Alfa menor durante o áudio e maior na última pós-sessão é uma descrição do sinal; não se transforma automaticamente em medida de tensão, ansiedade, relaxamento ou efeito residual. A razão alfa/beta baixa também fica abaixo da base nos quatro blocos do áudio.

Para reduzir a falsa precisão produzida por épocas próximas, a análise adicional usa médias de blocos cronológicos de 30 s. Só entram blocos com ao menos metade das épocas candidatas, mantendo os índices de aquisição. Os ICs abaixo reamostram blocos separadamente na base e em cada condição, 5.000 vezes:

| Condição | Blocos válidos / 10 | Delta alfa (p.p.) [IC 95%] |
| --- | --- | --- |
| BB 1 | 5 | -2.39 [-3.86; -1.02] |
| BB 2 | 6 | -2.49 [-4.78; 0.45] |
| BB 3 | 8 | -1.60 [-3.77; 0.75] |
| BB 4 | 8 | -4.78 [-5.95; -3.54] |
| Pós 1 | 10 | -1.89 [-3.09; -0.69] |
| Pós 2 | 9 | 1.37 [-1.23; 3.93] |

Estes deltas são diferenças de médias de blocos, não subtrações das medianas da tabela anterior. Os ICs são exploratórios, sem ajuste simultâneo, condicionados a um registro e não entre participantes. Blocos adjacentes podem continuar correlacionados. Não foram produzidos p-valores locais confirmatórios. A tabela completa contém todos os quatro descritores e a sensibilidade de 60 s por bloco. A baixa retenção em BB 1/BB 2 e o descarte de blocos podem introduzir seleção desigual.

## 5. Teste dos indicadores em 36 participantes externos

Diferença positiva significa repouso maior que cálculo mental. Alfa/theta estão em fração da potência total; 0,01 equivale a um ponto percentual.

| Indicador | Repouso - tarefa [IC 95%] | p Holm | Pessoas na direção positiva |
| --- | --- | --- | --- |
| Alfa relativo | 0.0569 [0.0320; 0.0808] | 0.0004 | 75.0% |
| Theta relativo | -0.0216 [-0.0362; -0.0080] | 0.0180 | 27.8% |
| ln(alfa / beta baixa) | 0.0707 [-0.0494; 0.1899] | 0.4527 | 58.3% |
| FAA relativa | 0.0314 [-0.0343; 0.1011] | 0.5497 | 55.6% |

Alfa e theta distinguiram as condições em média, após Holm. Alfa foi maior no repouso em 75% dos participantes, com diferença média de aproximadamente 5,69 pontos percentuais. Theta apresentou a direção inversa. Já a razão alfa/beta baixa e a FAA não estabeleceram diferença após correção. Portanto, não é justificável tratar todos esses descritores como marcadores intercambiáveis de relaxamento.

No piloto, a queda do alfa durante BB não acompanha a direção do repouso na EEGMAT. Theta local não apresenta uma mudança uniforme nos blocos do áudio. Essa comparação oferece evidência contra a interpretação simples “o áudio elevou um padrão de repouso”, mas não constitui teste causal de relaxamento: o áudio também pode mudar atenção, e o estado dos olhos e a qualidade de contato não estão controlados.

![Comparação do piloto, EEGMAT e questionários externos](../results/validation_2026_09_28/comparacao_relaxamento.png)

## 6. O modelo externo valida um detector de relaxamento?

Não. A regressão logística usou os quatro descritores, sem seleção pelos resultados locais, com padronização aprendida no treino. A validação deixou uma pessoa inteira de fora em cada iteração (36 folds). Cada pessoa contribuiu com uma observação de repouso e uma de tarefa. Resultado:

- Acurácia balanceada: **58.33%**; IC bootstrap por participante [51.39%; 65.28%].
- ROC-AUC: **0.614**.
- Controle: 199 permutações com troca dos rótulos dentro de cada pessoa, repetindo toda a validação. p empírico **0.075**; média nula 49.87%.

O IC descritivo da acurácia e o teste de permutação usam construções diferentes; não devem ser interpretados como testes equivalentes. O teste de permutação não atingiu 0,05. A amostra é pequena para certificar desempenho individual e não há calibração externa como escala de relaxamento.

A aplicação diagnóstica do modelo às janelas locais foi salva em `local_external_model_transfer.csv`. Entre 1 e 3 dos 4 descritores locais ficam fora do intervalo marginal de 95% dos dados externos, dependendo da condição. Isso sinaliza mudança de domínio; não é um teste multivariado formal. Por esse motivo, e pelo desempenho externo modesto, suas probabilidades não são apresentadas como “percentual de relaxamento”.

## 7. O que diz a base com áudios e questionário

Reanalisaram-se três escalas, PA, CT e PT, com dados pareados C1-C0, separadamente por grupo. Na codificação descrita pelos autores, valores maiores indicam maior relaxamento, apesar dos nomes CT/PT conterem “tensão”. Holm cobre os seis contrastes. Os ICs são bootstrap da diferença média por participante; p-valores vêm do teste pareado Wilcoxon.

| Grupo / escala | N | Binaural - monaural [IC 95%] | p Holm |
| --- | --- | --- | --- |
| G6 / PA | 22 | -1.64 [-5.64; 2.27] | 1.000 |
| G6 / CT | 22 | 0.73 [-0.86; 2.41] | 1.000 |
| G6 / PT | 22 | -0.55 [-2.64; 1.27] | 1.000 |
| G40 / PA | 19 | -3.21 [-10.74; 3.63] | 1.000 |
| G40 / CT | 19 | 2.11 [-1.21; 5.26] | 0.774 |
| G40 / PT | 19 | -0.58 [-2.32; 1.00] | 1.000 |

Nenhum contraste mostrou vantagem conclusiva do binaural sobre o monaural nesta amostra disponível. Isso não prova que os dois sejam equivalentes ou que ouvir áudio seja ineficaz. Também não é o mesmo contraste que comparar binaural com a avaliação anterior ao experimento: uma melhora em relação à base pode ocorrer sem vantagem específica sobre outro áudio. Não se atribuiu uma avaliação “pré” distinta a cada sessão, pois o protocolo não oferece essa estrutura.

Para EEG, calculou-se primeiro a média por canal, depois por região frontal F3/Fz/F4 e por bloco, com peso igual. Subtraiu-se a base de ruído inicial de cada condição e comparou-se C1-C0. Foram exigidos todos os canais e todos os três blocos. As bandas publicadas (alfa 9-12 Hz) e o denominador diferem do cálculo local; por isso se compara direção dos contrastes, não valores absolutos como se fossem a mesma métrica.

| Grupo | Pares completos | Delta alfa C1 - C0 [IC 95%] | p Holm |
| --- | --- | --- | --- |
| G6 | 7 | -0.0209 [-0.0371; -0.0031] | 0.219 |
| G40 | 5 | -0.0078 [-0.0288; 0.0158] | 0.625 |

Os testes não permaneceram abaixo de 0,05, com cobertura de apenas 7/5 participantes completos. A correlação exploratória entre a diferença de alfa e a diferença de CT também não foi conclusiva: G6, rho=0.364, N=7; G40, rho=-0.700, N=5; ambos p de Holm acima de 0,05. Não foi demonstrado que esse alfa frontal prediga o relaxamento subjetivo.

## 8. Conclusão para o seu experimento

**O que avançou:** a pergunta foi reformulada para relaxamento; os indicadores locais foram recalculados com procedimento comum a EEG externo; houve comparação com 36 pessoas e reanálise independente de uma base com áudio e questionário. Houve suporte externo para alfa como diferença média entre repouso e tarefa, mas não para convertê-lo isoladamente em detector de relaxamento.

**O que os dados locais sustentam:** o áudio coincidiu com mudanças no EEG. A queda de alfa durante os quatro blocos persiste nas sensibilidades testadas e não corresponde à direção de repouso encontrada na EEGMAT. Os indicadores escolhidos não confirmam aumento de relaxamento durante BB. Não se pode concluir o contrário, isto é, aumento de ansiedade ou ausência de relaxamento percebido.

**O que não foi validado:** efeito causal ou vantagem específica do binaural no participante; recuperação fisiológica nas pós-sessões; detector individual de relaxamento; equivalência entre estímulos. Permanecem N=1, condição confundida com tempo/arquivo, ausência de controle acústico local, ausência de escala subjetiva e diferenças de aquisição entre bases.

Uma nova coleta deve medir relaxamento antes/depois, documentar olhos e contato, registrar áudio e marcadores, repetir sessões e alternar aleatoriamente binaural e áudio controle pareado. Os dados atuais continuam úteis como piloto e para escolher critérios de qualidade, não como confirmação de eficácia.

## 9. Reproduzir e inspecionar

Os comandos abaixo usam o ambiente isolado de validação. Os dados externos têm aproximadamente 200 MB e não são incluídos automaticamente em commits. Não é necessário treinar novamente a EEGNet da auditoria anterior.

```bash
python3 -m venv --system-site-packages .venv-validation
.venv-validation/bin/python -m pip install -r requirements-validation.txt
.venv-validation/bin/python -m src.pipelines.fetch_validation_data
.venv-validation/bin/python -m src.pipelines.run_external_validation
.venv-validation/bin/python -m unittest discover -s tests -v
.venv-validation/bin/python -m reports.build_external_report
.venv-validation/bin/python reports/generate_pdf.py --markdown reports/estudo_relaxamento_validacao_externa.md --output reports/Estudo_Relaxamento_Validacao_Externa.pdf
```

Resultados ficam em `results/validation_2026_09_28/`: tabelas por participante, épocas locais, contrastes, cobertura/exclusões, previsões fora do treino, distribuição nula, figura e manifestos de arquivos/versões/parâmetros. O plano e as licenças estão em `reports/plano_validacao_externa.md` e `data/EXTERNAL_DATA.md`. Os dados brutos locais e o relatório da auditoria anterior foram preservados.

## Referências verificadas

- Zyma I et al. Electroencephalograms during Mental Arithmetic Task Performance. Data 2019, 4(1), 14. [DOI 10.3390/data4010014](https://doi.org/10.3390/data4010014). Base [EEGMAT 1.0.0](https://doi.org/10.13026/C2JQ1P), licença Open Data Commons Attribution v1.0.
- Sudre S et al. A new perspective on binaural beats: Investigating the effects of spatially moving sounds on human mental states. PLOS ONE 2024, 19(7), e0306427. [DOI 10.1371/journal.pone.0306427](https://doi.org/10.1371/journal.pone.0306427). Dados [OSF stzfd](https://osf.io/stzfd/), CC BY 4.0. Esta reanálise modifica agregação, seleção e testes; não é endossada pelos autores originais.
