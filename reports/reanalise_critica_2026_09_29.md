# Reanálise crítica do projeto EEG — 29/09/2026

## Pergunta e conclusão

A pergunta científica é se o áudio binaural teve efeito no **relaxamento**. Os registros locais não permitem responder causalmente: há um participante, uma sequência fixa de arquivos, nenhuma condição acústica de controle, nenhum marcador independente de início/fim do áudio e nenhuma medida subjetiva de relaxamento. É possível descrever o EEG e testar se os procedimentos de análise são estáveis. É incorreto transformar classificação de arquivos, potência alfa, rótulo de “Estado” ou índice de defasagem em prova de relaxamento. Também é incorreto afirmar que foi demonstrada ausência de efeito do áudio.

A nova reanálise confirma uma queda descritiva do alfa frontal nos blocos de áudio, sobretudo BB 3 e BB 4. A retenção de épocas muda muito entre condições. A EEGMAT reproduz uma diferença entre repouso e cálculo mental em 36 pessoas, inclusive com outra janela de repouso, mas essa tarefa não calibra um detector de relaxamento. A base de [Sudre et al.](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0306427) oferece um contraste acústico mais pertinente: os questionários C1 (binaural) versus C0 (monaural) não mostram vantagem conclusiva nesta reanálise, e o resultado frontal de EEG depende da cobertura dos eletrodos.

Todos os testes abaixo são **exploratórios e posteriores à primeira análise**. Os valores numéricos e os hashes dos arquivos de entrada estão em [results/critical_2026_09_29](../results/critical_2026_09_29/).

## Adendo: o “Estado” exibido pelo programa anterior

O [script histórico](../legacy/compute_frontal_alpha_asymmetry.py) chama `FAA = ln(alfa F4) - ln(alfa F3)` de “Estado”: FAA positiva recebe “Predomínio Esquerdo (Aproximação/Relaxamento)”; FAA não positiva recebe “Predomínio Direito (Evitação/Alerta)”. Reproduzi a [tabela histórica](../results/tables/tabela_assimetria_alfa_frontal.csv) diretamente dos arquivos brutos. A coluna de janelas usa a mesma fórmula da curva antiga, com janelas de 8 s e passo de 2 s inteiramente contidas em cada condição. Resultados e hashes estão em [legacy_state_audit.csv](../results/critical_2026_09_29/legacy_state_audit.csv) e [manifesto](../results/critical_2026_09_29/legacy_state_audit_manifest.json).

| Condição | FAA da condição | “Estado” histórico | Janelas de 8 s com FAA positiva |
| --- | ---: | --- | ---: |
| Base | +1,956 | Aproximação/Relaxamento | 47/171 (27,5%) |
| BB 1 | +0,465 | Aproximação/Relaxamento | 106/146 (72,6%) |
| BB 2 | +3,221 | Aproximação/Relaxamento | 109/146 (74,7%) |
| BB 3 | +1,620 | Aproximação/Relaxamento | 101/146 (69,2%) |
| BB 4 | +1,748 | Aproximação/Relaxamento | 74/146 (50,7%) |
| Pós 1 | +4,259 | Aproximação/Relaxamento | 20/152 (13,2%) |
| Pós 2 | +1,977 | Aproximação/Relaxamento | 54/156 (34,6%) |

O mesmo rótulo em Base e BB não mostra transição para “relaxado” com o som. A proporção de janelas positivas foi descritivamente maior em BB que na Base, mas janelas com sobreposição não são réplicas e a sequência fixa não identifica causa. FAA da condição e contagem de sinais das janelas são operações matemáticas diferentes; o contraste marcante no Pós 1 demonstra que o rótulo agregado não descreve cada instante. O programa antigo não rejeita épocas por amplitude. Na reanálise com potência relativa e controle de qualidade, a FAA mediana foi negativa nas sete condições, evidenciando sensibilidade ao processamento sem identificar qual etapa produziu a diferença. Nenhum desses sinais é escala validada de calma individual. Em [Harmon-Jones e Sigelman (2001)](https://pubmed.ncbi.nlm.nih.gov/11374750/), maior atividade frontal esquerda também acompanhou raiva com motivação de aproximação; em [Barry et al. (2007)](https://pubmed.ncbi.nlm.nih.gov/17911042/), abrir os olhos alterou o alfa. Essas são razões científicas para não converter a classe em diagnóstico de relaxamento.

Uma sensibilidade adicional manteve o filtro e a fórmula do programa antigo e excluiu janelas de 8 s cujo pico a pico máximo em F3 ou F4 excedia 150 ou 500 µV. Com 150 µV, a Base ficou com 22/138 janelas FAA positivas (15,9%), BB 2 com 4/24 (16,7%) e BB 4 com 4/63 (6,3%); BB 1 reteve apenas 12 de 146 janelas. Com 500 µV, os quatro blocos BB voltaram a ter proporções positivas acima da Base, mas em contagens diferentes. Ao trocar o passo de 2 s por 8 s, sem sobreposição, o corte de 150 µV deixou BB 2 com 0/6 e BB 4 com 1/18 janelas positivas. O limite aplicado a uma janela inteira de 8 s é mais restritivo que o mesmo número aplicado a épocas de 2 s, e nenhum corte aqui identifica artefatos com certeza. A instabilidade descritiva e a seleção desigual estão em [legacy_state_amplitude_sensitivity.csv](../results/critical_2026_09_29/legacy_state_amplitude_sensitivity.csv); não há teste causal ou p-valor adequado para essas janelas de uma pessoa.

## 1. Cronologia verificada nos timestamps

| Intervalo sem registro | Duração |
| --- | ---: |
| Fim da Base → início do arquivo BB | 125,4 s |
| Fim do arquivo BB → início do Pós 1 | **642,0 s (10 min 42 s)** |
| Fim do Pós 1 → início do Pós 2 | 181,2 s |

São intervalos **entre arquivos de EEG**, calculados dos timestamps numéricos, não dos nomes dos arquivos. Os timestamps repetidos são esperados em parte pelo registro de pares de amostras, mas não representam uma observação independente por época. Sem um marcador de desligamento do áudio, não se conhece o intervalo exato entre o áudio e o Pós 1. Portanto, “recuperação imediata” e “consolidação tardia” não são interpretações verificadas. A primeira análise posterior começa mais de dez minutos depois do fim do arquivo de EEG rotulado BB.

## 2. Reavaliação dos testes científicos existentes

| Teste/estimativa | Unidade válida e pergunta | Juízo crítico |
| --- | --- | --- |
| Bootstrap da FAA e inversão de sinais do índice de defasagem local | Épocas do mesmo arquivo | A autocorrelação e a seleção de épocas violam a interpretação simples de réplicas independentes. IC e p descrevem o registro sob pressupostos fortes; não são IC ou p para pessoas ou para efeito do áudio. O índice implementado é uma variante por época, não wPLI espectral debiased. |
| Comparação alfa local por blocos | Blocos de 30 s do mesmo participante | Reduz a falsa precisão de 2 s por época, mas blocos adjacentes ainda dependem uns dos outros. A diferença é descritiva e condicional ao registro; não há teste causal local. |
| Random Forest e EEGNet Base versus BB | Blocos retidos dos dois arquivos | O gap temporal impede vizinhança imediata entre treino e teste, mas condição, arquivo e tempo continuam confundidos. A filtragem da sessão inteira também usa informação global. Acurácia mede discriminação entre **estes arquivos**, não relaxamento nem efeito do binaural. O rótulo único embaralhado por fold é diagnóstico, não distribuição nula inferencial. |
| EEGMAT: Wilcoxon pareado e Holm | 36 pessoas, repouso versus cálculo mental | A unidade “pessoa” é correta; há diferença média de indicadores, sobretudo alfa. O alvo é repouso versus tarefa aritmética. A própria [descrição oficial da EEGMAT](https://physionet.org/content/eegmat/1.0.0/) não oferece escore subjetivo de relaxamento; a limpeza por ICA, eletrodos e tarefa também diferem do piloto. |
| EEGMAT: modelo com uma pessoa fora do treino e permutação por pessoa | 36 pessoas | Melhor separação de treino e teste do que embaralhar épocas, mas 58,3% de acurácia balanceada e p=0,075 não validam um detector individual. A aplicação ao piloto sofre mudança de domínio. |
| Sudre: C1–C0 em questionários e EEG | Pares por pessoa, em grupos separados | É o contraste mais próximo da **especificidade binaural**, pois C0 também contém áudio. Falta de significância não demonstra equivalência. O questionário tem mais pares que o EEG frontal completo; ordem e eletrodo merecem sensibilidade. |
| Sudre: correlação alfa–questionário | 7 pessoas em G6, 5 em G40 no recorte frontal completo | Amostra pequena e selecionada pela disponibilidade de todos os canais; não valida alfa como preditor do relaxamento percebido. |

Um detalhe estatístico merece cuidado: no EEG frontal G6 do relatório anterior, o IC bootstrap da **média** não inclui zero, mas o Wilcoxon pareado não é significativo. Eles resumem aspectos diferentes da distribuição, com apenas sete pares, e nenhum autoriza uma conclusão confirmatória.

### Testes pareados reescritos

Recalculei os quatro contrastes EEGMAT, seis contrastes de questionário Sudre e quatro variantes de EEG frontal com um teste bilateral de **troca do sinal da diferença de cada pessoa**. A estatística é a média das diferenças; as trocas são exatas até 17 pares e usam 99.999 sorteios com semente fixa nas amostras maiores. Holm é aplicado separadamente aos quatro indicadores EEGMAT, aos seis questionários, aos dois contrastes F3/Fz/F4 e aos dois contrastes secundários Fz. A hipótese de intercambialidade dos rótulos continua sendo uma suposição científica: em EEGMAT a ordem repouso→tarefa foi fixa, e em Sudre podem existir efeitos de período/carryover.

| Família | Resultado após Holm da própria família |
| --- | --- |
| EEGMAT, 36 pessoas | Alfa p=0,00032; theta p=0,01467; razão alfa/beta e FAA p≥0,5068. |
| Sudre, seis escalas C1–C0 | Todos p=1,000. |
| Sudre, EEG F3/Fz/F4 | G6 p=0,15625; G40 p=0,56250. |
| Sudre, EEG Fz (sensibilidade) | G6 p=0,31348; G40 p=0,64867. |

A associação alfa–CT de Sudre também foi refeita com **permutação exata das pessoas**: G6, rho=0,364, N=7, p ajustado=0,467; G40, rho=−0,700, N=5, p ajustado=0,467. A estimativa grande em G40 vem de cinco pares e não é evidência de predição. Estes testes foram escolhidos depois da análise inicial e não tornam a auditoria confirmatória. Valores individuais, p brutos e método estão em [participant_signflip.csv](../results/critical_2026_09_29/participant_signflip.csv) e [sudre_exact_correlations.csv](../results/critical_2026_09_29/sudre_exact_correlations.csv).

## 3. Seleção por qualidade no registro local

Reprocessei as sete janelas de cinco minutos diretamente dos arquivos OpenBCI, com F3/F4, filtro 1–30 Hz, sem notch de 20 Hz, limiar frontal de 150 µV e bordas de 2 s removidas. As **772 épocas retidas** e os quatro indicadores reproduzem exatamente a tabela primária de 28/09. Para esta sensibilidade, calculei alfa médio por bloco cronológico de 30 s e mantive blocos com pelo menos 50% ou 80% das épocas candidatas.

| Condição | Blocos ≥50% | Δ alfa versus Base | Blocos ≥80% | Δ alfa versus Base |
| --- | ---: | ---: | ---: | ---: |
| BB 1 | 5 | −2,39 p.p. | **1** | −4,61 p.p. |
| BB 2 | 7 | −2,47 p.p. | **1** | −1,28 p.p. |
| BB 3 | 8 | −1,60 p.p. | 5 | −2,14 p.p. |
| BB 4 | 8 | −4,78 p.p. | 4 | −4,59 p.p. |

No corte de 80%, a Base retém 9 blocos. Os números em BB 1 e BB 2 provêm de **um bloco cada**; não sustentam inferência sobre toda a condição. BB 3 e BB 4 conservam quedas descritivas com mais blocos. O filtro de qualidade também muda quais trechos são analisados, de modo que a persistência da direção não identifica a causa. Os valores da tabela são diferenças entre médias de blocos, diferentes das medianas por época do relatório anterior.

## 4. Outra janela da EEGMAT

Comparei os primeiros e os últimos 60 s do arquivo de repouso com os mesmos 60 s iniciais da tarefa em cada pessoa. O processamento e a regra de qualidade são os mesmos do estudo externo anterior.

| Janela de repouso | Pessoas | Alfa repouso − tarefa, média | Direção positiva | p Wilcoxon exploratório |
| --- | ---: | ---: | ---: | ---: |
| Primeiros 60 s | 36 | +6,22 p.p. | 27/36 | 0,00076 |
| Últimos 60 s | 36 | +5,69 p.p. | 27/36 | 0,00011 |

A diferença repouso–aritmética não depende exclusivamente da janela final escolhida antes. Isso melhora a robustez desse **contraste específico**, sem converter alfa em escala de relaxamento e sem transferir a estimativa para um participante com outro equipamento. Os dois p-valores são análise de sensibilidade posterior, não dois novos testes confirmatórios independentes.

## 5. Outra cobertura frontal em Sudre

O recorte original exigia F3, Fz e F4 em todas as fases de C0 e C1, deixando sete pares em G6 e cinco em G40. A análise alternativa usa somente Fz; ela aumenta a cobertura, mas deixa de ser comparável diretamente à média F3/F4 do piloto.

| Grupo | Canais | Pares completos | Média do contraste alfa C1–C0 | p Wilcoxon exploratório |
| --- | --- | ---: | ---: | ---: |
| G6 | F3/Fz/F4 | 7 | −0,02095 | 0,109 |
| G6 | Fz | 13 | −0,01016 | 0,146 |
| G40 | F3/Fz/F4 | 5 | −0,00782 | 0,625 |
| G40 | Fz | 17 | +0,00233 | 0,712 |

As médias estão em fração da potência relativa (0,01 = 1 ponto percentual). Os grupos de participantes mudam junto com o sensor, então a mudança de direção em G40 **não** pode ser atribuída só a Fz. Nenhum desses contrastes oferece evidência confirmatória de vantagem binaural. A cobertura completa anterior era pequena demais para generalizações sobre EEG frontal.

## 6. Ordem dos áudios em Sudre

A [pesquisa original](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0306427) randomizou a ordem de quatro estímulos por participante. Separei descritivamente a diferença C1–C0 conforme C1 veio antes ou depois de C0. Em G40, o escore CT médio foi **−1,0** quando C1 veio antes (11 pessoas) e **+6,38** quando veio depois (8 pessoas). PA e PT também mudam bastante entre esses subconjuntos; a tabela completa está em [sudre_order_sensitivity.csv](../results/critical_2026_09_29/sudre_order_sensitivity.csv).

Esses subconjuntos são pequenos e a estratificação foi escolhida após observar os dados. A divergência pode refletir acaso, ordem, período, experiência anterior ou composição dos grupos; não é demonstração de interação. Ela mostra por que um novo experimento deve contrabalancear e registrar a ordem, além de analisar cada pessoa como unidade.

## Perspectivas para a próxima coleta

1. Definir **um desfecho primário de relaxamento percebido**, com escala aplicada antes e depois de cada sessão. EEG, frequência cardíaca/variabilidade e comportamento podem ser desfechos secundários, cada um com interpretação própria.
2. Comparar binaural com um áudio controle de duração, volume, expectativa e conteúdo semelhantes; randomizar a ordem e repetir sessões em dias diferentes. A distinção C1–C0 de Sudre é um exemplo de pergunta mais específica que “antes versus durante”.
3. Registrar por hardware/software os instantes de início e fim do áudio dentro do arquivo EEG, além de olhos abertos/fechados, contato dos eletrodos, movimento e intervalos entre sessões. Sem isso, o “Pós” não é uma medida com latência conhecida.
4. Pré-especificar janelas, canais, limiares de artefato e família de hipóteses. Analisar diferenças **por pessoa ou sessão replicada**, manter os trechos rejeitados identificáveis e relatar a sensibilidade à perda de dados.
5. Testar separadamente mecanismos possíveis: resposta auditiva à frequência do batimento, mudança geral de atenção, relaxamento subjetivo e efeitos da apresentação espacial. A [revisão sistemática de Ingendoh et al.](https://pubmed.ncbi.nlm.nih.gov/37205669/) encontrou estudos de entrainment heterogêneos e resultados inconsistentes; uma mudança em uma banda de EEG não resolve as demais perguntas.

## Reproduzir

A execução usa os dados externos descritos em [data/EXTERNAL_DATA.md](../data/EXTERNAL_DATA.md) e preserva os arquivos originais:

```bash
python3 -m unittest discover -s tests -v
.venv-validation/bin/python -m src.pipelines.run_critical_reanalysis
```

O teste automatizado usa dados sintéticos para verificar pareamento, lacunas, filtragem, seleção e invariância. A reanálise real grava tabelas por sessão/bloco/pessoa e um manifesto SHA-256 em [results/critical_2026_09_29](../results/critical_2026_09_29/). A primeira rotina roda no Python do sistema sem depender de MNE; a segunda precisa do ambiente de validação com MNE e dos dados externos baixados. Ambos os comandos foram executados nesta revisão.
