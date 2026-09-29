# Áudio binaural e relaxamento: relatório científico consolidado

**Revisão de 29/09/2026.** Reanálise exploratória de quatro registros OpenBCI atribuídos a um participante, com duas bases públicas independentes. Consulte os [resultados da auditoria](../results/audit_2026_09_26/), da [validação externa](../results/validation_2026_09_28/) e da [revisão crítica](../results/critical_2026_09_29/). O [estudo externo de 28/09](estudo_relaxamento_validacao_externa.md) registra uma etapa anterior.

## Resumo

**O efeito do áudio binaural sobre o relaxamento desta pessoa permanece indeterminado.** O programa antigo exibiu “Aproximação/Relaxamento” nos quatro blocos com áudio, mas também na linha de base, antes do áudio. Reproduzimos esse resultado dos arquivos brutos. Há diferenças de EEG ao longo da sessão; elas não permitem afirmar que o som causou relaxamento, nem que não teve efeito. Condição, ordem temporal e arquivo coincidem, e faltam áudio controle, medida direta de relaxamento e marcador independente do som. Épocas e blocos de uma pessoa não são réplicas independentes.

O alfa frontal relativo foi menor nos quatro blocos rotulados BB que na Base. Ao exigir 80% de épocas válidas por bloco de 30 s, BB 1 e BB 2 ficam com apenas um bloco cada; BB 3 e BB 4 mantêm quedas descritivas. Na EEGMAT, alfa distingue repouso de cálculo mental em 36 pessoas, mas essa base não mede relaxamento percebido. Nos dados de [Sudre et al. (2024)](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0306427), o contraste binaural-monaural não forneceu vantagem conclusiva nesta reanálise. Ausência de significância não demonstra equivalência.

## 1. Pergunta e desenho

A pergunta é **se o binaural altera o relaxamento**. Três resultados precisam ser separados: (1) os registros de EEG têm diferenças descritivas; (2) o programa antigo atribuiu um rótulo de “Estado”; (3) a sensação de relaxamento e a causa de uma eventual mudança não foram medidas de modo que possam ser identificadas. Para avaliar um efeito causal específico, seria preciso comparar com áudio controle em sessões repetidas e medir relaxamento diretamente. Os quatro CSVs locais são rotulados Base, BB, Pós 1 e Pós 2. A atribuição a uma pessoa, a montagem e os nomes das condições vêm da documentação; não há marcador de estímulo no sinal que os confirme. A posição dos canais posteriores P3/O1 e P4/O2 é ambígua.

A [EEGMAT 1.0.0](https://physionet.org/content/eegmat/1.0.0/) contém registros antes e durante cálculo mental: testa descritores entre tarefas, sem validar uma escala de relaxamento. [Sudre et al.](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0306427) disponibilizaram questionários e potências de EEG derivadas para áudio monaural C0 e binaural C1, em grupos de 6 e 40 Hz. É um contraste acústico pertinente, com protocolo, participantes e canais diferentes dos do piloto. O [plano externo](plano_validacao_externa.md) foi escrito antes dos resultados de 28/09; as sensibilidades de 29/09 são posteriores e exploratórias.

## 2. Integridade e cronologia local

{{ACQUISITION_TABLE}}

Os intervalos vêm dos timestamps Unix dos CSVs. Existem **{{GAP_BB_POST_TEXT}} sem EEG entre o fim do arquivo BB e o início do Pós 1**. Como não se sabe quando o áudio foi desligado, “recuperação imediata” não é uma interpretação verificável. Mesmo supondo reprodução contínua durante todo o arquivo BB, a ausência de marcador impede verificar o instante real do estímulo e a ausência de uma condição acústica controle impede separar o efeito binaural de outros efeitos de ouvir um som. Timestamps repetem-se em aproximadamente metade das linhas, sem regressão observada na auditoria; as análises usam a ordem das amostras e 200 Hz nominais. Não há registro independente de olhos abertos/fechados, impedância, intensidade sonora ou movimento ocular.

## 3. Métodos e unidades de análise

A comparação com EEGMAT usa F3/F4, filtro 1-30 Hz por janela, reamostragem a 200 Hz quando necessária, remoção de 2 s nas bordas e épocas sem sobreposição de 2 s. O critério principal rejeita amplitude pico a pico frontal abaixo de 2 µV ou acima de 150 µV. Cada janela local de cinco minutos gera 148 épocas candidatas. Alfa relativo usa potência 8-13 Hz sobre 1-30 Hz; theta relativo usa 4-8 Hz no mesmo denominador. Também se calcularam ln(alfa/beta baixa) e assimetria alfa frontal. **Nenhum deles, isoladamente, é escala validada de relaxamento.**

A auditoria eletrofisiológica anterior usa outro processamento: notch de 60 Hz, opção de notches de 20/80 Hz, filtro 1-45 Hz e 150 épocas candidatas por bloco de cinco minutos. As duas análises não devem ter contagens e valores misturados. O notch de 20 Hz testa sensibilidade a possível interferência sem identificar sua origem. A rejeição por amplitude usa somente os canais frontais.

Nos dados externos, primeiro se agrega por pessoa e condição. Os testes reescritos usam a média das diferenças pareadas e inversões bilaterais do sinal **dentro de cada pessoa**: exatas até 17 pares e por 99.999 sorteios fixos nas amostras maiores. Holm corrige, separadamente, quatro indicadores EEGMAT, seis questionários Sudre, dois contrastes de EEG F3/Fz/F4 e dois de Fz. A interpretação pressupõe intercambialidade dos rótulos. A ordem repouso→tarefa foi fixa na EEGMAT; em Sudre podem existir efeitos de período e da sessão anterior. Esses testes posteriores não são confirmação pré-registrada.

## 4. Estados exibidos pelo programa antigo e EEG local

Em [compute_frontal_alpha_asymmetry.py](../legacy/compute_frontal_alpha_asymmetry.py), o “Estado” é uma regra aplicada à assimetria alfa frontal (FAA) de **potência absoluta**: `FAA = ln(alfa F4) - ln(alfa F3)`. Para `FAA > 0`, o programa escreve **“Predomínio Esquerdo (Aproximação/Relaxamento)”**; nos demais casos, **“Predomínio Direito (Evitação/Alerta)”**. Ele não recebe dados do reprodutor de áudio, não foi calibrado contra relatos de calma e não calcula a mudança de estado causada pelo som. A tabela abaixo reproduz a [saída histórica](../results/tables/tabela_assimetria_alfa_frontal.csv) diretamente dos quatro arquivos brutos; a última coluna usa a mesma regra da curva temporal do programa e conta janelas móveis de 8 s, com passo de 2 s, **inteiramente contidas em cada condição**. Consulte a [auditoria reproduzível](../results/critical_2026_09_29/legacy_state_audit.csv).

{{LEGACY_STATE_TABLE}}

**Leitura do rótulo:** todas as sete condições receberam “Aproximação/Relaxamento”, inclusive a Base. Portanto, a tabela do programa **não registra uma passagem de “não relaxado” para “relaxado” quando o áudio começa**. Ela também não significa que a pessoa permaneceu nesse estado ao longo de toda a gravação: na Base, 47 de 171 janelas móveis tiveram FAA positiva; no Pós 1, apenas 20 de 152, apesar do rótulo positivo da condição inteira. Isso ocorre porque o FAA da condição usa o logaritmo da razão de potências agregadas, enquanto a coluna de janelas conta sinais calculados separadamente. Janelas sobrepostas não são observações independentes, e o programa antigo não rejeita épocas por amplitude. Esses números descrevem a saída do algoritmo, sem teste válido do efeito do áudio.

O nome “Relaxamento” no rótulo excede o que a FAA identifica. Em um [experimento com indução de raiva](https://pubmed.ncbi.nlm.nih.gov/11374750/), maior atividade frontal esquerda acompanhou raiva e comportamento de aproximação; logo, a direção da assimetria não é exclusiva de calma. O [estudo de Barry et al.](https://pubmed.ncbi.nlm.nih.gov/17911042/) mostra que abrir ou fechar os olhos altera a potência alfa; o estado dos olhos não foi registrado neste piloto. Essas fontes sustentam a necessidade de medir diretamente a experiência da pessoa e controlar o protocolo, sem indicar qual fator concreto produziu cada trecho do nosso EEG.

### Reanálise com potência relativa e controle de qualidade

A tabela seguinte usa outro processamento: filtro 1-30 Hz, potência **relativa**, épocas de 2 s e rejeição por amplitude frontal. A FAA é a mediana das épocas aceitas, não a razão entre potências agregadas do programa antigo. A direção dessa FAA relativa foi negativa em todas as condições. Aplicar mecanicamente os nomes de “Estado” do código antigo a essa outra medida produziria rótulos diferentes; **nenhuma das duas versões foi validada como detector individual de relaxamento**. Filtro, normalização, rejeição de épocas e agregação mudaram em conjunto, de modo que esta comparação não isola qual etapa explica a divergência.

{{LOCAL_TABLE}}

São **medianas de épocas retidas**, não medidas de relaxamento. A retenção difere muito entre condições. Para verificar essa seleção, calculei alfa médio por bloco cronológico de 30 s e apliquei dois cortes mínimos. Δ é a diferença entre a média dos blocos da condição e da Base, em pontos percentuais:

{{QUALITY_TABLE}}

Com 80% de retenção, BB 1 e BB 2 têm **um bloco aproveitável cada**; o valor não representa a condição inteira. BB 3 e BB 4 mantêm a direção negativa com cinco e quatro blocos. Blocos adjacentes continuam correlacionados, e a regra de seleção altera os trechos avaliados. Com uma gravação por condição não cabe p local confirmatório. A descrição do Pós também não demonstra retorno ao repouso, pois a latência desde o desligamento do áudio é desconhecida.

![Comparação do piloto, EEGMAT e questionários externos](../results/validation_2026_09_28/comparacao_relaxamento.png)

## 5. Validação dos indicadores na EEGMAT

Diferença positiva significa repouso maior que cálculo mental. Alfa e theta são mostrados em pontos percentuais de potência relativa; razão alfa/beta e FAA ficam em unidades logarítmicas.

{{EEGMAT_TABLE}}

O alfa foi maior no repouso em 27 das 36 pessoas. A direção persiste ao trocar a janela de repouso:

{{WINDOW_TABLE}}

Essas janelas são sensibilidades escolhidas após a primeira análise. A limpeza anterior pelos autores, os eletrodos e a tarefa diferem do piloto. A FAA relativa não apresentou diferença conclusiva entre repouso e cálculo mental nessa base (p Holm = 0,50680), outro motivo para não usar o limiar zero como classificador de calma. A regressão logística com uma pessoa inteira fora do treino por fold obteve **{{MODEL_BALANCED}}% de acurácia balanceada** (ROC-AUC {{MODEL_AUC}}; p de {{MODEL_PERMUTATIONS}} permutações por pessoa = {{MODEL_P}}). Isso não valida um detector individual de relaxamento. As probabilidades transferidas ao piloto não são “percentual de relaxamento”.

## 6. Binaural versus monaural em Sudre

Diferença positiva significa escore C1 maior que C0. Na codificação do artigo, valores maiores de PA, CT e PT indicam mais relaxamento, embora CT/PT mencionem “tensão”. As contagens vêm dos arquivos disponíveis; diferem parcialmente da amostra publicada, cujas exclusões não puderam ser reconstruídas integralmente.

{{QUESTIONNAIRE_TABLE}}

Nenhuma das seis escalas apresentou contraste conclusivo após Holm **nesta comparação direta C1-C0**. O [artigo original de Sudre et al.](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0306427) relata melhora em alguns questionários após C1 em comparação com a avaliação pré-experimental; esse é outro contraste. Comparar C1 com C0 pergunta pela **vantagem específica do binaural sobre outro áudio**, enquanto comparar C1 com a avaliação anterior também pode captar passagem do tempo, expectativa e ouvir um som em geral. Não encontrar vantagem C1-C0 na nossa reanálise não prova ineficácia geral nem equivalência. Tampouco transfere automaticamente os resultados de outro protocolo e outras pessoas para este piloto.

No EEG, subtraiu-se a base de ruído de cada condição e comparou-se a média dos três blocos de estímulo C1-C0. A seleção estrita exige F3, Fz e F4 em todas as fases; a sensibilidade usa só Fz:

{{SENSOR_TABLE}}

A alternativa Fz muda simultaneamente participantes e canais. Em G40 a média muda de sinal, sem permitir atribuir a mudança somente ao sensor. A associação alfa-CT, com permutação exata por pessoa, não foi conclusiva: G6, N={{CORR_G6_N}}, rho={{CORR_G6_RHO}}, p Holm={{CORR_G6_P}}; G40, N={{CORR_G40_N}}, rho={{CORR_G40_RHO}}, p Holm={{CORR_G40_P}}.

A diferença CT C1-C0 conforme a ordem dos sons foi:

{{ORDER_TABLE}}

A ordem foi randomizada no estudo original, mas estes subconjuntos são pequenos e a estratificação é posterior. A divergência descritiva não prova efeito de ordem; justifica medir e contrabalancear a sequência na próxima coleta.

## 7. O que os outros testes locais não demonstram

O índice de defasagem F3-F4 da auditoria é uma variante por época da parte imaginária de Hilbert, **não o wPLI espectral debiased**. Suas inversões de sinal por época pressupõem intercambialidade apesar da dependência temporal:

{{CONNECTIVITY_TABLE}}

O menor p ajustado é {{CONNECTIVITY_MIN_P}}. Esses testes ocorrem **dentro de cada condição**, sem comparação direta BB-Base. O limiar de rejeição de 150 ou 500 µV altera os índices, inclusive o da Base. Diferença entre uma linha significativa e outra não significativa também não testaria a diferença entre condições.

Os classificadores temporais diagnosticam discriminação entre dois arquivos:

{{CLASSIFIER_TABLE}}

Os seis folds compartilham conjuntos de treino; o desvio é entre folds, não intervalo populacional. Cada classe corresponde a um arquivo e uma posição na sequência. A filtragem bidirecional da sessão completa também impede afirmar independência completa do pré-processamento. Acurácia acima do acaso não prova efeito binaural ou relaxamento.

## 8. Conclusão e próxima coleta

**Sustentado:** o rótulo histórico “Aproximação/Relaxamento” apareceu na Base, nos quatro blocos BB e nas duas sessões Pós; a curva de janelas contém ambos os sinais de FAA em cada fase. Os registros locais diferem em alguns descritores de EEG, e o alfa frontal relativo durante BB foi menor que na Base nas configurações examinadas, com perda desigual de épocas. Isso demonstra mudanças **observadas durante a sessão**, compatíveis com várias explicações; não demonstra ausência de resposta ao áudio. Alfa e theta diferenciam repouso de cálculo mental em média na EEGMAT. Nos dados de áudio reanalisados, C1-C0 não forneceu vantagem conclusiva nas escalas ou no EEG frontal após correção.

**Não identificado:** se o áudio, em comparação com outro som ou silêncio, causou alguma dessas mudanças; se a pessoa se sentiu mais calma; a resposta imediatamente após desligar o áudio; equivalência entre estímulos; arrastamento neural; ou benefício clínico. **Ausência de prova de benefício não é prova de efeito nulo.** A inferência causal falha porque há um participante, uma ordem fixa, nenhuma sessão controle comparável e nenhuma medida direta de relaxamento. Além disso, olhos, atenção e qualidade de contato não foram acompanhados. A [revisão sistemática de Ingendoh et al.](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0286023) encontrou resultados heterogêneos sobre arrastamento em EEG; uma banda isolada não determina estado psicológico.

Uma nova coleta deve definir escala de relaxamento **antes e depois de cada sessão** como desfecho primário; comparar binaural com áudio controle pareado; randomizar e contrabalancear; repetir sessões e pessoas; registrar no EEG marcadores do áudio, olhos, contato e intensidade. Janelas, canais e rejeição de artefatos devem ser definidos antes de ver os resultados. A unidade inferencial será a pessoa ou sessão replicada conforme o desenho, não cada época de 2 s.

## 9. Reproduzir e auditar

Os arquivos locais estão em data/raw/. Os dados externos ficam fora do Git; [fontes e licenças](../data/EXTERNAL_DATA.md), hashes e resultados derivados estão documentados. Na raiz do projeto:

~~~bash
python3 -m unittest discover -s tests -v
.venv-validation/bin/python -m src.pipelines.fetch_validation_data
.venv-validation/bin/python -m src.pipelines.run_external_validation
.venv-validation/bin/python -m src.pipelines.run_critical_reanalysis
python3 -m src.pipelines.audit_legacy_states
python3 reports/update_report.py
python3 reports/generate_pdf.py
~~~

A auditoria eletrofisiológica e os classificadores têm comandos próprios no [README](../README.md). O gerador lê CSVs auditados e monta este texto a partir do template, evitando transcrição manual. O PDF usa Chromium. A validação externa não reproduz o protocolo local.

## Fontes científicas e dados

**Dados:** [Zyma et al., Data (2019)](https://doi.org/10.3390/data4010014) e [EEGMAT no PhysioNet](https://physionet.org/content/eegmat/1.0.0/); [Sudre et al., PLOS ONE (2024)](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0306427) e [dados OSF](https://osf.io/stzfd/).

**Métodos e contexto:** [Harmon-Jones e Sigelman, Journal of Personality and Social Psychology (2001)](https://pubmed.ncbi.nlm.nih.gov/11374750/); [Barry et al., Clinical Neurophysiology (2007)](https://pubmed.ncbi.nlm.nih.gov/17911042/); [Ingendoh et al., PLOS ONE (2023)](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0286023); [Vinck et al., NeuroImage (2011)](https://doi.org/10.1016/j.neuroimage.2011.01.055); [Lawhern et al., Journal of Neural Engineering (2018)](https://doi.org/10.1088/1741-2552/aace8c). A EEGNet do projeto é uma variante.
