# Projeto EEG: piloto de estimulação binaural

[Relatório científico consolidado](reports/relatorio_cientifico.md) · [PDF](reports/Relatorio_Cientifico_EEG_Binaural.pdf) · [Reanálise crítica detalhada](reports/reanalise_critica_2026_09_29.md).

**Conclusão atual:** o programa histórico exibiu “Aproximação/Relaxamento” antes, durante e depois do trecho rotulado como áudio binaural. A auditoria reproduziu esses rótulos dos arquivos brutos; eles não demonstram que o áudio causou relaxamento nem que não teve efeito. Há mudanças descritivas no EEG, mas o efeito específico do som e a sensação de calma da pessoa permanecem indeterminados. O relatório principal integra o piloto, 36 participantes da EEGMAT e os contrastes binaural-monaural publicados por Sudre et al.

## Reanálise crítica - 29/09/2026

[Relatório crítico](reports/reanalise_critica_2026_09_29.md) · [Resultados e manifesto](results/critical_2026_09_29/) · [Testes reescritos](tests/).

Os timestamps mostram **642,0 s sem registro** entre o fim do arquivo BB e o início do Pós 1; não há marcador que informe quando o áudio foi desligado. A nova análise recalcula diretamente dos arquivos brutos as 772 épocas locais da configuração principal e examina retenção por bloco, escolha da janela EEGMAT, cobertura dos eletrodos e ordem dos áudios na base Sudre. A [auditoria dos estados históricos](results/critical_2026_09_29/legacy_state_audit.csv) reproduz a tabela antiga e verifica também as janelas de 8 s da curva: o rótulo por condição não significa um estado contínuo. Para reproduzir: `python3 -m unittest discover -s tests -v`, `python3 -m src.pipelines.audit_legacy_states` e `.venv-validation/bin/python -m src.pipelines.run_critical_reanalysis` após baixar os dados externos.

## Novo estudo de relaxamento com validação externa - 28/09/2026

[Relatório completo](reports/estudo_relaxamento_validacao_externa.md) · [PDF](reports/Estudo_Relaxamento_Validacao_Externa.pdf) · [Plano](reports/plano_validacao_externa.md) · [Fontes e licenças](data/EXTERNAL_DATA.md).

Reanálise de 72 EDF da EEGMAT (36 participantes) e dos CSVs públicos do estudo de Sudre et al. (2024), que incluem EEG e questionários para áudio binaural/monaural. O alfa frontal foi maior no repouso externo, mas caiu durante os quatro blocos de áudio no piloto, nas quatro sensibilidades testadas. Isso não confirma relaxamento pelos indicadores escolhidos e não prova piora subjetiva.

O classificador externo, validado deixando pessoas inteiras fora do treino, obteve 58,33% de acurácia balanceada (p de permutação=0,075). Não foi validado como detector individual de relaxamento. A comparação binaural-monaural dos questionários não mostrou vantagem conclusiva após Holm. A cobertura de EEG frontal da segunda base foi pequena e há diferenças de amostra em relação ao artigo, explicitadas no relatório.

Os dados externos ficam fora do Git; resultados derivados, hashes e scripts estão em `results/validation_2026_09_28/` e `src/pipelines/`. Para reproduzir, use `requirements-validation.txt` e os comandos do novo relatório. A auditoria de 26/09 abaixo permanece como histórico metodológico.

Análise exploratória de quatro gravações OpenBCI Ganglion (quatro canais, 200 Hz), atribuídas a um participante. A revisão de 26/09/2026 corrige o pipeline e substitui as conclusões anteriores. Condição, arquivo e ordem temporal estão confundidos; os dados não demonstram causalidade nem benefício clínico.

## Resultado da auditoria

- As contagens do relatório estavam misturadas: a tabela principal usava 150 µV e janelas de cinco minutos, enquanto o texto descrevia 500 µV e registros completos. A base correta é 144/150 épocas nessa tabela.
- Com 150 µV e notch, nenhum dos sete testes de defasagem permanece significativo após Holm (menor p ajustado: 0,0909). Com 500 µV, os resultados mudam, inclusive para a base. São testes dentro de condições, não comparações de BB versus base.
- Os blocos de treino/teste agora preservam índices originais e um intervalo nominal de 120 s entre bordas. O treino/teste anterior não garantia esse intervalo.
- Nova EEGNet, 45 épocas por fold, seis folds: acurácia balanceada 64,90% ± 21,02 pontos percentuais; controle diagnóstico 50,88% ± 4,93. Random Forest: 69,85% ± 14,26. Desvios são amostrais entre folds, não ICs populacionais.
- A filtragem ainda opera na sessão completa; não se afirma “zero vazamento”. O notch de 20 Hz é uma hipótese de tratamento de interferência, não prova de aliasing elétrico. Métricas posteriores requerem cautela: a rejeição frontal não controla todos os canais.

O [relatório científico consolidado](reports/relatorio_cientifico.md) substitui o texto principal desta auditoria. A cópia PDF da raiz é sincronizada pelo gerador. As tabelas e figuras antigas foram preservadas para comparação, mas suas legendas e conclusões não representam a análise revisada.

## Execução

Python 3.11 ou superior. Use um ambiente com as dependências de `requirements.txt` e uma instalação PyTorch compatível com CPU ou sua GPU. O ambiente `.venv` encontrado no projeto estava incompleto; a revisão foi executada com `python3` do host. As versões efetivamente usadas estão nos manifestos dos resultados.

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python3 -m unittest discover -s tests -v
python3 -m src.pipelines.run_electrophysiology
python3 -m src.pipelines.run_anti_bias_audit
python3 -m src.pipelines.run_eegnet_gpu --epochs 45 --device auto
python3 reports/update_report.py
python3 reports/generate_pdf.py
```

O PDF precisa de Chromium/Chrome local. A geração funciona offline com figuras incorporadas e falha explicitamente quando falta uma figura ou o navegador não produz o arquivo. Não há dependência de fórmulas LaTeX não renderizadas.

## Organização

- `data/raw/`: registros originais, preservados. `data/README.md`: metadados e ressalvas.
- `src/preprocessing.py`, `asymmetry.py`, `connectivity.py`: funções compartilhadas com validação de entrada e sementes fixas.
- `src/validation.py`: partições temporais e ajuste de Holm.
- `src/features.py`, `decoding.py`, `models/eegnet.py`: descritores, dados e variante EEGNet.
- `src/pipelines/`: eletrofisiologia com sensibilidade, Random Forest temporal/transferência e EEGNet temporal.
- `results/audit_2026_09_26/`: resultados revisados, figura e manifestos SHA-256 de dados/código, parâmetros e versões.
- `results/validation_2026_09_28/` e `results/critical_2026_09_29/`: validação externa, sensibilidades, testes pareados e manifestos.
- `results/tables/`, `results/figures/`, `legacy/`: material histórico; a tabela dos estados FAA foi reproduzida nesta revisão, os demais itens não foram revalidados integralmente.
- `reports/relatorio_cientifico_template.md` e `reports/update_report.py`: montam o relatório principal a partir dos CSVs; `generate_pdf.py`: renderiza e sincroniza os PDFs.
- `tests/`: regressões de leitura, rejeição, cronologia, estatística, filtragem e dimensões da rede.

## Limites e continuidade

Os próximos passos científicos exigem marcadores de áudio e qualidade de contato, sessões replicadas e contrabalanceadas com controle acústico, hipóteses pré-especificadas e inferência que respeite a dependência temporal. Uma única permutação por fold é apenas diagnóstico. Previsões nas pós-sessões não medem “retorno ao repouso”. O teste A/A histórico não foi repetido e não prova ausência de deriva.

As referências verificadas e o detalhamento dos resultados estão no relatório. Licença MIT para o código: ver `LICENSE`.
