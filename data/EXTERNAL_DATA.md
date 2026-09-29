# Bases externas da análise de 28/09/2026

Os arquivos baixados ficam em `data/external/`, ignorados pelo Git. Nenhum arquivo do participante local foi enviado para serviços externos. `results/validation_2026_09_28/download_manifest.json` contém URLs e SHA-256 de cada arquivo.

## EEGMAT 1.0.0

Zyma I, Tukaev S, Seleznov I, Kiyono K, Popov A, Chernykh M, Shpenkov O. Electroencephalograms during Mental Arithmetic Task Performance. Data 2019, 4(1), 14. DOI: https://doi.org/10.3390/data4010014.

Base: https://physionet.org/content/eegmat/1.0.0/ ; DOI: https://doi.org/10.13026/C2JQ1P.

Licença: Open Data Commons Attribution License v1.0. Foram baixados os 72 EDF de 36 participantes e metadados. Downloads pelo espelho público AWS indicado pelo PhysioNet, verificados contra `SHA256SUMS.txt` do repositório oficial. Os sufixos `_1` e `_2` significam repouso e cálculo mental, respectivamente; não são rótulos clínicos de relaxamento e estresse.

As durações e o filtro descritos resumidamente na página contêm ambiguidades; a análise registra as durações reais dos cabeçalhos EDF e aplica seu próprio filtro documentado. A limpeza prévia pelos autores não pode ser desfeita.

## Sudre et al. 2024 - OSF stzfd

Sudre S, Kronland-Martinet R, Petit L, Rozé J, Ystad S, Aramaki M. A new perspective on binaural beats: Investigating the effects of spatially moving sounds on human mental states. PLOS ONE 19(7), e0306427. DOI: https://doi.org/10.1371/journal.pone.0306427.

Dados: https://osf.io/stzfd/ ; licença do projeto consultada na API OSF: CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/). Foram baixados dois CSVs de potências espectrais por época, dois de questionários e o dicionário. Não são sinais brutos.

C0 é controle monaural e C1 é binaural. Usam-se apenas os contrastes planejados, separadamente para G6 e G40. O artigo informa 19 participantes por grupo; os CSVs incluem 22 questionários e 20 IDs de EEG em G6, e 19 em G40. Não se inventou uma lista de exclusões para reproduzir N=19. A reanálise inclui pares disponíveis segundo critérios explícitos. Em G40, 1.092 duplicatas exatas de EEG foram removidas na memória; todas pertencem à condição C2, fora do contraste principal. Os arquivos baixados permanecem intactos.

Os dados derivados e resultados desta reanálise devem manter estas atribuições. O código do projeto mantém sua licença MIT; isso não altera as licenças dos dados externos.
