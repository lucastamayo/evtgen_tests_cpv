# Validação local

Executada com as bibliotecas já instaladas em `~/local`. Detalhes do ambiente
e contagens estão em `provenance/validation.json`.

| Variante | Compilação dos fontes copiados | Geração SVS + SSD_CP | Análise |
|---|---|---|---|
| EvtGen 2.0.0 local | Passou | 1000 tentativas/modelo, seed 12345 | Passou |
| EvtGen 2.0.0 label patch | Passou | 1000 tentativas/modelo, seed 12345 | HepMC e truth passaram |
| EvtGen 2.2.3 local | Passou | 1000 tentativas/modelo, seed 12345 | Passou |

No teste com patch, os 836 candidatos SSD_CP têm 836 linhas de truth, em 672
eventos. A classificação de cadeias retornou zero casos `weird`; o baseline
2.0.0 retornou 33. Esses números verificam a execução e a exportação dos registros
em baixa estatística, sem estabelecer validação física de precisão.

Duas gerações independentes de 100 tentativas SVS, semente 67890 e EvtGen 2.0.0,
produziram HepMC idênticos byte a byte neste ambiente. O controle SSD_CP com
mixing=0 também gerou e analisou 100 tentativas.

Os hashes dos 3657 arquivos históricos copiados foram conferidos. Os scripts
shell passaram por `bash -n`; os utilitários Python passaram por compilação de
sintaxe. A receita de reconstrução completa das dependências a partir dos
pacotes originais ainda não foi executada em ambiente limpo. O FWLite/CMSSW
histórico não foi executado nesta validação.

O arquivo compactado dos 24 artefatos históricos foi lido integralmente:
a integridade gzip e todos os nomes e tamanhos internos foram conferidos contra
o inventário. As quatro partes têm SHA256 registrado em
`provenance/release-assets.json`.
