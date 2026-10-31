# EvtGen: testes de CPV em B0 → J/ψ K_S (SVS × SSD_CP)

Códigos, configurações de decaimento, análises e histórico recuperados de
`/home/darthlucas/programs/test` e `programs/test/cp_study`.

O objetivo é comparar o controle SVS com SSD_CP e investigar os rótulos de
flavour no registro HepMC. Inclui uma variante **experimental** do EvtGen 2.0.0
que preserva os rótulos exportados e registra separadamente o truth tag interno.

## Reproduzir nesta máquina

As dependências existentes em `$HOME/local` são usadas por padrão. Na raiz deste
repositório, com Python 3.11 ou superior, CMake, g++ e make:

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -r requirements.txt
bash scripts/build.sh 2.0.0
bash scripts/run.sh 2.0.0 1000 12345
```

Isso compila os fontes preservados, gera **novos** eventos SVS e SSD_CP e analisa
os arquivos em `runs/2.0.0_1000_12345/`. O número é de tentativas de geração;
eventos aceitos e candidatos de sinal podem ser menores. Para obter estatística
comparável ao estudo histórico, use 50000 em lugar de 1000. As saídas são
protegidas contra sobrescrita: escolha outra semente ou mova a pasta anterior.

```bash
bash scripts/build.sh 2.0.0-label-patch
bash scripts/run.sh 2.0.0-label-patch 50000 12345
bash scripts/build.sh 2.2.3
bash scripts/run.sh 2.2.3 50000 12345
```

A mesma semente permite repetir a configuração em um ambiente compatível;
versões distintas podem consumir números aleatórios em ordens diferentes.
Resultados históricos não têm semente registrada e não são prometidos bit a bit.
Consulte [o tutorial completo](docs/REPRODUZIR.md) para instalar dependências,
rodar mixing=0 e analisar os arquivos históricos.

## Conteúdo

| Pasta | Conteúdo |
|---|---|
| `src/`, `scripts/` | Gerador com modelo, semente, mixing e saídas explícitos; compilação e execução |
| `tests/` | Fontes originais, decay cards, tabelas, scripts e resultados tabulados/gráficos |
| `tests/cp_study/logs/` | Relatórios históricos das investigações |
| `vendor/evtgen-2.0.0/` | Snapshot local de R02-00-00, incluindo instrumentação preexistente |
| `vendor/evtgen-2.2.3/` | Snapshot local de R02-02-03 |
| `vendor/evtgen-2.0.0-label-patch/` | Snapshot experimental completo |
| `patches/` | Diferenças entre o snapshot 2.0.0 local e a variante experimental |
| `archives/` | Pacotes originais de EvtGen, HepMC, Pythia, PHOTOS e TAUOLA |
| `provenance/` | Origens, hashes SHA256 e inventário dos artefatos grandes |

## Cuidados na interpretação

- `tests/generator3.cc` aponta para SSD_CP, apesar de escrever `eventsSVS.hepmc`.
  O arquivo foi preservado como evidência. Não é possível provar o modelo usado
  nos eventos históricos apenas pelo nome ou pela versão atual do fonte.
  O novo gerador seleciona SVS e SSD_CP explicitamente.
- Os decay cards originais foram preservados: o card SVS define `MyJ/psi`, mas
  o decaimento B usa `J/psi`, sem o alias. Portanto não se deve assumir uma
  seleção exclusiva de `J/psi → mu+ mu-` a partir do comentário do card.
- Debug escrito dentro de `decay()` inclui tentativas rejeitadas e não equivale
  à amostra de eventos aceitos.
- Com o patch, o tag interno SSD_CP está no CSV de truth, separado dos PDG IDs.
  Os rótulos do HepMC, sozinhos, não recuperam esse tag. Consulte
  [as conclusões históricas](tests/cp_study/CP_violation_investigation_summary.md)
  e [a descrição das modificações](tests/cp_study/SSD_CP_modifications_explained.md).
- Os números citados nos relatórios são resultados históricos. Não representam
  validação geral do patch ou dos novos eventos gerados.

## Eventos brutos e publicação

Os HepMC e dumps grandes continuam intactos nas pastas originais. Seus caminhos
e tamanhos estão em `provenance/local-artifacts.json`; não entram diretamente no
Git. [O guia de publicação](docs/PUBLICAR.md) explica como preservar **todos**
esses artefatos em partes de uma release, com checksums para recuperação.

Os componentes de terceiros mantêm seus avisos e licenças originais; veja
`vendor/*/COPYING` e os arquivos de licença dentro dos pacotes. O script histórico
`historical/setupEvtGen.sh` é preservado para consulta e usa versões/URLs antigas;
o tutorial atual usa os pacotes locais, sem depender desses endereços.
