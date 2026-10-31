# Tutorial de reprodução

## 1. Ambiente

Linux, Python ≥3.11, CMake, make, g++, gfortran. Em Debian/Ubuntu, os pacotes
básicos são `build-essential cmake gfortran python3-venv`.
As versões dos pacotes preservados são: Pythia 8.306, HepMC 2.06.10,
PHOTOS 3.64, TAUOLA 1.1.8-LHC e EvtGen 2.0.0/2.2.3. O pacote HepMC 2.06.11
é mantido apenas como arquivo adicional; o estudo usa 2.06.10 e interface HepMC2.

Para usar as dependências existentes nesta máquina:

```bash
cd ~/evtgen-cpv-tests
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -r requirements.txt
```

Para instalar em outra máquina a partir dos pacotes incluídos:

```bash
export DEPS_PREFIX="$HOME/evtgen-cpv-deps"
bash scripts/build-dependencies.sh
```

O prefixo deve estar vazio. O script se recusa a substituir instalações existentes.
Essa instalação completa de dependências é uma receita separada da validação
com as bibliotecas já instaladas. Os pacotes antigos podem exigir ajustes com
compiladores mais recentes; consulte os READMEs dentro dos arquivos originais.
Não execute `historical/setupEvtGen.sh` para seguir este tutorial.

Os caminhos individuais podem ser configurados antes de compilar e executar:

```bash
export PYTHIA_PREFIX=/caminho/pythia8306
export HEPMC_PREFIX=/caminho/hepmc2.06.10
export PHOTOS_PREFIX=/caminho/photos364
export TAUOLA_PREFIX=/caminho/tauola118
```

## 2. Comparar SVS e SSD_CP

Execute da raiz do repositório, com o ambiente Python ativado:

```bash
bash scripts/build.sh 2.0.0
bash scripts/run.sh 2.0.0 1000 12345
```

Parâmetros de `run.sh`: variante, tentativas de eventos e semente Pythia
(inteiro de 1 a 900000000). O gerador usa pp a 13.6 TeV, produção bbbar,
`pTHatMin=8 GeV`, MPI desligado, mixing Pythia desligado e mixing EvtGen=1.
O diretório `runs/2.0.0_1000_12345/` contém os HepMC, logs, CSV e figuras.
No gerador sem patch o CSV de truth contém apenas o cabeçalho: os tags internos
não são exportados nessa variante.

A leitura dos eventos usa o parser HepMC2 incluído, sem precisar de pyhepmc.
Figuras de baixa estatística não validam uma curva de CPV. Use 50000 ou mais
para comparação com os relatórios históricos, conforme espaço e tempo disponíveis.

## 3. Variante com labels e truth tag

```bash
bash scripts/build.sh 2.0.0-label-patch
bash scripts/run.sh 2.0.0-label-patch 50000 12345
```

Além da análise HepMC, grava `SSD_CP_truth.csv` e `truth_analysis/`.
A variante inclui as alterações em `EvtCPUtil.cpp`, `EvtSSDCP.cpp` e a interface
Pythia modificada `tests/cp_study/EvtGenCleanLabels.h`. A atribuição de CPV deve
usar o tag interno apropriado; PDG IDs do registro não substituem esse tag.

O snapshot baseline local já tinha instrumentação `ssdcp_debug.txt` em SSD_CP;
por isso é identificado como snapshot local, não como release upstream intocada.
Os tarballs originais estão separados em `archives/`.

## 4. Comparar com EvtGen 2.2.3

```bash
bash scripts/build.sh 2.2.3
bash scripts/run.sh 2.2.3 50000 12345
```

Cada variante tem uma pasta de build e bibliotecas explícitas. O script de
execução coloca a biblioteca da variante primeiro no `LD_LIBRARY_PATH`.
O patch é específico do snapshot 2.0.0 e não é aplicado automaticamente ao 2.2.3.

## 5. Controle mixing=0

Depois de compilar 2.0.0, o gerador pode ser chamado diretamente:

```bash
source scripts/env.sh
export LD_LIBRARY_PATH="$REPO_DIR/build/2.0.0/evtgen/lib:$LD_LIBRARY_PATH"
mkdir -p runs/mixing0
cd runs/mixing0
"$REPO_DIR/build/2.0.0/generate" SSD_CP 1000 sample.hepmc truth.csv 12345 0 "$REPO_DIR/tests"
python3 "$REPO_DIR/tests/cp_study/hepmc_cp_asymmetry.py" \
  --sample SSD_CP=sample.hepmc --outdir analysis
```

A opção controla o argumento de mixing da interface EvtGen. Não constitui,
por si só, prova de que toda lógica interna do modelo CP está desligada.

## 6. Analisar os dados históricos

Nesta máquina, os dados originais continuam em `~/programs/test`:

```bash
python3 tests/cp_study/hepmc_cp_asymmetry.py \
  --sample "SVS=$HOME/programs/test/eventsSVS.hepmc" \
  --sample "SSD_CP=$HOME/programs/test/events_SSDCP2.hepmc" \
  --max-events 1000 --outdir runs/historical-smoke
```

A identificação SVS acima repete o nome histórico; a inconsistência do fonte
`generator3.cc` impede certificá-la apenas pelo nome. Os wrappers originais
em `tests/cp_study/` foram preservados e esperam caminhos/execuções locais da época.
Para novos estudos, use os scripts de `scripts/`.

O FWLite em `tests/cp_study/fwlite_cp_asymmetry.py` requer ambiente CMSSW e
arquivos EDM compatíveis. Não é necessário para a análise HepMC acima e não é
instalado por `requirements.txt`.

## 7. Conferir integridade

`provenance/copied-files.json` registra origem, destino, tamanho e SHA256 de cada
arquivo recuperado. Os novos scripts e documentos não pertencem a esse conjunto
de arquivos históricos. Para conferir:

```bash
python3 scripts/verify-provenance.py
```
