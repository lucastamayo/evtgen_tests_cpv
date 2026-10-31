# SSD_CP: o que foi modificado, por que, e como interpretar os outputs

Este documento explica a cadeia completa do estudo `Bd -> J/psi K_S` com `SVS` e
`SSD_CP`: o que era gerado antes, qual bug apareceu, o que foi modificado, e
qual output deve ser usado agora.

## Resumo curto

Antes do patch, o `SSD_CP` usava os rotulos PDG do proprio B no event record
para duas coisas ao mesmo tempo:

1. representar a cadeia fisica/exportada de mixing no HepMC;
2. carregar, indiretamente, a tag interna de sabor usada pelo modelo `SSD_CP`.

Isso contaminava o record HepMC. Em alguns eventos, cadeias que deveriam ser
interpretadas como oscillation-like apareciam como cadeias de mesmo sabor:

```text
B -> B -> X
```

ou, mais concretamente:

```text
anti-B0 -> anti-B0 -> J/psi K_S
```

Essas cadeias foram classificadas como `weird`.

Depois do patch, as duas ideias foram separadas:

1. os rotulos PDG exportados no HepMC descrevem a cadeia B de forma limpa;
2. a tag interna do `SSD_CP` e o tempo interno usados na amplitude CP sao
   exportados separadamente em `ssdcp_truth_tags.csv`.

Com isso:

```text
original SSD_CP: weird = 1235
patched SSD_CP:  weird = 0
```

Mas isso nao significa que `record_mixing_asymmetry.png` virou a assimetria CP.
Ele continua sendo um diagnostico de mixing do record HepMC. A assimetria CP do
`SSD_CP` deve ser lida de `ssdcp_truth_tag_asymmetry.png`.

## O que a teoria pede para este canal

O decaimento estudado e:

```text
Bd -> J/psi K_S
```

Esse final e um estado de CP, nao um final flavour-specific. Portanto o estado
final nao identifica se o B no decaimento era `B0` ou `anti-B0`. Ambos podem
chegar a `J/psi K_S`.

Para mixing flavour-specific, a curva natural e:

```text
A_mix(t) = (N_unmixed - N_mixed) / (N_unmixed + N_mixed)
         ~= cos(delta_m t)
```

Para CP dependente do tempo em um final de CP, a curva natural e:

```text
A_CP(t) = (N(tag inicial B0 -> f_CP) - N(tag inicial anti-B0 -> f_CP)) / sum
        ~= S sin(delta_m t) - C cos(delta_m t)
```

No estudo atual, a curva de referencia usada para `SSD_CP` foi:

```text
0.700 sin(delta_m t)
delta_m = 0.502 ps^-1
```

Por isso, contar `mixed` versus `unmixed` na arvore responde uma pergunta de
mixing. Para a curva CP, precisa-se de uma tag de sabor inicial.

## O que era gerado antes

Os samples existentes eram analisados principalmente com:

```bash
python3 cp_study/hepmc_cp_asymmetry.py \
  --sample SVS=eventsSVS.hepmc \
  --sample SSD_CP=events_SSDCP2.hepmc \
  --outdir cp_study/results_existing \
  --max-events 50000
```

O script lia os eventos HepMC aceitos, procurava candidatos:

```text
B0/anti-B0 -> J/psi K_S
```

e calculava quantidades a partir da arvore exportada no HepMC.

No `SSD_CP` original, a variavel interna:

```cpp
EvtId other_b;
```

existia dentro de `EvtSSDCP::decay()`, e era usada para calcular a amplitude CP.
Mas ela nao era exportada como uma truth label limpa para os eventos aceitos.

Havia tambem um arquivo de debug:

```text
ssdcp_debug.txt
```

mas ele nao e uma fonte confiavel de eventos aceitos. O metodo `decay()` pode
ser chamado em tentativas internas de accept/reject do EvtGen, entao esse debug
pode conter trials rejeitados.

## O bug das cadeias `B -> B -> X`

A causa principal estava na interacao entre:

1. `EvtSSDCP::decay()`;
2. `EvtCPUtil::OtherIncoherentB()`;
3. a interface Pythia8/EvtGen que exporta a arvore para Pythia/HepMC.

No codigo original, `EvtSSDCP::decay()` chamava:

```cpp
EvtCPUtil::getInstance()->OtherIncoherentB(p, t, other_b, 0.5);
```

O `OtherIncoherentB()` original fazia duas coisas problematicas:

1. se o B ja vinha de uma cadeia de mixing generica, ele reamostrava o lifetime
   do parent;
2. com probabilidade `probB`, ele fazia flip dos IDs do B na propria arvore:

```cpp
p->getParent()->setId(chargeConj(...));
p->setId(chargeConj(...));
```

Esse flip era usado para escolher a tag interna de CP, mas tambem alterava os
PDG labels exportados. A interface Pythia8/EvtGen, por sua vez, reconhece
oscilacao quando encontra uma particula com uma unica filha e ID diferente da
mae. Ao misturar essas convencoes, alguns eventos ficavam como cadeias B-para-B
de mesmo sabor.

O parser passou a classificar as cadeias assim:

```text
unmixed: root B e signal B sao a mesma particula
mixed:   existe cadeia B -> B -> X e root/signal tem sabores opostos
weird:   existe cadeia B -> B -> X e root/signal tem o mesmo sabor
```

No sample original de 50k:

```text
SVS:
  signals = 43362
  unmixed = 35210
  mixed   = 8152
  weird   = 0

SSD_CP original:
  signals = 43173
  unmixed = 35050
  mixed   = 6888
  weird   = 1235
```

O numero `mixed + weird = 8123` ficava perto do numero de mixed do `SVS`
(`8152`), mostrando que `weird` nao era uma nova populacao fisica. Era uma
populacao oscillation-like com label contaminada.

## O que foi modificado

### `cp_study/hepmc_cp_asymmetry.py`

Esse script foi expandido para fazer uma analise mais completa do HepMC aceito.

Ele agora:

1. le eventos HepMC2 em texto;
2. encontra candidatos `B -> J/psi K_S`;
3. sobe a arvore ate o root B;
4. calcula o tempo proprio a partir dos vertices;
5. classifica cada cadeia como `unmixed`, `mixed` ou `weird`;
6. escreve CSVs de candidatos e bins;
7. gera plots de CP tag, initial/final flavour, record mixing, lifetime e
   componentes de mixing.

Principais outputs:

```text
*_signal_events.csv
asymmetry_bins.csv
mixing_counts.csv
mixing_closure.csv
cp_tag_asymmetry.png
initial_final_flavor_asymmetry.png
record_mixing_asymmetry.png
mixing_components.png
proper_time_distribution.png
summary.json
```

O observavel `record_mixing` foi adicionado para contar:

```text
mixed-like = mixed + weird
A_record(t) = (N_unmixed - N_mixed_like) / sum
```

Isso corrigia o counting integrado do sample original, mas nao fazia o
`SSD_CP` original fechar contra `cos(delta_m t)` bin-a-bin. Esse resultado foi
importante: mostrou que os rotulos de cadeia B do `SSD_CP` nao eram uma fonte
confiavel para uma medida de mixing diferencial em tempo.

### `mixing_from_hepmc.py`

Esse script foi usado como cross-check independente para histogramar as classes:

```text
unmixed
mixed
weird
```

Ele ajudou a confirmar que o `SVS` nao tinha `weird`, enquanto o `SSD_CP`
original tinha uma fracao nao desprezivel.

### `cp_study/evtgen_patch/R02-00-00`

Foi criada uma copia local e experimental do EvtGen em:

```text
cp_study/evtgen_patch/R02-00-00
```

Importante: essa copia nao substitui o EvtGen instalado em:

```text
/home/darthlucas/local/evtgen200
```

Ela foi usada apenas para testar uma correcao isolada do `SSD_CP`.

### `EvtCPUtil.cpp` no EvtGen patch

O metodo modificado foi:

```cpp
EvtCPUtil::OtherIncoherentB(EvtParticle* p, double& t, EvtId& otherb, double probB)
```

Antes:

1. reamostrava lifetime em casos ja misturados;
2. fazia flip dos IDs de B na arvore exportada;
3. calculava `otherb` depois desses flips.

Depois:

1. se o B ja tem parent de mixing generico, usa o lifetime ja amostrado pelo
   passo de mixing;
2. escolhe a tag efetiva com probabilidade `probB`;
3. nao reescreve os IDs de B na arvore;
4. calcula:

```cpp
otherb = chargeConj(effectiveTagSource);
```

Motivo:

```text
usar o flip apenas como tag interna do SSD_CP, sem contaminar os labels HepMC.
```

### `EvtSSDCP.cpp` no EvtGen patch

Depois de chamar `OtherIncoherentB()`, o modelo agora salva atributos no
`EvtParticle`:

```cpp
p->setAttribute("ssdcp_truth", 1);
p->setAttribute("ssdcp_parent_pdg", EvtPDL::getStdHep(p->getId()));
p->setAttribute("ssdcp_other_b_pdg", EvtPDL::getStdHep(other_b));
p->setAttribute("ssdcp_flip", flip ? 1 : 0);
p->setAttributeDouble("ssdcp_t_mm", t);
```

Motivo:

```text
exportar a tag interna aceita pelo SSD_CP sem depender dos PDG labels do HepMC.
```

### `cp_study/EvtGenCleanLabels.h`

Foi criada uma copia local da interface Pythia8/EvtGen.

Ela faz duas coisas relevantes:

1. copia os atributos `ssdcp_*` do `EvtParticle` para um vetor C++;
2. mantem os labels de B exportados limpos.

A estrutura adicionada foi:

```cpp
struct SSDCPTruth {
  int pyIndex, parentPdg, otherBPdg, flip;
  double tMm;
};

vector<SSDCPTruth> ssdcpTruth;
```

Quando encontra um `EvtParticle` com `ssdcp_truth`, ela guarda:

```text
py_index
parent_pdg
other_b_pdg
flip
t_mm
```

### `cp_study/generator_ssdcp_label_patch.cc`

Esse gerador roda o sample patched.

Ele:

1. configura Pythia para produzir `bbbar`;
2. desliga o mixing/decay de B do proprio Pythia;
3. usa a interface local `EvtGenCleanLabels.h`;
4. carrega o decay file `Bd_JpsiKs_SSDCP.dec`;
5. escreve o HepMC patched;
6. escreve o CSV de truth tag.

Outputs principais:

```text
cp_study/results_ssdcp_label_patch/events_SSDCP_label_patch.hepmc
cp_study/results_ssdcp_label_patch/ssdcp_truth_tags.csv
```

O CSV tem:

```text
hepmc_event,py_index,parent_pdg,other_b_pdg,flip,t_mm,t_ps
```

### `cp_study/run_ssdcp_label_patch.sh`

Script de reproducao do fluxo patched.

Ele faz:

```bash
LD_LIBRARY_PATH=cp_study/evtgen_patch/build/lib:$LD_LIBRARY_PATH \
  cp_study/generator_ssdcp_label_patch 50000

python3 cp_study/hepmc_cp_asymmetry.py \
  --sample SVS=eventsSVS.hepmc \
  --sample SSD_CP_LABEL_PATCH=cp_study/results_ssdcp_label_patch/events_SSDCP_label_patch.hepmc \
  --outdir cp_study/results_ssdcp_label_patch/analysis_50000 \
  --max-events 50000

python3 cp_study/ssdcp_truth_asymmetry.py \
  cp_study/results_ssdcp_label_patch/ssdcp_truth_tags.csv \
  --outdir cp_study/results_ssdcp_label_patch/analysis_50000
```

### `cp_study/ssdcp_truth_asymmetry.py`

Esse e o script novo para a curva correta de CP truth tag do `SSD_CP`.

Ele le:

```text
cp_study/results_ssdcp_label_patch/ssdcp_truth_tags.csv
```

e define:

```text
truth_signal_initial_pdg = -other_b_pdg
```

Convencao:

```text
other_b_pdg = -511  => signal_initial =  511  => B0
other_b_pdg =  511  => signal_initial = -511  => anti-B0
```

Depois calcula, em bins de `t_ps`:

```text
A_truth(t) =
  (N(truth initial B0) - N(truth initial anti-B0))
  /
  (N(truth initial B0) + N(truth initial anti-B0))
```

e plota contra:

```text
0.700 sin(delta_m t)
```

Outputs:

```text
ssdcp_truth_tag_asymmetry.png
ssdcp_truth_asymmetry_bins.csv
```

## O que esta sendo gerado agora

O diretorio principal patched e:

```text
cp_study/results_ssdcp_label_patch/
```

Arquivos no nivel superior:

```text
events_SSDCP_label_patch.hepmc
ssdcp_truth_tags.csv
generation_5000.log
generation_50000.log
analysis_5000/
analysis_50000/
```

O diretorio mais importante atualmente e:

```text
cp_study/results_ssdcp_label_patch/analysis_50000/
```

Ele contem:

```text
summary.json
mixing_counts.csv
mixing_closure.csv
SVS_signal_events.csv
SSD_CP_LABEL_PATCH_signal_events.csv
asymmetry_bins.csv
cp_tag_asymmetry.png
initial_final_flavor_asymmetry.png
record_mixing_asymmetry.png
mixing_components.png
proper_time_distribution.png
ssdcp_truth_tag_asymmetry.png
ssdcp_truth_asymmetry_bins.csv
```

## Como interpretar cada plot/CSV agora

### `ssdcp_truth_tag_asymmetry.png`

Este e o plot correto para ver a assimetria CP truth do `SSD_CP` patched.

Ele usa:

```text
other_b_pdg
t_ps
```

vindos de `ssdcp_truth_tags.csv`.

Ele deve seguir aproximadamente:

```text
0.700 sin(delta_m t)
```

Para o sample patched de 50k:

```text
truth rows = 43426
RMS residual vs 0.700 sin(delta_m t) ~= 0.032
```

### `ssdcp_truth_asymmetry_bins.csv`

Tabela numerica por bin usada para o plot acima.

Colunas principais:

```text
bin_low_ps
bin_high_ps
bin_center_ps
count
truth_initial_B0
truth_initial_antiB0
asymmetry
error
```

### `record_mixing_asymmetry.png`

Este plot nao e a assimetria CP.

Ele mede o diagnostico de mixing do record HepMC:

```text
A_record(t) = (N_unmixed - N_mixed_like) / sum
mixed_like = mixed + weird
```

Para um sample flavour-specific limpo, isso deveria seguir:

```text
cos(delta_m t)
```

No `SVS`, fecha bem:

```text
SVS record_mixing RMS vs cos = 0.026
```

No `SSD_CP_LABEL_PATCH`, nao fecha bem:

```text
SSD_CP_LABEL_PATCH record_mixing RMS vs cos = 0.561
```

Isso significa que, mesmo com labels limpos, a cadeia HepMC do `SSD_CP` nao deve
ser usada como observavel fisico de mixing diferencial em tempo. Para a
assimetria CP, use a tag interna exportada.

### `initial_final_flavor_asymmetry.png`

Tambem e diagnostico de mixing, nao CP.

Ele calcula:

```text
(N(initial_pdg == final_pdg) - N(initial_pdg != final_pdg)) / sum
```

Isso e util para samples flavour-specific ou para diagnostico do record. Para
`Bd -> J/psi K_S`, nao e a curva CP.

### `cp_tag_asymmetry.png`

No sample original, esse plot tentava usar o root B flavour como tag CP. Depois
do patch, os labels HepMC foram limpos e a tag interna do SSD_CP foi separada.

Portanto:

```text
Nao use cp_tag_asymmetry.png do SSD_CP_LABEL_PATCH como a assimetria CP final.
Use ssdcp_truth_tag_asymmetry.png.
```

### `mixing_counts.csv`

Mostra a contagem integrada de classes:

```text
model,mixing_class,count,fraction
SVS,unmixed,35210,0.8120
SVS,mixed,8152,0.1880
SVS,weird,0,0.0
SSD_CP_LABEL_PATCH,unmixed,35085,0.8079
SSD_CP_LABEL_PATCH,mixed,8341,0.1921
SSD_CP_LABEL_PATCH,weird,0,0.0
```

Esse arquivo confirma que o bug visivel das cadeias `weird` foi removido no
sample patched.

### `mixing_closure.csv`

Compara os diagnosticos de mixing contra:

```text
cos(delta_m t)
```

Resultado atual:

```text
SVS initial_final RMS vs cos              = 0.026
SVS record_mixing RMS vs cos              = 0.026
SSD_CP_LABEL_PATCH initial_final RMS      = 0.561
SSD_CP_LABEL_PATCH record_mixing RMS      = 0.561
```

Interpretacao:

```text
SVS fecha como diagnostico de mixing.
SSD_CP_LABEL_PATCH nao deve ser usado para extrair mixing via B-chain labels.
```

## Antes versus agora

### Antes

O `SSD_CP` original gerava:

```text
events_SSDCP2.hepmc
```

e a unica forma pratica de inferir a tag CP era olhar para os labels de B no
HepMC. Isso era ambiguo porque os mesmos labels tambem estavam sendo usados para
representar mixing.

Consequencias:

```text
initial_B0_fraction ~ 0.616
weird = 1235 em 43173 signals
mixed sozinho subestimava a populacao oscillation-like
mixed + weird corrigia o counting integrado
mas nao corrigia a curva bin-a-bin contra cos(delta_m t)
```

### Agora

O `SSD_CP_LABEL_PATCH` gera:

```text
events_SSDCP_label_patch.hepmc
ssdcp_truth_tags.csv
```

Agora:

```text
HepMC B labels: descrevem a cadeia exportada limpa
ssdcp_truth_tags.csv: guarda a tag interna usada pelo SSD_CP
```

Consequencias:

```text
initial_B0_fraction ~ 0.502
weird = 0 em 43426 signals
truth-tag CP asymmetry fecha contra 0.700 sin(delta_m t)
record_mixing continua sendo apenas diagnostico de record, nao CP
```

## Como isso se conecta com reconstrucao real

Em uma analise reconstruida de `Bd -> J/psi K_S`, voce nao extrairia a
assimetria CP contando `mixed` e `unmixed` do proprio B sinal.

O procedimento conceitual e:

1. reconstruir `B0/anti-B0 -> J/psi K_S`;
2. medir o tempo proprio;
3. obter uma tag de sabor inicial com um flavour tagger;
4. formar a assimetria entre eventos taggeados como `B0` e `anti-B0`;
5. corrigir por mistag, resolucao temporal, aceitancia e assimetrias de detector.

No MC truth deste estudo, a tag de sabor inicial e:

```text
truth_signal_initial_pdg = -other_b_pdg
```

Em dados reais, o analogo seria a saida calibrada de um flavour tagger.

## Regra de uso daqui para frente

Use:

```text
ssdcp_truth_tags.csv
ssdcp_truth_tag_asymmetry.png
ssdcp_truth_asymmetry_bins.csv
```

para estudar a CP truth tag do `SSD_CP`.

Use:

```text
mixing_counts.csv
mixing_components.png
record_mixing_asymmetry.png
mixing_closure.csv
```

como diagnosticos do record HepMC e da limpeza dos labels.

Nao use:

```text
initial_pdg vs final_pdg
mixed vs unmixed
cp_tag_asymmetry.png do sample patched
```

como a medida final da assimetria CP do `SSD_CP`.

## Comandos de reproducao

Rodar o fluxo patched completo:

```bash
cd /home/darthlucas/programs/test
bash cp_study/run_ssdcp_label_patch.sh 50000
```

Rodar apenas a analise truth-tag:

```bash
cd /home/darthlucas/programs/test
python3 cp_study/ssdcp_truth_asymmetry.py \
  cp_study/results_ssdcp_label_patch/ssdcp_truth_tags.csv \
  --outdir cp_study/results_ssdcp_label_patch/analysis_50000
```

Rodar a analise HepMC dos samples existentes:

```bash
cd /home/darthlucas/programs/test
bash cp_study/run_existing_samples.sh \
  --max-events 1000000 \
  --outdir cp_study/results_1M_mixing
```

## Arquivos principais modificados/adicionados

Codigo de analise:

```text
cp_study/hepmc_cp_asymmetry.py
cp_study/ssdcp_truth_asymmetry.py
mixing_from_hepmc.py
cp_study/run_existing_samples.sh
cp_study/run_ssdcp_label_patch.sh
```

Patch experimental de geracao:

```text
cp_study/evtgen_patch/R02-00-00/src/EvtGenBase/EvtCPUtil.cpp
cp_study/evtgen_patch/R02-00-00/src/EvtGenModels/EvtSSDCP.cpp
cp_study/EvtGenCleanLabels.h
cp_study/generator_ssdcp_label_patch.cc
```

Documentacao/logs:

```text
cp_study/README.md
cp_study/logs/2026-05-11_mixing_chain_log.md
cp_study/logs/2026-05-11_weird_chain_investigation.md
cp_study/logs/2026-05-11_complete_modification_report.md
cp_study/SSD_CP_modifications_explained.md
```

## Frase final para nao se perder

O patch nao transforma a cadeia HepMC em uma curva CP. Ele faz algo mais
importante: limpa os labels HepMC e exporta a tag interna correta do `SSD_CP`.
Para `Bd -> J/psi K_S`, a assimetria CP vem da tag de sabor inicial versus
tempo, nao da contagem `mixed/unmixed` da arvore do B sinal.
