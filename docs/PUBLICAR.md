# Publicar os testes e preservar os eventos

Destino informado: https://github.com/lucastamayo/evtgen_tests_cpv

## Código e tutoriais

Depois de autenticar GitHub CLI (`gh auth login` e `gh auth setup-git`) ou configurar
outro método Git autorizado, envie o commit local:

```bash
cd ~/evtgen-cpv-tests
git remote -v
git push -u origin main
```

Confira os arquivos e a visibilidade no repositório de destino.

## Dados grandes como arquivos de release

Nesta preparação, os 24 artefatos (16.571.475.904 bytes originais) já foram
compactados e conferidos em `release-assets/historical-data/`: quatro partes,
6.626.720.034 bytes no total. O inventário e os checksums estão também em
`provenance/release-assets.json`. Esses arquivos estão locais; ainda não foram
enviados ao GitHub. O comando abaixo permite recriar o pacote em outro diretório.

Os eventos e dumps originais ultrapassam o tamanho adequado para Git. Para
preservar os artefatos inventariados em `provenance/local-artifacts.json`, rode
nesta máquina, com espaço livre suficiente para uma cópia integral dos dados:

```bash
python3 scripts/package-local-artifacts.py --output release-assets/historical-data
```

O script lê os originais sem modificá-los, cria um tar comprimido com gzip e dividido em partes de
1900 MiB e calcula SHA256 de cada parte. Também inclui binários e dumps históricos
inventariados. As partes não são arquivos tar independentes: junte-as na ordem.

Com o código já publicado e o GitHub CLI autenticado:

```bash
gh release create historical-data-v1 --repo lucastamayo/evtgen_tests_cpv \
  --title 'Eventos e artefatos históricos de CPV' \
  --notes 'Dados históricos; consulte proveniência e limitações no README.'
gh release upload historical-data-v1 release-assets/historical-data/* \
  --repo lucastamayo/evtgen_tests_cpv
```

Essa publicação só estará completa quando todos os arquivos forem enviados e
conferidos. O inventário no Git não é um backup dos eventos.

## Recuperar os eventos

```bash
mkdir recovered-data
cd recovered-data
gh release download historical-data-v1 --repo lucastamayo/evtgen_tests_cpv
sha256sum -c SHA256SUMS
cat evtgen-local.tar.gz.part-* | tar -xzf -
```

Os arquivos recuperados têm caminhos relativos à antiga pasta `programs/test`.
Os binários preservados dependem das bibliotecas da máquina original; prefira
recompilar os fontes usando o tutorial.
