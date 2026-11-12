# Passo a passo para um painel novo

Tempo típico: 1 dia para ~50 alvos (a maior parte é esperar o NCBI).

## 0. Preparar

```bash
cd massarray-designer
python -m unittest discover -s tests     # tem que dar OK
export NCBI_EMAIL="voce@exemplo.com"
```

## 1. Fazer a lista de alvos

Para cada item da lista da chefe, decida (use a tabela da `AULA.md`, seção 4):

- **Detecção ou resistência?**
- **Referência**: um acesso do GenBank (de preferência RefSeq `NC_...` ou um genoma completo) e
  o trecho (`start`/`end`) com o gene escolhido. Dica: procure o gene em
  https://www.ncbi.nlm.nih.gov/gene e use as coordenadas dele, com ~300 nt de folga de cada lado.
- **Inclusividade**: uma busca do NCBI que traga ~40 sequências do alvo
  (ex.: `Chlamydia trachomatis[orgn] AND ompA[gene]`). Teste a busca no site antes.
- **Exclusividade**: buscas dos **parentes próximos** e do que mais aparece na mesma amostra
  (para vaginal: outras Candida, outros HPV, outras bactérias da flora; para respiratório: os outros
  vírus respiratórios).
- Para **resistência**: a posição exata (1-based dentro do trecho baixado) e a base mutante (`alt`).

Escreva tudo em `config/<painel>.json` (modelo: `config/exemplo_painel.json`;
explicação: `config/COMO_PREENCHER.md`).

## 2. Baixar (`fetch`)

```bash
python -m massarray fetch config/<painel>.json
```
Confira `work/<painel>/cache/`: um FASTA por conjunto. Abra alguns. Se uma busca trouxe pouca coisa
(ou coisa errada), ajuste o termo no config, apague o arquivo do cache e rode de novo.

## 3. Escolher a base e montar a região (`markers`)

```bash
python -m massarray markers config/<painel>.json
```
Para cada alvo aparece: base-alvo, conservação (tem que ser ≥0,95), quantas amostras alinharam e
quantas posições variáveis. Se **"SEM POSIÇÃO BOA"**: o gene é parecido demais com os parentes ou
as amostras não alinharam → troque de gene/trecho.

## 4. Primers (`primers`)

```bash
python -m massarray primers config/<painel>.json
```
O campo `nivel` diz o quanto foi preciso relaxar as regras. "padrão" é o ideal. Se chegou em
"sem filtro de exclusão", **desconfie**: provavelmente vai ter reação cruzada → troque a região.

## 5. Testes (`check`)

```bash
python -m massarray check config/<painel>.json
```
Leia linha por linha:
- `inclusividade 38/40` — bom. Abaixo de ~85%: ver se as que falharam são sequências parciais
  (comum) ou cepas reais diferentes.
- `não-alvo com sinal 0/N` — tem que ser 0. Se não for: redesenhe (mude `target.scan` para outro
  trecho) ou explique na observação.
- `painel {...}` — reação cruzada com outro ensaio do painel. Mesma regra.
- `posicionamento OK` — obrigatório.

## 6. BLAST (`blast`)

```bash
python -m massarray blast config/<painel>.json
```
Um amplicon por vez. Pode interromper (Ctrl+C) e rodar de novo: continua de onde parou.
Para cada resultado com "outro mais parecido" ≥95% de identidade, investigue:
erro de cadastro? nome novo? parente real? (ver `LICOES_APRENDIDAS.md`).

## 7. Exportar (`export`)

```bash
python -m massarray export config/<painel>.json
```
Em `work/<painel>/saida/`:
- `<painel>_regioes.fasta` — regiões com IUPAC e `[X/Y]` (é o que a chefe chama de "fasta das regiões")
- `<painel>_ADS_input_IUPAC.txt` — para carregar direto no ADS
- `<painel>_ADS_input_simples.txt` — igual, sem IUPAC (só se o ADS reclamar)
- `<painel>_amplicons.fasta` — os produtos de PCR
- `<painel>_lista_final.csv` — planilha (abrir no Excel/Sheets)

## 8. Conferência final

Siga `docs/CHECKLIST.md` antes de enviar.
