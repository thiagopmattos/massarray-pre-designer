# massarray-designer

Ferramenta em Python para **desenhar e conferir ensaios de MassARRAY (Agena, química iPLEX)**
que detectam patógenos e mutações de resistência em amostras clínicas.

Você entra com uma lista de alvos (por exemplo "HPV16", "Chlamydia trachomatis",
"gyrA S91F do gonococo"). A ferramenta:

1. baixa as sequências do NCBI;
2. escolhe a melhor base para ler;
3. monta a região no formato da Agena;
4. desenha os primers;
5. testa tudo no computador e no BLAST;
6. entrega os arquivos que vão para o equipamento.

Ela foi construída a partir de dois painéis feitos no laboratório. Os dois estão em `examples/`:

| Painel | Ensaios | Pasta |
|---|---|---|
| Respiratório (swab) | 67 desenhados + 26 do artigo Zhao 2021 (só na planilha) | `examples/painel_respiratorio_67/` |
| Microbiota vaginal / IST | 70: 58 de detecção, 10 de resistência e 2 controles humanos | `examples/painel_vaginal_IST_70/` |

> Se você nunca trabalhou com MassARRAY, leia antes [`docs/AULA.md`](docs/AULA.md) (conceitos, 15 min).
> Para usar a ferramenta passo a passo, veja [`docs/PASSO_A_PASSO.md`](docs/PASSO_A_PASSO.md).
> Antes de enviar qualquer resultado, siga [`docs/CHECKLIST.md`](docs/CHECKLIST.md).
> O relato completo do que foi feito nos dois painéis, com glossário, está em [`docs/GUIA_COMPLETO.md`](docs/GUIA_COMPLETO.md).

---

## Sumário

1. [Instalação](#1-instalação)
2. [Uso rápido](#2-uso-rápido)
3. [O que sai no final](#3-o-que-sai-no-final)
4. [Como funciona, etapa por etapa (e por quê)](#4-como-funciona-etapa-por-etapa-e-por-quê)
5. [As regras de desenho e de onde elas vêm](#5-as-regras-de-desenho-e-de-onde-elas-vêm)
6. [Os testes e como interpretar](#6-os-testes-e-como-interpretar)
7. [Formatos da Agena: entrada do ADS e tabela TSV](#7-formatos-da-agena-entrada-do-ads-e-tabela-tsv)
8. [Configuração do painel](#8-configuração-do-painel)
9. [Comandos](#9-comandos)
10. [Estrutura do código](#10-estrutura-do-código)
11. [Testes automáticos](#11-testes-automáticos)
12. [Limitações conhecidas](#12-limitações-conhecidas)
13. [Histórico de mudanças](#13-histórico-de-mudanças)

---

## 1. Instalação

Só precisa de **Python 3.9 ou mais novo**. A ferramenta não usa nenhuma biblioteca externa:
tudo foi escrito com a biblioteca padrão, para rodar em qualquer computador do laboratório.

```bash
git clone https://github.com/<seu-usuario>/massarray-designer.git
cd massarray-designer
python -m unittest discover -s tests        # tem que terminar com "OK" (roda sem internet)
```

Para baixar dados do NCBI, informe seu e-mail. É uma exigência do NCBI.
A API key é opcional e deixa o download cerca de 3× mais rápido.

```bash
export NCBI_EMAIL="voce@exemplo.com"
export NCBI_API_KEY="..."            # opcional: https://www.ncbi.nlm.nih.gov/account/settings/
```
No Windows (PowerShell), use `$env:NCBI_EMAIL="voce@exemplo.com"`.

---

## 2. Uso rápido

```bash
cp config/exemplo_painel.json config/meu_painel.json    # edite com seus alvos (seção 8)
python -m massarray all config/meu_painel.json           # faz tudo
```

Na primeira vez, prefira rodar etapa por etapa para olhar cada resultado (seção 9).
Os resultados ficam em `work/<nome_do_painel>/saida/`.

Para um painel que já está pronto (só gerar os arquivos da Agena ou conferir):

```bash
python -m massarray verify examples/painel_vaginal_IST_70/markers.json examples/painel_vaginal_IST_70/primers.json
python -m massarray agena  examples/painel_vaginal_IST_70/markers.json examples/painel_vaginal_IST_70/primers.json saida/ painel_vaginal
```

---

## 3. O que sai no final

Na pasta `work/<painel>/saida/`:

| Arquivo | Para que serve | Quem usa |
|---|---|---|
| `<painel>_ADS_entrada.txt` | Entrada do **Assay Design Suite** (ADS) da Agena: região de cada ensaio com os **primers de PCR em MAIÚSCULO** e o resto em minúsculo | o ADS (programa da Agena) |
| `<painel>_ensaios_agena.tsv` | **Tabela de ensaios** no layout da Agena (colunas WELL, SNP_ID, primers, UEP, massas...), com poços sugeridos | o software do equipamento / conferência |
| `<painel>_lista_final.csv` | Planilha com primers, status dos testes e observações | pessoas (chefe, laboratório) |
| `<painel>_regioes.fasta` | As regiões com o alvo `[X/Y]` e as posições variáveis em IUPAC | visualização e registro |
| `<painel>_amplicons.fasta` | O trecho copiado por cada par de primers | conferência e BLAST |

Arquivos intermediários (ficam em `work/<painel>/`, servem para continuar de onde parou):

| Arquivo | Conteúdo |
|---|---|
| `cache/*.fasta` | tudo o que foi baixado do NCBI (o trabalho fica reprodutível) |
| `markers.json` | as regiões: consenso (`plain`), versão com IUPAC (`iu`), posição do alvo (`F`) |
| `primers.json` | os primers, posições, Tm, massa e bases estendidas |
| `check.json` | resultado dos testes in silico |
| `blast.json` | resumo do BLAST de cada amplicon |

---

## 4. Como funciona, etapa por etapa (e por quê)

```
 config.json ──► fetch ──► markers ──► primers ──► check ──► blast ──► export
                 (NCBI)    (alvo +     (F, R,      (PCR no   (banco    (ADS, TSV,
                           região)     UEP)        computador) inteiro)  planilha)
```

### 4.1 `fetch`: baixar as sequências

**O que faz:** para cada alvo, baixa do NCBI três conjuntos de sequências:

- **referência**: um trecho do genoma em volta do gene (ex.: `NC_000117.1`, posições 779500–780900);
- **inclusividade**: até ~40 sequências de cepas diferentes do próprio alvo;
- **exclusividade**: sequências de parentes próximos e de organismos que aparecem na mesma amostra.

**Como faz:** usa as E-utilities do NCBI (`esearch` para buscar e `efetch` para baixar), respeitando
o limite de pedidos por segundo. Tudo é salvo em `cache/`.

**Por quê:**

- O alvo só é bom se for **igual em todas as cepas**. Por isso precisamos de várias cepas: com uma só,
  não dá para saber o que varia.
- O alvo também precisa ser **diferente dos parentes**. Por isso precisamos dos parentes.
- Registros marcados como **UNVERIFIED** são descartados. No painel vaginal, vários "HPV68"
  desse tipo eram, na verdade, HPV39.
- Sequências com mais de 2% de `N` (bases desconhecidas) também são descartadas.
- O cache deixa o trabalho **reprodutível**: dá para refazer as contas com os mesmos dados e
  mostrar exatamente o que foi usado.

### 4.2 `markers`: escolher a base-alvo e montar a região

**O que faz:**

1. **Recorta** cada sequência baixada em volta da região da referência e coloca todas na mesma fita
   (`align.orient_and_trim`). Isso usa três "âncoras" de 30 nt: início, meio e fim.
2. **Escolhe a posição-alvo** (só para detecção; em resistência a posição vem do config).
   O programa percorre a região posição por posição (`marker.scan_positions`) e olha um trecho de
   31 nt em volta de cada uma:
   - **found**: em quantas cepas do alvo esse trecho aparece com no máximo 5 diferenças.
     Precisa de ≥85% das cepas.
   - **base**: em quantas a base do meio é igual à da referência. Precisa de ≥95%.
   - **mm**: média de diferenças. Precisa de ≤1,2.
   - **close**: fração dos não-alvos em que o trecho aparece com ≤3 diferenças.
     Quanto menor, menos reação cruzada.

   Ganha a posição com menos semelhança com não-alvos e, em caso de empate, a mais conservada.
3. **Monta a região** de 241 nt (120 de cada lado do alvo) com `marker.build_marker`:
   - cada cepa é alinhada na região (`align.align_window`): primeiro o núcleo de 31 nt em volta do
     alvo, depois pedaços de 20 nt para cada lado, aceitando deslocamentos de até ±6 nt
     (pequenas inserções e deleções);
   - em cada posição contamos as bases: a mais comum forma o **consenso** (`plain`);
   - se a segunda base aparece em **≥10% das cepas**, a posição vira **código IUPAC** (`iu`),
     por exemplo R = A ou G, Y = C ou T.

**Por quê:**

- **241 nt** é o tamanho que o ADS espera (cerca de 100–150 nt de cada lado). Dá espaço para
  posicionar os primers sem pegar regiões longe demais.
- O **IUPAC** avisa o desenho de primers (o nosso e o do ADS) para **fugir** das posições que variam.
  Um primer em cima de uma posição variável deixaria de grudar em parte das cepas, o que causa
  falso negativo.
- O limite de **10%** evita que um erro de sequenciamento isolado vire "variação", sem deixar de
  marcar variantes reais e frequentes.

### 4.3 `primers`: desenhar F, R e UEP

**O que faz:** `design.design` testa todas as combinações possíveis dentro da região e fica com a de
melhor pontuação:

- **primers de PCR (F e R)**: 18–25 nt, Tm 55–64 °C (ideal 60 °C), GC 30–70%, sem 5 bases iguais
  seguidas, no máximo 3 G/C nas 5 últimas bases;
- **UEP (primer de extensão)**: 15–28 nt, encostado na base-alvo (antes dela na fita de cima, ou
  depois dela na fita de baixo), Tm de extensão ≥47 °C, **nunca com base ambígua**;
- **amplicon** de 70–130 pb (ideal ~100);
- nenhum dímero de ponta 3' entre os três primers;
- pontuação = distância da Tm ideal + diferença de Tm entre F e R + distância do amplicon ideal +
  Tm do UEP.

Se não acha nada com as regras padrão, `design_with_fallbacks` relaxa as regras aos poucos e anota
qual foi o nível usado (`nivel`): k-mer de 9, amplicon de 60–140, primer degenerado permitido e,
por último, sem filtro de exclusão. Quanto mais relaxado, mais atenção nos testes.

As regras e a origem de cada uma estão na seção 5.

### 4.4 `check`: testes no computador

Para cada ensaio, uma PCR in silico (`insilico.amplifies`) contra:

- as cepas do alvo, que dá a **inclusividade**;
- os não-alvos, que dá a **especificidade**. Se algum amplifica, o programa vê se o UEP também gruda
  e qual base ele estenderia (`extension_base`); só é problema se gerar sinal;
- as cepas de **todos os outros alvos do painel**, que dá a **reação cruzada no painel**.

Também roda `verify_placement`, que confere a geometria de cada ensaio (seção 6).

### 4.5 `blast`: comparação com o banco inteiro

Cada amplicon é comparado com todo o GenBank (`core_nt`) usando megablast, **um amplicon por pedido**.
O resumo diz quantos acertos são do organismo esperado e qual é o "outro" organismo mais parecido
(identidade e cobertura). Registros sem espécie (não cultivados, vetores, metagenomas) são ignorados.
O resultado é salvo a cada ensaio, então dá para interromper e continuar depois.

**Por que um por vez:** no painel vaginal, um pedido com 10 amplicons juntos (incluindo HPV16 e
HPV18, que têm milhares de registros) ficou mais de 1 hora na fila do NCBI. Um por um, cada pedido
levou cerca de 30 segundos.

### 4.6 `export`: arquivos finais

Gera a entrada do ADS, a tabela TSV com poços sugeridos, a planilha CSV e os FASTA (seções 3 e 7).

---

## 5. As regras de desenho e de onde elas vêm

| Regra | Valor | Por quê |
|---|---|---|
| Etiqueta (tag) nos primers de PCR | `ACGTTGGATG` (10 nt) na ponta 5' | Padrão do iPLEX: os primers de PCR que sobram ficam com massa fora da janela de leitura e não atrapalham o espectro |
| Tamanho do amplicon | 70–130 pb (ideal 80–120) | DNA de swab costuma estar fragmentado; produtos curtos amplificam melhor e de forma mais parecida em multiplex |
| Tm dos primers de PCR | 55–64 °C, F e R parecidas | Todos os ensaios de um poço rodam com o mesmo programa de temperatura |
| ≤3 G/C nas 5 últimas bases | – | Uma ponta 3' "grudenta" demais liga em lugares errados e forma dímeros |
| **Base ambígua só na metade 5'** do primer de PCR (no máximo 2) | regra do laboratório (Hema) | A enzima começa a copiar pela **ponta 3'**: uma base errada ali derruba a PCR. Perto do 5' o primer tolera a diferença |
| **UEP nunca degenerado** | só A/C/G/T | No iPLEX o UEP é pesado: um UEP com base ambígua é uma mistura de oligos de massas diferentes e gera picos duplicados |
| UEP encostado no alvo | termina 1 base antes | Na extensão entra só **uma** base (terminador), a base-alvo |
| Filtro de exclusão (k-mer) | as 10 últimas bases de cada primer de PCR **não podem existir** nos não-alvos | Evita primers "universais" (comuns em 16S, 18S e ITS) que amplificam várias espécies. Resolveu 12 ensaios do painel vaginal |
| Primers só em posições sem variação | IUPAC proibido nos primers, salvo o caso acima | Garante que o primer grude em todas as cepas |
| Controle humano em íntron | GAPDH íntron 1, ACTB íntron 3 | Pseudogenes processados (cópias sem íntron) não geram produto falso |
| Resistência no 23S | numeração da *E. coli* convertida para a do organismo | A literatura cita A2058/A2059 na numeração da *E. coli*; em *M. genitalium* é A2071/A2072 |
| Mutações vizinhas | a vizinha vira IUPAC (`target.mask`) | O UEP de A2058G não pode ficar em cima de A2059 e vice-versa |

---

## 6. Os testes e como interpretar

### PCR in silico (`insilico.amplifies`)

Um par "amplifica" uma sequência se:

- as **4 últimas bases** de cada primer batem perfeitamente;
- há **no máximo 3 diferenças** no resto de cada primer;
- os primers ficam em fitas opostas, virados um para o outro;
- o produto tem 40–800 pb.

Essa regra é propositalmente **generosa**: é melhor acusar uma possível reação cruzada que não
acontece do que deixar passar uma que acontece.

### Sinal do UEP (`insilico.extension_base`)

Se um não-alvo amplifica, o programa vê se o UEP gruda nesse produto (com no máximo 1 diferença) e
qual base ele estenderia.

- Amplifica, mas o UEP **não** gruda: não há sinal. Exemplo: os primers do HPV56 amplificam 3 de 24
  HPV66, mas o UEP não gruda neles.
- Amplifica e o UEP **gruda**: falso positivo. Foi o caso do primeiro desenho de *C. parapsilosis*,
  que pegava *C. tropicalis* e foi redesenhado.

### Conferência de posicionamento (`insilico.verify_placement`)

Para cada ensaio, confere:

- F e R batem na posição indicada;
- o UEP está encostado na base-alvo;
- a base estendida é a correta (complemento quando o UEP está na fita de baixo);
- nenhum primer de PCR se sobrepõe ao UEP;
- o tamanho do amplicon confere;
- o UEP não é degenerado;
- as bases ambíguas dos primers de PCR estão só na metade 5'.

Também avisa quando há posição variável embaixo do UEP. O comando `verify` roda isso em um
painel pronto.

### BLAST: como ler

| Resultado | Interpretação |
|---|---|
| "100/100 acertos são X" | ótimo |
| "236/246 são X; raros com outro nome" | quase sempre **erro de cadastro** no GenBank: 1–3 registros "de outra espécie" idênticos no meio de centenas corretos |
| "outro mais parecido: Y 95%" | ~5 diferenças em 100 pb; em geral o UEP ou os primers não grudam, mas confira |
| "outro: Y 100%" com muitos acertos de Y | reação cruzada real: redesenhe ou anote (ex.: 23S do gonococo é igual ao de outras *Neisseria*) |
| nome diferente, mesma espécie | **nome novo**: *Haemophilus ducreyi* = *Pseudactinobacillus ducreyi*; BVAB2 = *Amygdalobacter*; *C. glabrata* = *Nakaseomyces glabratus*; *C. krusei* = *Pichia kudriavzevii*; *Atopobium vaginae* = *Fannyhessea vaginae*. Coloque os dois nomes no campo `expect` |

---

## 7. Formatos da Agena: entrada do ADS e tabela TSV

### 7.1 Entrada do ADS: `<painel>_ADS_entrada.txt`

Uma linha por ensaio, sem cabeçalho: `nome<TAB>sequência`.

```
HPV6_E6E7	tgcaagacg...gtgtttgcaagaatgcaCTGACCACTGCAGAGATTTATTCAtatgcatataaacasctaaaggtcctgtttcga[G/A]gcggctatccatatgcagcctgcgcGTGCTGCCTAGAATTTCATGGAaaaathaacc...
```

- **MAIÚSCULO = primers de PCR** (F na fita de cima e o complemento reverso do R). O ADS **mantém**
  esses primers.
- **minúsculo = todo o resto**, inclusive o lugar onde fica o UEP. Assim o ADS **desenha o UEP**
  escolhendo massas que não se sobrepõem dentro de cada poço.
- `[X/Y]` = base-alvo / base alternativa.
- Letras minúsculas como `r`, `y`, `s`, `h` = posições que variam entre cepas.

**Por que esse formato:** orientação do Hema. O ADS só monta o multiplex direito se puder escolher os
UEPs; os primers de PCR já foram testados e devem ser mantidos. `agena.check_ads_sequence` garante
que as maiúsculas são exatamente os primers F e R.

### 7.2 Tabela de ensaios: `<painel>_ensaios_agena.tsv`

Arquivo separado por tabulação, com o layout de colunas da saída do ADS que o software do
equipamento importa:

| Coluna | Conteúdo |
|---|---|
| `WELL` | poço sugerido (W1, W2...) |
| `TERM` | `iPLEX` |
| `SNP_ID` | nome do ensaio |
| `2nd-PCRP` | primer **reverse**, com tag |
| `1st-PCRP` | primer **forward**, com tag |
| `AMP_LEN` | tamanho do amplicon (sem as tags) |
| `UP_CONF`, `MP_CONF` | em branco (são notas de confiança calculadas pelo ADS) |
| `Tm(NN)` | Tm do UEP (vizinhos mais próximos) |
| `PcGC` | % de GC do UEP |
| `PWARN` | avisos: primer degenerado, amplicon >120, UEP fora da janela de massa, UEP a definir no ADS |
| `UEP_DIR` | `F` = UEP na fita da sequência de entrada; `R` = fita oposta |
| `UEP_MASS` | massa do UEP (Da) |
| `UEP_SEQ` | sequência do UEP |
| `EXT1_CALL` / `EXT1_MASS` / `EXT1_SEQ` | alelo normal: base, massa do produto, UEP + base |
| `EXT2_CALL` / `EXT2_MASS` / `EXT2_SEQ` | alelo alternativo (mutante) |

**Como as massas são calculadas:**

- UEP = soma das bases (A 313,21; C 289,18; G 329,21; T 304,2) − 61,96 Da, para um oligo de DNA sem
  fosfato 5';
- produto = UEP + terminador iPLEX (A 271,2; C 247,2; G 287,2; T 327,1 Da).

Os terminadores da química iPLEX são modificados justamente para aumentar a diferença de massa entre
as bases (a menor diferença, entre A e G, é de 16 Da).

**Poços sugeridos (`agena.plex`):** agrupamento simples, ensaio por ensaio, em ordem de massa. Um
ensaio entra no primeiro poço onde:

- todas as massas (UEP e os dois produtos) ficam a **≥30 Da** de todas as massas já presentes,
  porque o aparelho precisa separar os picos;
- não há outro ensaio da **mesma região do genoma**. Ensaios do mesmo lugar competem e geram
  produtos extras: ERG11 Y132F e K143R, MG 23S A2058G e A2059G e os três rpoB da TB ficam em
  poços diferentes;
- não há **dímero** de ponta 3' entre os primers;
- há vaga (máximo de 36 ensaios por poço).

O ADS faz isso melhor: ele também alonga UEPs com bases extras na ponta 5' para espalhar as massas.
Nos exemplos, o agrupamento simples deu 16–17 poços; o ADS costuma conseguir bem menos. Use os poços
sugeridos como ponto de partida ou para conferência.

> **Importante:** confira o cabeçalho do TSV com um arquivo de exemplo do seu equipamento ou da sua
> versão do software. Se for diferente, mude só a lista `COLUMNS` em `massarray/agena.py`.

Ensaios cujo UEP seria degenerado aparecem com as colunas do UEP em branco e
`PWARN = UEP a definir no ADS`.

---

## 8. Configuração do painel

Arquivo JSON com um bloco por ensaio. Modelo em `config/exemplo_painel.json` e explicação de cada
campo em [`config/COMO_PREENCHER.md`](config/COMO_PREENCHER.md). Resumo:

```json
{
  "panel": "meu_painel",
  "markers": [
    {
      "id": "CT_ompA",
      "organism": "Chlamydia trachomatis",
      "gene": "ompA",
      "type": "Detecção",
      "expect": "Chlamydia trachomatis",
      "reference": {"acc": "NC_000117.1", "start": 779500, "end": 780900},
      "inclusivity": {"term": "Chlamydia trachomatis[orgn] AND ompA[gene]", "max": 40},
      "exclusivity": [{"term": "Chlamydia muridarum[orgn] AND ompA", "max": 40}],
      "target": {"scan": [200, 1200]}
    },
    {
      "id": "NG_gyrA_S91F",
      "organism": "Neisseria gonorrhoeae",
      "type": "Resistência",
      "reference": {"acc": "NC_002946.2", "start": 620300, "end": 621500},
      "inclusivity": {"term": "Neisseria gonorrhoeae[orgn] AND gyrA", "max": 40},
      "target": {"position": 619, "alt": "T"}
    }
  ]
}
```

- **Detecção**: dê a região (`reference`) e deixe o programa escolher a base (`target.scan`).
- **Resistência**: dê a posição exata (`target.position`, 1-based dentro do trecho baixado) e a base
  mutante (`target.alt`). Use `target.codon` para conferir o aminoácido e `target.mask` para
  mutações vizinhas.

---

## 9. Comandos

| Comando | O que faz |
|---|---|
| `python -m massarray fetch   config.json` | baixa referência, cepas do alvo e não-alvos (cache) |
| `python -m massarray markers config.json` | escolhe a base-alvo e monta as regiões → `markers.json` |
| `python -m massarray primers config.json` | desenha F, R e UEP → `primers.json` |
| `python -m massarray check   config.json` | PCR in silico, sinal do UEP, reação cruzada no painel, posicionamento → `check.json` |
| `python -m massarray blast   config.json` | BLAST de cada amplicon (lento; continua de onde parou) → `blast.json` |
| `python -m massarray export  config.json` | entrada do ADS, TSV da Agena, planilha e FASTAs |
| `python -m massarray all     config.json` | todas as etapas em sequência |
| `python -m massarray verify  markers.json primers.json` | confere um painel pronto (offline) |
| `python -m massarray agena   markers.json primers.json pasta/ nome` | gera entrada do ADS e TSV de um painel pronto (offline) |

---

## 10. Estrutura do código

```
massarray/
  seq.py       complemento reverso, IUPAC, Tm (SantaLucia + sal), massa, regras de primer,
               regra "ambiguidade só no 5'", tradução de códons
  ncbi.py      busca/download no NCBI com limite de pedidos e cache
  align.py     busca aproximada rápida (princípio da casa dos pombos) e alinhamento na janela
  marker.py    varredura de posições (scan_positions) e região de consenso com IUPAC (build_marker)
  design.py    desenho dos 3 primers, filtro de exclusão por k-mer e níveis de relaxamento
  insilico.py  PCR in silico, base estendida pelo UEP, avaliação e conferência de posicionamento
  blast.py     BLAST pela URL API do NCBI e resumo
  agena.py     entrada do ADS (maiúsculo/minúsculo), tabela TSV, massas e poços sugeridos
  export.py    escreve todos os arquivos finais
  cli.py       linha de comando
config/        modelo de configuração + explicação dos campos
docs/          aula, passo a passo, lições aprendidas, checklist
examples/      os dois painéis feitos (entradas, primers e arquivos entregues)
tests/         testes automáticos
```

Detalhes de implementação que vale saber:

- **Busca aproximada rápida (`align.find_approx`)**: para achar um trecho com até N diferenças, o
  trecho é dividido em N+1 pedaços. Pelo menos um deles aparece sem nenhuma diferença, e esse pedaço
  é encontrado rapidamente com `str.find`; só depois a sequência inteira é conferida. Assim dá para
  buscar em genomas inteiros em Python puro.
- **Tm**: modelo de vizinhos mais próximos (SantaLucia 1998) com correção de sal. Padrão para PCR:
  primer 100 nM, Na⁺ 50 mM, Mg²⁺ 2 mM, dNTP 0,5 mM. O UEP usa as condições da extensão.
- **Bases ambíguas** contam como "combina" nas buscas (`seq.compatible`): R combina com A e com G.

---

## 11. Testes automáticos

```bash
python -m unittest discover -s tests -v
```

São 18 testes e todos rodam sem internet:

- busca aproximada comparada com a força bruta;
- alinhamento e recorte;
- região de consenso;
- escolha de posição fugindo de um trecho compartilhado com um "parente";
- desenho e avaliação completos com dados sintéticos;
- UEP nunca degenerado e regra da ambiguidade no 5';
- formato da entrada do ADS (maiúsculas = exatamente F e R) e da tabela TSV;
- poços sugeridos: massas ≥30 Da e ensaios da mesma região separados;
- conferência dos 70 ensaios do painel vaginal;
- o fluxo inteiro pela linha de comando, com um cache falso.

---

## 12. Limitações conhecidas

- **Primers provisórios**: o ADS ainda monta o multiplex e pode ajustar primers por massa e dímeros
  entre ensaios.
- **Testes in silico** dependem do GenBank: são uma amostra (até ~40 cepas por alvo), não a população.
  Alvos com poucas sequências públicas (*Ureaplasma*, BVAB2, *Megasphaera* tipo 1, HPV raros) precisam
  de validação em bancada.
- **Poços sugeridos** são uma aproximação: o ADS otimiza melhor.
- **Cabeçalho do TSV**: confira com um arquivo do seu equipamento (seção 7.2).
- **Detecção = presença/ausência**: o MassARRAY não quantifica.
- **Resistência em genes compartilhados** (23S de *Neisseria*/*Mycoplasma*): só interpretar com a
  detecção do organismo positiva.
- **BLAST e NCBI** dependem de internet e podem ficar lentos. O resultado é salvo a cada ensaio.

---

## 13. Histórico de mudanças

- **v1.1**:
  - regra "ambiguidade só na metade 5'" no desenho e na conferência;
  - entrada do ADS com primers de PCR em MAIÚSCULO e o resto em minúsculo (o ADS desenha o UEP);
  - tabela TSV no layout da Agena, com massas e poços sugeridos;
  - comando `agena`;
  - exemplos atualizados: 5 primers do vaginal e 5 ensaios do respiratório revisados;
  - os 3 UEPs degenerados do respiratório ficam "a definir no ADS".
- **v1.0**: pipeline completo (fetch → markers → primers → check → blast → export) e os dois painéis.

Licença: MIT (veja `LICENSE`).
