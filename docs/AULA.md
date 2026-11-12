# Aula: o que é um ensaio de MassARRAY e por que fizemos cada coisa

Texto para quem está começando. Sem pressa: cada ideia depende da anterior.

## 1. A ideia geral

O MassARRAY (Agena, química **iPLEX**) responde a uma pergunta simples para muitos alvos ao mesmo
tempo: **"nesta posição do DNA, qual base está presente?"**. Ele pesa moléculas num espectrômetro
de massa (MALDI-TOF). Como cada base (A, C, G, T) tem um peso diferente, dá para saber qual base foi
colocada.

Usamos isso de dois jeitos:

- **Detecção**: escolhemos uma posição que só existe no DNA do patógeno. Se o produto aparece,
  o patógeno está na amostra.
- **Resistência / genotipagem (SNP)**: a posição é uma mutação conhecida (ex.: `gyrA S91F` no
  gonococo). A base lida diz se a cepa é normal ou resistente.

## 2. As três etapas da reação (e os três primers)

```
DNA da amostra
   │  1) PCR com primer F e primer R  → copia um pedaço de ~100 pb (o AMPLICON)
   │  2) SAP                          → limpa os nucleotídeos que sobraram
   │  3) extensão com o UEP           → o UEP gruda colado na base-alvo e ganha UMA base
   ▼
espectrômetro: mede a massa do UEP + 1 base → descobre qual base foi
```

- **Primer F e R (PCR)**: 18–25 nt, Tm ~60 °C. Recebem na ponta 5' uma etiqueta de 10 nt,
  `ACGTTGGATG`, que a Agena usa para padronizar (por isso nos arquivos eles aparecem com `ACGTTGGATG...`).
- **Amplicon**: o pedaço copiado. No iPLEX o ideal é **80–120 pb** (fragmentos curtos funcionam
  melhor em DNA de swab, que costuma estar quebrado).
- **UEP (primer de extensão)**: 15–28 nt, termina **exatamente uma base antes** do alvo. Na extensão
  só entram nucleotídeos "terminadores": ele ganha uma base e para. O UEP pode estar na fita de cima
  (`+`, antes do alvo) ou de baixo (`-`, depois do alvo); se estiver na `-`, a base lida é o
  **complemento** (por isso na planilha "A/G" às vezes aparece como "T/C").

## 3. A "região-alvo" e o formato do ADS

O programa da Agena (ADS) recebe, para cada ensaio, um texto assim:

```
CT_ompA    ...ACCAAGCCTTATGATCGACGGAATTCTGTGGGAAGGTTTCGGCGGAGATCCT[T/C]GCGATCCTTGCACC...
```

- ~120 bases de cada lado da base-alvo, que fica **entre colchetes** com as duas opções.
- Letras como **R, Y, S, W, K, M, N** (código **IUPAC**) marcam posições que **variam entre cepas**
  (R = A ou G, Y = C ou T, ...). Elas avisam o ADS: "não ponha primer aqui".

Como achamos essas posições variáveis: baixamos dezenas de sequências do mesmo patógeno do GenBank,
alinhamos cada uma na região e contamos as bases. Se a base minoritária aparece em ≥10% das cepas,
a posição vira IUPAC. A base mais comum forma a sequência "consenso".

## 4. Como escolher a base-alvo de DETECÇÃO

Queremos uma posição que seja, ao mesmo tempo:

1. **Igual em todas as cepas do alvo** (senão perdemos cepas → falso negativo);
2. **Diferente em parentes próximos** (senão detectamos o parente → falso positivo);
3. Em um gene/região com **sequências suficientes** no banco para conferir.

Por isso usamos genes clássicos de cada grupo:

| Grupo | Região usada | Por quê |
|---|---|---|
| HPV | E6/E7 | são os oncogenes; ficam no DNA mesmo quando o vírus se integra ao genoma humano (L1 pode ser perdido) |
| Bactérias | 16S rRNA ou gene específico (ompA, porA, cfb, ureC...) | 16S tem muitas sequências; gene específico quando o 16S é parecido demais com parentes |
| Fungos (Candida) | ITS / 5.8S | é o "código de barras" oficial de fungos |
| Herpesvírus | gene tipo-específico (US4 = gG separa HSV-1 de HSV-2) | |
| Controle humano | íntron de GAPDH / ACTB | íntron não existe nos pseudogenes processados, então não há cópias "falsas" |

## 5. Como escolher a posição de RESISTÊNCIA

Aqui a posição é dada pela literatura (ex.: 23S rRNA **A2058G** e **A2059G** em macrolídeos).
Cuidados:

- **Numeração**: mutações do 23S são citadas na numeração da *E. coli*. Em outra bactéria, a mesma
  posição tem outro número (A2058 da E. coli = A2071 em *M. genitalium*). É preciso alinhar.
- **Códon**: para mutações em proteína (ex.: ERG11 **Y132F**), traduzimos o gene para confirmar que
  o códon 132 é mesmo Tirosina (TAT) e que a troca A→T gera Fenilalanina (TTT).
- **Mutações vizinhas**: quando duas mutações estão lado a lado (2058/2059), cada ensaio marca a
  outra como IUPAC, para o UEP não ficar em cima dela.
- **Cópias múltiplas**: o gonococo tem 4 cópias do 23S; uma amostra pode ter parte delas mutadas →
  pico misto. Isso é esperado.

## 6. Como sabemos que está certo (testes)

- **PCR in silico**: para cada sequência (alvo e não-alvo), procuramos onde F e R grudam
  (ponta 3' — últimas 4 bases — perfeita; até 3 erros no resto) e se formam um produto de 40–800 pb.
- **Sinal do UEP**: se um não-alvo amplifica, verificamos se o UEP também gruda nele. Se não gruda,
  não há sinal (foi o caso do HPV56 x HPV66).
- **Reação cruzada no painel**: cada par de primers contra as sequências de todos os outros alvos.
- **BLAST**: o amplicon é comparado com todo o GenBank. Esperamos que quase todos os acertos sejam do
  próprio alvo. Acertos "estranhos" muitas vezes são **erros de cadastro** no banco (ex.: uma
  "C. albicans" que na verdade é C. dubliniensis) ou **nomes novos** (ex.: *Haemophilus ducreyi* hoje
  também é *Pseudactinobacillus ducreyi*; BVAB2 virou *Amygdalobacter*).

## 7. O que a planilha final contém

| Coluna | Significado |
|---|---|
| Ensaio | nome curto: alvo_gene ou alvo_mutação |
| Patógeno | organismo |
| Tipo | Detecção / Resistência / Controle |
| Gene / região | onde está a base-alvo |
| Primer PCR F / R | 5'→3', já com a etiqueta ACGTTGGATG |
| Primer de extensão | UEP, 5'→3' (nunca degenerado) |
| Base estendida | normal / alternativa (na fita do UEP) |
| Amplicon | tamanho do produto de PCR, sem as etiquetas |
| Status do teste | resumo do BLAST + inclusividade + reação cruzada |
| Observação | limitações e cuidados de interpretação |

## 8. Aplicações

- Painéis sindrômicos (respiratório, IST, gastrointestinal, sepse) com dezenas de alvos por poço.
- Genotipagem de resistência barata em larga escala.
- Vigilância (ex.: tipagem de HPV em estudos populacionais).
- Farmacogenética e SNPs humanos (mesma química, outro tipo de alvo).
