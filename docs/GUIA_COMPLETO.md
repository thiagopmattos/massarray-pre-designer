# Guia completo: o que fizemos nos painéis de MassARRAY

*Thiago Mattos — novembro de 2026*

Este guia explica, do zero, tudo o que foi feito nos dois painéis (respiratório e vaginal/IST): os conceitos, cada decisão, os testes, os problemas que apareceram e como foram resolvidos. Está na ordem em que as coisas acontecem. No fim tem um glossário e respostas para as perguntas mais prováveis da chefe.

---

## 1. Visão geral

**Objetivo:** desenhar ensaios para o MassARRAY (Agena, química iPLEX) que detectam patógenos e mutações de resistência em amostras clínicas (swab).

| | Painel respiratório | Painel vaginal / IST |
|---|---|---|
| Ensaios desenhados por nós | 67 | 70 |
| Ensaios de outra fonte | 26 do artigo Zhao 2021 | — (os HPV da chefe ficaram como referência) |
| Total na planilha | 93 | 70 |
| Detecção | vírus, bactérias, fungos respiratórios | 28 HPV, ISTs, Candida, bactérias da vaginose |
| Resistência | influenza (oseltamivir, baloxavir), TB (rifampicina, isoniazida), Aspergillus (azóis), macrolídeos, carbapenemases | MG, NG, TP (azitromicina/quinolonas) e C. albicans (azóis) |
| Controles | — | GAPDH e ACTB humanos |

**O que foi entregue:** planilha com os primers e testes, regiões em FASTA, arquivo de entrada do programa da Agena (ADS) e amplicons.

---

## 2. Conceitos básicos

**DNA e bases.** O DNA é uma fita com 4 letras (A, C, G, T). A fita tem duas pontas: 5' (início) e 3' (fim). A outra fita é o "complemento reverso": A↔T, C↔G, lida ao contrário.

**PCR.** Reação que copia um pedaço do DNA milhões de vezes. Precisa de dois **primers** (pedacinhos de 18–25 letras) que grudam nas pontas do trecho: o **forward** (F) numa fita, o **reverse** (R) na outra. A enzima começa a copiar sempre a partir da **ponta 3' do primer** — por isso a ponta 3' precisa grudar perfeitamente.

**Amplicon.** O pedaço copiado pela PCR (do início do F ao fim do R). No MassARRAY o ideal é 80–120 pares de base (pb).

**Tm (temperatura de melting).** Temperatura em que metade do primer se solta do DNA. Os primers de uma PCR devem ter Tm parecidas (~60 °C). Calculamos pelo método dos "vizinhos mais próximos" (SantaLucia), que considera a sequência exata.

**SNP / posição-alvo.** Uma posição do DNA onde a base pode variar. No nosso caso, a base que o MassARRAY vai ler.

**Código IUPAC.** Letras que representam "mais de uma base possível" numa posição: R = A ou G; Y = C ou T; S = C ou G; W = A ou T; K = G ou T; M = A ou C; N = qualquer uma. Usamos para marcar posições que variam entre cepas.

**Primer degenerado.** Primer com letra IUPAC. Na prática é uma mistura de primers (ex.: um com A e outro com G naquela posição), para pegar cepas diferentes.

---

## 3. Como funciona o MassARRAY (iPLEX)

O aparelho é um espectrômetro de massa: ele **pesa** moléculas. Como cada base tem um peso diferente, dá para saber qual base foi adicionada.

```
1) PCR (primers F e R)      → copia o amplicon (~100 pb)
2) SAP                      → limpa as bases que sobraram
3) Extensão (primer UEP)    → o UEP gruda colado na base-alvo e ganha UMA base
4) Espectrômetro            → mede o peso do UEP + 1 base → descobre a base
```

**Os três primers de cada ensaio:**
- **F e R (PCR):** recebem na ponta 5' a etiqueta `ACGTTGGATG` (10 nt). Ela padroniza o peso dos primers que sobram e evita que atrapalhem a leitura. Por isso na planilha os primers começam com ACGTTGGATG.
- **UEP (primer de extensão):** termina exatamente **uma base antes** do alvo. Pode estar na fita de cima (antes do alvo) ou na de baixo (depois do alvo, lido ao contrário). Na fita de baixo, a base lida é o complemento — por isso "A/G" às vezes aparece como "T/C" na planilha.

**Regras do UEP:** nunca pode ter base ambígua (seria uma mistura de pesos → picos duplicados). Cada UEP de um mesmo poço precisa ter peso diferente dos outros — quem organiza isso é o programa ADS.

**Multiplex.** Vários ensaios (até ~40) rodam juntos no mesmo poço. O ADS decide quais ensaios vão juntos.

**Detecção x genotipagem.** Para detectar um patógeno, escolhemos uma posição que só existe no DNA dele: se aparece produto, ele está na amostra. Para resistência, a posição é a mutação: a base lida diz se a cepa é normal ou resistente.

---

## 4. Escolha dos genes (onde procurar)

| Grupo | Gene/região | Por quê |
|---|---|---|
| HPV | E6/E7 | São os oncogenes. Ficam no DNA mesmo quando o vírus se integra ao genoma humano (o L1, alvo de muitos testes, pode ser perdido na integração). |
| Bactérias da flora | 16S rRNA | Tem milhares de sequências no banco para comparar. |
| Bactérias específicas | ompA (CT), porA (NG), MgPa (MG), ureC (Ureaplasma), polA (TP), cfb (S. agalactiae), lytA (pneumococo)… | Quando o 16S é parecido demais com parentes, usamos um gene que só aquela espécie tem. |
| Candida/fungos | ITS / 5.8S | É o "código de barras" oficial dos fungos. |
| Herpesvírus | US4/gG (HSV), ORF29 (VZV), UL83 (CMV) | gG é diferente entre HSV-1 e HSV-2, então separa os dois. |
| Vírus respiratórios | genes conservados de cada vírus (N, M, NS, HA…) | |
| Controle humano | íntron do GAPDH e do ACTB | Os pseudogenes (cópias "falsas" desses genes espalhadas no genoma) não têm íntron, então não há produto falso. |

---

## 5. Os dados: de onde vieram as sequências

Todas as sequências vieram do **GenBank (NCBI)**, baixadas automaticamente.

Para cada alvo usamos dois conjuntos:
- **Inclusividade** (o próprio alvo): até 40 sequências de cepas diferentes. Serve para ver o que é conservado e para testar se o ensaio pega todas.
- **Exclusividade** (parentes e vizinhos): organismos parecidos ou que aparecem na mesma amostra. Serve para garantir que o ensaio **não** pega o que não deve.

**Quantos genomas por alvo (vaginal):**
- HPV: ~25 por tipo (raros: HPV70 4, HPV26 5, HPV40/44 8; HPV68 37).
- Cada Candida: 40. Bactérias da vaginose: 33–40 (Sneathia 15, F. vaginae 11).
- CT 35–40; MG 39; NG, TP, VZV 8; H. ducreyi 7; HSV e CMV 10; S. agalactiae 31.
- Ureaplasma, BVAB2, Megasphaera tipo 1: só 1–2 sequências públicas → validação fraca.
- Não-alvos: HPV contra 109 outros tipos; bactérias contra 150–350 sequências; Candida contra 130.

**Limite honesto:** é uma amostra do que existe no banco, não a população inteira. Os genes escolhidos são clássicos e o BLAST comparou cada amplicon com o GenBank inteiro, mas o ideal é confirmar com amostras positivas no laboratório.

**Cuidado com o banco:** registros "UNVERIFIED" foram descartados (alguns "HPV68" eram HPV39).

---

## 6. Escolha da posição-alvo e montagem da região

**Passo 1 — alinhar.** Cada sequência do alvo é colocada lado a lado com a referência na região do gene.

**Passo 2 — consenso.** Em cada posição contamos as bases. A mais comum vira o "consenso". Se a segunda base aparece em **≥10% das cepas**, a posição vira uma letra IUPAC (ex.: R). Assim, nenhum primer cai em posição variável.

**Passo 3 — procurar a posição.** Percorremos o gene posição por posição, olhando um trecho de 31 nt em volta de cada uma. Uma posição só serve se:
- o trecho aparece em ≥85–90% das cepas do alvo, com ~1 diferença no máximo em média;
- a base do meio é igual em ≥95–97% das cepas;
- o trecho é **diferente** dos parentes e dos outros alvos do painel.

Entre as aprovadas, escolhemos a com menos semelhança com não-alvos.

**Passo 4 — a região.** Recortamos 241 nt: 120 de cada lado da base-alvo, que fica no meio entre colchetes:

```
...ACCAAGCCTTATGATCGACGGA[T/C]GCGATCCTTGCACC...
```

Esse é o formato que o ADS lê.

---

## 7. Desenho dos primers

**Regras dos primers de PCR (F e R):**
- 18–25 nt, Tm 55–64 °C (alvo 60 °C), F e R com Tm parecidas.
- GC entre 30% e 70%; nada de 5 letras iguais seguidas (AAAAA).
- No máximo 3 G/C nas 5 últimas bases (ponta 3' não pode ser "grudenta" demais).
- Sem formar dímero (um primer grudando no outro pela ponta 3').
- **Ambiguidade só na metade 5'** (regra do Hema): a ponta 3' tem que ser exata, porque é dali que a enzima começa.
- Amplicon de 70–130 pb (ideal ~100).

**Regras do UEP:** 15–28 nt, sem base ambígua, Tm de extensão ≥47 °C, encostado na base-alvo.

**Filtro anti-reação-cruzada (importante).** A ponta 3' (últimas 10 bases) de cada primer de PCR **não pode existir** em nenhuma sequência não-alvo. Isso evita primers "universais" (que amplificam várias espécies, comum no 16S) — foi o que resolveu 12 ensaios do vaginal.

**Formato pedido pelo Hema para o ADS:** primers em MAIÚSCULO, resto da sequência em minúsculo. O ADS mantém nossos primers de PCR e desenha o UEP:

```
...gtgtttgcaagaatgcaCTGACCACTGCAGAGATTTATTCAtatgca...[G/A]...gcgGTGCTGCCTAGAATTTCATGGAaaaat...
```

---

## 8. Ensaios de resistência

Aqui a posição não é escolhida: vem da literatura (mutação conhecida). Cuidados:

- **Numeração do 23S.** Mutações de macrolídeos (A2058G, A2059G) são citadas na numeração da *E. coli*. Em outra bactéria a mesma posição tem outro número (A2058 da E. coli = A2071 em *M. genitalium*). Alinhamos para achar a posição certa.
- **Códon.** Para mutações em proteína (ex.: ERG11 **Y132F** = tirosina 132 → fenilalanina), traduzimos o gene para confirmar: códon TAT (Y) → TTT (F).
- **Mutações vizinhas.** 2058 e 2059 são lado a lado: cada ensaio marca a outra como variável, para o UEP não ficar em cima dela.
- **Cópias múltiplas.** O gonococo tem 4 cópias do 23S: uma amostra pode ter parte delas mutadas → pico misto. É esperado.
- **Interpretação.** O 23S do gonococo é igual ao de outras *Neisseria*; o de *M. genitalium* é 95% igual ao de *M. pneumoniae*. Por isso **só interpretar resistência se a detecção do organismo for positiva**.

Ensaios do vaginal: MG 23S A2058G/A2059G e parC S83I; NG 23S A2059G/C2611T e gyrA S91F; TP 23S A2058G/A2059G; C. albicans ERG11 Y132F e K143R.

---

## 9. Os testes (como sabemos que está certo)

**PCR in silico (no computador).** Para cada sequência, procuramos onde F e R grudam: últimas 4 bases da ponta 3' perfeitas, até 3 erros no resto, em fitas opostas, produto de 40–800 pb.
- Nas cepas do alvo → **inclusividade** (ex.: "detecta 38 de 40").
- Nos parentes → **especificidade**.
- Nos outros alvos do painel → **reação cruzada no painel**.

**Sinal do UEP.** Se um não-alvo amplifica, vemos se o UEP também gruda nele. Se não gruda, não há sinal (ex.: HPV56 amplifica 3 de 24 HPV66, mas o UEP não gruda → sem falso positivo).

**BLAST.** Cada amplicon é comparado com todo o GenBank. Esperamos que quase todos os acertos sejam do alvo. Acertos "estranhos" geralmente são:
- **erro de cadastro no banco** (ex.: uma "C. albicans" que é C. dubliniensis);
- **nome novo** da mesma espécie (*H. ducreyi* = *Pseudactinobacillus ducreyi*; BVAB2 = *Amygdalobacter*; *C. glabrata* = *Nakaseomyces glabratus*; *C. krusei* = *Pichia kudriavzevii*; *Atopobium vaginae* = *Fannyhessea vaginae*);
- **vetores de laboratório** ou isolados sem nome de tipo.

Dica do BLAST: mandar **um amplicon por vez** (10 juntos ficaram mais de 1 h na fila; um por um, ~30 s cada).

**Conferência de posicionamento.** Um script confere, em todos os ensaios, que F e R estão no lugar, que o UEP encosta no alvo, que a base estendida está certa e que nenhum primer se sobrepõe ao UEP.

---

## 10. Problemas encontrados e como foram resolvidos

| Problema | Solução |
|---|---|
| Alguns "HPV68" do banco eram HPV39 | Descartar registros UNVERIFIED e refazer o HPV68 |
| Primers que amplificavam várias bactérias/fungos (16S, ITS) | Filtro da ponta 3' contra não-alvos; 12 ensaios redesenhados (CT, HD, Lactobacillus, Mobiluncus, NG, C. glabrata, Megasphaera, HPV51/56/58/68) |
| C. parapsilosis amplificava C. tropicalis | Redesenhado em outra região do ITS2 |
| Primer reverse do CT desenhado numa base rara da referência | Corrigido para o consenso das 35 cepas |
| UEP com base ambígua (3 no respiratório: H5, PIV4, MCA) | Não funciona no iPLEX → no formato novo o ADS desenha o UEP |
| Ambiguidade perto da ponta 3' (10 primers) | Trocados por primers com ambiguidade só na metade 5'; testados de novo (cobertura igual ou melhor; rinovírus subiu de 45 para 56 de 60) |
| Sequências com nome errado no banco (Candida) | Contar quantos são: 1–3 "estranhos" no meio de centenas corretos = erro de cadastro |
| BLAST lento | Um amplicon por pedido |

---

## 11. Os arquivos entregues e a planilha

**Arquivos:**
- **Planilha** — resumo de tudo (primers, testes, observações).
- **ADS_input_..._primers_maiusculos.txt** — o que a chefe carrega no ADS (primers em maiúsculo).
- **Regiões em FASTA** — as mesmas regiões com o alvo entre colchetes (para visualizar).
- **Amplicons em FASTA** — o trecho copiado por cada par (para conferência/BLAST).

**Colunas da planilha (vaginal):**

| Coluna | O que é |
|---|---|
| Ensaio | Nome curto: organismo_gene ou organismo_mutação |
| Patógeno | Organismo |
| Tipo | Detecção / Resistência / Controle interno |
| Gene / região | Onde está a base-alvo |
| Origem | Desenho novo; "redundância" = HPV que a chefe já tinha |
| Primer PCR forward / reverse | 5'→3', já com a etiqueta ACGTTGGATG |
| Primer de extensão | UEP 5'→3' ("a definir no ADS" quando o nosso não servia) |
| Base estendida | normal / alternativa, na fita do UEP |
| Amplicon | tamanho do produto, sem as etiquetas |
| Status do teste | BLAST + inclusividade + reação cruzada |
| Observação | limitações e cuidados de interpretação |

**Por que 93 ensaios e 67 FASTAs no respiratório:** 67 desenhados por nós (com região) + 26 do artigo Zhao 2021 (primers do artigo, sem região).

---

## 12. Limitações e próximos passos

- **Primers provisórios:** o ADS ainda monta os poços (multiplex) e pode ajustar primers por peso e dímeros entre ensaios.
- **Testes in silico** dependem do que existe no banco; alvos com poucas sequências (Ureaplasma, BVAB2, Megasphaera) precisam de validação com amostra positiva.
- **TV** também amplifica outros *Trichomonas* (oral e de aves): risco baixo em amostra vaginal.
- **Amplicons >120 pb** (HSV-2, CMV, rinovírus, OXA-48…): funcionam, mas o ADS pode encurtar.
- **ERG11 Y132F e K143R** em poços diferentes (primers sobrepostos).
- **MassARRAY dá presença/ausência**, não quantidade.
- **Próximo passo de bancada:** testar cada ensaio com controle positivo (cepa ou plasmídeo), um negativo, e amostras reais conhecidas.

---

## 13. Como repetir numa próxima vez

Tudo virou uma ferramenta no seu GitHub (**massarray-designer**, privado):

1. Escrever a lista de alvos em `config/meu_painel.json` (modelo e explicação na pasta `config/`).
2. Rodar `python -m massarray all config/meu_painel.json`.
3. A ferramenta baixa as sequências, escolhe a posição, desenha os primers, testa, faz BLAST e gera os arquivos.
4. Conferir com `docs/CHECKLIST.md` antes de enviar.

Os documentos `docs/AULA.md`, `docs/PASSO_A_PASSO.md` e `docs/LICOES_APRENDIDAS.md` complementam este guia. Atenção: a regra nova do Hema (ambiguidade só no 5' e primers em maiúsculo no arquivo do ADS) ainda precisa ser incluída no repositório.

---

## 14. Glossário rápido

- **Amplicon:** pedaço de DNA copiado pela PCR.
- **ADS (Assay Design Suite):** programa da Agena que monta os multiplex.
- **BLAST:** busca de uma sequência contra o banco inteiro do NCBI.
- **Consenso:** sequência com a base mais comum em cada posição.
- **Exclusividade / especificidade:** o ensaio não pega o que não deve.
- **Inclusividade:** o ensaio pega todas as cepas do alvo.
- **In silico:** feito no computador.
- **IUPAC:** letras para posições com mais de uma base possível.
- **Multiplex:** vários ensaios no mesmo poço.
- **Primer:** pedacinho de DNA que inicia a cópia.
- **Ponta 3':** fim do primer, onde a enzima começa; tem que grudar perfeitamente.
- **Tag (ACGTTGGATG):** etiqueta do iPLEX nos primers de PCR.
- **Tm:** temperatura em que o primer se solta.
- **UEP:** primer de extensão, que ganha a base lida pelo aparelho.

## Perguntas que a chefe pode fazer (e respostas curtas)

- **Como as regiões foram escolhidas?** Gene clássico de cada grupo → alinhamento de até 40 sequências → posição conservada em ≥95% das cepas e diferente dos parentes → PCR in silico e BLAST.
- **Quantos genomas?** Até 40 por alvo (mediana ~25); menos nos raros. Comparação com 100–350 sequências de parentes.
- **Os primers estão desenhados?** Sim, testados; provisórios porque o ADS ainda monta os poços.
- **Por que 93 e 67?** 67 nossos + 26 do artigo Zhao.
- **E a ambiguidade?** Só na metade 5' dos primers de PCR; nenhuma no UEP.
- **Está validado?** In silico sim; falta a validação de bancada com amostras positivas.
