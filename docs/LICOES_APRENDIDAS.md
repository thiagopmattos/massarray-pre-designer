# Lições aprendidas nos dois painéis

Coisas que deram errado (ou quase) e como foram resolvidas. A ferramenta já trata a maioria delas;
as outras precisam de olho humano.

## Dados do GenBank

1. **Registros "UNVERIFIED" podem estar com o nome errado.** Vários "HPV68" eram, na verdade, HPV39.
   → `ncbi.fetch_set` descarta registros UNVERIFIED.
2. **Espécies com nome trocado no banco.** Em Candida, várias sequências de "C. albicans" eram
   C. dubliniensis, C. parapsilosis etc. (achamos até "Malassezia" e "Trichophyton" idênticos a
   C. parapsilosis, o que é impossível). Sinal típico: 0 erros contra 1–3 sequências "de outra
   espécie" no meio de dezenas corretas. → Não deixe 1 registro estranho derrubar um bom ensaio;
   conte quantos são (no BLAST, "236 de 246 são C. parapsilosis" é ótimo).
3. **Nomes novos de organismos** confundem o BLAST:
   *Haemophilus ducreyi* → *Pseudactinobacillus ducreyi*; BVAB2 → *Amygdalobacter*;
   Megasphaera tipo 1 → *M. lornae*; *C. glabrata* → *Nakaseomyces glabratus*;
   *C. krusei* → *Pichia kudriavzevii*; *C. lusitaniae* → *Clavispora lusitaniae*;
   *C. parapsilosis* aparece como *Lodderomyces parapsilosis*; *Atopobium vaginae* → *Fannyhessea vaginae*;
   *U. parvum* era "U. urealyticum biovar 1" antes de 2002. → Coloque os nomes antigos e novos no
   campo `expect` (é uma expressão regular: `"glabrat|Nakaseomyces"`).
4. **Sequências parciais** derrubam a "inclusividade" artificialmente (ex.: CAND_pan 12/70; MG MgPa
   com amostras de 231 pb). Olhe o tamanho das sequências que falharam antes de concluir.

## Desenho

5. **O UEP nunca pode ter base degenerada (IUPAC).** No iPLEX, um UEP degenerado é uma mistura de
   oligos de massas diferentes → picos duplicados. No painel respiratório saíram 3 assim
   (H5_HA, PIV4, MCA_copB) e tiveram que ser revistos; no vaginal o filtro já existia.
   → `design.py` só aceita UEP com A/C/G/T; `verify` acusa se houver.
6. **Primers "universais" escondidos.** Em 16S e 18S/ITS, o primer pode cair num trecho idêntico em
   muitas bactérias/fungos. O par amplifica tudo e só o UEP dá a especificidade — arriscado.
   → filtro: a ponta 3' (10 nt) de cada primer de PCR não pode existir nas sequências não-alvo.
   Resolveu CT, HD, LACT, MOB, NG, CGLA, MEGA1 e 4 HPV.
7. **Conferir posição contra o consenso, não só contra a referência.** O primer reverse do CT foi
   desenhado numa referência que tinha uma base rara; o consenso das 35 cepas tinha outra. Corrigido
   e estendido em 2 nt para manter a Tm. → `verify_placement` compara com `plain` e `iu`.
8. **Janela na ponta da referência.** Se o alvo está a <120 nt do fim da sequência baixada, a janela
   fica assimétrica e a base-alvo sai do meio. → `build_marker` encolhe a janela simetricamente.
   Melhor ainda: baixe com folga (`start`/`end` com ±300 nt).
9. **Parentes que amplificam mas não dão sinal** são aceitáveis (HPV56 amplifica 3/24 HPV66 mas o UEP
   não gruda). Anote na observação.
10. **Primers de detecção de resistência só valem com a detecção positiva.** O 23S de *N. gonorrhoeae*
    é idêntico ao de *N. meningitidis* e Neisseria comensais; o 23S de *M. genitalium* é 95% igual
    ao de *M. pneumoniae*.
11. **Amplicons >120 pb** funcionam, mas o ADS pode encurtar. Anote.

12b. **Base ambígua só na metade 5' do primer de PCR** (orientação do Hema). A enzima começa pela
    ponta 3'; uma base errada ali derruba a PCR. Dez primers dos dois painéis tinham ambiguidade perto
    do 3' e foram trocados. Em todos, a cobertura das cepas ficou igual ou melhor (rinovírus: 45 → 56
    de 60). → `seq.ambiguity_only_5prime` no desenho e na conferência.
12c. **Formato do ADS**: primers de PCR em MAIÚSCULO e o resto em minúsculo; o ADS mantém os primers
    e desenha o UEP. → `agena.write_ads_input`.
12d. **Ensaios da mesma região não vão no mesmo poço** (ex.: ERG11 Y132F e K143R têm primers
    sobrepostos). → `agena.plex` separa automaticamente.

## NCBI

12. **BLAST: um amplicon por pedido.** Dez juntos (com HPV16/18) ficaram >1 h na fila; um por um,
    ~30 s cada. Use megablast.
13. **Limite de acesso**: ≤3 pedidos/s sem API key. A ferramenta já espera entre pedidos.
14. **Guarde tudo o que baixou** (`work/<painel>/cache/`). Assim o trabalho é reprodutível e dá para
    refazer as contas sem baixar de novo.

## Comunicação

15. Mande só quando **tudo** estiver conferido (BLAST pendente = não terminado).
16. Primers da chefe: confira também. No painel vaginal, o primer de F. vaginae dela termina com
    "a" minúsculo onde as sequências têm G — erro na última base 3' costuma impedir a PCR.
17. Explique as limitações na coluna Observação: é o que evita resultado mal interpretado depois.
