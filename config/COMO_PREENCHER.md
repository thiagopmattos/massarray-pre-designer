# Como preencher o config do painel

Arquivo JSON. Um bloco por ensaio dentro de `"markers"`. Modelo: `exemplo_painel.json`.

```json
{
  "panel": "nome_curto_do_painel",
  "markers": [
    {
      "id": "CT_ompA",
      "organism": "Chlamydia trachomatis",
      "gene": "ompA (região 5')",
      "type": "Detecção",
      "expect": "Chlamydia trachomatis",
      "reference": {"acc": "NC_000117.1", "start": 779500, "end": 780900},
      "inclusivity": {"term": "Chlamydia trachomatis[orgn] AND ompA[gene]", "max": 40},
      "exclusivity": [{"term": "Chlamydia muridarum[orgn] AND ompA", "max": 40}],
      "target": {"scan": [200, 1200]},
      "note": "texto que vai para a coluna Observação"
    }
  ]
}
```

| Campo | Obrigatório | O que é |
|---|---|---|
| `id` | sim | Nome do ensaio. Sem espaços. Padrão: `ORGANISMO_gene` ou `ORGANISMO_gene_MUTAÇÃO` |
| `organism` | sim | Nome do patógeno (vai para a planilha) |
| `gene` | não | Gene/região (vai para a planilha) |
| `type` | não | `Detecção` (padrão), `Resistência` ou `Controle interno` |
| `expect` | sim p/ BLAST | Expressão regular com o(s) nome(s) do organismo no GenBank. Use `\|` para nomes alternativos: `"glabrat\|Nakaseomyces"` |
| `reference.acc` | sim | Acesso do GenBank da referência |
| `reference.start/end` | recomendado | Trecho a baixar (1-based). Deixe ~300 nt de folga em volta do alvo |
| `inclusivity.term` | recomendado | Busca no NCBI que traz cepas do alvo. `max` = quantas baixar (40 é bom) |
| `exclusivity` | recomendado | Lista de buscas de parentes/organismos que não podem dar sinal |
| `target.scan` | detecção | `[início, fim]` (1-based no trecho) onde procurar a base-alvo. Sem isso, procura no trecho todo |
| `target.position` | resistência | Posição exata da base (1-based no trecho baixado) |
| `target.alt` | resistência | Base mutante (na fita +) |
| `target.mask` | opcional | Posições vizinhas a marcar como variáveis, relativas ao alvo: `{"1": "R"}` marca a base seguinte como A/G (usado em 23S 2058/2059) |
| `target.codon` | opcional | Confere o aminoácido: `{"cds_start": 101, "aa": 132, "expect": "Y"}` (cds_start = posição 1-based do A do ATG no trecho) |
| `note` | não | Observação para a planilha |

## Dicas de busca no NCBI

- `Organismo[orgn]` restringe ao organismo (inclui sub-táxons).
- `gene[gene]` funciona quando o gene está anotado; senão, use o nome no título: `ompA[title]`.
- `AND 200:5000[slen]` limita o tamanho (evita genomas gigantes ou fragmentos minúsculos).
- `NOT UNVERIFIED[title]` (a ferramenta já remove esses).
- Teste a busca em https://www.ncbi.nlm.nih.gov/nuccore antes de colocar no config.

## Como achar a posição de uma mutação

1. Pegue a sequência do gene na referência (ex.: `gyrA` de NC_002946.2).
2. Ache o códon: aminoácido 91 → nucleotídeos 271–273 do CDS (3×91−2 até 3×91).
3. A posição no genoma = início do CDS + 270 (+ 0, 1 ou 2 conforme a base do códon que muda).
4. Se o gene está na fita negativa, a base no config é o **complemento**. Use `target.codon` para conferir.
5. Para 23S/16S: as mutações são citadas na numeração da *E. coli*. Alinhe o 23S do seu organismo com
   o da *E. coli* (J01695) e ache a posição correspondente.
