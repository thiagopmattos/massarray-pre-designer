# Exemplos: os dois painéis feitos

Cada pasta tem:
- `markers.json` — as regiões-alvo (campo `plain` = consenso, `iu` = com IUPAC, `F` = posição da base-alvo)
- `primers.json` — os primers (F, R, UEP), posições, Tm, massa e bases estendidas
- os arquivos que foram entregues (FASTA das regiões, entrada do ADS, amplicons, planilha CSV)

Para conferir qualquer um:
```bash
python -m massarray verify examples/painel_vaginal_IST_70/markers.json examples/painel_vaginal_IST_70/primers.json
```

## Painel vaginal / IST (70 ensaios) — versão final enviada
Conferido: 0 problemas (1 aviso: HPV68 tem uma posição variável dentro do UEP, anotada na observação).

## Painel respiratório (67 ensaios desenhados)
O `verify` acusa 3 ensaios com **UEP degenerado** (H5_HA, PIV4, MCA_copB) — isso não funciona no
iPLEX e precisa ser corrigido (ver `docs/LICOES_APRENDIDAS.md`, item 5). E 1 ensaio sem primers
(SPN_sorogrupo18_wzy, região rica em A/T: deixado para o ADS). Os 27 ensaios do artigo Zhao 2021
estão só na planilha (`lista_final_ensaios.csv`), pois os primers são do artigo.
