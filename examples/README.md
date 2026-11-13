# Exemplos: os dois painéis feitos

Cada pasta tem:
- `markers.json`: as regiões-alvo (`plain` = consenso, `iu` = com IUPAC, `F` = posição da base-alvo)
- `primers.json`: os primers (F, R, UEP), posições, Tm, massa e bases estendidas (versão final, já com a regra do 5')
- `*_ADS_entrada.txt`: **o que foi enviado para o ADS** (primers em MAIÚSCULO, resto em minúsculo)
- `*_ensaios_agena.tsv`: tabela de ensaios no layout da Agena, com poços sugeridos
- `lista_final*.csv`: a planilha entregue
- FASTA das regiões e dos amplicons

Conferir / regenerar:
```bash
python -m massarray verify examples/painel_vaginal_IST_70/markers.json examples/painel_vaginal_IST_70/primers.json
python -m massarray agena  examples/painel_vaginal_IST_70/markers.json examples/painel_vaginal_IST_70/primers.json saida/ painel_vaginal_70
```

## Painel vaginal / IST (70 ensaios)
`verify`: 0 problemas (1 aviso: o HPV68 tem uma posição variável embaixo do UEP, anotada na observação).
Revisão final: HPV82, HSV-2, Candida spp., NG gyrA e ERG11 Y132F com primers novos (ambiguidade só no 5').

## Painel respiratório (67 ensaios desenhados)
`verify` acusa 3 ensaios com **UEP degenerado** (H5_HA, PIV4, MCA_copB): no arquivo do ADS o UEP fica
em minúsculo e **o ADS desenha um novo**; no TSV aparecem como "UEP a definir no ADS".
1 ensaio sem primers (SPN_sorogrupo18_wzy, região rica em A/T: deixado para o ADS).
Revisão final: H5, rinovírus, PIV4 e OXA-48 com primers novos.
Os 26 ensaios do artigo Zhao 2021 estão só na planilha (`lista_final_ensaios.csv`): os primers são do artigo.
