# Checklist antes de enviar

Marque cada item. Se algum falhar, não envie ainda.

## Arquivos
- [ ] `python -m massarray verify markers.json primers.json` → "0 com problema"
- [ ] FASTA de regiões: 1 cabeçalho por ensaio, cada sequência com exatamente um `[X/Y]`
- [ ] Entrada do ADS: 1 linha por ensaio, `nome<TAB>sequência`, só os primers F e R em MAIÚSCULO
- [ ] Tabela TSV da Agena: cabeçalho conferido com um arquivo do equipamento; nenhum PWARN sem explicação
- [ ] Amplicons: tamanho igual ao da planilha
- [ ] Planilha: todas as linhas com F, R, UEP, base e amplicon

## Primers
- [ ] Nenhum UEP com letra diferente de A/C/G/T
- [ ] F e R começam com `ACGTTGGATG`
- [ ] Primers degenerados (só PCR) com a base ambígua SÓ na metade 5' (máx. 2) e anotados na observação
- [ ] Amplicons fora de 80–120 pb anotados

## Testes
- [ ] Inclusividade ≥85% (ou explicada: sequências parciais, poucos dados)
- [ ] Não-alvo com sinal = 0 (ou explicado)
- [ ] Sem reação cruzada no painel (ou explicada: ex.: Candida spp. pega as espécies, de propósito)
- [ ] BLAST concluído para **todos** os ensaios
- [ ] "Outro mais parecido" ≥95% investigado (erro de cadastro? nome novo? parente real?)

## Resistência
- [ ] Posição confirmada pela numeração certa (E. coli para 23S)
- [ ] Códon traduzido e conferido (aminoácido esperado)
- [ ] Observação: "interpretar só se a detecção do organismo for positiva" quando o gene é compartilhado

## Mensagem
- [ ] Lista do que foi feito, do que ficou de fora e por quê
- [ ] Números conferidos (ex.: quantos ensaios de cada tipo)
- [ ] Primers provisórios: o ADS ainda monta o multiplex
