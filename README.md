# massarray-designer

Ferramenta para **desenhar e conferir ensaios de MassARRAY (Agena iPLEX)** que detectam
patógenos (vírus, bactérias, fungos, parasitas) e mutações de resistência em amostras clínicas.

Ela faz, de forma automática e reprodutível, o que foi feito à mão nos dois painéis do laboratório:

| Painel | Ensaios | Pasta |
|---|---|---|
| Respiratório (swab) | 67 desenhados + 27 do artigo Zhao 2021 | `examples/painel_respiratorio_67/` |
| Microbiota vaginal / IST | 70 (58 detecção, 10 resistência, 2 controles humanos) | `examples/painel_vaginal_IST_70/` |

> **Se você nunca viu MassARRAY:** leia primeiro [`docs/AULA.md`](docs/AULA.md) (15 min). Depois
> [`docs/PASSO_A_PASSO.md`](docs/PASSO_A_PASSO.md) para usar a ferramenta, e
> [`docs/LICOES_APRENDIDAS.md`](docs/LICOES_APRENDIDAS.md) antes de mandar qualquer coisa para a chefe.

---

## O que a ferramenta entrega

Para cada patógeno/mutação da sua lista:

1. **Região-alvo** (~241 nt) com a base que será lida no meio, entre colchetes: `...ACGT[G/A]TTCA...`
   — em FASTA e no formato de entrada do **Assay Design Suite (ADS)** da Agena.
2. **Três primers**: PCR forward e reverse (já com a etiqueta `ACGTTGGATG`) e o **primer de extensão (UEP)**.
3. **Testes in silico**:
   - *Inclusividade*: o ensaio pega as várias cepas/variantes do alvo?
   - *Especificidade*: amplifica parentes próximos? Se amplifica, o UEP gera sinal neles?
   - *Reação cruzada no painel*: um ensaio pega o alvo de outro ensaio?
   - *Posicionamento*: confere que cada primer está no lugar certo e que o UEP não é degenerado.
4. **BLAST** de cada amplicon no NCBI (resumido: quantos acertos são do alvo e qual o "outro" mais parecido).
5. **Planilha (CSV)** com colunas no padrão que a chefe já conhece.

## Instalação

Só precisa de **Python 3.9+**. Não usa nenhuma biblioteca externa.

```bash
git clone https://github.com/<seu-usuario>/massarray-designer.git
cd massarray-designer
python -m unittest discover -s tests      # deve terminar com "OK" (não precisa de internet)
```

Para baixar dados do NCBI, informe seu e-mail (regra do NCBI). Opcional: uma API key deixa 3× mais rápido.

```bash
export NCBI_EMAIL="voce@exemplo.com"
export NCBI_API_KEY="..."        # opcional: https://www.ncbi.nlm.nih.gov/account/settings/
```
No Windows (PowerShell): `$env:NCBI_EMAIL="voce@exemplo.com"`.

## Uso rápido

1. Copie `config/exemplo_painel.json` e edite com os seus alvos (explicação de cada campo em
   [`config/COMO_PREENCHER.md`](config/COMO_PREENCHER.md)).
2. Rode tudo:

```bash
python -m massarray all config/meu_painel.json
```

ou passo a passo (recomendado na primeira vez, para olhar cada resultado):

```bash
python -m massarray fetch   config/meu_painel.json   # 1. baixa referência + cepas do alvo + parentes
python -m massarray markers config/meu_painel.json   # 2. escolhe a base-alvo e monta a região
python -m massarray primers config/meu_painel.json   # 3. desenha F, R e UEP
python -m massarray check   config/meu_painel.json   # 4. testes in silico
python -m massarray blast   config/meu_painel.json   # 5. BLAST (lento: ~30 s a 2 min por ensaio)
python -m massarray export  config/meu_painel.json   # 6. FASTA, entrada do ADS e planilha
```

Resultados em `work/<nome_do_painel>/saida/`. Os arquivos intermediários
(`markers.json`, `primers.json`, `check.json`, `blast.json`) e todas as sequências baixadas
(`cache/`) ficam salvos: dá para parar, corrigir o config e continuar.

Para **conferir um painel pronto** (sem internet):

```bash
python -m massarray verify examples/painel_vaginal_IST_70/markers.json examples/painel_vaginal_IST_70/primers.json
```

## Estrutura

```
massarray/            código
  seq.py              complemento reverso, IUPAC, Tm (vizinhos mais próximos), massa, regras de primer
  ncbi.py             busca e download no NCBI (E-utilities) com cache
  align.py            busca aproximada rápida e alinhamento da amostra na região
  marker.py           escolha da base-alvo e região de consenso com IUPAC
  design.py           desenho dos 3 primers (com filtro anti-reação-cruzada)
  insilico.py         PCR in silico, sinal do UEP, conferência de posicionamento
  blast.py            BLAST pela internet e resumo
  export.py           FASTA / ADS / CSV
  cli.py              linha de comando
config/               modelo de configuração do painel
docs/                 aula, passo a passo, lições aprendidas, checklist
examples/             os dois painéis feitos (entradas, primers e arquivos entregues)
tests/                testes automáticos (rodam offline)
```

## Limites (leia!)

- Os primers são **provisórios**: o ADS da Agena ainda precisa montar os multiplex (quais ensaios
  vão juntos em cada poço) e pode ajustar primers por causa de massa e dímeros entre ensaios.
- Os testes são **in silico**: dependem do que existe no GenBank. Alvo com poucas sequências
  públicas (ex.: Ureaplasma, BVAB2) tem validação fraca — teste com amostra positiva.
- MassARRAY dá **presença/ausência**, não quantidade.
- A ferramenta não substitui a revisão de alguém da área (ver [`docs/CHECKLIST.md`](docs/CHECKLIST.md)).

## Licença

MIT — veja `LICENSE`.
