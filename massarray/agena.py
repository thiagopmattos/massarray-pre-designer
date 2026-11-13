"""Arquivos no formato da Agena (MassARRAY / iPLEX).

Dois formatos:

1) ENTRADA DO ADS (Assay Design Suite) — o que se carrega no programa da Agena:
       nome<TAB>sequência
   - sequência de ~241 nt com o alvo entre colchetes: ...acgt[A/G]tgca...
   - PRIMERS DE PCR EM MAIÚSCULO e TODO O RESTO EM MINÚSCULO (inclusive onde fica
     o primer de extensão). Assim o ADS sabe que deve MANTER os nossos primers de
     PCR e desenha ele mesmo o primer de extensão (UEP), já escolhendo massas que
     não se sobrepõem dentro de cada poço. (Orientação do Hema.)
   - letras minúsculas r, y, s... marcam posições variáveis entre cepas.

2) TABELA DE ENSAIOS (TSV) — uma linha por ensaio, no layout de colunas da saída
   do ADS que o software do equipamento (Typer) importa:
       WELL TERM SNP_ID 2nd-PCRP 1st-PCRP AMP_LEN UP_CONF MP_CONF Tm(NN) PcGC PWARN
       UEP_DIR UEP_MASS UEP_SEQ EXT1_CALL EXT1_MASS EXT1_SEQ EXT2_CALL EXT2_MASS EXT2_SEQ
   - 1st-PCRP = primer forward, 2nd-PCRP = primer reverse (os dois COM a tag de 10 nt).
   - UEP_DIR: F = UEP na mesma fita da sequência de entrada; R = fita oposta.
   - EXT1 = alelo normal/alvo, EXT2 = alelo alternativo/mutante; *_SEQ = UEP + base.
   - Massas em Da: UEP = oligo de DNA sem fosfato 5'; produto = UEP + terminador
     iPLEX (massas dos terminadores modificados da química iPLEX Gold).
   - WELL: poço sugerido por um agrupamento simples (ver plex()). O ADS faz isso
     melhor; use como ponto de partida ou para conferir.
   CONFIRA o cabeçalho com um arquivo de exemplo do seu equipamento/versão do
   software: se for diferente, mude só a lista COLUMNS abaixo.
"""
import csv

from .seq import TAG, gc, mass, rc, three_prime_dimer, tm_ext

# massa acrescentada ao UEP por cada terminador iPLEX (Da)
TERMINATOR_MASS = {"A": 271.2, "C": 247.2, "G": 287.2, "T": 327.1}

COLUMNS = ["WELL", "TERM", "SNP_ID", "2nd-PCRP", "1st-PCRP", "AMP_LEN", "UP_CONF", "MP_CONF",
           "Tm(NN)", "PcGC", "PWARN", "UEP_DIR", "UEP_MASS", "UEP_SEQ", "EXT1_CALL", "EXT1_MASS",
           "EXT1_SEQ", "EXT2_CALL", "EXT2_MASS", "EXT2_SEQ"]


# ---------------------------------------------------------------- 1) entrada do ADS
def ads_sequence(marker, assay=None):
    """Região com primers de PCR em MAIÚSCULO, resto em minúsculo, alvo [X/Y]."""
    s = list(marker["iu"].lower())
    F = marker["F"]
    if assay:
        fs, fe = assay["Fpos"]
        rs, re_ = assay["Rpos"]
        s[fs:fe + 1] = list(assay["F"].upper())
        s[rs:re_ + 1] = list(rc(assay["R"]).upper())
    ref = marker.get("target") or marker["plain"][F]
    alt = marker.get("alt") or marker.get("a2") or {"A": "G", "G": "A", "C": "T", "T": "C"}[ref]
    return "".join(s[:F]) + f"[{ref.upper()}/{alt.upper()}]" + "".join(s[F + 1:])


def check_ads_sequence(seq, assay=None):
    """Confere o formato: 1 par de colchetes; maiúsculas = exatamente F + rc(R)."""
    probs = []
    if seq.count("[") != 1 or seq.count("]") != 1:
        probs.append("precisa de exatamente um [X/Y]")
        return probs
    outside = seq.split("[")[0] + seq.split("]")[1]
    up = "".join(c for c in outside if c.isupper())
    if assay and up != assay["F"].upper() + rc(assay["R"]).upper():
        probs.append("letras maiúsculas não são exatamente os primers F e R")
    return probs


def write_ads_input(path, markers, primers):
    """Arquivo de entrada do ADS (sem cabeçalho): nome<TAB>sequência."""
    lines = []
    for m in markers:
        a = primers.get(m["id"])
        seq = ads_sequence(m, a)
        probs = check_ads_sequence(seq, a)
        if probs:
            raise ValueError(f"{m['id']}: {'; '.join(probs)}")
        lines.append(f"{m['id']}\t{seq}")
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(lines) + "\n")
    return path


# ---------------------------------------------------------------- 2) tabela de ensaios
def assay_row(marker, a, well=""):
    E = a["E"]
    if not E or any(c not in "ACGT" for c in E):
        # UEP degenerado não funciona no iPLEX: fica em branco para o ADS desenhar
        return {**{c: "" for c in COLUMNS}, "WELL": well, "TERM": "iPLEX", "SNP_ID": marker["id"],
                "2nd-PCRP": TAG + a["R"], "1st-PCRP": TAG + a["F"], "AMP_LEN": a["amp"],
                "PWARN": "UEP a definir no ADS"}
    d = "F" if a["Estrand"] == "+" else "R"
    b1, b2 = a["ext_wt"], a["ext_mut"]
    um = mass(E)
    warn = []
    if any(c not in "ACGT" for c in a["F"] + a["R"]):
        warn.append("primer PCR degenerado")
    if a["amp"] > 120:
        warn.append("amplicon>120")
    return {
        "WELL": well, "TERM": "iPLEX", "SNP_ID": marker["id"],
        "2nd-PCRP": TAG + a["R"], "1st-PCRP": TAG + a["F"], "AMP_LEN": a["amp"],
        "UP_CONF": "", "MP_CONF": "", "Tm(NN)": round(tm_ext(E), 1), "PcGC": round(100 * gc(E), 1),
        "PWARN": ";".join(warn), "UEP_DIR": d, "UEP_MASS": round(um, 1), "UEP_SEQ": E,
        "EXT1_CALL": b1, "EXT1_MASS": round(um + TERMINATOR_MASS[b1], 1), "EXT1_SEQ": E + b1,
        "EXT2_CALL": b2, "EXT2_MASS": round(um + TERMINATOR_MASS[b2], 1), "EXT2_SEQ": E + b2,
    }


def _masses(row):
    return [row["UEP_MASS"], row["EXT1_MASS"], row["EXT2_MASS"]]


def _amplicon(marker, a):
    return marker["plain"][a["Fpos"][0]:a["Rpos"][1] + 1]


def _overlap(m1, a1, m2, a2, k=12):
    """True se um ensaio pode se ligar no amplicon do outro (mesmo lugar do genoma):
    aí os dois não podem ficar no mesmo poço (competem e geram produtos extras)."""
    A2 = _amplicon(m2, a2)
    A1 = _amplicon(m1, a1)
    for p in (a1["F"], a1["R"], a1["E"]):
        t = p[-k:]
        if t in A2 or rc(t) in A2:
            return True
    for p in (a2["F"], a2["R"], a2["E"]):
        t = p[-k:]
        if t in A1 or rc(t) in A1:
            return True
    return False


def _dimer(a1, a2):
    ps1 = [a1["F"], a1["R"], a1["E"]]
    ps2 = [a2["F"], a2["R"], a2["E"]]
    return any(three_prime_dimer(x, y) for x in ps1 for y in ps2)


def plex(markers, primers, min_sep=30.0, max_per_well=36):
    """Agrupamento simples em poços (multiplex), do jeito "guloso":
    ensaios em ordem de massa do UEP; cada um entra no primeiro poço onde
      - todas as massas (UEP, produto normal, produto alternativo) ficam a >= min_sep Da
        de todas as massas dos ensaios já presentes (o aparelho precisa separar os picos);
      - não há sobreposição de região (_overlap) nem dímero de ponta 3' (_dimer);
      - o poço ainda tem vaga (max_per_well).
    Retorna {id: 'W1', ...}."""
    M = {m["id"]: m for m in markers}
    rows = {i: assay_row(M[i], a) for i, a in primers.items() if a and i in M}
    rows = {i: r for i, r in rows.items() if r["UEP_MASS"] != ""}   # sem UEP: o ADS decide
    order = sorted(rows, key=lambda i: rows[i]["UEP_MASS"])
    wells = []
    where = {}
    for i in order:
        r = rows[i]
        placed = False
        for w in wells:
            if len(w) >= max_per_well:
                continue
            ok = True
            for j in w:
                if any(abs(x - y) < min_sep for x in _masses(r) for y in _masses(rows[j])):
                    ok = False; break
                if _overlap(M[i], primers[i], M[j], primers[j]) or _dimer(primers[i], primers[j]):
                    ok = False; break
            if ok:
                w.append(i); placed = True
                break
        if not placed:
            wells.append([i])
        where[i] = f"W{[k for k, w in enumerate(wells) if i in w][0] + 1}"
    return where


def write_assay_tsv(path, markers, primers, wells=None):
    """Tabela de ensaios (TSV) no layout de colunas do ADS/Typer."""
    M = {m["id"]: m for m in markers}
    wells = wells or {}
    rows = []
    for m in markers:
        a = primers.get(m["id"])
        if not a:
            continue
        r = assay_row(m, a, wells.get(m["id"], ""))
        if r["UEP_MASS"] != "" and not (4300 <= r["UEP_MASS"] <= 9000):
            r["PWARN"] = ";".join(x for x in [r["PWARN"], "UEP fora da janela de massa"] if x)
        rows.append(r)
    rows.sort(key=lambda r: (int(r["WELL"][1:]) if r["WELL"] else 10**6, r["UEP_MASS"] or 0))
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    return path
