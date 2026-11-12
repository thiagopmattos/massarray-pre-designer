"""Funções básicas de sequência: complemento reverso, IUPAC, Tm, massa, checagens de primer.

Tudo aqui é Python puro (sem dependências), para rodar em qualquer computador.
"""
import math

TAG = "ACGTTGGATG"  # etiqueta de 10 nt que o iPLEX coloca na ponta 5' dos primers de PCR

IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
}
IUPAC_REV = {"".join(sorted(v)): k for k, v in IUPAC.items()}
_COMP = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")


def rc(s):
    """Complemento reverso (aceita códigos IUPAC)."""
    return s.upper().translate(_COMP)[::-1]


def comp(base):
    return base.translate(_COMP)


def iupac_code(bases):
    """Conjunto de bases -> código IUPAC ('AG' -> 'R')."""
    key = "".join(sorted(set(bases)))
    return IUPAC_REV.get(key, "N")


def compatible(a, b):
    """True se as duas bases (IUPAC) podem ser iguais."""
    return bool(set(IUPAC.get(a, "N")) & set(IUPAC.get(b, "N")))


def is_clean(s):
    return all(c in "ACGT" for c in s)


def gc(s):
    s = s.upper()
    return (s.count("G") + s.count("C")) / max(len(s), 1)


def has_run(s, n=5):
    """Repetição de n bases iguais seguidas (ex.: AAAAA)."""
    return any(c * n in s for c in "ACGT")


def djb(s):
    """Checksum simples usado para conferir que uma sequência foi copiada sem erro."""
    h = 5381
    for ch in s:
        h = ((h * 33) & 0xFFFFFFFF) ^ ord(ch)
    return format(h, "x")


# --- Tm por vizinhos mais próximos (SantaLucia 1998) com correção de sal (von Ahsen) ---
NN = {
    "AA": (-7.9, -22.2), "TT": (-7.9, -22.2), "AT": (-7.2, -20.4), "TA": (-7.2, -21.3),
    "CA": (-8.5, -22.7), "TG": (-8.5, -22.7), "GT": (-8.4, -22.4), "AC": (-8.4, -22.4),
    "CT": (-7.8, -21.0), "AG": (-7.8, -21.0), "GA": (-8.2, -22.2), "TC": (-8.2, -22.2),
    "CG": (-10.6, -27.2), "GC": (-9.8, -24.4), "GG": (-8.0, -19.9), "CC": (-8.0, -19.9),
}


def tm(s, oligo=1e-7, na=0.05, mg=0.002, dntp=0.0005):
    """Tm (°C). Bases degeneradas são tratadas como A (estimativa conservadora).

    Padrões: primer de PCR 100 nM, Na+ 50 mM, Mg2+ 2 mM, dNTP 0,5 mM.
    Para o primer de extensão use tm(s, oligo=6e-7, mg=0).
    """
    s = "".join(c if c in "ACGT" else "A" for c in s.upper())
    dH, dS = 0.2, -5.7
    if s[0] in "AT":
        dH += 2.2; dS += 6.9
    if s[-1] in "AT":
        dH += 2.2; dS += 6.9
    for i in range(len(s) - 1):
        h, e = NN[s[i:i + 2]]
        dH += h; dS += e
    naeq = na + 3.3 * math.sqrt(max(mg - dntp, 0))
    dS += 0.368 * (len(s) - 1) * math.log(naeq)
    return dH * 1000 / (dS + 1.987 * math.log(oligo / 4)) - 273.15


def tm_ext(s):
    """Tm do primer de extensão (UEP), condições da reação iPLEX."""
    return tm(s, oligo=6e-7, na=0.05, mg=0.0, dntp=0.0)


def mass(s):
    """Massa aproximada (Da) de um oligo de DNA sem fosfato 5'."""
    s = s.upper()
    return round(s.count("A") * 313.21 + s.count("T") * 304.2 + s.count("C") * 289.18
                 + s.count("G") * 329.21 - 61.96, 1)


def ok_pcr_primer(p, gc_min=0.30, gc_max=0.70):
    """Regras simples de primer de PCR: sem IUPAC, GC 30–70%, sem 5 bases iguais,
    no máximo 3 G/C nas 5 últimas bases (ponta 3' não muito 'grudenta')."""
    p = p.upper()
    if not is_clean(p):
        return False
    if not (gc_min <= gc(p) <= gc_max):
        return False
    if has_run(p):
        return False
    if sum(c in "GC" for c in p[-5:]) > 3:
        return False
    return True


def three_prime_dimer(a, b, n=5):
    """True se as últimas n bases de 'a' são complementares a algum trecho de 'b'."""
    return a[-n:] in rc(b)


CODON = {}
_b = "TCAG"
_aa = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for i, x in enumerate(_b):
    for j, y in enumerate(_b):
        for k, z in enumerate(_b):
            CODON[x + y + z] = _aa[16 * i + 4 * j + k]


def translate(s):
    s = s.upper()
    return "".join(CODON.get(s[i:i + 3], "X") for i in range(0, len(s) - 2, 3))
