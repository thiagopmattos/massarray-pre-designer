"""Desenho dos 3 primers de um ensaio iPLEX a partir da janela do marcador.

  primer F  ->  ..........[UEP]X..........  <- primer R
  O UEP (primer de extensão) termina exatamente uma base antes do alvo X,
  na fita + (antes do alvo) ou na fita - (depois do alvo, lido ao contrário).

Regras (as mesmas usadas nos painéis respiratório e vaginal):
  - PCR: 18–25 nt, Tm 55–64 °C (alvo 60), GC 30–70%, sem 5 bases iguais,
    <=3 G/C nas 5 últimas bases. De preferência SEM base variável (IUPAC);
    se for inevitável (vírus muito variáveis), até 2 bases ambíguas e SÓ na
    metade 5' do primer (nunca perto da ponta 3').
  - Amplicon: 70–130 pb (ideal ~100; o iPLEX trabalha bem com 80–120).
  - UEP: 15–28 nt, SEM base variável, Tm >=47 °C (condição da extensão).
  - Sem dímero de ponta 3' entre os três primers.
  - Opcional: 'exc' (sequências não-alvo). Se dado, a ponta 3' (10 nt) de cada
    primer de PCR não pode existir em nenhuma delas. Foi isso que resolveu os
    primers 'universais' (ex.: 16S que amplificava várias bactérias).
  - Opcional: 'inc' (sequências do alvo). Se dado, cada primer de PCR tem que
    se ligar (<=1 erro, ponta 3' exata) em >= inc_frac delas.
"""
from .insilico import primer_hits
from .seq import (ambiguity_only_5prime, comp, has_run, is_clean, mass, ok_pcr_primer, rc,
                  three_prime_dimer, tm, tm_ext)

DEFAULTS = dict(len_min=18, len_max=25, tm_min=55.0, tm_max=64.0, tm_opt=60.0,
                amp_min=70, amp_max=130, amp_opt=100, uep_min=15, uep_max=28,
                uep_tm_min=47.0, uep_tm_max=62.0, kmer=10, inc_frac=0.9,
                allow_degenerate=False)


def kmer_set(seqs, k):
    S = set()
    for s0 in seqs:
        for s in (s0, rc(s0)):
            for i in range(len(s) - k + 1):
                S.add(s[i:i + k])
    return S


def _inc_ok(p, inc, frac):
    if not inc:
        return True
    ok = sum(1 for s in inc if primer_hits(p, s, maxmm=1))
    return ok >= frac * len(inc)


def _pcr_ok(p, allow_deg):
    if ok_pcr_primer(p):
        return True
    if not allow_deg:
        return False
    if not ambiguity_only_5prime(p):
        return False
    q = "".join(c if c in "ACGT" else "A" for c in p)
    return ok_pcr_primer(q, 0.25, 0.72)


def design(marker, inc=None, exc=None, **opts):
    o = dict(DEFAULTS, **opts)
    s, Fc = marker["iu"], marker["F"]
    L = len(s)
    excK = kmer_set(exc, o["kmer"]) if exc else None
    inc_long = [x for x in (inc or []) if len(x) >= 250][:25]

    def good(p):
        if not _pcr_ok(p, o["allow_degenerate"]):
            return None
        t = tm(p)
        if not (o["tm_min"] <= t <= o["tm_max"]):
            return None
        if excK is not None and p[-o["kmer"]:] in excK:
            return None
        return t

    Fc_list, Rc_list = [], []
    for ln in range(o["len_min"], o["len_max"] + 1):
        for st in range(0, Fc - ln + 1):
            p = s[st:st + ln]
            t = good(p)
            if t is not None:
                Fc_list.append((st, st + ln - 1, p, t))
        for st in range(Fc + 1, L - ln + 1):
            p = rc(s[st:st + ln])
            t = good(p)
            if t is not None:
                Rc_list.append((st, st + ln - 1, p, t))
    if inc_long:
        Fc_list = [f for f in Fc_list if _inc_ok(f[2], inc_long, o["inc_frac"])]
        Rc_list = [r for r in Rc_list if _inc_ok(r[2], inc_long, o["inc_frac"])]

    best = None
    for strand in "+-":
        for el in range(o["uep_min"], o["uep_max"] + 1):
            if strand == "+":
                st = Fc - el
                if st < 0:
                    break
                e = s[st:Fc]; reg = (st, Fc - 1)
            else:
                en = Fc + 1 + el
                if en > L:
                    break
                e = rc(s[Fc + 1:en]); reg = (Fc + 1, en - 1)
            if not is_clean(e) or has_run(e):
                continue
            te = tm_ext(e)
            if te < o["uep_tm_min"]:
                continue
            if te > o["uep_tm_max"]:
                break
            fset = [f for f in Fc_list if f[1] < min(reg[0], Fc)]
            rset = [r for r in Rc_list if r[0] > max(reg[1], Fc)]
            for fs, fe, fp, tf in fset:
                for rs, re_, rp, tr in rset:
                    amp = re_ - fs + 1
                    if amp < o["amp_min"] or amp > o["amp_max"]:
                        continue
                    score = (abs(tf - o["tm_opt"]) + abs(tr - o["tm_opt"]) + abs(tf - tr) * 1.5
                             + abs(amp - o["amp_opt"]) * 0.08 + abs(te - 53) * 0.3)
                    if best is not None and score >= best["score"]:
                        continue
                    if (three_prime_dimer(fp, rp) or three_prime_dimer(rp, fp)
                            or three_prime_dimer(e, fp) or three_prime_dimer(e, rp)):
                        continue
                    best = dict(score=round(score, 2), F=fp, R=rp, Fpos=[fs, fe], Rpos=[rs, re_],
                                amp=amp, tF=round(tf, 1), tR=round(tr, 1), E=e, Estrand=strand,
                                tE=round(te, 1), Emass=mass(e))
            break  # para cada fita, usa o menor UEP que passou
    if best:
        tgt, alt = marker["target"], marker.get("alt") or {"A": "G", "G": "A", "C": "T", "T": "C"}[marker["target"]]
        best["ext_wt"] = tgt if best["Estrand"] == "+" else comp(tgt)
        best["ext_mut"] = alt if best["Estrand"] == "+" else comp(alt)
    return best


def design_with_fallbacks(marker, inc=None, exc=None, **opts):
    """Tenta do mais rígido ao mais permissivo e diz qual nível funcionou."""
    tries = [
        ("padrão", {}),
        ("k-mer 9", {"kmer": 9}),
        ("amplicon 60–140", {"amp_min": 60, "amp_max": 140}),
        ("primer degenerado permitido", {"allow_degenerate": True, "amp_min": 60, "amp_max": 140}),
        ("sem filtro de exclusão", {"amp_min": 60, "amp_max": 140, "allow_degenerate": True, "_noexc": True}),
    ]
    for name, extra in tries:
        o = dict(opts, **{k: v for k, v in extra.items() if not k.startswith("_")})
        d = design(marker, inc, None if extra.get("_noexc") else exc, **o)
        if d:
            d["nivel"] = name
            return d
    return None
