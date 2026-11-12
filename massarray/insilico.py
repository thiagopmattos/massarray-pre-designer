"""PCR 'no computador' (in silico): o par de primers amplificaria esta sequência?

Regra usada (a mesma dos dois painéis):
  - as 4 últimas bases da ponta 3' do primer têm que bater perfeitamente;
  - no resto do primer, aceita até 3 erros;
  - os dois primers têm que estar em fitas opostas, virados um para o outro,
    gerando um produto de 40 a 800 pb.
Depois, para os produtos que 'amplificam', conferimos se o primer de extensão
(UEP) também se liga e qual base ele estenderia: é isso que gera o sinal.
"""
from .seq import IUPAC, comp, rc


def _ok(p, s):
    return s in IUPAC.get(p, "N") or s == p


def _occ(k, s):
    i = s.find(k)
    while i != -1:
        yield i
        i = s.find(k, i + 1)


def primer_hits(p, s, maxmm=3, exact3=4):
    """Sítios de ligação do primer p em s: lista de (pos, fita, erros).
    Rápido: primeiro acha a ponta 3' exata com str.find, depois conta os erros."""
    L = len(p)
    pr = rc(p)
    out = []
    tail, head = p[-exact3:], pr[:exact3]
    clean = all(c in "ACGT" for c in tail)
    plus = (o - (L - exact3) for o in _occ(tail, s)) if clean else range(len(s) - L + 1)
    for i in plus:
        if i < 0 or i + L > len(s):
            continue
        if not clean and not all(_ok(p[j], s[i + j]) for j in range(L - exact3, L)):
            continue
        m = 0
        for j in range(L - exact3):
            if not _ok(p[j], s[i + j]):
                m += 1
                if m > maxmm:
                    break
        if m <= maxmm:
            out.append((i, "+", m))
    minus = _occ(head, s) if clean else range(len(s) - L + 1)
    for i in minus:
        if i + L > len(s):
            continue
        if not clean and not all(_ok(pr[j], s[i + j]) for j in range(exact3)):
            continue
        m = 0
        for j in range(exact3, L):
            if not _ok(pr[j], s[i + j]):
                m += 1
                if m > maxmm:
                    break
        if m <= maxmm:
            out.append((i, "-", m))
    return out


def amplifies(F, R, s, maxmm=3, min_len=40, max_len=800):
    """Retorna (tamanho, erros_totais, inicio, fim) do melhor produto, ou None."""
    fh = primer_hits(F, s, maxmm)
    if not fh:
        return None
    rh = primer_hits(R, s, maxmm)
    if not rh:
        return None
    best = None
    for (fi, fs, fm) in fh:
        for (ri, rs, rm) in rh:
            if fs == "+" and rs == "-" and ri > fi:
                st, en = fi, ri + len(R)
            elif fs == "-" and rs == "+" and fi > ri:
                st, en = ri, fi + len(F)
            else:
                continue
            ln = en - st
            if min_len <= ln <= max_len and (best is None or fm + rm < best[1]):
                best = (ln, fm + rm, st, en)
    return best


def extension_base(E, s, start=0, end=None, maxmm=2):
    """Onde o UEP se liga dentro de s[start:end] e qual base ele estenderia.
    Retorna lista de (erros, base_estendida)."""
    end = len(s) if end is None else end
    sub = s[start:end]
    out = []
    for (i, st, m) in primer_hits(E, sub, maxmm, exact3=4):
        if st == "+":
            k = i + len(E)
            if k < len(sub):
                out.append((m, sub[k]))
        else:
            k = i - 1
            if k >= 0:
                out.append((m, comp(sub[k])))
    return out


def evaluate(assay, inc=(), exc=(), panel=None):
    """Avalia um ensaio: inclusividade, especificidade e reação cruzada no painel.

    assay: dict com F, R, E, ext_wt (base esperada).
    panel: dict {outro_id: [sequências]} para checar reação cruzada entre ensaios.
    Um produto 'perigoso' é o que amplifica E em que o UEP se liga (<=1 erro).
    """
    F, R, E = assay["F"], assay["R"], assay["E"]
    res = {"inc_ok": 0, "inc_n": len(inc), "exc_amp": 0, "exc_signal": 0, "exc_n": len(exc), "panel": {}}
    for s in inc:
        if amplifies(F, R, s):
            res["inc_ok"] += 1
    for s in exc:
        a = amplifies(F, R, s)
        if a:
            res["exc_amp"] += 1
            if any(m <= 1 for m, b in extension_base(E, s, a[2], a[3])):
                res["exc_signal"] += 1
    for oid, seqs in (panel or {}).items():
        amp = sig = 0
        for s in seqs:
            a = amplifies(F, R, s)
            if a:
                amp += 1
                if any(m <= 1 for m, b in extension_base(E, s, a[2], a[3])):
                    sig += 1
        if amp:
            res["panel"][oid] = {"amplifica": amp, "com_sinal": sig, "de": len(seqs)}
    return res


def verify_placement(marker, assay):
    """Confere, na própria janela do marcador, que tudo está no lugar certo:
    F e R no lugar, UEP encostado na base-alvo, base estendida correta e
    nenhum primer sobreposto ao UEP. Retorna lista de problemas (vazia = OK)."""
    from .seq import compatible
    pl, iu, Fc = marker["plain"], marker["iu"], marker["F"]
    E = assay["E"]
    probs = []

    def match(p, w):
        return len(p) == len(w) and all(compatible(a, b) for a, b in zip(p, w))

    fs, fe = assay["Fpos"]
    rs, re_ = assay["Rpos"]
    if not match(assay["F"], iu[fs:fe + 1]) and not match(assay["F"], pl[fs:fe + 1]):
        probs.append("primer F não bate na posição indicada")
    if not match(rc(assay["R"]), iu[rs:re_ + 1]) and not match(rc(assay["R"]), pl[rs:re_ + 1]):
        probs.append("primer R não bate na posição indicada")
    if assay["Estrand"] == "+":
        e_lo, e_hi = Fc - len(E), Fc - 1
        okE = pl[e_lo:Fc] == E or match(E, iu[e_lo:Fc])
        base = pl[Fc]
    else:
        e_lo, e_hi = Fc + 1, Fc + len(E)
        okE = rc(E) == pl[Fc + 1:Fc + 1 + len(E)] or match(rc(E), iu[Fc + 1:Fc + 1 + len(E)])
        base = comp(pl[Fc])
    if not okE:
        probs.append("primer de extensão não está encostado na base-alvo")
    if assay.get("ext_wt") and assay["ext_wt"] != base:
        probs.append(f"base estendida deveria ser {base}, está {assay['ext_wt']}")
    if not (fe < e_lo and rs > e_hi):
        probs.append("primer de PCR sobreposto ao primer de extensão")
    if any(c not in "ACGT" for c in E):
        probs.append("primer de extensão degenerado (não pode no iPLEX)")
    var_e = [i for i in range(e_lo, e_hi + 1) if iu[i] not in "ACGT"]
    if var_e:
        probs.append(f"aviso: posição variável embaixo do primer de extensão ({len(var_e)})")
    if assay["amp"] != re_ - fs + 1:
        probs.append("tamanho do amplicon não confere")
    return probs
