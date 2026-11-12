"""Escolha da posição-alvo e montagem da 'janela' (região) de cada marcador.

Um marcador = uma janela de 2*flank+1 nt (padrão 241) com a base-alvo no meio.
- 'plain': sequência de consenso (base mais comum em cada posição)
- 'iu'   : a mesma, mas com código IUPAC nas posições que variam entre cepas
           (>=10% das amostras). Isso avisa o desenho de primers para fugir delas.
"""
from .align import align_window, best_hit, find_approx, min_mismatch
from .seq import comp, iupac_code


def build_marker(mid, ref_seq, p, samples, flank=120, minor=0.10):
    """p = posição do alvo (0-based) em ref_seq. samples = lista de sequências (str)."""
    ref_win = ref_seq[p - flank:p + flank + 1]
    if len(ref_win) != 2 * flank + 1:
        # perto da ponta da referência: encolhe a janela de forma simétrica
        f2 = min(p, len(ref_seq) - 1 - p)
        ref_win = ref_seq[p - f2:p + f2 + 1]
        flank = f2
    F = flank
    al = [a for a in (align_window(ref_win, s, F) for s in samples) if a]
    plain, iu, n_var = [], [], 0
    for i, rb in enumerate(ref_win):
        cnt = {}
        for a in al:
            c = a[i]
            if c in "ACGT":
                cnt[c] = cnt.get(c, 0) + 1
        t = sum(cnt.values())
        if t < max(3, len(al) * 0.5):
            plain.append(rb); iu.append(rb)
            continue
        ks = sorted(cnt, key=lambda k: -cnt[k])
        maj = ks[0]
        plain.append(maj)
        keep = [k for k in ks if cnt[k] / t >= minor]
        if len(keep) > 1 and i != F:
            iu.append(iupac_code(keep)); n_var += 1
        else:
            iu.append(maj)
    tc = {}
    for a in al:
        c = a[F]
        if c in "ACGT":
            tc[c] = tc.get(c, 0) + 1
    tt = sum(tc.values())
    plain = "".join(plain)
    return {
        "id": mid, "F": F, "plain": plain, "iu": "".join(iu), "target": plain[F],
        "n": len(samples), "aligned": len(al),
        "target_cons": round(tc.get(plain[F], 0) / tt, 3) if tt else None,
        "target_counts": tc, "n_var": n_var,
    }


def scan_positions(ref, start, end, inc, exc, step=3, half=15, top=5,
                   min_found=0.85, min_base=0.95, max_mean_mm=1.2):
    """Procura boas posições-alvo para DETECÇÃO entre start e end (0-based).

    Para cada posição p, pega o trecho de 31 nt em volta e mede:
      - found: em quantas amostras do alvo (inc) ele aparece com <=5 erros
      - base:  fração em que a base do meio é igual à referência (conservação)
      - mm:    média de erros nas amostras (quanto menor, mais conservado)
      - close: fração das sequências NÃO-alvo (exc) com <=3 erros (reação cruzada)
    Ordena: menos reação cruzada primeiro, depois mais conservado.
    """
    out = []
    lo = max(start, 130)
    hi = min(end, len(ref) - 131)
    for p in range(lo, hi + 1, step):
        q = ref[p - half:p + half + 1]
        if any(c not in "ACGT" for c in q):
            continue
        found = same = mms = 0
        for x in inc:
            h = best_hit(q, x, 5)
            if not h:
                continue
            pos, strand, mm = h
            found += 1; mms += mm
            b = x[pos + half] if strand == "+" else comp(x[pos + half])
            if b == ref[p]:
                same += 1
        n = len(inc)
        if n and (found < min_found * n or same < min_base * found or mms / found > max_mean_mm):
            continue
        close = sum(1 for o in exc if find_approx(q, o, 3))
        out.append({"p": p, "found": found, "n": n, "mean_mm": round(mms / found, 2) if found else 0,
                    "close_frac": round(close / len(exc), 3) if exc else 0.0,
                    "cross_min_mm": min_mismatch(q, exc) if exc else None})
    out.sort(key=lambda r: (r["close_frac"] * 100 + r["mean_mm"] * 2, -(r["cross_min_mm"] or 0)))
    return out[:top]
