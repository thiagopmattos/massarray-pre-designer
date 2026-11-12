"""Busca aproximada (com até N erros) rápida, em Python puro.

Usa o 'princípio da casa dos pombos': se um trecho de 31 nt aparece com no
máximo 3 erros, então pelo menos um de 4 pedaços dele aparece SEM erro.
Procuramos esses pedaços com str.find (muito rápido) e só depois conferimos
a sequência inteira.
"""
from .seq import rc, compatible


def hamming(a, b, limit=None):
    m = 0
    for x, y in zip(a, b):
        if x != y and not compatible(x, y):
            m += 1
            if limit is not None and m > limit:
                return m
    return m


def _chunks(q, k):
    n = len(q)
    size = n // k
    out = []
    for i in range(k):
        st = i * size
        en = n if i == k - 1 else st + size
        out.append((st, q[st:en]))
    return out


def find_approx(q, s, maxmm=3, strands="+-"):
    """Todas as ocorrências de q em s com até maxmm erros.
    Retorna lista de (pos, fita, erros); pos é na fita + de s."""
    hits = {}
    L = len(q)
    for strand in strands:
        qq = q if strand == "+" else rc(q)
        k = min(maxmm + 1, max(1, L // 4))
        for off, ch in _chunks(qq, k):
            if any(c not in "ACGT" for c in ch):
                continue
            i = s.find(ch)
            while i != -1:
                st = i - off
                if 0 <= st <= len(s) - L and (st, strand) not in hits:
                    m = hamming(qq, s[st:st + L], maxmm)
                    if m <= maxmm:
                        hits[(st, strand)] = m
                i = s.find(ch, i + 1)
    return sorted([(p, st, m) for (p, st), m in hits.items()], key=lambda x: (x[2], x[0]))


def best_hit(q, s, maxmm=6):
    h = find_approx(q, s, maxmm)
    return h[0] if h else None


def min_mismatch(q, seqs, limit=6):
    """Menor número de erros de q contra um conjunto de sequências.
    Valores acima de 'limit' (padrão 6) são devolvidos como limit+1 ('mais que 6')."""
    best = limit + 1
    for s in seqs:
        for mm in range(0, min(best, limit + 1)):
            if find_approx(q, s, mm):
                best = mm
                break
        if best == 0:
            return 0
    return best


def orient_and_trim(region, seq, pad=300, anchor_len=30):
    """Coloca 'seq' na mesma fita da região de referência e recorta em volta dela.
    Usa 3 âncoras (início, meio, fim) para achar a região mesmo em genomas inteiros."""
    n = len(region)
    anchors = [region[n // 2 - 15:n // 2 + 15], region[min(150, n - anchor_len):min(150, n - anchor_len) + anchor_len],
               region[max(0, n - 180):max(0, n - 180) + anchor_len]]
    for a in anchors:
        h = best_hit(a, seq, 4)
        if h:
            pos, strand, _ = h
            y = seq if strand == "+" else rc(seq)
            p = pos if strand == "+" else len(seq) - pos - len(a)
            ra = region.find(a)
            return y[max(0, p - ra - pad): p - ra + n + pad]
    return None


def align_window(ref_win, seq, F, core=15, chunk=20, slack=6):
    """Alinha uma amostra à janela de referência (sem gaps longos).

    1) acha o núcleo de 31 nt em volta do alvo; 2) estende para os lados em
    pedaços de 20 nt, permitindo deslocamento de até ±6 nt (pequenas inserções/
    deleções). Retorna a lista de bases alinhadas ('-' = não alinhado)."""
    cq = ref_win[F - core:F + core + 1]
    h = best_hit(cq, seq, 6)
    if not h:
        return None
    pos, strand, mm = h
    s = seq if strand == "+" else rc(seq)
    if strand == "-":
        pos = len(seq) - pos - len(cq)
    out = ["-"] * len(ref_win)
    for j in range(len(cq)):
        out[F - core + j] = s[pos + j]
    base0 = pos - (F - core)
    # direita
    off = base0
    c = F + core + 1
    while c < len(ref_win):
        ch = ref_win[c:min(c + chunk, len(ref_win))]
        best = None
        for d in range(-slack, slack + 1):
            i = off + c + d
            if i < 0 or i + len(ch) > len(s):
                continue
            m = hamming(ch, s[i:i + len(ch)])
            if best is None or m < best[0]:
                best = (m, d)
        if best is None or best[0] > max(5, len(ch) * 0.35):
            break
        off += best[1]
        for j in range(len(ch)):
            out[c + j] = s[off + c + j]
        c += chunk
    # esquerda
    off = base0
    c = F - core
    while c > 0:
        st = max(0, c - chunk)
        ch = ref_win[st:c]
        best = None
        for d in range(-slack, slack + 1):
            i = off + st + d
            if i < 0 or i + len(ch) > len(s):
                continue
            m = hamming(ch, s[i:i + len(ch)])
            if best is None or m < best[0]:
                best = (m, d)
        if best is None or best[0] > max(5, len(ch) * 0.35):
            break
        off += best[1]
        for j in range(len(ch)):
            out[st + j] = s[off + st + j]
        c -= chunk
    return out
