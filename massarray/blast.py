"""BLAST no NCBI pela URL API (sem instalar nada).

Lições aprendidas:
  - Mande UMA sequência por pedido. Vários amplicons juntos (sobretudo de HPV16/18,
    que têm milhares de registros) ficaram presos na fila por mais de 1 hora;
    um por um, cada pedido levou ~30 s.
  - Use megablast para amplicons (~100 pb): é rápido e acha tudo com >~85% de identidade,
    que é o que importa para reação cruzada.
  - Não mande pedidos em paralelo: o NCBI pede no máximo 1 pedido a cada 10 s
    e checagem de status no máximo 1 vez por minuto por RID (usamos 15–20 s com cuidado).
"""
import json
import re
import time
import urllib.parse
import urllib.request

URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"

# nomes que NÃO contam como "outro organismo" no resumo (registros sem espécie)
GENERIC = re.compile(r"uncultured|unidentified|MAG:|vector|synthetic|construct|metagenome|"
                     r"bacterium clone|environmental|isolate .* genomic sequ", re.I)


def put(fasta, program="blastn", megablast=True, database="core_nt", hitlist=250, entrez_query=None):
    p = {"CMD": "Put", "PROGRAM": program, "DATABASE": database, "QUERY": fasta,
         "HITLIST_SIZE": str(hitlist), "FORMAT_TYPE": "JSON2_S", "MEGABLAST": "on" if megablast else "off"}
    if entrez_query:
        p["ENTREZ_QUERY"] = entrez_query
    data = urllib.parse.urlencode(p).encode()
    with urllib.request.urlopen(urllib.request.Request(URL, data=data), timeout=120) as r:
        t = r.read().decode()
    rid = re.search(r"RID = (\S+)", t)
    if not rid:
        raise RuntimeError("BLAST não devolveu RID")
    return rid.group(1)


def get(rid, poll=20, timeout=1800):
    t0 = time.time()
    while time.time() - t0 < timeout:
        time.sleep(poll)
        with urllib.request.urlopen(f"{URL}?CMD=Get&FORMAT_OBJECT=SearchInfo&RID={rid}", timeout=120) as r:
            info = r.read().decode()
        st = re.search(r"Status=(\w+)", info)
        st = st.group(1) if st else "?"
        if st == "READY":
            with urllib.request.urlopen(f"{URL}?CMD=Get&FORMAT_TYPE=JSON2_S&RID={rid}", timeout=300) as r:
                return json.loads(r.read().decode())
        if st in ("FAILED", "UNKNOWN"):
            raise RuntimeError(f"BLAST {rid}: {st}")
    raise TimeoutError(f"BLAST {rid} não terminou em {timeout}s")


def summarize(report, expect_regex):
    """Resume um resultado: quantos acertos são do organismo esperado e qual o
    'outro' organismo mais parecido (identidade % / cobertura %)."""
    s = report["results"]["search"]
    qlen = s["query_len"]
    hits = s.get("hits", [])
    exp = re.compile(expect_regex, re.I)
    match, others = 0, {}
    for h in hits:
        d = h["description"][0]
        title = d.get("title", "")
        sci = d.get("sciname") or " ".join(title.split()[:2])
        hs = h["hsps"][0]
        idt = round(100 * hs["identity"] / hs["align_len"], 1)
        cov = round(100 * hs["align_len"] / qlen)
        if exp.search(title + " " + sci):
            match += 1
            continue
        if GENERIC.search(title + " " + sci):
            continue
        k = sci
        if k not in others or others[k]["score"] < idt * cov:
            others[k] = {"identidade": idt, "cobertura": cov, "score": idt * cov, "n": 0}
        others[k]["n"] += 1
    top = sorted(others.items(), key=lambda kv: -kv[1]["score"])[:5]
    return {"acertos": len(hits), "do_alvo": match,
            "outros": [{"organismo": k, "identidade": v["identidade"], "cobertura": v["cobertura"], "n": v["n"]}
                       for k, v in top]}


def blast_one(name, seq, expect_regex, **kw):
    rid = put(f">{name}\n{seq}", **{k: v for k, v in kw.items() if k in ("program", "megablast", "database", "hitlist", "entrez_query")})
    js = get(rid)
    return summarize(js["BlastOutput2"][0]["report"], expect_regex)


def status_text(summary, alvo):
    """Frase curta para a planilha."""
    n, m = summary["acertos"], summary["do_alvo"]
    if n == 0:
        return "BLAST: nenhum acerto"
    txt = f"BLAST: {m}/{n} acertos são {alvo}"
    if summary["outros"]:
        o = summary["outros"][0]
        txt += f"; outro mais parecido: {o['organismo']} ({o['identidade']}% / cobertura {o['cobertura']}%)"
    return txt
