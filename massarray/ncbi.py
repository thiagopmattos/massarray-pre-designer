"""Acesso ao NCBI (E-utilities) com cache em disco.

Regras do NCBI: no máximo 3 pedidos/segundo sem chave (10 com API key) e
informar um e-mail. Configure com variáveis de ambiente:
    NCBI_EMAIL=voce@exemplo.com
    NCBI_API_KEY=xxxxxxxx   (opcional, crie em https://www.ncbi.nlm.nih.gov/account/)
"""
import json
import os
import re
import time
import urllib.parse
import urllib.request

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
_last = [0.0]


def _params(extra):
    p = dict(extra)
    if os.environ.get("NCBI_EMAIL"):
        p["email"] = os.environ["NCBI_EMAIL"]
    if os.environ.get("NCBI_API_KEY"):
        p["api_key"] = os.environ["NCBI_API_KEY"]
    p["tool"] = "massarray-designer"
    return p


def _get(endpoint, params, tries=4):
    wait = 0.12 if os.environ.get("NCBI_API_KEY") else 0.4
    url = EUTILS + endpoint + "?" + urllib.parse.urlencode(_params(params))
    for i in range(tries):
        dt = time.time() - _last[0]
        if dt < wait:
            time.sleep(wait - dt)
        _last[0] = time.time()
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                return r.read().decode()
        except Exception as e:  # rede instável / limite do NCBI
            if i == tries - 1:
                raise RuntimeError(f"falha ao acessar NCBI: {url}\n{e}")
            time.sleep(2 * (i + 1))


def parse_fasta(text):
    recs = []
    for block in text.split(">"):
        block = block.strip()
        if not block:
            continue
        lines = block.splitlines()
        recs.append({"title": lines[0], "seq": "".join(lines[1:]).upper()})
    return recs


def esearch(term, retmax=40, db="nuccore"):
    t = _get("esearch.fcgi", {"db": db, "term": term, "retmax": retmax, "retmode": "json"})
    return json.loads(t)["esearchresult"]["idlist"]


def efetch_fasta(ids, start=None, stop=None, db="nuccore"):
    """Baixa sequências em FASTA. start/stop (1-based) recortam a região."""
    if isinstance(ids, str):
        ids = [ids]
    out = []
    for i in range(0, len(ids), 50):
        p = {"db": db, "id": ",".join(ids[i:i + 50]), "rettype": "fasta", "retmode": "text"}
        if start:
            p.update(seq_start=start, seq_stop=stop)
        out += parse_fasta(_get("efetch.fcgi", p))
    return out


def esummary(ids, db="nuccore"):
    t = json.loads(_get("esummary.fcgi", {"db": db, "id": ",".join(ids), "retmode": "json"}))["result"]
    return [{"id": i, "acc": t[i].get("accessionversion"), "title": t[i].get("title"),
             "len": t[i].get("slen")} for i in ids if i in t]


def fetch_set(term, retmax=40, max_n_frac=0.02, drop_unverified=True):
    """Busca + baixa um conjunto de sequências e limpa as ruins.

    - Remove registros 'UNVERIFIED' (lição aprendida: alguns 'HPV68' eram HPV39).
    - Remove sequências com mais de 2% de N.
    """
    ids = esearch(term, retmax)
    if not ids:
        return []
    recs = efetch_fasta(ids)
    good = []
    for r in recs:
        if drop_unverified and re.search(r"UNVERIFIED", r["title"], re.I):
            continue
        if not r["seq"] or r["seq"].count("N") > max_n_frac * len(r["seq"]):
            continue
        good.append(r)
    return good


class Cache:
    """Guarda tudo o que foi baixado em uma pasta, para não baixar de novo
    e para o trabalho ser reprodutível (os FASTA ficam salvos)."""

    def __init__(self, folder):
        self.folder = folder
        os.makedirs(folder, exist_ok=True)

    def path(self, name):
        return os.path.join(self.folder, re.sub(r"[^\w.-]", "_", name) + ".fasta")

    def get(self, name, loader):
        p = self.path(name)
        if os.path.exists(p):
            with open(p) as f:
                return parse_fasta(f.read())
        recs = loader()
        with open(p, "w") as f:
            for r in recs:
                f.write(f">{r['title']}\n{r['seq']}\n")
        return recs
