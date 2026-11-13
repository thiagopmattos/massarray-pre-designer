"""Linha de comando.

Uso (dentro da pasta do repositório):
    python -m massarray all     config/meu_painel.json         # faz tudo
    python -m massarray fetch   config/meu_painel.json         # 1. baixa sequências do NCBI
    python -m massarray markers config/meu_painel.json         # 2. escolhe o alvo e monta as janelas
    python -m massarray primers config/meu_painel.json         # 3. desenha os primers
    python -m massarray check   config/meu_painel.json         # 4. testes in silico
    python -m massarray blast   config/meu_painel.json         # 5. BLAST dos amplicons (lento)
    python -m massarray export  config/meu_painel.json         # 6. FASTA, ADS e planilha
    python -m massarray verify  markers.json primers.json      # confere um painel já pronto (offline)
    python -m massarray agena   markers.json primers.json saida/   # gera os arquivos da Agena de um painel pronto

Tudo fica em work/<nome_do_painel>/ (cache do NCBI, JSONs intermediários e saídas).
"""
import json
import os
import sys

from . import blast as bl
from . import ncbi
from .align import orient_and_trim
from .design import design_with_fallbacks
from .export import write_all
from .insilico import evaluate, verify_placement
from .marker import build_marker, scan_positions
from .seq import translate


def load_cfg(path):
    with open(path, encoding="utf-8") as f:
        cfg = json.load(f)
    cfg.setdefault("panel", os.path.splitext(os.path.basename(path))[0])
    work = os.path.join("work", cfg["panel"])
    os.makedirs(work, exist_ok=True)
    return cfg, work


def _J(work, name, data=None):
    p = os.path.join(work, name)
    if data is None:
        if not os.path.exists(p):
            return {}
        with open(p, encoding="utf-8") as f:
            return json.load(f)
    with open(p, "w", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False, indent=1)
    return data


def log(*a):
    print(*a, flush=True)


# ---------------------------------------------------------------- 1. fetch
def cmd_fetch(cfg, work):
    cache = ncbi.Cache(os.path.join(work, "cache"))
    for m in cfg["markers"]:
        r = m["reference"]
        cache.get(m["id"] + "__ref", lambda: ncbi.efetch_fasta(r["acc"], r.get("start"), r.get("end")))
        if m.get("inclusivity"):
            q = m["inclusivity"]
            cache.get(m["id"] + "__inc", lambda: ncbi.fetch_set(q["term"], q.get("max", 40)))
        for k, q in enumerate(m.get("exclusivity", [])):
            cache.get(f"{m['id']}__exc{k}", lambda: ncbi.fetch_set(q["term"], q.get("max", 40)))
        log("baixado:", m["id"])


def _seqs(work, m):
    cache = ncbi.Cache(os.path.join(work, "cache"))
    ref = cache.get(m["id"] + "__ref", lambda: [])[0]["seq"]
    inc = [r["seq"] for r in cache.get(m["id"] + "__inc", lambda: [])] if m.get("inclusivity") else []
    exc = []
    for k in range(len(m.get("exclusivity", []))):
        exc += [r["seq"] for r in cache.get(f"{m['id']}__exc{k}", lambda: [])]
    # recorta tudo em volta da região (genomas inteiros ficariam lentos)
    inc_t = [x for x in (orient_and_trim(ref, s) for s in inc) if x]
    exc_t = [x for x in (orient_and_trim(ref, s) for s in exc) if x] or exc
    return ref, inc_t, exc_t


# ---------------------------------------------------------------- 2. markers
def cmd_markers(cfg, work):
    out = _J(work, "markers.json")
    for m in cfg["markers"]:
        ref, inc, exc = _seqs(work, m)
        t = m.get("target", {})
        if t.get("position"):
            p = int(t["position"]) - 1           # 1-based dentro da região baixada
        else:
            a, b = t.get("scan", [1, len(ref)])
            cand = scan_positions(ref, a - 1, b - 1, inc, exc)
            if not cand:
                log("SEM POSIÇÃO BOA:", m["id"], "(tente outra região/gene)")
                continue
            p = cand[0]["p"]
            log(m["id"], "candidatas:", cand)
        mk = build_marker(m["id"], ref, p, inc)
        for off, code in (t.get("mask") or {}).items():
            k = mk["F"] + int(off)
            mk["iu"] = mk["iu"][:k] + code + mk["iu"][k + 1:]
        mk.update(organism=m.get("organism", ""), gene=m.get("gene", ""), type=m.get("type", "Detecção"),
                  acc=m["reference"]["acc"], pos=(m["reference"].get("start") or 1) + p, alt=t.get("alt"))
        if t.get("codon"):
            c = t["codon"]   # {"cds_start": 1-based na região, "strand": "+", "aa": 132, "expect": "Y"}
            from .seq import rc
            cds = ref[c["cds_start"] - 1:] if c.get("strand", "+") == "+" else rc(ref[:c["cds_start"]])
            aa = translate(cds)[c["aa"] - 1] if len(cds) >= 3 * c["aa"] else "?"
            mk["codon_check"] = f"aa {c['aa']} = {aa} (esperado {c.get('expect')})"
            log(m["id"], mk["codon_check"])
        out[m["id"]] = mk
        log(f"{m['id']}: alvo {mk['target']} | conservação {mk['target_cons']} | "
            f"{mk['aligned']}/{mk['n']} alinhadas | {mk['n_var']} posições variáveis")
    _J(work, "markers.json", out)


# ---------------------------------------------------------------- 3. primers
def cmd_primers(cfg, work):
    mk = _J(work, "markers.json")
    out = _J(work, "primers.json")
    for m in cfg["markers"]:
        if m["id"] not in mk:
            continue
        ref, inc, exc = _seqs(work, m)
        d = design_with_fallbacks(mk[m["id"]], inc, exc)
        out[m["id"]] = d
        log(m["id"], "->", f"F {d['F']} R {d['R']} E {d['E']} ({d['amp']} pb, nível: {d['nivel']})" if d else "FALHOU")
    _J(work, "primers.json", out)


# ---------------------------------------------------------------- 4. check
def cmd_check(cfg, work):
    mk, pr = _J(work, "markers.json"), _J(work, "primers.json")
    sets = {}
    for m in cfg["markers"]:
        if m["id"] in mk:
            sets[m["id"]] = _seqs(work, m)
    res = {}
    for m in cfg["markers"]:
        i = m["id"]
        if not pr.get(i):
            continue
        ref, inc, exc = sets[i]
        panel = {o: s[1] for o, s in sets.items() if o != i}
        ev = evaluate(pr[i], inc, exc, panel)
        ev["posicionamento"] = verify_placement(mk[i], pr[i]) or ["OK"]
        res[i] = ev
        log(i, f"inclusividade {ev['inc_ok']}/{ev['inc_n']}",
            f"| não-alvo com sinal {ev['exc_signal']}/{ev['exc_n']}",
            f"| painel {ev['panel'] or '-'} | {'; '.join(ev['posicionamento'])}")
    _J(work, "check.json", res)


# ---------------------------------------------------------------- 5. blast
def cmd_blast(cfg, work):
    mk, pr = _J(work, "markers.json"), _J(work, "primers.json")
    res = _J(work, "blast.json")
    for m in cfg["markers"]:
        i = m["id"]
        if i in res or not pr.get(i):
            continue
        p = pr[i]
        amp = mk[i]["plain"][p["Fpos"][0]:p["Rpos"][1] + 1]
        log("BLAST", i, "...")
        try:
            res[i] = bl.blast_one(i, amp, m.get("expect", m.get("organism", ".")))
        except Exception as e:
            log("  erro:", e)
            continue
        _J(work, "blast.json", res)   # salva a cada um: pode parar e continuar depois
        log(" ", bl.status_text(res[i], m.get("organism", "")))


# ---------------------------------------------------------------- 6. export
def cmd_export(cfg, work):
    mk, pr = _J(work, "markers.json"), _J(work, "primers.json")
    ck, bj = _J(work, "check.json"), _J(work, "blast.json")
    status, info = {}, {}
    for m in cfg["markers"]:
        i = m["id"]
        st = []
        if i in bj:
            st.append(bl.status_text(bj[i], m.get("organism", "")))
        if i in ck:
            c = ck[i]
            if c["inc_n"]:
                st.append(f"in silico: detecta {c['inc_ok']} de {c['inc_n']} sequências")
            if c["exc_n"]:
                st.append(f"sinal em {c['exc_signal']} de {c['exc_n']} não-alvo")
            st.append("sem reação cruzada no painel" if not c["panel"] else f"amplifica também: {c['panel']}")
        status[i] = "; ".join(st)
        info[i] = m.get("note", "")
    ms = [dict(mk[m["id"]]) for m in cfg["markers"] if m["id"] in mk]
    out = write_all(ms, pr, os.path.join(work, "saida"), cfg["panel"], info, status)
    log("arquivos em", out)


def cmd_verify(markers_path, primers_path):
    with open(markers_path, encoding="utf-8") as f:
        ms = json.load(f)
    with open(primers_path, encoding="utf-8") as f:
        pr = json.load(f)
    if isinstance(ms, dict):
        ms = list(ms.values())
    bad = 0
    for m in ms:
        p = pr.get(m["id"])
        if not p:
            log(m["id"], ": sem primers"); bad += 1; continue
        probs = [x for x in verify_placement(m, p) if not x.startswith("aviso")]
        warn = [x for x in verify_placement(m, p) if x.startswith("aviso")]
        if probs:
            bad += 1
            log(m["id"], ": PROBLEMA ->", "; ".join(probs))
        elif warn:
            log(m["id"], ": ok (", "; ".join(warn), ")")
    log(f"{len(ms)} ensaios conferidos, {bad} com problema")
    return bad


def cmd_agena(markers_path, primers_path, outdir, panel="painel"):
    """Gera, de um painel já pronto, a entrada do ADS (primers em maiúsculo) e a
    tabela de ensaios TSV com poços sugeridos."""
    from .agena import plex, write_ads_input, write_assay_tsv
    with open(markers_path, encoding="utf-8") as f:
        ms = json.load(f)
    with open(primers_path, encoding="utf-8") as f:
        pr = json.load(f)
    if isinstance(ms, dict):
        ms = list(ms.values())
    os.makedirs(outdir, exist_ok=True)
    good = {k: v for k, v in pr.items() if v}
    wells = plex(ms, good)
    a = write_ads_input(os.path.join(outdir, f"{panel}_ADS_entrada.txt"), ms, pr)
    b = write_assay_tsv(os.path.join(outdir, f"{panel}_ensaios_agena.tsv"), ms, good, wells)
    n = len(set(wells.values()))
    log(f"{a}\n{b}\n{len(good)} ensaios, {n} poços sugeridos")


def main(argv=None):
    argv = argv or sys.argv[1:]
    if not argv or argv[0] in ("-h", "--help"):
        print(__doc__); return 0
    cmd = argv[0]
    if cmd == "verify":
        return 1 if cmd_verify(argv[1], argv[2]) else 0
    if cmd == "agena":
        cmd_agena(argv[1], argv[2], argv[3], argv[4] if len(argv) > 4 else "painel")
        return 0
    cfg, work = load_cfg(argv[1])
    steps = {"fetch": cmd_fetch, "markers": cmd_markers, "primers": cmd_primers,
             "check": cmd_check, "blast": cmd_blast, "export": cmd_export}
    if cmd == "all":
        for k in ("fetch", "markers", "primers", "check", "blast", "export"):
            log(f"\n=== {k} ===")
            steps[k](cfg, work)
    elif cmd in steps:
        steps[cmd](cfg, work)
    else:
        print("comando desconhecido:", cmd); return 2
    return 0
