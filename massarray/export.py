"""Arquivos de saída: entrada do ADS (primers em maiúsculo), tabela TSV da Agena,
FASTA das regiões, amplicons e planilha (CSV)."""
import csv
import os

from .agena import plex, write_ads_input, write_assay_tsv
from .seq import TAG


def _write(path, text):
    with open(path, "w", encoding="utf-8") as f:
        f.write(text)


def region_with_brackets(m, use_iupac=True):
    s = m["iu"] if use_iupac else m["plain"]
    F = m["F"]
    alt = m.get("alt") or {"A": "G", "G": "A", "C": "T", "T": "C"}[m["target"]]
    return s[:F] + f"[{m['target']}/{alt}]" + s[F + 1:]


def write_all(markers, primers, outdir, panel="painel", info=None, status=None):
    """markers: lista de dicts (build_marker + metadados id/organism/gene/acc/pos/type)
    primers: {id: assay}; info: {id: observação}; status: {id: texto do teste}."""
    os.makedirs(outdir, exist_ok=True)
    info = info or {}
    status = status or {}
    fa, amps = [], []
    rows = []
    for m in markers:
        i = m["id"]
        p = primers.get(i)
        head = f">{i} | {m.get('organism', '')} | {m.get('gene', '')} | {m.get('acc', '')}:{m.get('pos', '')}"
        fa.append(head + "\n" + region_with_brackets(m))
        if p:
            a = m["plain"][p["Fpos"][0]:p["Rpos"][1] + 1]
            amps.append(f">{i} | {m.get('organism', '')} | {len(a)} pb\n{a}")
        rows.append([i, m.get("organism", ""), m.get("type", "Detecção"), m.get("gene", ""),
                     TAG + p["F"] if p else "", TAG + p["R"] if p else "", p["E"] if p else "",
                     f"{p['ext_wt']} / {p['ext_mut']}" if p else "", p["amp"] if p else "",
                     status.get(i, ""), info.get(i, "" if p else "Sem primers: desenhar no ADS")])
    _write(os.path.join(outdir, f"{panel}_regioes.fasta"), "\n".join(fa) + "\n")
    # entrada do ADS: primers de PCR em MAIÚSCULO, resto (inclusive o lugar do UEP) em minúsculo
    write_ads_input(os.path.join(outdir, f"{panel}_ADS_entrada.txt"), markers, primers)
    # tabela de ensaios no layout da Agena (TSV), com poços sugeridos
    good = {k: v for k, v in primers.items() if v}
    write_assay_tsv(os.path.join(outdir, f"{panel}_ensaios_agena.tsv"), markers, good, plex(markers, good))
    _write(os.path.join(outdir, f"{panel}_amplicons.fasta"), "\n".join(amps) + "\n")
    with open(os.path.join(outdir, f"{panel}_lista_final.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["Ensaio", "Patógeno", "Tipo", "Gene / região", "Primer PCR forward (5→3, com tag)",
                    "Primer PCR reverse (5→3, com tag)", "Primer de extensão (5→3)",
                    "Base estendida (normal / alternativa)", "Amplicon (pb)", "Status do teste", "Observação"])
        w.writerows(rows)
    return outdir
