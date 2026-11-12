"""Arquivos de saída: FASTA das regiões, entrada do ADS, amplicons e planilha (CSV)."""
import csv
import os

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
    fa, ads, ads_pl, amps = [], [], [], []
    rows = []
    for m in markers:
        i = m["id"]
        p = primers.get(i)
        head = f">{i} | {m.get('organism', '')} | {m.get('gene', '')} | {m.get('acc', '')}:{m.get('pos', '')}"
        fa.append(head + "\n" + region_with_brackets(m))
        ads.append(f"{i}\t{region_with_brackets(m)}")
        ads_pl.append(f"{i}\t{region_with_brackets(m, False)}")
        if p:
            a = m["plain"][p["Fpos"][0]:p["Rpos"][1] + 1]
            amps.append(f">{i} | {m.get('organism', '')} | {len(a)} pb\n{a}")
        rows.append([i, m.get("organism", ""), m.get("type", "Detecção"), m.get("gene", ""),
                     TAG + p["F"] if p else "", TAG + p["R"] if p else "", p["E"] if p else "",
                     f"{p['ext_wt']} / {p['ext_mut']}" if p else "", p["amp"] if p else "",
                     status.get(i, ""), info.get(i, "" if p else "Sem primers: desenhar no ADS")])
    _write(os.path.join(outdir, f"{panel}_regioes.fasta"), "\n".join(fa) + "\n")
    _write(os.path.join(outdir, f"{panel}_ADS_input_IUPAC.txt"), "\n".join(ads) + "\n")
    _write(os.path.join(outdir, f"{panel}_ADS_input_simples.txt"), "\n".join(ads_pl) + "\n")
    _write(os.path.join(outdir, f"{panel}_amplicons.fasta"), "\n".join(amps) + "\n")
    with open(os.path.join(outdir, f"{panel}_lista_final.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["Ensaio", "Patógeno", "Tipo", "Gene / região", "Primer PCR forward (5→3, com tag)",
                    "Primer PCR reverse (5→3, com tag)", "Primer de extensão (5→3)",
                    "Base estendida (normal / alternativa)", "Amplicon (pb)", "Status do teste", "Observação"])
        w.writerows(rows)
    return outdir
