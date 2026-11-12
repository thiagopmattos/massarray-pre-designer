"""Testes (rodam sem internet):  python -m unittest discover -s tests -v"""
import json
import os
import random
import sys
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from massarray.align import find_approx, hamming, orient_and_trim  # noqa: E402
from massarray.design import design, design_with_fallbacks  # noqa: E402
from massarray.insilico import amplifies, extension_base, verify_placement, evaluate  # noqa: E402
from massarray.marker import build_marker, scan_positions  # noqa: E402
from massarray.seq import rc, tm, tm_ext, iupac_code, translate, djb  # noqa: E402

EX = os.path.join(ROOT, "examples")
rnd = random.Random(42)


def _load(p):
    with open(p, encoding="utf-8") as f:
        return json.load(f)


def rand_seq(n):
    return "".join(rnd.choice("ACGT") for _ in range(n))


def mutate(s, rate):
    return "".join(rnd.choice([b for b in "ACGT" if b != c]) if rnd.random() < rate else c for c in s)


class TestSeq(unittest.TestCase):
    def test_rc_iupac(self):
        self.assertEqual(rc("ACGTRY"), "RYACGT")

    def test_tm_range(self):
        # primer conhecido do painel (Tm ~60 °C nas condições padrão)
        self.assertAlmostEqual(tm("TGCTGAACCAAGCCTTATGATC"), 60.0, delta=0.5)
        self.assertGreater(tm_ext("TTCGGCGGAGATCCT"), 47)

    def test_iupac(self):
        self.assertEqual(iupac_code("AG"), "R")
        self.assertEqual(iupac_code("CT"), "Y")

    def test_translate(self):
        self.assertEqual(translate("ATGTATAAATAA"), "MYK*")

    def test_checksum(self):
        self.assertEqual(djb("ACGT"), djb("ACGT"))
        self.assertNotEqual(djb("ACGT"), djb("ACGA"))


class TestApproxSearch(unittest.TestCase):
    def test_find_approx_matches_bruteforce(self):
        for _ in range(30):
            s = rand_seq(600)
            q = mutate(s[200:231], 0.06)
            if rnd.random() < 0.5:
                s = rc(s)
            got = {(p, st) for p, st, m in find_approx(q, s, 3)}
            exp = set()
            for st, qq in (("+", q), ("-", rc(q))):
                for i in range(len(s) - 31 + 1):
                    if hamming(qq, s[i:i + 31]) <= 3:
                        exp.add((i, st))
            self.assertEqual(got, exp)

    def test_orient_and_trim(self):
        region = rand_seq(800)
        genome = rand_seq(3000) + rc(mutate(region, 0.02)) + rand_seq(3000)
        t = orient_and_trim(region, genome, pad=50)
        self.assertIsNotNone(t)
        self.assertLess(hamming(region, t[50:850]), 40)


class TestMarkerAndDesign(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # "espécie" sintética: referência + 30 cepas com 1% de variação
        cls.ref = rand_seq(1500)
        cls.inc = [rand_seq(200) + mutate(cls.ref, 0.01) + rand_seq(200) for _ in range(30)]
        # "espécie parente": 15% diferente, mas com um trecho idêntico (300–420)
        rel = list(mutate(cls.ref, 0.15))
        rel[300:420] = cls.ref[300:420]
        cls.exc = ["".join(rel)]

    def test_build_marker(self):
        m = build_marker("X", self.ref, 700, self.inc)
        self.assertEqual(len(m["plain"]), 241)
        self.assertEqual(m["F"], 120)
        self.assertEqual(m["aligned"], 30)
        self.assertGreaterEqual(m["target_cons"], 0.9)

    def test_scan_avoids_shared_region(self):
        cands = scan_positions(self.ref, 200, 1300, self.inc, self.exc, step=10)
        self.assertTrue(cands)
        best = cands[0]["p"]
        self.assertFalse(300 <= best <= 420, "escolheu o trecho igual ao parente")
        self.assertEqual(cands[0]["close_frac"], 0.0)

    def test_design_and_checks(self):
        m = build_marker("X", self.ref, 700, self.inc)
        d = design_with_fallbacks(m, self.inc, self.exc)
        self.assertIsNotNone(d)
        self.assertEqual(verify_placement(m, d), [])
        a = amplifies(d["F"], d["R"], self.ref)
        self.assertIsNotNone(a)
        self.assertEqual(a[0], d["amp"])
        # o UEP estende a base-alvo esperada
        bases = [b for mm, b in extension_base(d["E"], self.ref, a[2], a[3]) if mm == 0]
        self.assertIn(d["ext_wt"], bases)
        ev = evaluate(d, self.inc, self.exc)
        self.assertGreaterEqual(ev["inc_ok"], 27)
        self.assertEqual(ev["exc_signal"], 0)

    def test_uep_never_degenerate(self):
        m = build_marker("X", self.ref, 700, self.inc)
        F = m["F"]
        iu = list(m["iu"])
        for k in range(F - 30, F):
            iu[k] = "R" if k % 3 == 0 else iu[k]
        m["iu"] = "".join(iu)
        d = design(m)
        if d:
            self.assertTrue(all(c in "ACGT" for c in d["E"]))


class TestExamplePanels(unittest.TestCase):
    def _check(self, folder):
        ms = _load(os.path.join(EX, folder, "markers.json"))
        pr = _load(os.path.join(EX, folder, "primers.json"))
        probs = {}
        for m in ms:
            p = pr.get(m["id"])
            if p:
                e = [x for x in verify_placement(m, p) if not x.startswith("aviso")]
                if e:
                    probs[m["id"]] = e
        return probs

    def test_vaginal_panel_ok(self):
        self.assertEqual(self._check("painel_vaginal_IST_70"), {})

    def test_vaginal_amplicons_match_primers(self):
        ms = {m["id"]: m for m in _load(os.path.join(EX, "painel_vaginal_IST_70", "markers.json"))}
        pr = _load(os.path.join(EX, "painel_vaginal_IST_70", "primers.json"))
        for i, p in pr.items():
            amp = ms[i]["plain"][p["Fpos"][0]:p["Rpos"][1] + 1]
            self.assertEqual(len(amp), p["amp"], i)


class TestCliOffline(unittest.TestCase):
    """Roda os passos markers -> primers -> check -> export com um 'cache' falso
    (sem internet), para garantir que a linha de comando funciona de ponta a ponta."""

    def test_cli(self):
        import tempfile
        from massarray import cli
        tmp = tempfile.mkdtemp()
        cwd = os.getcwd()
        os.chdir(tmp)
        try:
            refA, refB = rand_seq(1500), rand_seq(1500)
            cfg = {"panel": "teste", "markers": [
                {"id": "A_gene", "organism": "Bicho A", "gene": "gA", "expect": "Bicho A",
                 "reference": {"acc": "FAKE_A"}, "inclusivity": {"term": "x"}, "exclusivity": [{"term": "y"}]},
                {"id": "B_snp", "organism": "Bicho B", "gene": "gB", "type": "Resistência",
                 "reference": {"acc": "FAKE_B"}, "inclusivity": {"term": "x"},
                 "target": {"position": 801, "alt": "T"}}]}
            os.makedirs("config")
            with open("config/teste.json", "w") as f:
                json.dump(cfg, f)
            cache = os.path.join("work", "teste", "cache")
            os.makedirs(cache)

            def w(name, seqs):
                with open(os.path.join(cache, name + ".fasta"), "w") as f:
                    for k, x in enumerate(seqs):
                        f.write(f">s{k}\n{x}\n")
            w("A_gene__ref", [refA])
            w("A_gene__inc", [mutate(refA, 0.01) for _ in range(12)])
            w("A_gene__exc0", [mutate(refA, 0.2) for _ in range(3)])
            w("B_snp__ref", [refB])
            w("B_snp__inc", [mutate(refB, 0.01) for _ in range(12)])
            c, work = cli.load_cfg("config/teste.json")
            cli.cmd_markers(c, work)
            cli.cmd_primers(c, work)
            cli.cmd_check(c, work)
            cli.cmd_export(c, work)
            out = os.path.join(work, "saida")
            self.assertTrue(os.path.exists(os.path.join(out, "teste_ADS_input_IUPAC.txt")))
            with open(os.path.join(work, "check.json")) as f:
                ck = json.load(f)
            for i in ("A_gene", "B_snp"):
                self.assertEqual(ck[i]["posicionamento"], ["OK"], i)
        finally:
            os.chdir(cwd)


if __name__ == "__main__":
    unittest.main()
