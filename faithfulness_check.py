#!/usr/bin/env python3
"""
Deterministic checks that an agent's final answer is faithful to ClawBio skill outputs.

No LLM is involved: every check is plain code, so the verdict itself is reproducible.

Subcommands
-----------
check      Verify one run directory.
stability  Compare several runs of the same input (label agreement, pass rate, citation overlap).

Expected run directory layout
-----------------------------
run_dir/
    outputs/*.json        skill outputs (must include classification.json)
    run_log.jsonl         one JSON record per line: skill calls + one agent metadata record
    final_answer.json     the agent's structured answer

See README.md for the schemas.
"""
import argparse
import hashlib
import itertools
import json
import math
import re
import sys
from collections import Counter
from pathlib import Path

PHENOTYPES = ("inflamed", "excluded", "desert")

DEFAULT_CONFIG = {
    # Skills that must appear in run_log.jsonl, in this order (other skills may be interleaved).
    "expected_order": ["load_qc", "score_cells", "spatial_stats", "classify", "report"],
    # Skills that use randomness (permutation tests) and must log a seed.
    "stochastic_skills": ["spatial_stats"],
    # Numbers allowed in free text even if no skill produced them (e.g. [95] for "95% CI").
    "allowed_constants": [],
    # Text fields of final_answer.json scanned for numbers.
    "text_fields": ["narrative", "conclusion"],
    # Terms the caveat must contain (case-insensitive, substring match).
    "required_caveat_terms": ["not", "validated", "predict"],
    # If one member of a pair is reported, the other must be too (effect size with its p-value).
    "paired_stats": [
        ["classification.immune_phenotype.nhood_z", "classification.immune_phenotype.nhood_p"],
        ["classification.biomarker.morans_i", "classification.biomarker.morans_p"],
    ],
    "ref_id_pattern": r"\bREF-\d{3}\b",
}

NUM_RE = re.compile(r"(?<![\w.\-])-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?%?")


# ----------------------------------------------------------------------------- helpers

def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def flatten(obj, prefix=""):
    """Yield (dotted_path, number) for every numeric leaf."""
    if isinstance(obj, dict):
        for k, v in obj.items():
            yield from flatten(v, f"{prefix}.{k}" if prefix else k)
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            yield from flatten(v, f"{prefix}[{i}]")
    elif isinstance(obj, bool):
        return
    elif isinstance(obj, (int, float)) and not (isinstance(obj, float) and math.isnan(obj)):
        yield prefix, float(obj)


def _sig_round(x, sig):
    if x == 0:
        return 0.0
    return round(x, sig - 1 - int(math.floor(math.log10(abs(x)))))


def number_matches(token, value):
    """True if `token` (as written in text) is a faithful rendering of `value`."""
    pct = token.endswith("%")
    s = token.rstrip("%")
    target = float(s)
    x = value * 100 if pct else value
    if "e" in s.lower():
        mant = s.lower().split("e")[0].replace("-", "").replace(".", "").lstrip("0")
        sig = max(len(mant), 1)
        return _sig_round(x, sig) == _sig_round(target, sig)
    decimals = len(s.split(".")[1]) if "." in s else 0
    return abs(round(x, decimals) - target) < 10 ** -(decimals + 6)


def extract_numbers(text, ref_pattern):
    """Return list of (token, comparator) where comparator is '<', '<=' or '='."""
    text = re.sub(ref_pattern, " ", text.replace("\u2212", "-"))
    found = []
    for m in NUM_RE.finditer(text):
        before = text[max(0, m.start() - 3):m.start()].replace(" ", "")
        comp = "<=" if before.endswith(("<=", "\u2264")) else "<" if before.endswith("<") else "="
        found.append((m.group(), comp))
    return found


def _result(name, status, detail):
    return {"check": name, "status": status, "detail": detail}


def load_run(run_dir):
    run_dir = Path(run_dir)
    outputs = {}
    for p in sorted((run_dir / "outputs").glob("*.json")):
        outputs[p.stem] = json.loads(p.read_text())
    log = []
    log_path = run_dir / "run_log.jsonl"
    if log_path.exists():
        log = [json.loads(line) for line in log_path.read_text().splitlines() if line.strip()]
    final = json.loads((run_dir / "final_answer.json").read_text())
    return outputs, log, final


# ----------------------------------------------------------------------------- checks

def check_numbers(final, outputs, cfg):
    leaves = [(f"{stem}.{path}", v) for stem, obj in outputs.items() for path, v in flatten(obj)]
    matched_paths, results = set(), []
    for field in cfg["text_fields"]:
        for token, comp in extract_numbers(final.get(field, "") or "", cfg["ref_id_pattern"]):
            if comp == "=":
                hits = [p for p, v in leaves if number_matches(token, v)]
            else:  # threshold such as "p < 0.001": some p-value-like leaf must satisfy it
                thr = float(token.rstrip("%"))
                hits = [p for p, v in leaves if p.rsplit(".", 1)[-1].endswith(("_p", "pval", "p_value"))
                        and (v < thr if comp == "<" else v <= thr)]
            if hits:
                matched_paths.update(hits)
            elif float(token.rstrip("%")) in cfg["allowed_constants"]:
                continue
            else:
                results.append(_result("numeric_faithfulness", "fail",
                                       f"{field}: '{comp if comp != '=' else ''}{token}' matches no skill output"))
    if not results:
        results.append(_result("numeric_faithfulness", "pass",
                               f"all numbers traced to skill outputs ({len(matched_paths)} leaves referenced)"))
    for a, b in cfg["paired_stats"]:
        if (a in matched_paths) != (b in matched_paths):
            present, missing = (a, b) if a in matched_paths else (b, a)
            results.append(_result("paired_statistics", "fail",
                                   f"reported {present} without {missing}"))
    return results


def check_labels(final, outputs):
    res = []
    cls = outputs.get("classification")
    if cls is None:
        return [_result("label_faithfulness", "fail", "outputs/classification.json missing")]
    pairs = [
        ("immune_phenotype", cls.get("immune_phenotype", {}).get("label")),
        ("biomarker_pattern", cls.get("biomarker", {}).get("label")),
        ("biomarker_gene", cls.get("biomarker", {}).get("gene")),
    ]
    for key, skill_val in pairs:
        if final.get(key) != skill_val:
            res.append(_result("label_faithfulness", "fail",
                                f"final_answer.{key}='{final.get(key)}' but skill says '{skill_val}'"))
    abstained = [k for k, v in pairs[:2] if v == "unclassifiable"]
    if abstained and not (final.get("abstention_reason") or "").strip():
        res.append(_result("abstention", "fail",
                           f"skill abstained on {abstained} but final answer gives no abstention_reason"))
    label = cls.get("immune_phenotype", {}).get("label")
    conclusion = (final.get("conclusion") or "").lower()
    others = [p for p in PHENOTYPES if p != label and re.search(rf"\b{p}\b", conclusion)]
    if others:
        res.append(_result("label_consistency_in_text", "warn",
                            f"conclusion mentions {others} although label is '{label}'"))
    if not any(r["status"] == "fail" for r in res):
        res.insert(0, _result("label_faithfulness", "pass", "labels match classification skill output"))
    return res


def check_citations(final, refs, cfg):
    res = []
    ref_by_id = {r["id"]: r for r in refs}
    cited = list(final.get("citations", []))
    text = " ".join(final.get(f, "") or "" for f in cfg["text_fields"])
    in_text = set(re.findall(cfg["ref_id_pattern"], text))
    for rid in sorted(in_text - set(cited)):
        res.append(_result("citation_grounding", "fail", f"{rid} appears in text but not in citations list"))
    labels = {final.get("immune_phenotype"), final.get("biomarker_pattern")}
    subtype = final.get("sample_subtype", "unknown")
    for rid in sorted(set(cited) | in_text):
        ref = ref_by_id.get(rid)
        if ref is None:
            res.append(_result("citation_grounding", "fail", f"{rid} is not in the curated reference file"))
            continue
        if not ref.get("verified", False):
            res.append(_result("citation_verified", "warn", f"{rid} is marked verified=false in the reference file"))
        if not labels & set(ref.get("phenotypes", [])):
            res.append(_result("citation_relevance", "warn",
                                f"{rid} covers {ref.get('phenotypes')} but reported labels are {sorted(l for l in labels if l)}"))
        subtypes = ref.get("subtypes", [])
        if "any" not in subtypes and subtype not in subtypes and not final.get("subtype_mismatch_acknowledged"):
            res.append(_result("subtype_match", "fail",
                                f"{rid} is from {subtypes}, sample is '{subtype}', mismatch not acknowledged"))
    if not any(r["status"] == "fail" for r in res):
        res.insert(0, _result("citation_grounding", "pass", f"{len(set(cited))} citation(s) grounded in curated file"))
    return res


def check_required(final, outputs, cfg):
    res = []
    caveat = (final.get("caveat") or "").lower()
    missing = [t for t in cfg["required_caveat_terms"] if t.lower() not in caveat]
    res.append(_result("caveat", "fail" if missing else "pass",
                       f"caveat missing terms {missing}" if missing else "caveat present"))
    qc = outputs.get("classification", {}).get("qc", {})
    n_in, n_pass, n_drop = qc.get("n_cells_input"), qc.get("n_cells_passed"), qc.get("n_dropped")
    if None in (n_in, n_pass, n_drop):
        res.append(_result("qc_reporting", "fail", "classification.qc lacks n_cells_input/n_cells_passed/n_dropped"))
    else:
        if n_in - n_pass != n_drop:
            res.append(_result("qc_internal_consistency", "fail",
                                f"skill QC inconsistent: {n_in} - {n_pass} != {n_drop}"))
        if final.get("dropped_cells_reported") != n_drop:
            res.append(_result("qc_reporting", "fail",
                                f"final answer reports {final.get('dropped_cells_reported')} dropped cells, skill says {n_drop}"))
        else:
            res.append(_result("qc_reporting", "pass", f"{n_drop} dropped cells reported"))
    return res


def check_trajectory(log, run_dir, cfg):
    res = []
    skills = [r for r in log if r.get("type") == "skill"]
    names = [r.get("skill") for r in skills]
    pos = -1
    for name in cfg["expected_order"]:
        idx = [i for i, n in enumerate(names) if n == name]
        if not idx:
            res.append(_result("trajectory_order", "fail", f"required skill '{name}' never called"))
        elif idx[0] < pos:
            res.append(_result("trajectory_order", "fail", f"'{name}' first called before its predecessor"))
        else:
            pos = idx[0]
    produced = {}
    for r in skills:
        name = r.get("skill")
        if not r.get("version"):
            res.append(_result("provenance", "fail", f"{name}: no skill version logged"))
        if name in cfg["stochastic_skills"] and r.get("seed") is None:
            res.append(_result("provenance", "fail", f"{name}: stochastic skill ran without a logged seed"))
        for path, digest in (r.get("inputs") or {}).items():
            if path in produced and produced[path] != digest:
                res.append(_result("chaining", "fail",
                                    f"{name} consumed {path} with a different checksum than was produced"))
        for path, digest in (r.get("outputs") or {}).items():
            produced[path] = digest
            f = Path(run_dir) / path
            if not f.exists():
                res.append(_result("output_checksums", "warn", f"{path} logged but not present"))
            elif sha256(f) != digest:
                res.append(_result("output_checksums", "fail", f"{path} changed after the skill wrote it"))
    agent = [r for r in log if r.get("type") == "agent"]
    if not agent or not all(agent[0].get(k) is not None for k in ("model", "prompt_version", "temperature")):
        res.append(_result("agent_metadata", "warn", "agent record lacks model / prompt_version / temperature"))
    if not any(r["status"] == "fail" for r in res):
        res.insert(0, _result("trajectory", "pass", f"{len(skills)} skill calls in expected order with provenance"))
    return res


def check_truth(final, outputs, truth):
    cls = outputs.get("classification", {})
    got = {"immune": cls.get("immune_phenotype", {}).get("label"),
           "pattern": cls.get("biomarker", {}).get("label")}
    exp = {"immune": truth["expected_immune_phenotype"], "pattern": truth["expected_biomarker_pattern"]}
    res = []
    for k in got:
        ok = got[k] == exp[k]
        res.append(_result(f"ground_truth_{k}", "pass" if ok else "fail",
                            f"expected '{exp[k]}', skill returned '{got[k]}'"))
    return res


def run_checks(run_dir, refs, cfg, truth=None):
    outputs, log, final = load_run(run_dir)
    results = []
    results += check_labels(final, outputs)
    results += check_numbers(final, outputs, cfg)
    results += check_citations(final, refs, cfg)
    results += check_required(final, outputs, cfg)
    results += check_trajectory(log, run_dir, cfg)
    if truth is not None:
        results += check_truth(final, outputs, truth)
    status = "fail" if any(r["status"] == "fail" for r in results) else "pass"
    return {"run_dir": str(run_dir), "status": status,
            "n_fail": sum(r["status"] == "fail" for r in results),
            "n_warn": sum(r["status"] == "warn" for r in results),
            "results": results}


# ----------------------------------------------------------------------------- stability

def stability(run_dirs, refs, cfg):
    reports, finals = [], []
    for d in run_dirs:
        reports.append(run_checks(d, refs, cfg))
        finals.append(load_run(d)[2])

    def agreement(key):
        vals = [f.get(key) for f in finals]
        mode, count = Counter(vals).most_common(1)[0]
        return {"mode": mode, "agreement": count / len(vals), "values": vals}

    cites = [set(f.get("citations", [])) for f in finals]
    jacc = [len(a & b) / len(a | b) if (a | b) else 1.0 for a, b in itertools.combinations(cites, 2)]
    return {
        "n_runs": len(run_dirs),
        "immune_phenotype": agreement("immune_phenotype"),
        "biomarker_pattern": agreement("biomarker_pattern"),
        "faithfulness_pass_rate": sum(r["status"] == "pass" for r in reports) / len(reports),
        "citation_jaccard_mean": sum(jacc) / len(jacc) if jacc else 1.0,
        "citation_jaccard_min": min(jacc) if jacc else 1.0,
        "per_run_status": [r["status"] for r in reports],
    }


# ----------------------------------------------------------------------------- CLI

def _print_report(rep):
    icon = {"pass": "PASS", "fail": "FAIL", "warn": "WARN"}
    print(f"\n{rep['run_dir']}: {rep['status'].upper()}  ({rep['n_fail']} fail, {rep['n_warn']} warn)")
    for r in rep["results"]:
        print(f"  [{icon[r['status']]}] {r['check']}: {r['detail']}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("check")
    c.add_argument("--run-dir", required=True)
    c.add_argument("--refs", required=True, help="curated reference JSON")
    c.add_argument("--config", help="JSON overriding DEFAULT_CONFIG keys")
    c.add_argument("--truth", help="ground_truth.json from synthetic_tissue.py")
    c.add_argument("--out", help="write report JSON here")
    s = sub.add_parser("stability")
    s.add_argument("--runs", nargs="+", required=True)
    s.add_argument("--refs", required=True)
    s.add_argument("--config")
    s.add_argument("--out")
    args = ap.parse_args()

    cfg = dict(DEFAULT_CONFIG)
    if args.config:
        cfg.update(json.loads(Path(args.config).read_text()))
    refs = json.loads(Path(args.refs).read_text())

    if args.cmd == "check":
        truth = json.loads(Path(args.truth).read_text()) if args.truth else None
        rep = run_checks(args.run_dir, refs, cfg, truth)
        _print_report(rep)
        if args.out:
            Path(args.out).write_text(json.dumps(rep, indent=2))
        sys.exit(0 if rep["status"] == "pass" else 1)
    else:
        rep = stability(args.runs, refs, cfg)
        print(json.dumps(rep, indent=2))
        if args.out:
            Path(args.out).write_text(json.dumps(rep, indent=2))
        ok = (rep["immune_phenotype"]["agreement"] == 1.0 and rep["biomarker_pattern"]["agreement"] == 1.0
              and rep["faithfulness_pass_rate"] == 1.0)
        sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
