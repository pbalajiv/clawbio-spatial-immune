#!/usr/bin/env python3
"""
T8, part 1: run the full synthetic ground-truth suite and confirm every label matches.

This is the standalone, human-readable counterpart to test_classify_reference.py's
test_ground_truth_labels (which asserts the same thing inside pytest). Produces the
same table shape as CLASSIFICATION_SPEC.md section 10.

Usage:
    python scripts/validate_synthetic_suite.py --out runs/t8_synthetic_suite.json
"""
import argparse
import json
import sys
from pathlib import Path

HERE = Path(__file__).parent.parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import classify_reference as cr
import synthetic_tissue as st


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-cells", type=int, default=20000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    cfg = cr.load_config()
    rows = []
    all_pass = True

    print(f"{'tissue':<32} {'expected':<28} {'actual':<28} {'match'}")
    print("-" * 100)
    for scenario, her2, extra in st.SUITE:
        name = f"{scenario}_{her2}" + ("_lowqc" if extra.get("low_quality") else "")
        tmp_dir = HERE / "runs" / "t8_synthetic_suite_tissues" / name
        truth = st.generate(scenario, her2, tmp_dir, n_cells=args.n_cells, seed=args.seed, **extra)
        X, genes, xy, prov = cr.load_input(tmp_dir)
        result = cr.classify(X, genes, xy, cfg, seed=args.seed, provenance=prov)

        exp_imm, exp_bio = truth["expected_immune_phenotype"], truth["expected_biomarker_pattern"]
        act_imm, act_bio = result["immune_phenotype"]["label"], result["biomarker"]["label"]
        match = (exp_imm == act_imm) and (exp_bio == act_bio)
        all_pass &= match
        print(f"{name:<32} {exp_imm + '/' + exp_bio:<28} {act_imm + '/' + act_bio:<28} {'OK' if match else 'MISMATCH'}")
        rows.append({"tissue": name, "expected_immune_phenotype": exp_imm, "expected_biomarker_pattern": exp_bio,
                    "actual_immune_phenotype": act_imm, "actual_biomarker_pattern": act_bio, "match": match})

    print("-" * 100)
    print(f"{sum(r['match'] for r in rows)}/{len(rows)} tissues match ground truth"
          f" ({'ALL PASS' if all_pass else 'SOME MISMATCHES'})")

    if args.out:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out).write_text(json.dumps({"all_pass": all_pass, "rows": rows}, indent=2))
        print(f"wrote {args.out}")

    raise SystemExit(0 if all_pass else 1)


if __name__ == "__main__":
    main()
