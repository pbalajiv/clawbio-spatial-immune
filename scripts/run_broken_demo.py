#!/usr/bin/env python3
"""
T8, part 3: a deliberately broken run, to show faithfulness_check.py actually catches it.

Starts from a real, passing run (runs/t6_real) and plants six independent errors, each
targeting a different check:
  1. immune_phenotype label silently "upgraded"           -> label_faithfulness
  2. an invented reference (REF-999) cited                -> citation_grounding
  3. REF-005 (DCIS-only) cited without acknowledging the
     subtype mismatch this time                           -> subtype_match
  4. dropped_cells_reported doesn't match qc.n_dropped     -> qc_reporting
  5. caveat weakened, missing required terms               -> caveat
  6. a fabricated z-score not matching any skill output    -> numeric_faithfulness
  7. outputs/classification.json edited after run_log.jsonl
     recorded its checksum                                 -> output_checksums (chaining)

Usage:
    python scripts/run_broken_demo.py --source runs/t6_real --out runs/t8_broken --refs curated_references.json
"""
import argparse
import json
import shutil
import sys
from pathlib import Path

HERE = Path(__file__).parent.parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", default=str(HERE / "runs" / "t6_real"))
    ap.add_argument("--out", default=str(HERE / "runs" / "t8_broken"))
    ap.add_argument("--refs", default=str(HERE / "curated_references.json"))
    args = ap.parse_args()

    src, out = Path(args.source), Path(args.out)
    if out.exists():
        shutil.rmtree(out)
    shutil.copytree(src, out)

    final = json.loads((out / "final_answer.json").read_text())

    # 1. label silently upgraded
    real_label = final["immune_phenotype"]
    final["immune_phenotype"] = "inflamed" if real_label != "inflamed" else "excluded"

    # 2 & 3. an invented reference, and a real one cited without acknowledging its subtype mismatch
    final["citations"] = ["REF-005", "REF-999"]
    final["subtype_mismatch_acknowledged"] = False

    # 4. wrong dropped-cell count
    final["dropped_cells_reported"] = 500

    # 5. weakened caveat (drops "validated" and "predict")
    final["caveat"] = "Interpret with care."

    # 6. fabricated statistic (does not match any skill output)
    final["narrative"] = (
        "Tumour-immune neighbourhood enrichment z = -3.87 suggests strong infiltration (REF-999). "
        f"{final['biomarker_gene']} pattern p = 0.001."
    )
    final["conclusion"] = f"Immune-{final['immune_phenotype']} tumour with notable {final['biomarker_gene']} signal."

    (out / "final_answer.json").write_text(json.dumps(final, indent=2))

    # 7. tamper the classification.json output after the fact (run_log.jsonl's checksum is now stale)
    classification_path = out / "outputs" / "classification.json"
    result = json.loads(classification_path.read_text())
    result["biomarker"]["morans_i"] = 0.999  # edited post-hoc, log still has the old hash
    classification_path.write_text(json.dumps(result, indent=2))

    print(f"Planted 7 errors in {out}. Running the gate...\n")
    import orchestrator as orch
    report = orch.gate(out, args.refs)

    expected_checks = {
        "label_faithfulness", "citation_grounding", "subtype_match", "qc_reporting",
        "caveat", "numeric_faithfulness", "output_checksums",
    }
    caught = {r["check"] for r in report["results"] if r["status"] == "fail"}
    print(f"\nExpected failing check categories: {sorted(expected_checks)}")
    print(f"Actually caught:                    {sorted(caught)}")
    print(f"All planted errors caught: {expected_checks <= caught}")
    raise SystemExit(0 if report["status"] == "fail" else 1)  # exit 0 means "the gate correctly failed it"


if __name__ == "__main__":
    main()
