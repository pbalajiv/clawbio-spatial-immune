#!/usr/bin/env python3
"""
T8, part 2: 5-10 repeated pipeline runs on real data, then faithfulness_check.py stability.

Varies the permutation-test seed across runs (0..N-1) -- the point is to check whether the
*label* (immune phenotype, biomarker pattern) is robust to different random permutation draws,
even though the exact z/p numbers will wobble slightly run to run. Label agreement should be
1.0: the decision rule's thresholds are supposed to be robust to permutation noise, not to
happen to land on the right side of a boundary by luck of the seed.

final_answer.json for each run is drafted mechanically (orchestrator.draft_final_answer_template),
not hand-authored per run -- see that function's docstring for why that's fine here specifically.

Usage:
    python scripts/run_stability_demo.py --n-runs 10 --out-dir runs/t8_stability
"""
import argparse
import json
import sys
from pathlib import Path

HERE = Path(__file__).parent.parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import orchestrator as orch


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", default=str(HERE / "data" / "tissue.h5ad"))
    ap.add_argument("--config", default=str(HERE / "data" / "real_data_config.json"))
    ap.add_argument("--n-runs", type=int, default=10)
    ap.add_argument("--out-dir", default=str(HERE / "runs" / "t8_stability"))
    ap.add_argument("--refs", default=str(HERE / "curated_references.json"))
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    run_dirs = []
    for seed in range(args.n_runs):
        run_dir = out_dir / f"run_{seed}"
        print(f"\n=== run {seed} ===")
        result = orch.run(args.input, run_dir, config_path=args.config, seed=seed)
        final = orch.draft_final_answer_template(result)
        (run_dir / "final_answer.json").write_text(json.dumps(final, indent=2))
        run_dirs.append(str(run_dir))

    import faithfulness_check as fc
    refs = json.loads(Path(args.refs).read_text())
    report = fc.stability(run_dirs, refs, fc.DEFAULT_CONFIG)
    print("\n" + "=" * 60)
    print("STABILITY REPORT")
    print("=" * 60)
    print(json.dumps(report, indent=2))

    stable = (report["immune_phenotype"]["agreement"] == 1.0
             and report["biomarker_pattern"]["agreement"] == 1.0
             and report["faithfulness_pass_rate"] == 1.0)
    (out_dir / "stability_report.json").write_text(json.dumps(report, indent=2))
    print(f"\nwrote {out_dir / 'stability_report.json'}")
    raise SystemExit(0 if stable else 1)


if __name__ == "__main__":
    main()
