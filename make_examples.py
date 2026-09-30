#!/usr/bin/env python3
"""Build example runs for faithfulness_check.py: one faithful, one with planted errors.

Also writes curated_references.template.json. Its entries are starting points only:
every entry is verified=false until someone on the team checks the paper and the claim.
"""
import json
from pathlib import Path

from faithfulness_check import sha256

HERE = Path(__file__).parent / "examples"

REFS = [
    {"id": "REF-001", "verified": False, "subtypes": ["any"],
     "phenotypes": ["inflamed", "excluded", "desert"],
     "citation": "Chen DS, Mellman I. Elements of cancer immunity and the cancer-immune set point. Nature 2017.",
     "claim": "Defines inflamed, excluded and desert immune phenotypes; inflamed tumours are generally "
              "more likely to respond to PD-1/PD-L1 blockade."},
    {"id": "REF-002", "verified": False, "subtypes": ["TNBC"], "phenotypes": ["inflamed", "excluded"],
     "citation": "Keren L et al. A structured tumor-immune microenvironment in triple negative breast cancer "
                 "revealed by multiplexed ion beam imaging. Cell 2018.",
     "claim": "Compartmentalised vs mixed tumour-immune architecture associated with survival in TNBC."},
    {"id": "REF-003", "verified": False, "subtypes": ["TNBC", "HER2+"], "phenotypes": ["inflamed", "desert"],
     "citation": "Salgado R et al. The evaluation of tumor-infiltrating lymphocytes (TILs) in breast cancer: "
                 "recommendations by an International TILs Working Group 2014. Ann Oncol 2015.",
     "claim": "Standardised stromal TIL scoring; higher TILs associated with better prognosis and "
              "neoadjuvant response in TNBC and HER2+ disease."},
    {"id": "REF-004", "verified": False, "subtypes": ["HER2+"], "phenotypes": ["heterogeneous"],
     "citation": "Filho OM et al. Impact of HER2 heterogeneity on treatment response of early-stage "
                 "HER2-positive breast cancer: phase II neoadjuvant clinical trial of T-DM1 combined "
                 "with pertuzumab. Cancer Discov 2021.",
     "claim": "HER2 heterogeneity associated with lower pCR to neoadjuvant T-DM1 plus pertuzumab."},
]

CLASSIFICATION = {
    "immune_phenotype": {"label": "excluded", "reason": None, "nhood_z": -4.21, "nhood_p": 0.001,
                         "t_cell_fraction": 0.083, "infiltration_fraction": 0.012, "infiltration_ratio": 0.108,
                         "n_permutations": 1000},
    "biomarker": {"gene": "ERBB2", "label": "heterogeneous", "morans_i": 0.412, "morans_p": 0.001,
                  "n_tumor_cells": 8123},
    "qc": {"n_cells_input": 25000, "n_cells_passed": 24312, "n_dropped": 688,
           "drop_reasons": {"low_umi": 601, "low_genes": 87}},
}
QC = {"n_cells_input": 25000, "n_cells_passed": 24312, "n_dropped": 688}
SCORES = {"marker_concordance_spearman": {"EPCAM": 0.81, "KRT8": 0.77, "CD8A": 0.64, "PTPRC": 0.72},
          "concordance_ok": True}
SPATIAL = {"nhood_z": -4.21, "nhood_p": 0.001, "morans_i": 0.412, "morans_p": 0.001}

CAVEAT = ("This is pattern classification of one tissue section against known categories. It is not a "
          "validated per-patient prediction of treatment response.")

GOOD = {
    "immune_phenotype": "excluded", "biomarker_pattern": "heterogeneous", "biomarker_gene": "ERBB2",
    "sample_subtype": "HER2+", "subtype_mismatch_acknowledged": False,
    "citations": ["REF-001", "REF-004"], "dropped_cells_reported": 688, "caveat": CAVEAT,
    "abstention_reason": None,
    "conclusion": "Immune-excluded tumour with spatially heterogeneous ERBB2 expression.",
    "narrative": ("The section is classified as immune-excluded: tumour-immune neighbourhood enrichment "
                  "z = -4.21 (p = 0.001, 1000 permutations); immune cells make up 8.3% of cells, of which "
                  "1.2% lie within tumour regions. ERBB2 is spatially heterogeneous among 8123 tumour cells "
                  "(Moran's I = 0.412, p = 0.001). 688 of 25000 cells were removed at QC. The immune "
                  "phenotype framework is described in REF-001; ERBB2 heterogeneity has been associated "
                  "with lower pCR to HER2-targeted neoadjuvant therapy (REF-004)."),
}

BAD = dict(GOOD)
BAD.update({
    "immune_phenotype": "inflamed",                       # label upgraded by the agent
    "citations": ["REF-001", "REF-002", "REF-999"],        # TNBC ref unacknowledged + invented ref
    "dropped_cells_reported": 500,                         # wrong QC count
    "caveat": "Interpret with care.",                      # caveat missing required terms
    "conclusion": "Immune-inflamed tumour with heterogeneous ERBB2.",
    "narrative": ("Tumour-immune enrichment z = -3.9 suggests infiltration (REF-002). "   # fabricated z
                  "ERBB2 heterogeneity p = 0.001 (REF-999)."),                         # p without Moran's I
})


def write_run(name, final, drop_seed=False):
    d = HERE / name
    (d / "outputs").mkdir(parents=True, exist_ok=True)
    files = {"qc": QC, "scores": SCORES, "spatial_stats": SPATIAL, "classification": CLASSIFICATION}
    for stem, obj in files.items():
        (d / "outputs" / f"{stem}.json").write_text(json.dumps(obj, indent=2))
    (d / "outputs" / "report.md").write_text("# report placeholder\n")
    s = {k: sha256(d / "outputs" / f"{k}.json") for k in files}
    s["report"] = sha256(d / "outputs" / "report.md")
    log = [
        {"type": "agent", "model": "MODEL-ID-HERE", "prompt_version": "orchestrator-v1", "temperature": 0},
        {"type": "skill", "skill": "load_qc", "version": "0.1.0",
         "inputs": {"data/tissue.h5ad": "<sha256 of input>"}, "outputs": {"outputs/qc.json": s["qc"]}},
        {"type": "skill", "skill": "score_cells", "version": "0.1.0",
         "inputs": {"outputs/qc.json": s["qc"]}, "outputs": {"outputs/scores.json": s["scores"]}},
        {"type": "skill", "skill": "spatial_stats", "version": "0.1.0", "seed": None if drop_seed else 42,
         "inputs": {"outputs/scores.json": s["scores"]},
         "outputs": {"outputs/spatial_stats.json": s["spatial_stats"]}},
        {"type": "skill", "skill": "classify", "version": "0.1.0",
         "inputs": {"outputs/spatial_stats.json": s["spatial_stats"], "outputs/qc.json": s["qc"]},
         "outputs": {"outputs/classification.json": s["classification"]}},
        {"type": "skill", "skill": "report", "version": "0.1.0",
         "inputs": {"outputs/classification.json": s["classification"]},
         "outputs": {"outputs/report.md": s["report"]}},
    ]
    (d / "run_log.jsonl").write_text("\n".join(json.dumps(r) for r in log) + "\n")
    (d / "final_answer.json").write_text(json.dumps(final, indent=2))


if __name__ == "__main__":
    HERE.mkdir(exist_ok=True)
    (HERE / "curated_references.template.json").write_text(json.dumps(REFS, indent=2))
    write_run("run_good", GOOD)
    write_run("run_bad", BAD, drop_seed=True)
    print(f"examples written to {HERE}")
