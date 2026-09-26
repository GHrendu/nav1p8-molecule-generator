from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

try:
    from rdkit import Chem
except ImportError:  # pragma: no cover
    Chem = None

from .pipeline import run_generation_pipeline
from .reinvent import prepare_reinvent_config
from .config import load_config
from .protein_rules import ProteinWordScorer, load_protein_sequence
from .site_rules import load_site_profile, score_site_c
from .scoring import score_candidate


def _generate(args: argparse.Namespace) -> int:
    result = run_generation_pipeline(args.config, args.output)
    print(json.dumps({"accepted": len(result["accepted"]), "rejected": len(result["rejected"])}, indent=2))
    return 0


def _prepare_reinvent(args: argparse.Namespace) -> int:
    result = prepare_reinvent_config(
        args.output,
        args.prior,
        device=args.device,
        num_smiles=args.num_smiles,
        project_config=args.config,
        python_path=args.python,
    )
    print(json.dumps(result, indent=2))
    return 0


def _score_stdin(args: argparse.Namespace) -> int:
    cfg = load_config(args.config) if args.config else {"filter_config": None}
    filters = cfg["filter_config"]
    protein_scorer = None
    site_profile = None
    if args.config and cfg.get("protein_sequence_file"):
        protein_scorer = ProteinWordScorer(load_protein_sequence(cfg["protein_sequence_file"]))
    if args.config and cfg.get("site_profile_file"):
        site_profile = load_site_profile(cfg["site_profile_file"])
    lines = [line.strip() for line in sys.stdin if line.strip()]
    rewards = []
    for line in lines:
        smiles = line
        try:
            parsed = json.loads(line)
            if isinstance(parsed, dict):
                smiles = str(parsed.get("smiles", parsed.get("SMILES", "")))
        except json.JSONDecodeError:
            pass
        if not smiles:
            rewards.append(0.0)
            continue
        mol = Chem.MolFromSmiles(smiles) if Chem is not None else smiles
        if Chem is not None and mol is None:
            rewards.append(0.0)
            continue
        protein_metrics = protein_scorer.score(smiles) if protein_scorer else None
        site_metrics = score_site_c(mol, site_profile) if site_profile else None
        _, metrics = score_candidate(mol, None, filters or None, protein_metrics, site_metrics)
        rewards.append(float(metrics["reward"]))
    print(json.dumps({"version": 1, "payload": {"nav18_reward": rewards}}))
    return 0


def _score_reinvent_file(args: argparse.Namespace) -> int:
    import pandas as pd

    cfg = load_config(args.config)
    filters = cfg["filter_config"]
    protein_scorer = (
        ProteinWordScorer(load_protein_sequence(cfg["protein_sequence_file"]))
        if cfg.get("protein_sequence_file")
        else None
    )
    site_profile = load_site_profile(cfg["site_profile_file"]) if cfg.get("site_profile_file") else None
    frame = pd.read_csv(args.input)
    smiles_column = args.smiles_column
    if smiles_column not in frame.columns:
        raise ValueError(f"Input file does not contain SMILES column: {smiles_column}")
    rows = []
    for value in frame[smiles_column].fillna(""):
        smiles = str(value).strip()
        mol = Chem.MolFromSmiles(smiles) if Chem is not None and smiles else None
        if mol is None:
            rows.append({"smiles": smiles, "project_reward": 0.0, "valid": False})
            continue
        protein_metrics = protein_scorer.score(smiles) if protein_scorer else None
        site_metrics = score_site_c(mol, site_profile) if site_profile else None
        _, metrics = score_candidate(mol, None, filters, protein_metrics, site_metrics)
        rows.append({"smiles": smiles, "project_reward": metrics["reward"], "valid": True, **metrics})
    scored = pd.concat([frame.reset_index(drop=True), pd.DataFrame(rows)], axis=1)
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    scored.to_csv(args.output, index=False)
    print(json.dumps({"input_rows": len(frame), "output": str(args.output)}, indent=2))
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(prog="navgen")
    subparsers = parser.add_subparsers(dest="command", required=True)

    gen = subparsers.add_parser("generate", help="Generate and score candidate molecules")
    gen.add_argument("--config", required=True, help="JSON configuration path")
    gen.add_argument("--output", required=True, help="Directory for generated outputs")
    gen.set_defaults(func=_generate)

    prep = subparsers.add_parser("prepare-reinvent", help="Create REINVENT4 sampling and learning config stubs")
    prep.add_argument("--config", help="Project config file", default=None)
    prep.add_argument("--prior", required=True, help="Mol2Mol prior file path")
    prep.add_argument("--output", required=True, help="Output directory for reinforcement-learning config")
    prep.add_argument("--device", default="cpu")
    prep.add_argument("--num-smiles", type=int, default=1000)
    prep.add_argument("--python", default="python", help="Python executable in the scoring environment")
    prep.set_defaults(func=_prepare_reinvent)

    score = subparsers.add_parser("score-stdin", help="Score SMILES for REINVENT4 ExternalProcess")
    score.add_argument("--config", help="Project JSON configuration for the shared scorer")
    score.set_defaults(func=_score_stdin)

    score_file = subparsers.add_parser("score-reinvent-file", help="Apply the project scorer to a REINVENT4 CSV")
    score_file.add_argument("--config", required=True, help="Project JSON configuration")
    score_file.add_argument("--input", required=True, help="REINVENT4 CSV output")
    score_file.add_argument("--output", required=True, help="Scored CSV destination")
    score_file.add_argument("--smiles-column", default="SMILES")
    score_file.set_defaults(func=_score_reinvent_file)

    args = parser.parse_args()
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
