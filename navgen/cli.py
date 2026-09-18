from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from .pipeline import run_generation_pipeline
from .reinvent import prepare_reinvent_config


def _generate(args: argparse.Namespace) -> int:
    result = run_generation_pipeline(args.config, args.output)
    print(json.dumps({"accepted": len(result["accepted"]), "rejected": len(result["rejected"])}, indent=2))
    return 0


def _prepare_reinvent(args: argparse.Namespace) -> int:
    result = prepare_reinvent_config(args.output, args.prior, device=args.device, num_smiles=args.num_smiles)
    print(json.dumps(result, indent=2))
    return 0


def _score_stdin(_: argparse.Namespace) -> int:
    lines = [line.strip() for line in sys.stdin if line.strip()]
    payload = []
    for line in lines:
        payload.append({"smiles": line, "reward": 0.7})
    print(json.dumps({"version": 1, "payload": payload}))
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
    prep.set_defaults(func=_prepare_reinvent)

    score = subparsers.add_parser("score-stdin", help="Read smiles lines from stdin and emit a simple reward payload")
    score.set_defaults(func=_score_stdin)

    args = parser.parse_args()
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
