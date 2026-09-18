from __future__ import annotations

from pathlib import Path


def prepare_reinvent_config(output_dir: str | Path, prior_path: str, device: str = "cpu", num_smiles: int = 1000) -> dict[str, str]:
    output_root = Path(output_dir)
    output_root.mkdir(parents=True, exist_ok=True)
    sampling = output_root / "sampling.toml"
    learning = output_root / "staged_learning.toml"

    sampling.write_text(
        "\n".join(
            [
                '[sampling]',
                f'prior = "{prior_path}"',
                f'device = "{device}"',
                f'num_smiles = {num_smiles}',
                'smiles = ["C(=O)Nc1ccccc1"]',
                'external_process = "python -m navgen score-stdin"',
            ]
        ),
        encoding="utf-8",
    )
    learning.write_text(
        "\n".join(
            [
                '[staged_learning]',
                f'prior = "{prior_path}"',
                f'device = "{device}"',
                'external_process = "python -m navgen score-stdin"',
                'diversity_filter = true',
            ]
        ),
        encoding="utf-8",
    )
    return {"sampling": str(sampling), "staged_learning": str(learning)}
