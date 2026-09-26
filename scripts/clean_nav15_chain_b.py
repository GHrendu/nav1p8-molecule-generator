from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path

from clean_nav18_chain_a import RESIDUE_REQUIRED


def clean_chain_b_pdb(input_path: Path, output_path: Path) -> tuple[int, int]:
    lines = input_path.read_text(encoding="utf-8", errors="ignore").splitlines()
    residues: dict[tuple[int, str, str], dict[str, list[str]]] = defaultdict(
        lambda: defaultdict(list)
    )
    for line in lines:
        if not line.startswith("ATOM  ") or line[21:22] != "B":
            continue
        residue_name = line[17:20].strip()
        if residue_name not in RESIDUE_REQUIRED:
            continue
        altloc = line[16:17]
        if altloc not in (" ", "A"):
            continue
        residue_key = (int(line[22:26]), line[26:27], residue_name)
        atom_name = line[12:16].strip()
        residues[residue_key][atom_name].append(line)

    kept_residues: list[tuple[tuple[int, str, str], dict[str, str]]] = []
    for residue_key, atom_lines in sorted(residues.items()):
        selected_atoms: dict[str, str] = {}
        for atom_name, alternatives in atom_lines.items():
            selected_atoms[atom_name] = max(
                alternatives,
                key=lambda line: (
                    line[16:17] == " ",
                    float(line[54:60].strip() or 0),
                ),
            )
        residue_name = residue_key[2]
        if RESIDUE_REQUIRED[residue_name].issubset(selected_atoms):
            kept_residues.append((residue_key, selected_atoms))

    output_lines: list[str] = []
    atom_index = 1
    for residue_index, (_, atom_lines) in enumerate(kept_residues, start=1):
        for atom_name, line in atom_lines.items():
            if atom_name == "OXT":
                continue
            output_lines.append(
                line[:6]
                + f"{atom_index:5d}"
                + line[11:16]
                + " "
                + line[17:20]
                + " A"
                + f"{residue_index:4d}"
                + " "
                + line[27:]
            )
            atom_index += 1
    output_lines.extend(("TER", "END"))
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(output_lines) + "\n", encoding="utf-8")
    return len(kept_residues), atom_index - 1


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Clean Nav1.5 chain B to complete standard amino-acid residues."
    )
    parser.add_argument("--input", default="Nav1.5_structure/6LQA.pdb")
    parser.add_argument("--output", default="outputs/nav15_chainB_6LQA_clean.pdb")
    args = parser.parse_args()
    residue_count, atom_count = clean_chain_b_pdb(Path(args.input), Path(args.output))
    print(
        f"Wrote {args.output}: {residue_count} complete standard residues, "
        f"{atom_count} ATOM records."
    )


if __name__ == "__main__":
    main()
