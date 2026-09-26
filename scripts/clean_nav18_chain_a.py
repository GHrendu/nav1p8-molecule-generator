from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path

RESIDUE_REQUIRED = {
    "ALA": {"N", "CA", "C", "O", "CB"},
    "ARG": {"N", "CA", "C", "O", "CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"},
    "ASN": {"N", "CA", "C", "O", "CB", "CG", "OD1", "ND2"},
    "ASP": {"N", "CA", "C", "O", "CB", "CG", "OD1", "OD2"},
    "CYS": {"N", "CA", "C", "O", "CB", "SG"},
    "GLN": {"N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "NE2"},
    "GLU": {"N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "OE2"},
    "GLY": {"N", "CA", "C", "O"},
    "HIS": {"N", "CA", "C", "O", "CB", "CG", "ND1", "CD2", "CE1", "NE2"},
    "ILE": {"N", "CA", "C", "O", "CB", "CG1", "CG2", "CD1"},
    "LEU": {"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2"},
    "LYS": {"N", "CA", "C", "O", "CB", "CG", "CD", "CE", "NZ"},
    "MET": {"N", "CA", "C", "O", "CB", "CG", "SD", "CE"},
    "PHE": {"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
    "PRO": {"N", "CA", "C", "O", "CB", "CG", "CD"},
    "SER": {"N", "CA", "C", "O", "CB", "OG"},
    "THR": {"N", "CA", "C", "O", "CB", "OG1", "CG2"},
    "TRP": {"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"},
    "TYR": {"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"},
    "VAL": {"N", "CA", "C", "O", "CB", "CG1", "CG2"},
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Clean a Nav1.8 receptor to a chain-A only, standard-residue-only PDB")
    parser.add_argument("--input", default="Nav1.8_structure/7WE4.pdb", help="Input PDB file")
    parser.add_argument("--output", default="outputs/nav18_chainA_complete.pdb", help="Output PDB file")
    return parser.parse_args()


def is_standard_residue(name: str) -> bool:
    return name in RESIDUE_REQUIRED


def clean_chain_a_pdb(input_path: Path, output_path: Path) -> None:
    lines = input_path.read_text(encoding='utf-8', errors='ignore').splitlines()
    residues: dict[int, dict[str, set[str]]] = defaultdict(lambda: defaultdict(set))
    for line in lines:
        if not line.startswith('ATOM'):
            continue
        if line[21:22] != 'A':
            continue
        res_name = line[17:20].strip()
        if not is_standard_residue(res_name):
            continue
        res_index = int(line[22:26])
        atom_name = line[12:16].strip()
        residues[res_index][res_name].add(atom_name)

    kept_resids = []
    for res_index, res_to_atoms in sorted(residues.items()):
        for res_name, atom_names in res_to_atoms.items():
            required = RESIDUE_REQUIRED[res_name]
            if required.issubset(atom_names):
                kept_resids.append((res_index, res_name))
                break

    kept_resids.sort(key=lambda item: item[0])
    old_to_new = {old_index: idx + 1 for idx, (old_index, _) in enumerate(kept_resids)}

    out_lines: list[str] = []
    atom_index = 1
    for line in lines:
        if not line.startswith('ATOM'):
            continue
        if line[21:22] != 'A':
            continue
        res_name = line[17:20].strip()
        if not is_standard_residue(res_name):
            continue
        res_index = int(line[22:26])
        if res_index not in old_to_new:
            continue
        new_res_index = old_to_new[res_index]
        atom_name = line[12:16].strip()
        if atom_name == 'OXT':
            continue
        # Keep the original coordinates and atom fields exactly as encoded, but replace the
        # residue sequence number with a continuous chain-A numbering so the cleaned structure is
        # accepted by downstream receptor preparation tools.
        new_line = line[:6] + f"{atom_index:5d}" + line[11:17] + f"{res_name:<3s}" + " A" + f"{new_res_index:4d}" + line[26:]
        out_lines.append(new_line.rstrip())
        atom_index += 1

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text('\n'.join(out_lines) + '\n', encoding='utf-8')
    print(f"Wrote {output_path} with {len(out_lines)} ATOM lines and {len(old_to_new)} residues.")


if __name__ == '__main__':
    args = parse_args()
    clean_chain_a_pdb(Path(args.input), Path(args.output))
