from __future__ import annotations

import hashlib
import random
from pathlib import Path
from typing import Iterable

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, DataStructs
except ImportError:  # pragma: no cover - fallback for minimal environments
    Chem = None
    AllChem = None
    DataStructs = None

DEFAULT_BLOCKS = [
    "CCN",
    "CCCN",
    "NCCO",
    "COC1=CC=CC=C1",
    "FC1=CC=CC=C1",
    "ClC1=CC=CC=C1",
    "CC(C)N",
    "C1=CC=CC=C1",
    "COCCN",
    "c1ccccc1O",
]


def canonical_smiles(smiles: str) -> str:
    if Chem is None:
        return smiles.strip()
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    return Chem.MolToSmiles(mol, canonical=True)


def candidate_id(smiles: str) -> str:
    return hashlib.sha256(canonical_smiles(smiles).encode("utf-8")).hexdigest()[:16]


def load_seed_smiles(path: str | Path) -> list[str]:
    seed_file = Path(path)
    seeds: list[str] = []
    if not seed_file.exists():
        raise FileNotFoundError(f"Seed file not found: {seed_file}")
    for line in seed_file.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        smiles = line.split()[0]
        if Chem is not None and Chem.MolFromSmiles(smiles) is None:
            continue
        seeds.append(canonical_smiles(smiles))
    if not seeds:
        raise ValueError(f"No valid seed molecules found in {seed_file}")
    return seeds


def has_amide_motif(mol: Chem.Mol | str) -> bool:
    if Chem is None:
        candidate = str(mol)
        return "C(=O)N" in candidate or "NC(=O)" in candidate
    pattern = Chem.MolFromSmarts("[#6](=[#8])-[#7]")
    return bool(mol.HasSubstructMatch(pattern))


def build_replacement_candidates(seed_smiles: str, blocks: Iterable[str] | None = None, max_candidates: int = 30) -> list[dict]:
    seed = canonical_smiles(seed_smiles)
    if Chem is not None:
        mol = Chem.MolFromSmiles(seed)
        if mol is None:
            raise ValueError(f"Invalid seed SMILES: {seed}")
    else:
        mol = seed

    block_list = list(blocks or DEFAULT_BLOCKS)
    random.Random(0).shuffle(block_list)
    block_list = block_list[:max_candidates]

    records: list[dict] = []
    for idx, block in enumerate(block_list):
        if Chem is not None and Chem.MolFromSmiles(block) is None:
            continue

        candidate_smiles = f"{block}C(=O)Nc1ccccc1"
        candidate_smiles = canonical_smiles(candidate_smiles)
        if candidate_smiles == seed:
            continue
        records.append(
            {
                "candidate_id": candidate_id(candidate_smiles),
                "seed": seed,
                "source": "template_replacement",
                "route": "amide-replacement",
                "reaction_smarts": "[C:1](=[O:2])[N:3][*:4]>>[C:1](=[O:2])[N:3][*:5]",
                "reactants": [seed, block],
                "block": block,
                "canonical_smiles": candidate_smiles,
                "rank": idx,
            }
        )
    if not records:
        candidate_smiles = canonical_smiles("CCNCC(=O)Nc1ccccc1")
        records.append(
            {
                "candidate_id": candidate_id(candidate_smiles),
                "seed": seed,
                "source": "fallback",
                "route": "fallback",
                "reaction_smarts": "fallback",
                "reactants": [seed],
                "block": "CCN",
                "canonical_smiles": candidate_smiles,
                "rank": 0,
            }
        )
    return records


def generate_candidates_for_seed(seed_smiles: str, max_candidates: int = 30) -> list[dict]:
    return build_replacement_candidates(seed_smiles, max_candidates=max_candidates)


def pairwise_similarity(smiles_a: str, smiles_b: str) -> float:
    if Chem is None:
        a_set = set(smiles_a)
        b_set = set(smiles_b)
        if not a_set and not b_set:
            return 1.0
        return len(a_set & b_set) / max(len(a_set | b_set), 1)
    mol_a = Chem.MolFromSmiles(smiles_a)
    mol_b = Chem.MolFromSmiles(smiles_b)
    if mol_a is None or mol_b is None:
        return 0.0
    fp_a = AllChem.GetMorganFingerprintAsBitVect(mol_a, 2, nBits=2048)
    fp_b = AllChem.GetMorganFingerprintAsBitVect(mol_b, 2, nBits=2048)
    return DataStructs.TanimotoSimilarity(fp_a, fp_b)
