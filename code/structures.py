"""Ligand-to-residue distance measurement and pocket assignment.

Given a structure with an inhibitor bound (a co-crystal, or an AlphaFold model
with the inhibitor aligned into the pocket), measure the minimum heavy-atom
distance from the ligand to every residue of the protein chain, then classify
screen hits as inside or outside the pocket.

The measuring half runs once, in ``build_inputs``, and its results are stored
in ``Inputs-PocketDistances.xlsx``; a notebook run reads that workbook through
:func:`load_pocket_distances` and never opens an mmCIF file.
"""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.PDB import MMCIFParser

from . import config as cfg

AA_THREE_TO_ONE = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
    'TER': '.',
}


def is_protein_residue(res) -> bool:
    """Biopython flags standard residues with a blank hetero field."""
    return res.id[0] == " "


def is_het_residue(res) -> bool:
    return res.id[0].strip() != ""


def heavy_atom_coords(res) -> np.ndarray:
    """Coordinates of every non-hydrogen atom in a residue."""
    coords = [atom.coord for atom in res.get_atoms() if atom.element != "H"]
    return np.array(coords, dtype=float)


def min_distance(coords_a: np.ndarray, coords_b: np.ndarray) -> float:
    """Smallest pairwise distance between two coordinate sets."""
    diffs = coords_a[:, None, :] - coords_b[None, :, :]
    return float(np.sqrt(np.min(np.sum(diffs * diffs, axis=2))))


def residue_ligand_distances(
    structure_file: str | Path,
    ligand_name: str,
    ligand_chain: str,
    protein_chain: str,
) -> pd.DataFrame:
    """Minimum ligand distance for every residue of the protein chain.

    Returns a frame of Distance / Chain / Residue / Position sorted by
    distance. The original notebooks wrote the chain id of whichever chain the
    ligand search loop happened to end on, rather than the protein chain being
    measured; this records ``protein_chain``.
    """
    structure = MMCIFParser(QUIET=True).get_structure("struct", str(structure_file))
    model = structure[0]

    ligands = [
        res for chain in model for res in chain
        if is_het_residue(res)
        and res.get_resname().strip() == ligand_name
        and chain.id == ligand_chain
    ]
    if not ligands:
        raise ValueError(
            f"No ligand {ligand_name!r} in chain {ligand_chain!r} of {structure_file}"
        )

    lig_coords = heavy_atom_coords(ligands[0])
    if lig_coords.size == 0:
        raise ValueError(f"Ligand {ligand_name!r} has no heavy atoms")

    rows, residues = [], []
    for res in model[protein_chain]:
        if not is_protein_residue(res):
            continue
        residues.append(res)
        res_coords = heavy_atom_coords(res)
        if res_coords.size == 0:
            continue
        rows.append({
            'Distance': round(min_distance(res_coords, lig_coords), 2),
            'Chain': protein_chain,
            'Residue': res.get_resname().strip(),
            'Position': res.id[1],
        })

    df = pd.DataFrame(rows).sort_values('Distance').reset_index(drop=True)
    df.attrs['sequence'] = format_residues(residues)
    return df


def format_residues(residues) -> str:
    """One-letter sequence with '-' at unresolved positions."""
    if not residues:
        return ''
    protein_len = residues[-1].get_full_id()[3][1]
    chars = ['-'] * protein_len
    for res in residues:
        index = res.get_full_id()[3][1]
        if 1 <= index <= protein_len:
            chars[index - 1] = AA_THREE_TO_ONE.get(res.resname, 'X')
    return ''.join(chars)


@lru_cache(maxsize=1)
def _distance_workbook(path: str) -> dict[str, pd.DataFrame]:
    """Every sheet of the distance workbook, parsed once per session.

    Both notebooks ask for the co-crystal and the AlphaFold tables separately,
    and the interactive notebook asks again per figure; the workbook is
    eighteen sheets, so it is read once and handed out from here.
    """
    sheets = pd.read_excel(path, sheet_name=None)
    return {name: table for name, table in sheets.items()
            if name != cfg.POCKET_STRUCTURE_SHEET}


def load_distance_tables(path=None) -> dict[str, pd.DataFrame]:
    """Residue-to-ligand distances for every measured structure."""
    return _distance_workbook(str(path or cfg.POCKET_DISTANCES_XLSX))


def load_pocket_distances(structure_specs, out_dir: Path) -> dict[str, pd.DataFrame]:
    """Distance table for each configured structure, also written to CSV.

    ``structure_specs`` is one of ``config.PDB_STRUCTURES`` /
    ``config.AF_STRUCTURES``. The tables are read from
    ``Inputs-PocketDistances.xlsx`` rather than re-measured: the mmCIF parse
    above produced that workbook once, in ``build_inputs``, and is no longer
    on the path of a notebook run.
    """
    tables = load_distance_tables()
    out_dir.mkdir(parents=True, exist_ok=True)
    selected = {}
    for gene, structure, _inhibitor, lig_name, lig_chain, prot_chain, _cond in structure_specs:
        df = tables[structure].copy()
        df.to_csv(out_dir / f"{structure}-distances.csv", index=False)
        selected[structure] = df
    return selected


def _near_any(hit_pos: float, positions, window: int = cfg.GUIDE_POSITION_WINDOW) -> bool:
    """Whether any residue position falls within ``window`` of a hit."""
    return any(hit_pos - window <= p <= hit_pos + window for p in positions)


def split_in_out_pocket(
    hits: pd.DataFrame,
    distance_df: pd.DataFrame,
    pos_col: str = 'pos',
    cutoff: float = cfg.POCKET_DISTANCE_CUTOFF,
    window: int = cfg.GUIDE_POSITION_WINDOW,
) -> dict[str, pd.DataFrame]:
    """Partition hits into in-pocket, outside-pocket and structurally unresolved.

    A hit is in-pocket when a residue within ``cutoff`` angstroms of the ligand
    lies within ``window`` residues of the edited position; outside-pocket when
    a resolved residue is nearby but none is close to the ligand; unresolved
    when the structure covers no residue near that position at all.
    """
    within = distance_df.loc[distance_df['Distance'] <= cutoff, 'Position'].tolist()
    without = distance_df.loc[distance_df['Distance'] > cutoff, 'Position'].tolist()

    labels = []
    for pos in hits[pos_col]:
        if _near_any(pos, within, window):
            labels.append('InPocket')
        elif _near_any(pos, without, window):
            labels.append('OutPocket')
        else:
            labels.append('Unresolved')

    labelled = hits.assign(pocket=labels)
    return {
        name: labelled[labelled['pocket'] == name].copy()
        for name in ('InPocket', 'OutPocket', 'Unresolved')
    }


def interpolate_positions(df: pd.DataFrame, xcol: str = 'Position') -> pd.DataFrame:
    """Reindex a distance table onto every integer position in its range.

    Gaps in a structure would otherwise draw as straight lines between the
    flanking resolved residues in the lollipop distance panel.
    """
    df = df.sort_values(xcol).set_index(xcol)
    df = df.reindex(np.arange(df.index.min(), df.index.max() + 1))
    df = df.reset_index()
    df[xcol] = df[xcol].astype(int)
    return df
