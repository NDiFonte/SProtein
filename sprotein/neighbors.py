"""Closest protein residues to every node of the paths found by ``sprotein run``."""

from __future__ import annotations

from pathlib import Path

import MDAnalysis as mda
import numpy as np
from scipy.spatial import cKDTree


def annotate_neighbors(
    U: mda.Universe,
    src: Path,
    dst: Path,
    n: int = 1,
    cutoff: float = 7.0,
    anchors: dict[str, str] | None = None,
    progress: bool = True,
) -> int:
    """Append the ``n`` closest protein residues (by centre of mass, within
    ``cutoff`` A) to each node line of a ``*_nodes.dat`` file.

    ``anchors`` maps a residue name to an atom name to use instead of that
    residue's centre of mass (e.g. ``{"PAL": "C1"}`` for a long ligand).
    Returns the number of frames processed.
    """
    anchors = anchors or {}
    protein = U.select_atoms("protein")
    residues = protein.residues
    labels = [f"{r.resname}{r.resid}" for r in residues]

    lines = Path(src).read_text().splitlines()
    n_frames = sum(1 for line in lines if line.startswith("FRAME"))
    bar = None
    if progress:
        try:
            from tqdm.auto import tqdm
            bar = tqdm(total=n_frames, unit="frame", desc="Neighbours")
        except ImportError:
            pass

    tree = None
    with open(dst, "w") as out:
        for line in lines:
            if line.startswith("FRAME"):
                U.trajectory[int(line.split()[1]) - 1]
                coms = protein.center_of_mass(compound="residues")
                tree = cKDTree(coms)
                out.write(line + "\n\n")
                if bar:
                    bar.update()
                continue
            parts = line.split()
            if len(parts) < 3:
                out.write("\n")
                continue
            atom = U.atoms[int(parts[0]) - 1]
            residue = atom.residue
            anchor = anchors.get(residue.resname)
            if anchor:
                sel = residue.atoms.select_atoms(f"name {anchor}")
                point = sel.positions[0] if len(sel) else residue.atoms.center_of_mass()
            else:
                point = residue.atoms.center_of_mass()
            idx = tree.query_ball_point(point, cutoff)
            if idx:
                d = np.linalg.norm(coms[idx] - point, axis=1)
                idx = [idx[k] for k in np.argsort(d, kind="stable")[:n]]
            out.write(line.strip() + "\t" + " ".join(labels[i] for i in idx) + "\n")
    if bar:
        bar.close()
    return n_frames
