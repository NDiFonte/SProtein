"""Helpers to build atom selections without interactive prompts."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass

import MDAnalysis as mda

COMMON_WATER_NAMES = ("SOL", "HOH", "WAT", "TIP3", "TIP3P", "TIP4", "TIP4P", "TIP5", "SPC", "SPCE", "T3P", "T4P", "H2O")


class SelectionError(ValueError):
    """Raised when a selection is empty or ambiguous; the message tells the user how to fix it."""


@dataclass
class WaterModel:
    resname: str
    oxygen: str
    hydrogens: tuple[str, ...]  # empty for oxygen-only water (e.g. crystal structures)


def detect_water(U: mda.Universe, resname: str = "auto") -> WaterModel:
    """Find the water residue and its O / H atom names.

    With ``resname="auto"`` the most abundant residue with a common water name
    (SOL, HOH, WAT, TIP3, ...) is used.
    """
    counts = Counter(U.residues.resnames)
    if resname == "auto":
        found = [r for r in counts if r.upper() in COMMON_WATER_NAMES]
        if not found:
            raise SelectionError(
                "Could not detect the water residue automatically. Set `water = \"<RESNAME>\"` in the "
                "config (run `sprotein inspect <topology>` to see residue names)."
            )
        resname = max(found, key=lambda r: counts[r])
    if resname not in counts:
        raise SelectionError(f"Water residue '{resname}' not found in the topology.")
    names = list(U.select_atoms(f"resname {resname}").residues[0].atoms.names)
    oxygens = [n for n in names if n.upper().startswith("O")]
    hydrogens = [n for n in names if n.upper().startswith("H")]
    if not oxygens:
        raise SelectionError(f"Residue '{resname}' does not look like water (atoms: {' '.join(names)}).")
    # Crystal structures usually have oxygen-only waters: fine for undirected graphs.
    return WaterModel(resname, oxygens[0], tuple(hydrogens[:2]) if len(hydrogens) >= 2 else ())


def select_one_atom(U: mda.Universe, selection: str, what: str = "reference atom") -> mda.AtomGroup:
    ag = U.select_atoms(selection)
    if len(ag) == 0:
        raise SelectionError(f"The {what} selection '{selection}' matches no atom.")
    if len(ag) > 1:
        listing = ", ".join(f"{a.resname}{a.resid}:{a.name}(index {a.index})" for a in ag[:5])
        raise SelectionError(
            f"The {what} selection '{selection}' matches {len(ag)} atoms ({listing}{', ...' if len(ag) > 5 else ''}). "
            "Make it unique, e.g. by adding `and index N`."
        )
    return ag


def hydrogens_of_selections(U: mda.Universe, selections: list[str], h_distance: float) -> list[str]:
    """Selection strings for hydrogens bonded (within ``h_distance``) to atoms in ``selections``."""
    out: set[str] = set()
    for sel in selections:
        atoms = U.select_atoms(sel)
        if len(atoms) == 0:
            raise SelectionError(f"Node selection '{sel}' matches no atom.")
        for atom in atoms:
            near = U.select_atoms(f"(around {h_distance} (index {atom.index})) and name H*")
            for h in near:
                if h.resid == atom.resid and h.resname == atom.resname:
                    out.add(f"resname {h.resname} and name {h.name}")
    return sorted(out)


def describe_topology(U: mda.Universe, resname: str | None = None) -> str:
    """Human readable summary used by ``sprotein inspect``."""
    lines = []
    if resname is None:
        counts = Counter(U.residues.resnames)
        prot = U.select_atoms("protein")
        lines.append(f"Atoms: {U.atoms.n_atoms}   Residues: {U.residues.n_residues}   Frames: {U.trajectory.n_frames}")
        if prot.n_atoms:
            lines.append(f"Protein: {prot.n_residues} residues (resid {prot.resids.min()}-{prot.resids.max()})")
        try:
            w = detect_water(U)
            h = ",".join(w.hydrogens) if w.hydrogens else "none (oxygen only: undirected graphs only)"
            lines.append(f"Water (auto-detected): {w.resname}  O={w.oxygen}  H={h}")
        except SelectionError as e:
            lines.append(f"Water: {e}")
        lines.append("")
        lines.append(f"{'resname':<8} {'count':>7}  atoms (first residue)")
        for name in sorted(counts):
            atoms = U.select_atoms(f"resname {name}").residues[0].atoms.names
            shown = " ".join(atoms[:14]) + (" ..." if len(atoms) > 14 else "")
            lines.append(f"{name:<8} {counts[name]:>7}  {shown}")
        lines.append("")
        lines.append("Use `sprotein inspect <topology> --resname NAME` to list atoms, resids and indices of a residue.")
    else:
        res = U.select_atoms(f"resname {resname}").residues
        if len(res) == 0:
            raise SelectionError(f"No residue named '{resname}'.")
        lines.append(f"{resname}: {len(res)} residue(s)")
        for r in res[:20]:
            lines.append(f"  resid {r.resid}:")
            for a in r.atoms:
                lines.append(f"    {a.name:<6} index {a.index:<7} -> \"resname {resname} and resid {r.resid} and name {a.name}\"")
        if len(res) > 20:
            lines.append(f"  ... {len(res) - 20} more")
    return "\n".join(lines)
