"""Main analysis: shortest water/protein paths from a reference atom to bulk water."""

from __future__ import annotations

import csv
import json
import random
import time
from contextlib import ExitStack
from dataclasses import dataclass

import MDAnalysis as mda
import numpy as np

from . import graph
from .config import Config, asdict, validate
from .selection import SelectionError, detect_water, hydrogens_of_selections, select_one_atom

OUTPUT_FILES = {
    "sp_xyz": "traj_sp.xyz",
    "sp_nodes": "traj_path_nodes.dat",
    "all_xyz": "traj_all_paths.xyz",
    "all_nodes": "traj_all_paths_nodes.dat",
    "summary": "summary.csv",
    "params": "parameters.json",
}


@dataclass
class FrameResult:
    frame: int
    time: float
    n_nodes: int
    n_edges: int
    n_bulk: int
    paths: list[list[int]]
    shortest: list[int]


def _or_join(selections: list[str]) -> str:
    return " or ".join(f"({s})" for s in selections)


class PathFinder:
    """Builds the node selections once and analyses frames one at a time."""

    def __init__(self, U: mda.Universe, cfg: Config):
        validate(cfg)
        self.U = U
        self.cfg = cfg
        self.protein = U.select_atoms("protein")
        if self.protein.n_atoms == 0:
            raise SelectionError("The topology contains no protein atoms (MDAnalysis 'protein' selection is empty).")
        self.water = detect_water(U, cfg.water)
        self.ref = select_one_atom(U, cfg.reference)
        self.ref_index = int(self.ref[0].index)
        O = self.water.oxygen

        parts = [f"({cfg.reference})"]
        if cfg.selections:
            parts.append(f"({_or_join(cfg.selections)})")
        parts.append(f"(name {O} and around {cfg.protein_cutoff} group gprotein)")
        self.nodes = U.select_atoms(" or ".join(parts), gprotein=self.protein, updating=True)

        self.hydrogens = None
        if cfg.directed:
            h1, h2 = self.water.hydrogens
            h_sels = hydrogens_of_selections(U, cfg.selections, cfg.h_distance)
            h_parts = ([f"({_or_join(h_sels)})"] if h_sels else []) + [
                f"((name {h1} {h2}) and around {cfg.protein_cutoff + 3} group gprotein)"
            ]
            self.hydrogens = U.select_atoms(" or ".join(h_parts), gprotein=self.protein, updating=True)
        self._protein_resix = np.unique(self.protein.resindices)
        self.rng = random.Random(cfg.seed)

    def analyse_frame(self, ts) -> FrameResult:
        cfg = self.cfg
        nodes = self.nodes
        coords = nodes.positions
        box = ts.dimensions
        if cfg.pbc and box is None:
            raise SelectionError("pbc = true but the trajectory has no box information; set pbc = false.")
        pbc_box = box if cfg.pbc else None

        src = int(np.searchsorted(nodes.indices, self.ref_index))
        A = graph.adjacency_matrix(coords, cfg.edge_distance, pbc_box)

        water_idx = np.flatnonzero(nodes.names == self.water.oxygen)
        A_ww = A[water_idx][:, water_idx]
        degree = np.asarray(A_ww.sum(axis=1)).ravel()
        bulk = water_idx[degree >= cfg.bulk_degree]

        if cfg.directed:
            H = self.hydrogens
            is_protein = np.isin(nodes.resindices, self._protein_resix)
            G = graph.directed_hbond_matrix(
                A, coords, nodes.resindices, is_protein, H.positions, H.resindices,
                box, cfg.pbc, cfg.hbond_angle,
            )
            if cfg.source_is_donor:
                G = G.T.tocsr()
                G.sort_indices()
        else:
            G = A

        pred = graph.bfs_predecessors(G, src)
        paths = graph.channel_paths(pred, src, bulk)
        shortest = graph.pick_shortest(paths, self.rng)
        return FrameResult(ts.frame, float(ts.time), len(nodes), int(G.nnz), len(bulk), paths, shortest)


def _node_line(atom) -> str:
    # 1-based atom index (as in .gro/.pdb files), residue name, residue id
    return f"{atom.index + 1} {atom.resname} {atom.resid}\n"


def _xyz_line(xyz) -> str:
    return f"{xyz[0]:.3f} {xyz[1]:.3f} {xyz[2]:.3f}\n"


def run(cfg: Config, progress: bool = True, log=print) -> dict:
    """Run the analysis described by ``cfg`` and write the output files."""
    validate(cfg)
    top, trj = cfg.path(cfg.topology), cfg.path(cfg.trajectory)
    for p in (top, trj):
        if not p.exists():
            raise FileNotFoundError(f"Input file '{p}' not found.")
    log(f"Loading {top.name} + {trj.name} ...")
    U = mda.Universe(str(top), str(trj))
    finder = PathFinder(U, cfg)
    w = finder.water
    log(f"Water: {w.resname} (O={w.oxygen})   Reference: {finder.ref[0].resname}{finder.ref[0].resid}:{finder.ref[0].name}"
        f"   Graph: {'directed' if cfg.directed else 'undirected'}, PBC {'on' if cfg.pbc else 'off'}")

    out = cfg.out
    out.mkdir(parents=True, exist_ok=True)
    params = {k: (str(v) if not isinstance(v, (int, float, str, bool, list, dict, type(None))) else v)
              for k, v in asdict(cfg).items() if k != "base_dir"}
    params["water_detected"] = {"resname": w.resname, "oxygen": w.oxygen, "hydrogens": list(w.hydrogens)}
    (out / OUTPUT_FILES["params"]).write_text(json.dumps(params, indent=2))

    frames = U.trajectory[cfg.start:cfg.stop:cfg.step]
    iterator = frames
    if progress:
        try:
            from tqdm.auto import tqdm
            iterator = tqdm(frames, total=len(frames), unit="frame", desc="Shortest paths")
        except ImportError:
            pass

    lengths, n_frames, t0 = [], 0, time.perf_counter()
    with ExitStack() as stack:
        f = {k: stack.enter_context(open(out / OUTPUT_FILES[k], "w", newline="" if k == "summary" else None))
             for k in ("sp_xyz", "sp_nodes", "summary")}
        if cfg.all_paths:
            for k in ("all_xyz", "all_nodes"):
                f[k] = stack.enter_context(open(out / OUTPUT_FILES[k], "w"))
        summary = csv.writer(f["summary"])
        summary.writerow(["frame", "time_ps", "nodes", "edges", "bulk_waters", "exits", "shortest_length", "shortest_path"])

        for ts in iterator:
            r = finder.analyse_frame(ts)
            atoms, coords = finder.nodes, finder.nodes.positions
            header = f"FRAME {r.frame + 1}\n"
            for k in f:
                if k != "summary":
                    f[k].write(header)
            if cfg.all_paths:
                for p in r.paths:
                    f["all_nodes"].writelines(_node_line(atoms[i]) for i in p)
                    f["all_xyz"].writelines(_xyz_line(coords[i]) for i in p)
                    f["all_nodes"].write("\n")
                    f["all_xyz"].write("\n")
            f["sp_nodes"].writelines(_node_line(atoms[i]) for i in r.shortest)
            f["sp_xyz"].writelines(_xyz_line(coords[i]) for i in r.shortest)
            for k in f:
                if k != "summary":
                    f[k].write("\n")
            n = len(r.shortest) - 1 if r.shortest else ""
            if r.shortest:
                lengths.append(n)
            route = " ".join(f"{atoms[i].resname}{atoms[i].resid}" for i in r.shortest)
            summary.writerow([r.frame + 1, f"{r.time:.3f}", r.n_nodes, r.n_edges, r.n_bulk, len(r.paths), n, route])
            n_frames += 1
            for fh in f.values():
                fh.flush()

    elapsed = time.perf_counter() - t0
    stats = {
        "frames": n_frames,
        "frames_with_path": len(lengths),
        "mean_length": float(np.mean(lengths)) if lengths else None,
        "seconds": elapsed,
        "output_dir": str(out),
    }
    log(f"\nDone: {n_frames} frames in {elapsed:.1f} s. Path found in {len(lengths)}/{n_frames} frames"
        + (f", mean length {stats['mean_length']:.2f} edges." if lengths else "."))
    log(f"Results in {out}/  (summary.csv gives one line per frame)")
    return stats
