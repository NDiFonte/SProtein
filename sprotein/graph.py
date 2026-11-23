"""Graph construction and path search.

All functions work on plain NumPy arrays and SciPy CSR matrices so they can be
tested without a trajectory. Node indices always refer to positions in the
per-frame node array.
"""

from __future__ import annotations

import random
from typing import Sequence

import numpy as np
from MDAnalysis.lib.distances import calc_bonds, self_capped_distance
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import breadth_first_order
from scipy.spatial import KDTree

NO_PREDECESSOR = -9999  # value used by scipy.sparse.csgraph


def adjacency_matrix(coords: np.ndarray, cutoff: float, box: np.ndarray | None = None) -> csr_matrix:
    """Symmetric adjacency matrix connecting nodes closer than ``cutoff`` (inclusive).

    Uses a cell-list neighbour search instead of a dense N x N distance matrix,
    then re-checks candidate pairs with exact (float64) distances so that edges
    sitting exactly on the cutoff are treated like a full distance matrix would.

    Parameters
    ----------
    coords : (N, 3) array of positions in Angstrom.
    cutoff : edge distance in Angstrom.
    box : MDAnalysis box ``[lx, ly, lz, alpha, beta, gamma]`` for periodic
        boundary conditions, or ``None`` for no PBC.
    """
    coords = np.ascontiguousarray(coords, dtype=np.float32)
    n = len(coords)
    if n < 2:
        return csr_matrix((n, n), dtype=np.int8)
    # Candidate search with a small margin, exact filter below. 'pkdtree' is used
    # because 'nsgrid' was found to miss a few pairs in triclinic boxes.
    pairs, _ = self_capped_distance(coords, cutoff + 1e-3, box=box, method="pkdtree")
    if len(pairs):
        d = calc_bonds(coords[pairs[:, 0]], coords[pairs[:, 1]], box=box)
        pairs = pairs[d <= cutoff]
    i = np.concatenate([pairs[:, 0], pairs[:, 1]])
    j = np.concatenate([pairs[:, 1], pairs[:, 0]])
    A = csr_matrix((np.ones(len(i), dtype=np.int8), (i, j)), shape=(n, n))
    A.sum_duplicates()
    A.sort_indices()
    return A


def _minimum_image_shift(ref: np.ndarray, other: np.ndarray, lengths: np.ndarray) -> np.ndarray:
    """Shift ``other`` by one box length along each axis where it is more than half
    a box away from ``ref`` (same convention as the original SProtein code)."""
    other = other.copy()
    far = np.abs(ref - other) > lengths / 2
    up = far & (ref >= other)
    down = far & ~(ref >= other)
    other[up] += np.broadcast_to(lengths, other.shape)[up]
    other[down] -= np.broadcast_to(lengths, other.shape)[down]
    return other


def directed_hbond_matrix(
    A: csr_matrix,
    coords: np.ndarray,
    node_resindex: np.ndarray,
    node_is_protein: np.ndarray,
    h_coords: np.ndarray,
    h_resindex: np.ndarray,
    box: np.ndarray | None,
    pbc: bool,
    angle_cutoff_deg: float = 30.0,
) -> csr_matrix:
    """Keep only edges i -> j where i can donate a hydrogen bond to j.

    Edge i -> j survives when one of the (up to 3) hydrogens closest to node i
    that belong to the same residue as i makes an angle H-i-j below
    ``angle_cutoff_deg``. Rules carried over from the original implementation:

    * non-protein nodes with no hydrogen of their own lose all outgoing edges;
    * edges between two protein nodes of the same residue are always kept.
    """
    n = A.shape[0]
    coords = np.asarray(coords, dtype=np.float32)
    h_coords = np.asarray(h_coords, dtype=np.float32)
    rows = np.repeat(np.arange(n), np.diff(A.indptr))
    cols = A.indices.copy()
    if n == 0 or len(rows) == 0:
        return csr_matrix((n, n), dtype=np.int8)

    k = min(3, len(h_coords))
    if k == 0:
        partner = np.zeros((n, 1), dtype=bool)
        partner_xyz = np.zeros((n, 1, 3), dtype=np.float32)
    else:
        _, idx = KDTree(h_coords).query(coords, k=k)
        idx = np.asarray(idx).reshape(n, k)
        partner = h_resindex[idx] == node_resindex[:, None]
        partner_xyz = h_coords[idx]  # (n, k, 3)

    lengths = np.asarray(box[:3], dtype=np.float32) if box is not None else None
    if pbc and lengths is not None:
        partner_xyz = _minimum_image_shift(coords[:, None, :], partner_xyz, lengths)

    o1 = coords[rows]
    o2 = coords[cols]
    if lengths is not None:
        o2 = _minimum_image_shift(o1, o2, lengths)

    v2 = o2 - o1
    v2 = v2 / np.linalg.norm(v2, axis=1, keepdims=True)
    v1 = partner_xyz[rows] - o1[:, None, :]  # (E, k, 3)
    with np.errstate(invalid="ignore", divide="ignore"):
        v1 = v1 / np.linalg.norm(v1, axis=2, keepdims=True)
        cos = np.clip(np.einsum("ekd,ed->ek", v1, v2), -1.0, 1.0)
        angle = np.arccos(cos)
    lower = angle_cutoff_deg * np.pi / 180
    angle_ok = (partner[rows] & (angle <= lower)).any(axis=1)

    same_res_protein = node_is_protein[rows] & node_is_protein[cols] & (node_resindex[rows] == node_resindex[cols])
    dead_row = ~partner.any(axis=1) & ~node_is_protein
    keep = ~dead_row[rows] & (same_res_protein | angle_ok)

    D = csr_matrix((np.ones(keep.sum(), dtype=np.int8), (rows[keep], cols[keep])), shape=(n, n))
    D.sort_indices()
    return D


def bfs_predecessors(A: csr_matrix, source: int) -> np.ndarray:
    """Breadth-first search from ``source``; returns the predecessor of every node
    (``NO_PREDECESSOR`` if unreachable). One search per frame replaces one search
    per bulk-water target."""
    _, pred = breadth_first_order(A, source, directed=True, return_predecessors=True)
    return pred


def path_to(pred: np.ndarray, source: int, target: int) -> list[int]:
    """Reconstruct the path source -> target from a predecessor array ([] if unreachable)."""
    if target == source:
        return [source]
    if pred[target] == NO_PREDECESSOR:
        return []
    path = [target]
    node = target
    while node != source:
        node = int(pred[node])
        path.append(node)
    path.reverse()
    return path


def channel_paths(pred: np.ndarray, source: int, targets: Sequence[int]) -> list[list[int]]:
    """Shortest paths from ``source`` to every reachable target, truncated at the
    first target encountered along the way and de-duplicated.

    This gives the set of distinct "exits" from the source into bulk water.
    Paths are returned sorted by length, then lexicographically.
    """
    paths = [p for p in (path_to(pred, source, int(t)) for t in targets) if p]
    reached = {p[-1] for p in paths}
    unique = set()
    for p in paths:
        for pos, node in enumerate(p):
            if node in reached:
                unique.add(tuple(p[: pos + 1]))
                break
    return [list(p) for p in sorted(unique, key=lambda p: (len(p), p))]


def pick_shortest(paths: Sequence[Sequence[int]], rng: random.Random) -> list[int]:
    """Shortest path among ``paths``; ties are broken at random with ``rng``."""
    if not paths:
        return []
    m = min(len(p) for p in paths)
    return list(rng.choice([p for p in paths if len(p) == m]))
