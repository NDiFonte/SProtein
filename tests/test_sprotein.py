import random
from pathlib import Path

import numpy as np
import pytest

from sprotein import graph
from sprotein.tools import join_paths, path_lengths, read_frames, vmd_trajectory

EXAMPLES = Path(__file__).resolve().parent.parent / "examples"


# ------------------------------------------------------------------ graph
def test_adjacency_chain_and_cutoff_inclusive():
    coords = np.array([[0, 0, 0], [3, 0, 0], [6, 0, 0], [20, 0, 0]], dtype=np.float32)
    A = graph.adjacency_matrix(coords, 3.0)
    assert A.nnz == 4  # 0-1, 1-2 in both directions
    assert A[0, 1] and A[1, 0] and A[1, 2] and not A[0, 2] and not A[2, 3]
    assert (A != A.T).nnz == 0 and A.diagonal().sum() == 0


def test_adjacency_pbc():
    box = np.array([10, 10, 10, 90, 90, 90], dtype=np.float32)
    coords = np.array([[0.5, 5, 5], [9.5, 5, 5]], dtype=np.float32)
    assert graph.adjacency_matrix(coords, 2.0, box).nnz == 2
    assert graph.adjacency_matrix(coords, 2.0, None).nnz == 0


def test_bfs_paths_truncate_and_dedupe():
    # 0 - 1 - 2 - 3,  1 - 4 ; targets 2, 3, 4  -> 3 is behind 2, so it is truncated
    coords = np.array([[0, 0, 0], [3, 0, 0], [6, 0, 0], [9, 0, 0], [3, 3, 0]], dtype=np.float32)
    A = graph.adjacency_matrix(coords, 3.0)
    pred = graph.bfs_predecessors(A, 0)
    paths = graph.channel_paths(pred, 0, [2, 3, 4])
    assert paths == [[0, 1, 2], [0, 1, 4]]
    assert graph.pick_shortest(paths, random.Random(0)) in paths
    assert graph.path_to(pred, 0, 0) == [0]


def test_unreachable_target():
    coords = np.array([[0, 0, 0], [50, 0, 0]], dtype=np.float32)
    pred = graph.bfs_predecessors(graph.adjacency_matrix(coords, 3.0), 0)
    assert graph.channel_paths(pred, 0, [1]) == []
    assert graph.pick_shortest([], random.Random(0)) == []


def test_directed_hbond_angle():
    # donor water O at origin with H pointing to +x; acceptors on +x (kept) and -x (dropped)
    coords = np.array([[0, 0, 0], [2.8, 0, 0], [-2.8, 0, 0]], dtype=np.float32)
    node_res = np.array([0, 1, 2])
    h = np.array([[0.96, 0, 0], [-0.3, 0.9, 0]], dtype=np.float32)
    h_res = np.array([0, 0])
    A = graph.adjacency_matrix(coords, 3.5)
    D = graph.directed_hbond_matrix(A, coords, node_res, np.zeros(3, bool), h, h_res, None, False)
    assert D[0, 1] == 1 and D[0, 2] == 0
    # nodes 1 and 2 have no hydrogens and are not protein: no outgoing edges
    assert D[1].nnz == 0 and D[2].nnz == 0


# ------------------------------------------------------------------ tools
SP = """FRAME 1

FRAME 2
1.000 2.000 3.000
4.000 5.000 6.000
7.000 8.000 9.000

FRAME 3
1.000 1.000 1.000
2.000 2.000 2.000

"""

ALL = """FRAME 1

FRAME 2
1 1 1
2 2 2

1 1 1
3 3 3


"""


def test_read_frames_and_lengths(tmp_path):
    f = tmp_path / "traj_sp.xyz"
    f.write_text(SP)
    frames = read_frames(f)
    assert [h for h, _ in frames] == ["FRAME 1", "FRAME 2", "FRAME 3"]
    s = path_lengths(f)
    assert s["frames"] == 3 and s["frames_with_path"] == 2 and s["lengths"] == [2, 1]


def test_vmd_trajectory(tmp_path):
    src, dst = tmp_path / "traj_sp.xyz", tmp_path / "vmd.xyz"
    src.write_text(SP)
    assert vmd_trajectory(src, dst) == 3
    lines = dst.read_text().splitlines()
    assert len(lines) == 4 * (3 + 2)  # frame 0 + 3 frames, each: count, title, 3 atoms
    assert lines[:2] == ["3", "FRAME 0"]
    assert lines[10:13] == ["3", "FRAME 2", "  Se  1.000 2.000 3.000"]
    assert lines[-1].startswith("  Se  0.000")


def test_join(tmp_path):
    src, dst = tmp_path / "all.xyz", tmp_path / "joint.xyz"
    src.write_text(ALL)
    join_paths(src, dst)
    assert dst.read_text() == "FRAME 1\n\nFRAME 2\n1 1 1\n2 2 2\n1 1 1\n3 3 3\n\n"


# ------------------------------------------------------ end-to-end example
@pytest.mark.skipif(not (EXAMPLES / "traj.xtc").exists(), reason="example data missing")
def test_example_run(tmp_path):
    from sprotein.config import load_config
    from sprotein.core import run

    cfg = load_config(EXAMPLES / "config.toml")
    cfg.output_dir = str(tmp_path)
    stats = run(cfg, progress=False, log=lambda *a: None)
    assert stats["frames"] == 11
    assert stats["frames_with_path"] == 1  # only frame 11 with edge_distance = 3.5
    nodes = read_frames(tmp_path / "traj_path_nodes.dat")
    path = nodes[10][1][0]
    assert path[0] == "8563 PAL 651" and len(path) == 7
    assert (tmp_path / "summary.csv").read_text().count("\n") == 12


# ------------------------------------------- crystal-like structure (water = O only)
def _write_crystal_pdb(path):
    atoms = [("N", "ALA", 1, (0, 0, 0)), ("CA", "ALA", 1, (1.5, 0, 0)),
             ("C", "ALA", 1, (3, 0, 0)), ("O", "ALA", 1, (3, 1.2, 0))]
    # chain of waters from the protein to a small "bulk" cluster around x = 21
    xyz = [(x, 1.2, 0) for x in (6, 9, 12, 15, 18, 21)] + [(21, 4.2, 0), (21, -1.8, 0), (24, 1.2, 0)]
    atoms += [("O", "HOH", 100 + k, p) for k, p in enumerate(xyz)]
    lines = [
        f"{'ATOM  ' if res == 'ALA' else 'HETATM'}{i:5d}  {name:<3} {res} A{resid:4d}    "
        f"{p[0]:8.3f}{p[1]:8.3f}{p[2]:8.3f}  1.00  0.00          {name[0]:>2}"
        for i, (name, res, resid, p) in enumerate(atoms, 1)
    ]
    path.write_text("\n".join(lines) + "\nEND\n")


def test_oxygen_only_water(tmp_path):
    from sprotein.config import Config
    from sprotein.core import PathFinder, run
    from sprotein.selection import SelectionError, detect_water
    import MDAnalysis as mda

    pdb = tmp_path / "crystal.pdb"
    _write_crystal_pdb(pdb)
    U = mda.Universe(str(pdb))
    w = detect_water(U)
    assert (w.resname, w.oxygen, w.hydrogens) == ("HOH", "O", ())

    cfg = Config(topology=str(pdb), trajectory="", reference="resname ALA and name O",
                 protein_cutoff=30.0, bulk_degree=4, pbc=False, output_dir=str(tmp_path / "out"))
    stats = run(cfg, progress=False, log=lambda *a: None)
    assert stats["frames"] == 1 and stats["frames_with_path"] == 1
    path = read_frames(tmp_path / "out" / "traj_path_nodes.dat")[0][1][0]
    assert len(path) == 7 and path[0].endswith("ALA 1") and path[-1].endswith("HOH 105")

    cfg.directed = True
    with pytest.raises(SelectionError, match="hydrogens"):
        PathFinder(U, cfg)
