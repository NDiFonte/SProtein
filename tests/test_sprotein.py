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
    assert stats["frames_with_path"] == 4  # frames 2, 9, 10, 11 (same as SProtein 1.x)
    nodes = read_frames(tmp_path / "traj_path_nodes.dat")
    first = nodes[1][1][0]
    assert first[0] == "8563 PAL 651" and len(first) == 10
    assert (tmp_path / "summary.csv").read_text().count("\n") == 12
