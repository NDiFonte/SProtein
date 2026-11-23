"""Post-processing tools working on the files written by ``sprotein run``."""

from __future__ import annotations

import statistics
from collections import Counter
from pathlib import Path

import numpy as np


def read_frames(path: Path) -> list[tuple[str, list[list[str]]]]:
    """Parse a SProtein output file into ``[(frame_header, [path_lines, ...]), ...]``.

    Files are made of blocks: a ``FRAME n`` line, then paths (one line per node)
    separated by blank lines.
    """
    frames: list[tuple[str, list[list[str]]]] = []
    current: list[str] = []
    for raw in Path(path).read_text().splitlines():
        line = raw.strip()
        if line.startswith("FRAME"):
            if frames and current:
                frames[-1][1].append(current)
            current = []
            frames.append((line, []))
        elif not line:
            if frames and current:
                frames[-1][1].append(current)
            current = []
        elif frames:
            current.append(line)
    if frames and current:
        frames[-1][1].append(current)
    return frames


# --------------------------------------------------------------------- VMD
def vmd_trajectory(src: Path, dst: Path, frame_zero: bool = True) -> int:
    """Convert ``traj_sp.xyz`` into a multi-frame XYZ that VMD can load.

    Every frame gets the same number of atoms (the longest path); path nodes are
    written as ``Se`` atoms and padding atoms sit at the origin. Frames without a
    path are filled with ``O`` dummies at the origin. Returns the atoms per frame.
    """
    frames = read_frames(src)
    paths = [p[0] if p else [] for _, p in frames]
    k = max((len(p) for p in paths), default=0)
    zero_o = "  O  0.000     0.000     0.000\n"
    zero_se = "  Se  0.000     0.000     0.000\n"
    with open(dst, "w") as out:
        if frame_zero:
            out.write(f"{k}\nFRAME 0\n" + zero_o * k)
        for (header, _), path in zip(frames, paths):
            out.write(f"{k}\n{header}\n")
            if path:
                out.writelines(f"  Se  {line}\n" for line in path)
                out.write(zero_se * (k - len(path)))
            else:
                out.write(zero_o * k)
    return k


# -------------------------------------------------------------------- join
def join_paths(src: Path, dst: Path) -> None:
    """Merge all paths of each frame into a single block (one point cloud per frame)."""
    with open(dst, "w") as out:
        for header, paths in read_frames(src):
            out.write(header + "\n")
            for p in paths:
                out.writelines(line + "\n" for line in p)
            out.write("\n")


# ----------------------------------------------------------------- lengths
def path_lengths(src: Path) -> dict:
    """Path lengths (edges = nodes - 1) for every path in ``src``."""
    frames = read_frames(src)
    lengths = [len(p) - 1 for _, paths in frames for p in paths]
    with_path = sum(1 for _, paths in frames if paths)
    return {
        "frames": len(frames),
        "frames_with_path": with_path,
        "paths": len(lengths),
        "lengths": lengths,
        "mean": float(np.mean(lengths)) if lengths else float("nan"),
        "std": statistics.stdev(lengths) if len(lengths) > 1 else float("nan"),
    }


def plot_lengths(lengths: list[int], save: Path | None = None, show: bool = True) -> None:
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.hist(lengths, bins=range(0, max(lengths) + 2), align="left", rwidth=0.8)
    ax.set_xlabel("Path length (edges)")
    ax.set_ylabel("Frequency")
    ax.set_title("Distribution of path lengths")
    fig.tight_layout()
    _finish(fig, save, show)


# ----------------------------------------------------------- frequencies
def residue_frequencies(src: Path) -> Counter:
    """Count residues listed in the neighbour columns (4th onward) of a neighbours file."""
    counter: Counter = Counter()
    for line in Path(src).read_text().splitlines():
        s = line.strip()
        if not s or s.startswith("FRAME"):
            continue
        counter.update(s.split()[3:])
    return counter


def plot_frequencies(counter: Counter, n: int = 10, save: Path | None = None, show: bool = True) -> None:
    import matplotlib.pyplot as plt

    top = counter.most_common(n)
    if not top:
        raise ValueError("No neighbour residues to plot.")
    labels, counts = zip(*top)
    fig, ax = plt.subplots(figsize=(max(7, 0.6 * len(labels) + 2), 6))
    ax.bar(labels, counts, color="tab:blue")
    ax.tick_params(axis="x", rotation=45, labelsize=14)
    ax.set_ylabel("Frequency", fontsize=16)
    ax.set_title(f"Top {len(labels)} residues along the paths", fontsize=18)
    fig.tight_layout()
    _finish(fig, save, show)


def _finish(fig, save, show):
    import matplotlib.pyplot as plt

    if save:
        fig.savefig(save, dpi=200)
    if show:
        plt.show()
    plt.close(fig)
