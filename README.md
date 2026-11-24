# SProtein

SProtein finds, frame by frame in a molecular dynamics trajectory, the **shortest
water/protein path** that connects a reference atom (e.g. a buried ligand or an
active site atom) to **bulk water**. It treats water oxygens near the protein, plus
any protein or ligand atoms you choose, as nodes of a graph. It then runs a
breadth-first search from the reference atom to every bulk water molecule.

The graph can be **undirected** (nodes closer than a distance cutoff are connected) or
**directed** (edges follow hydrogen bonds, donor H → acceptor).

> Version 2 is a rewrite of the original scripts (still available on the
> [`old`](https://github.com/NDiFonte/SProtein/tree/old) branch). It gives the same
> paths and is about 10x faster. All parameters live in a config file, and the
> interactive prompts are gone.

## Installation

```bash
git clone https://github.com/NDiFonte/SProtein.git
cd SProtein
pip install .            # or: pip install -e .  for development
```

This installs the `sprotein` command (you can also run `python -m sprotein`).
Requirements: Python ≥ 3.9, MDAnalysis ≥ 2.0, NumPy, SciPy, Matplotlib, tqdm.

## Quick start (example data included)

```bash
cd examples
sprotein run -c config.toml          # ~5 s on the 11-frame example
sprotein neighbors -c config.toml    # closest residues along each path
sprotein plot-freq -c config.toml -n 15
sprotein lengths -c config.toml
sprotein vmd -c config.toml          # then: vmd sprotein_output/VMD_path_traj.xyz
```

## Using it on your system

1. **Look at your topology** to find the names you need:

   ```bash
   sprotein inspect structure.gro              # residues, atom names, detected water model
   sprotein inspect structure.gro -r PAL       # atoms, resids and indices of one residue
   ```

2. **Create a config file** and edit it (every option is commented):

   ```bash
   sprotein init config.toml --reference "resname PAL and resid 651 and name O2"
   ```

3. **Run**:

   ```bash
   sprotein run -c config.toml
   ```

   Any option can also be given (or overridden) on the command line, e.g.
   `sprotein run -c config.toml --directed --stop 500 -o run_directed`.
   See `sprotein run --help`.

### Main parameters (`[graph]` section)

| option | default | meaning |
|---|---|---|
| `reference` | – | source atom of every path (MDAnalysis selection, must match one atom) |
| `selections` | `[]` | extra protein/ligand atoms used as nodes, e.g. `["resname GLU and name OE1 OE2"]` |
| `water` | `"auto"` | water residue name (auto-detects SOL, HOH, WAT, TIP3, ...) |
| `protein_cutoff` | 5.0 Å | water oxygens within this distance of the protein become nodes |
| `edge_distance` | 3.5 Å | nodes closer than this are connected |
| `bulk_degree` | 7 | waters with ≥ this many water neighbours are bulk (path targets) |
| `pbc` | `true` | periodic boundary conditions |
| `directed` | `false` | follow hydrogen bonds instead of plain distances |
| `source_is_donor` | `true` | directed graphs only |
| `hbond_angle` | 30° | directed graphs only: max H–donor–acceptor angle |

Frames can be restricted with `start`, `stop`, `step` in `[input]` (or `--start/--stop/--step`).

### Crystal structures (single PDB, waters without hydrogens)

SProtein also works on a single structure such as a crystal structure:

```toml
[input]
topology   = "3KS3.pdb"
trajectory = ""          # empty: single structure

[graph]
pbc         = false      # the CRYST1 cell is not a periodic simulation box
bulk_degree = 4          # resolved crystal waters are sparse: lower the bulk threshold
```

Waters with only the oxygen (e.g. `HOH` with atom `O`) are accepted for undirected
graphs. Directed graphs need hydrogens and give a clear error otherwise.
Crystal waters are far sparser than in a simulation, so `bulk_degree` (and possibly
`edge_distance`) usually need lower/higher values than for MD.

## Output (`sprotein_output/` by default)

| file | content |
|---|---|
| `summary.csv` | one row per frame: number of nodes/edges/bulk waters, number of exits, shortest path length and residues |
| `traj_sp.xyz` / `traj_path_nodes.dat` | shortest path of each frame: coordinates / `atom_index resname resid` (1-based atom index) |
| `traj_all_paths.xyz` / `traj_all_paths_nodes.dat` | every distinct path to bulk water (disable with `all_paths = false`) |
| `parameters.json` | the parameters used, for reproducibility |

Each `.xyz` / `.dat` file is made of `FRAME n` blocks, with paths separated by blank lines.
When several paths share the shortest length, one is picked at random using `seed`,
so runs are reproducible.

## Post-processing commands

| command | what it does |
|---|---|
| `sprotein neighbors` | adds the `n` closest protein residues (by centre of mass, within `cutoff`) to every node → `neighbors_traj_path_nodes.dat`. `--anchor PAL:C1` uses one atom instead of a residue's centre of mass |
| `sprotein plot-freq -n 15 [--save fig.png]` | bar plot of the residues that appear most often next to the paths |
| `sprotein lengths [--save hist.png]` | path length statistics and histogram (use `-i .../traj_all_paths.xyz` for all paths) |
| `sprotein vmd` | `VMD_path_traj.xyz`: a fixed-size XYZ trajectory of the shortest paths for VMD (`--no-frame-zero` to drop the dummy first frame) |
| `sprotein join` | merges all paths of a frame into one block (`joint_traj_all_paths.xyz`) |

They all read from `sprotein_output/` by default; use `-c config.toml` (to use the
config's output folder), `-d folder` or `-i file` to point elsewhere.

## Changes from version 1

* **Speed**: one breadth-first search per frame instead of one per bulk water
  molecule. One neighbour search per frame (KD-tree) instead of two dense N×N distance
  matrices. Vectorised hydrogen-bond angle filter. Selections are built once.
  On the example: ~5 s vs ~47 s, with about a third of the memory.
* **Same results**: checked against the original code on the example trajectory for
  undirected/directed graphs, donor/acceptor sources, with and without PBC: the
  adjacency matrices and the set of paths are identical in every frame.
* **Usability**: a single `sprotein` command with `--help`. TOML config instead of
  editing scripts. Automatic water detection instead of an interactive prompt.
  `inspect` helps you write selections. Clear error messages, progress bar,
  `summary.csv` and `parameters.json`. Frame slicing. Reproducible tie-breaking.
* **Fixes**: works with networkx ≥ 3 (removed dependency on networkx and dynetx).
  Histidine variants (HISD/HISE/HISH, ...) count as protein in directed graphs.
  `neighbors` labels residues correctly for non-contiguous residue numbering.
  The PAL→C1 special case is no longer hard-coded: use `anchors` / `--anchor`.

## Tests

```bash
pip install -e ".[test]"
pytest
```
