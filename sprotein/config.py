"""Run configuration: defaults, TOML loading and the template written by ``sprotein init``."""

from __future__ import annotations

import dataclasses
from dataclasses import dataclass, field, fields
from pathlib import Path

try:
    import tomllib
except ModuleNotFoundError:  # Python < 3.11
    try:
        import tomli as tomllib  # type: ignore
    except ModuleNotFoundError:  # pragma: no cover
        tomllib = None


@dataclass
class Config:
    # [input]
    topology: str = "structure.gro"
    trajectory: str = "traj.xtc"
    start: int | None = None
    stop: int | None = None
    step: int | None = None
    # [graph]
    reference: str = ""
    selections: list[str] = field(default_factory=list)
    water: str = "auto"
    protein_cutoff: float = 5.0
    edge_distance: float = 3.5
    bulk_degree: int = 7
    pbc: bool = True
    directed: bool = False
    source_is_donor: bool = True
    h_distance: float = 1.3
    hbond_angle: float = 30.0
    # [output]
    output_dir: str = "sprotein_output"
    all_paths: bool = True
    seed: int = 0
    # [neighbors]
    neighbors_n: int = 1
    neighbors_cutoff: float = 7.0
    anchors: dict[str, str] = field(default_factory=dict)
    # folder of the config file; relative paths are resolved against it
    base_dir: Path = field(default=Path("."), repr=False)

    def path(self, p: str) -> Path:
        """Resolve a path from the config relative to the config file's folder."""
        q = Path(p).expanduser()
        return q if q.is_absolute() else self.base_dir / q

    @property
    def out(self) -> Path:
        return self.path(self.output_dir)


SECTIONS = {
    "input": ["topology", "trajectory", "start", "stop", "step"],
    "graph": [
        "reference", "selections", "water", "protein_cutoff", "edge_distance", "bulk_degree",
        "pbc", "directed", "source_is_donor", "h_distance", "hbond_angle",
    ],
    "output": ["output_dir", "all_paths", "seed"],
    "neighbors": ["neighbors_n", "neighbors_cutoff", "anchors"],
}
# short names accepted inside [neighbors]
ALIASES = {("neighbors", "n"): "neighbors_n", ("neighbors", "cutoff"): "neighbors_cutoff", ("output", "directory"): "output_dir"}


class ConfigError(ValueError):
    pass


def load_config(path: str | Path | None) -> Config:
    cfg = Config()
    if path is None:
        return cfg
    path = Path(path)
    if not path.exists():
        raise ConfigError(f"Config file '{path}' not found. Create one with `sprotein init {path}`.")
    if tomllib is None:
        raise ConfigError("Reading TOML needs Python >= 3.11 or `pip install tomli`.")
    with open(path, "rb") as fh:
        try:
            data = tomllib.load(fh)
        except tomllib.TOMLDecodeError as e:
            raise ConfigError(f"{path}: {e}") from None
    valid = {f.name for f in fields(Config)}
    for section, values in data.items():
        if section not in SECTIONS or not isinstance(values, dict):
            raise ConfigError(f"{path}: unknown section [{section}] (expected one of: {', '.join(SECTIONS)}).")
        for key, value in values.items():
            name = ALIASES.get((section, key), key)
            if name not in valid or name not in SECTIONS[section]:
                raise ConfigError(f"{path}: unknown option '{key}' in [{section}].")
            setattr(cfg, name, value)
    cfg.base_dir = path.resolve().parent
    return cfg


def apply_overrides(cfg: Config, **overrides) -> Config:
    for key, value in overrides.items():
        if value is not None:
            setattr(cfg, key, value)
    return cfg


def validate(cfg: Config) -> None:
    if not cfg.reference:
        raise ConfigError(
            "No reference (source) atom given. Set `reference = \"...\"` in [graph] or pass --reference. "
            "Use `sprotein inspect <topology> --resname NAME` to find it."
        )
    for name in ("protein_cutoff", "edge_distance", "h_distance"):
        if getattr(cfg, name) <= 0:
            raise ConfigError(f"'{name}' must be positive.")
    if cfg.bulk_degree < 1:
        raise ConfigError("'bulk_degree' must be >= 1.")
    if isinstance(cfg.selections, str):
        cfg.selections = [cfg.selections]


TEMPLATE = '''\
# SProtein configuration file. Run with:  sprotein run -c {name}
# Paths are relative to this file. Selections use MDAnalysis syntax:
# https://docs.mdanalysis.org/stable/documentation_pages/selections.html

[input]
topology   = "structure.gro"     # structure / topology file (gro, pdb, tpr, psf, ...)
trajectory = "traj.xtc"          # trajectory file (xtc, trr, dcd, ...); "" for a single structure (e.g. a PDB)
# start = 0                      # optional frame slicing (0-based, stop excluded)
# stop  = 100
# step  = 1

[graph]
# Source node of every path (must select exactly one atom).
# Tip: `sprotein inspect structure.gro --resname PAL` lists atoms and indices.
reference = "{reference}"

# Extra protein/ligand atoms to include as graph nodes, e.g.
# selections = ["resname GLU and name OE1 OE2", "resname ASP and name OD1 OD2"]
selections = []

water          = "auto"   # water residue name, or "auto" (SOL, HOH, WAT, TIP3, ...)
protein_cutoff = 5.0      # A, water oxygens within this distance of the protein become nodes
edge_distance  = 3.5      # A, two nodes closer than this are connected
bulk_degree    = 7        # water nodes with at least this many water neighbours are "bulk" (path targets)
pbc            = true     # periodic boundary conditions

# Directed graphs follow hydrogen bonds (donor H -> acceptor) instead of plain distances.
directed        = false
source_is_donor = true    # only used when directed = true
hbond_angle     = 30.0    # degrees, max H-donor-acceptor angle (directed only)
h_distance      = 1.3     # A, H atoms this close to a selected atom are its hydrogens (directed only)

[output]
directory = "sprotein_output"
all_paths = true          # also write every distinct exit path, not only the shortest one
seed      = 0             # random seed to break ties between equally short paths

[neighbors]
n      = 1                # closest protein residues reported for each path node
cutoff = 7.0              # A, search radius (residue centres of mass)
# For these residues use one atom instead of the residue centre of mass:
# anchors = {{ PAL = "C1" }}
'''


def write_template(path: Path, reference: str = "") -> None:
    path.write_text(TEMPLATE.format(name=path.name, reference=reference))


def asdict(cfg: Config) -> dict:
    return dataclasses.asdict(cfg)
