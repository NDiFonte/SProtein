"""Command line interface: ``sprotein <command> ...`` (or ``python -m sprotein``)."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import __version__
from .config import Config, ConfigError, apply_overrides, load_config, write_template

EPILOG = """\
typical workflow:
  sprotein inspect structure.gro               see residues, water model, atom indices
  sprotein init                                write a commented config.toml
  sprotein run -c config.toml                  find the paths (writes sprotein_output/)
  sprotein neighbors -c config.toml            closest residues along the shortest paths
  sprotein plot-freq -n 15                     most frequent residues
  sprotein lengths                             path length statistics
  sprotein vmd                                 XYZ trajectory of the paths for VMD
"""


def _abs(p):
    return str(Path(p).resolve()) if p else None


def _cfg(args) -> Config:
    cfg = load_config(getattr(args, "config", None))
    return cfg


def _outdir(args, cfg: Config) -> Path:
    return Path(args.dir) if getattr(args, "dir", None) else cfg.out


# ---------------------------------------------------------------- commands
def cmd_init(args):
    path = Path(args.path)
    if path.exists() and not args.force:
        sys.exit(f"'{path}' already exists (use --force to overwrite).")
    write_template(path, reference=args.reference or "")
    print(f"Wrote {path}. Edit the [graph] reference atom, then run:  sprotein run -c {path}")


def cmd_inspect(args):
    import MDAnalysis as mda
    from .selection import describe_topology

    U = mda.Universe(args.topology, *([args.trajectory] if args.trajectory else []))
    print(describe_topology(U, args.resname))


def cmd_run(args):
    from .core import run

    cfg = _cfg(args)
    apply_overrides(
        cfg,
        topology=_abs(args.topology), trajectory=_abs(args.trajectory), reference=args.reference,
        selections=args.selection, water=args.water, protein_cutoff=args.protein_cutoff,
        edge_distance=args.edge_distance, bulk_degree=args.bulk_degree, directed=args.directed,
        source_is_donor=args.source_is_donor, pbc=args.pbc, start=args.start, stop=args.stop,
        step=args.step, output_dir=_abs(args.output), seed=args.seed,
    )
    run(cfg, progress=not args.quiet)


def cmd_neighbors(args):
    import MDAnalysis as mda
    from .neighbors import annotate_neighbors

    cfg = _cfg(args)
    out = _outdir(args, cfg)
    src = Path(args.input) if args.input else out / "traj_path_nodes.dat"
    dst = Path(args.output) if args.output else src.with_name("neighbors_" + src.name)
    if not src.exists():
        sys.exit(f"'{src}' not found. Run `sprotein run` first (or pass --input).")
    if dst.exists() and not args.force:
        sys.exit(f"'{dst}' already exists (use --force to overwrite).")
    anchors = dict(cfg.anchors)
    for a in args.anchor or []:
        if ":" not in a:
            sys.exit(f"--anchor expects RESNAME:ATOM, got '{a}'.")
        res, atom = a.split(":", 1)
        anchors[res] = atom
    top = _abs(args.topology) or str(cfg.path(cfg.topology))
    trj = _abs(args.trajectory) or str(cfg.path(cfg.trajectory))
    U = mda.Universe(top, trj)
    annotate_neighbors(
        U, src, dst,
        n=args.n if args.n is not None else cfg.neighbors_n,
        cutoff=args.cutoff if args.cutoff is not None else cfg.neighbors_cutoff,
        anchors=anchors, progress=not args.quiet,
    )
    print(f"Wrote {dst}")


def cmd_plot_freq(args):
    from .tools import plot_frequencies, residue_frequencies

    out = _outdir(args, _cfg(args))
    src = Path(args.input) if args.input else out / "neighbors_traj_path_nodes.dat"
    if not src.exists():
        sys.exit(f"'{src}' not found. Run `sprotein neighbors` first (or pass --input).")
    counter = residue_frequencies(src)
    for label, count in counter.most_common(args.n):
        print(f"{label:<10} {count}")
    plot_frequencies(counter, args.n, save=args.save, show=not args.no_show)
    if args.save:
        print(f"Saved {args.save}")


def cmd_lengths(args):
    from .tools import path_lengths, plot_lengths

    out = _outdir(args, _cfg(args))
    src = Path(args.input) if args.input else out / "traj_sp.xyz"
    if not src.exists():
        sys.exit(f"'{src}' not found. Run `sprotein run` first (or pass --input).")
    s = path_lengths(src)
    if not s["lengths"]:
        sys.exit("No paths in the file.")
    pct = 100 * s["frames_with_path"] / s["frames"] if s["frames"] else 0
    print(f"Total frames:        {s['frames']}")
    print(f"Frames with a path:  {s['frames_with_path']} ({pct:.2f}%)")
    print(f"Paths:               {s['paths']}")
    print(f"Average path length: {s['mean']:.2f}")
    print(f"Standard deviation:  {s['std']:.2f}")
    lengths_file = src.with_name("lengths_" + src.stem + ".dat")
    lengths_file.write_text("".join(f"{x}\n" for x in s["lengths"]))
    print(f"Lengths written to {lengths_file}")
    if not args.no_plot:
        plot_lengths(s["lengths"], save=args.save, show=not args.no_show)


def cmd_vmd(args):
    from .tools import vmd_trajectory

    out = _outdir(args, _cfg(args))
    src = Path(args.input) if args.input else out / "traj_sp.xyz"
    dst = Path(args.output) if args.output else src.with_name("VMD_path_traj.xyz")
    if not src.exists():
        sys.exit(f"'{src}' not found. Run `sprotein run` first (or pass --input).")
    k = vmd_trajectory(src, dst, frame_zero=not args.no_frame_zero)
    print(f"Wrote {dst} ({k} atoms per frame). Load it in VMD with:  vmd {dst}")


def cmd_join(args):
    from .tools import join_paths

    out = _outdir(args, _cfg(args))
    src = Path(args.input) if args.input else out / "traj_all_paths.xyz"
    dst = Path(args.output) if args.output else src.with_name("joint_" + src.name)
    if not src.exists():
        sys.exit(f"'{src}' not found.")
    join_paths(src, dst)
    print(f"Wrote {dst}")


# ------------------------------------------------------------------ parser
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="sprotein",
        description="Shortest water/protein paths from a reference atom to bulk water in MD trajectories.",
        epilog=EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--version", action="version", version=f"sprotein {__version__}")
    sub = p.add_subparsers(dest="command", metavar="<command>")

    def add(name, func, help_):
        sp = sub.add_parser(name, help=help_, description=help_)
        sp.set_defaults(func=func)
        return sp

    def common(sp, config=True, outdir=True):
        if config:
            sp.add_argument("-c", "--config", help="TOML config file (see `sprotein init`)")
        if outdir:
            sp.add_argument("-d", "--dir", help="results folder (default: output directory from the config, "
                                                "or sprotein_output)")

    sp = add("init", cmd_init, "write a commented configuration file")
    sp.add_argument("path", nargs="?", default="config.toml")
    sp.add_argument("--reference", help="reference atom selection to put in the file")
    sp.add_argument("-f", "--force", action="store_true", help="overwrite an existing file")

    sp = add("inspect", cmd_inspect, "list residues, atom names and indices to help write selections")
    sp.add_argument("topology")
    sp.add_argument("trajectory", nargs="?")
    sp.add_argument("-r", "--resname", help="show atoms, resids and indices of this residue")

    sp = add("run", cmd_run, "find the shortest paths in every frame")
    common(sp, outdir=False)
    g = sp.add_argument_group("overrides (take precedence over the config file)")
    g.add_argument("-s", "--topology")
    g.add_argument("-t", "--trajectory")
    g.add_argument("-r", "--reference", help="reference atom selection (source node)")
    g.add_argument("--selection", action="append", help="extra node selection (repeatable)")
    g.add_argument("--water", help="water residue name (default: auto)")
    g.add_argument("--protein-cutoff", type=float)
    g.add_argument("--edge-distance", type=float)
    g.add_argument("--bulk-degree", type=int)
    g.add_argument("--directed", action=argparse.BooleanOptionalAction, default=None)
    g.add_argument("--source-is-donor", action=argparse.BooleanOptionalAction, default=None)
    g.add_argument("--pbc", action=argparse.BooleanOptionalAction, default=None)
    g.add_argument("--start", type=int)
    g.add_argument("--stop", type=int)
    g.add_argument("--step", type=int)
    g.add_argument("--seed", type=int)
    g.add_argument("-o", "--output", help="output directory")
    sp.add_argument("-q", "--quiet", action="store_true", help="no progress bar")

    sp = add("neighbors", cmd_neighbors, "closest protein residues to each path node")
    common(sp)
    sp.add_argument("-i", "--input", help="nodes file (default: <dir>/traj_path_nodes.dat)")
    sp.add_argument("-o", "--output", help="default: neighbors_<input>")
    sp.add_argument("-s", "--topology")
    sp.add_argument("-t", "--trajectory")
    sp.add_argument("-n", type=int, help="neighbours per node (default 1)")
    sp.add_argument("--cutoff", type=float, help="search radius in A (default 7)")
    sp.add_argument("--anchor", action="append", metavar="RES:ATOM",
                    help="use this atom instead of the centre of mass for residue RES (repeatable)")
    sp.add_argument("-f", "--force", action="store_true", help="overwrite the output file")
    sp.add_argument("-q", "--quiet", action="store_true")

    sp = add("plot-freq", cmd_plot_freq, "bar plot of the residues most often next to the path")
    common(sp)
    sp.add_argument("-i", "--input", help="default: <dir>/neighbors_traj_path_nodes.dat")
    sp.add_argument("-n", type=int, default=10, help="how many residues (default 10)")
    sp.add_argument("--save", help="save the figure (png, pdf, svg, ...)")
    sp.add_argument("--no-show", action="store_true", help="do not open a window")

    sp = add("lengths", cmd_lengths, "path length statistics and histogram")
    common(sp)
    sp.add_argument("-i", "--input", help="default: <dir>/traj_sp.xyz")
    sp.add_argument("--save", help="save the histogram")
    sp.add_argument("--no-show", action="store_true", help="do not open a window")
    sp.add_argument("--no-plot", action="store_true", help="statistics only")

    sp = add("vmd", cmd_vmd, "XYZ trajectory of the shortest paths, loadable in VMD")
    common(sp)
    sp.add_argument("-i", "--input", help="default: <dir>/traj_sp.xyz")
    sp.add_argument("-o", "--output", help="default: VMD_path_traj.xyz next to the input")
    sp.add_argument("--no-frame-zero", action="store_true", help="omit the dummy frame 0")

    sp = add("join", cmd_join, "merge all paths of each frame into one block")
    common(sp)
    sp.add_argument("-i", "--input", help="default: <dir>/traj_all_paths.xyz")
    sp.add_argument("-o", "--output")
    return p


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    if not getattr(args, "func", None):
        parser.print_help()
        return 1
    from MDAnalysis.exceptions import SelectionError as MDASelectionError

    from .selection import SelectionError

    try:
        args.func(args)
    except (ConfigError, SelectionError, MDASelectionError, FileNotFoundError) as e:
        sys.exit(f"error: {e}")
    except KeyboardInterrupt:
        sys.exit("interrupted")
    return 0


if __name__ == "__main__":
    sys.exit(main())
