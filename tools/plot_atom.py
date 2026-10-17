#!/usr/bin/env python3
"""Plot atomx results (the <prefix>.json files written with output='json').

Usage:
    plot_atom.py RUN.json [RUN2.json ...] [-o figure.png]

With several files, potentials, densities and levels are overlaid, e.g. to
compare KLI with OEP or LDA for the same atom. Panels (Hartree atomic units):
  1. radial orbitals r*R(r) of the first run (small components dashed)
  2. radial density 4 pi r^2 n(r)
  3. r * V_xc(r): shell structure and the -1/r tail of exact exchange
  4. r * V(r) of the first run: nucleus, Hartree, xc and Kohn-Sham
  5. eigenvalue levels of every run
  6. SCF convergence
"""
import argparse
import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

RY_TO_HA = 0.5


def load(path):
    with open(path) as f:
        d = json.load(f)
    d["r"] = np.asarray(d["grid"]["r"])
    d["name"] = f"{d['element']} {d['dft']}" + (
        " Dirac" if d["relativity"] == "dirac" else "")
    return d


def plot_orbitals(ax, d, only=None):
    r = d["r"]
    for i, o in enumerate(d["orbitals"]):
        if o["occupation"] <= 0 or (only and o["label"] not in only):
            continue
        line, = ax.plot(r, o["large"], lw=1.2, label=o["label"])
        if "small" in o:
            ax.plot(r, o["small"], lw=0.9, ls="--", color=line.get_color())
    ax.set_xscale("log")
    ax.set_xlabel("r (bohr)")
    ax.set_ylabel("r R(r)")
    ax.set_title(f"Orbitals: {d['name']}")
    ax.axhline(0, color="k", lw=0.5)
    ax.legend(fontsize=7, ncol=2)


def plot_density(ax, runs):
    for d in runs:
        for s, rho in enumerate(d["density"]):
            label = d["name"] + ((" up", " down")[s] if d["nspin"] == 2 else "")
            ax.plot(d["r"], rho, lw=1.2, label=label)
    ax.set_xscale("log")
    ax.set_xlabel("r (bohr)")
    ax.set_ylabel(r"$4\pi r^2 n(r)$")
    ax.set_title("Radial density")
    ax.legend(fontsize=7)


def plot_xc(ax, runs):
    for d in runs:
        for s, v in enumerate(d["potentials_ry"]["xc"]):
            label = d["name"] + ((" up", " down")[s] if d["nspin"] == 2 else "")
            rv = d["r"] * np.asarray(v) * RY_TO_HA
            # Orbital-dependent potentials are only defined where the density
            # is non-zero (the solvers set far tails to exactly 0).
            rv[np.asarray(d["density"][s]) <= 0.0] = np.nan
            ax.plot(d["r"], rv, lw=1.2, label=label)
    ax.axhline(-1.0, color="k", lw=0.6, ls=":", label=r"$-1/r$ limit")
    ax.set_xscale("log")
    ax.set_xlabel("r (bohr)")
    ax.set_ylabel(r"$r\,V_{xc}(r)$ (Ha bohr)")
    ax.set_title("Exchange-correlation potential")
    ax.legend(fontsize=7)


def plot_potentials(ax, d):
    r, p = d["r"], d["potentials_ry"]
    empty = np.asarray(d["density"][0]) <= 0.0  # see plot_xc

    def masked(v):
        rv = r * np.asarray(v) * RY_TO_HA
        rv[empty] = np.nan
        return rv

    ax.plot(r, r * np.asarray(p["nucleus"]) * RY_TO_HA, label="nucleus")
    ax.plot(r, r * np.asarray(p["hartree"]) * RY_TO_HA, label="Hartree")
    ax.plot(r, masked(p["xc"][0]), label="xc")
    ax.plot(r, masked(p["kohn_sham"][0]), "k", lw=1.6, label="Kohn-Sham")
    ax.set_xscale("log")
    ax.set_xlabel("r (bohr)")
    ax.set_ylabel(r"$r\,V(r)$ (Ha bohr)")
    ax.set_title(f"Potentials: {d['name']}")
    ax.legend(fontsize=7)


def plot_levels(ax, runs):
    for x, d in enumerate(runs):
        levels = sorted((o["energy"], o["label"]) for o in d["orbitals"]
                        if o["occupation"] > 0)
        for i, (e, label) in enumerate(levels):
            e *= RY_TO_HA
            ax.hlines(e, x - 0.35, x + 0.35, lw=1.5)
            # alternate sides so that near-degenerate labels do not overlap
            ax.text(x + (0.38 if i % 2 else -0.38), e, label, fontsize=6,
                    va="center", ha="left" if i % 2 else "right")
    ax.set_yscale("symlog", linthresh=1.0)
    ax.set_xticks(range(len(runs)), [d["name"] for d in runs], fontsize=7)
    ax.set_xlim(-0.8, len(runs) - 0.2)
    ax.set_ylabel("eigenvalue (Ha)")
    ax.set_title("Levels")


def plot_scf(ax, runs):
    for d in runs:
        res = d["scf"]["residuals"]
        ax.semilogy(range(1, len(res) + 1), res, marker=".", label=d["name"])
    ax.set_xlabel("SCF iteration")
    ax.set_ylabel("mean squared potential change (Ry$^2$)")
    ax.set_title("SCF convergence")
    ax.legend(fontsize=7)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("runs", nargs="+", help="atomx JSON files")
    p.add_argument("-o", "--out", default="atomx.png", help="output image")
    p.add_argument("--orbitals", default="",
                   help="comma-separated labels to show in the orbital panel, "
                        "e.g. 6S+,5D-,5D+ (default: all)")
    args = p.parse_args()
    runs = [load(f) for f in args.runs]

    fig, axes = plt.subplots(2, 3, figsize=(16, 9), constrained_layout=True)
    only = {x.strip() for x in args.orbitals.split(",") if x.strip()}
    plot_orbitals(axes[0, 0], runs[0], only)
    plot_density(axes[0, 1], runs)
    plot_xc(axes[0, 2], runs)
    plot_potentials(axes[1, 0], runs[0])
    plot_levels(axes[1, 1], runs)
    plot_scf(axes[1, 2], runs)
    e = ", ".join(f"{d['name']}: {d['energies_ry']['total'] * RY_TO_HA:.6f} Ha"
                  for d in runs)
    fig.suptitle(e, fontsize=10)
    fig.savefig(args.out, dpi=130)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
