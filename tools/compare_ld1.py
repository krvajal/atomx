#!/usr/bin/env python3
"""Run atomx and ld1.x on the same &input files and compare the results.

Usage: compare_ld1.py ATOMX LD1_X WORKDIR CASE [CASE ...]

Each CASE is "name: namelist body", e.g.
    "ne_kli: zed=10., config='1s2 2s2 2p6', rel=0, dft='KLI'"
Prints the largest eigenvalue and total-energy differences per case.
"""
import os
import re
import subprocess
import sys


def run(exe, inp, cwd):
    os.makedirs(cwd, exist_ok=True)
    return subprocess.run([exe], input=inp, capture_output=True, text=True,
                          cwd=cwd).stdout


def parse(out):
    energy = {}
    for key in ("Etot", "Ekin", "Encl", "Eh", "Exc"):
        m = re.search(r"\b" + key + r"\s*=\s*(\S+) Ry", out)
        energy[key] = float(m.group(1)) if m else None
    # eigenvalue table before any "Averaged results" block
    head = out.split("Averaged results")[0]
    eig = [float(e) for e in re.findall(
        r"^\s+\d \d(?: \d\.\d)?\s+\d[SPDF] \d\(\s*[-\d.]+\)\s+(\S+)", head, re.M)]
    m = re.search(r"final scf error:\s+(\S+) reached in\s+(\d+)", out)
    return energy, eig, (m.group(1), int(m.group(2))) if m else None


def main():
    atomx, ld1, work = sys.argv[1:4]
    worst = 0.0
    for case in sys.argv[4:]:
        name, body = case.split(":", 1)
        inp = "&input\n " + body.strip() + "\n/\n"
        ea, va, sa = parse(run(atomx, inp, os.path.join(work, name, "atomx")))
        el, vl, sl = parse(run(ld1, inp, os.path.join(work, name, "ld1")))
        if ea["Etot"] is None or el["Etot"] is None or len(va) != len(vl):
            print(f"{name:14s} FAILED  atomx={ea['Etot']} ld1={el['Etot']} "
                  f"n_eig={len(va)}/{len(vl)}")
            worst = float("inf")
            continue
        de = max(abs(ea[k] - el[k]) for k in ea)
        dv = max((abs(a - b) for a, b in zip(va, vl)), default=0.0)
        worst = max(worst, de)
        print(f"{name:14s} Etot {ea['Etot']:16.6f} Ry  max|dE| {de:.1e} Ry  "
              f"max|de_i| {dv:.1e} Ry  iter atomx/ld1 {sa[1] if sa else '-'}/"
              f"{sl[1] if sl else '-'}")
    print(f"worst energy difference: {worst:.1e} Ry")


if __name__ == "__main__":
    main()
