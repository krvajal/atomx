# atomx

All-electron density-functional calculations for spherical atoms, in C++20.
It is a port of the all-electron part of Quantum ESPRESSO's `ld1.x`, built
from small composable objects, with a pluggable linear-algebra backend
(LAPACK, or SLATE on CPUs/GPUs/many nodes).

| | |
|---|---|
| Kinetic energy | Schroedinger (`rel=0`), Dirac (`rel=2`) |
| Exchange-correlation | LDA (Slater + Perdew-Zunger), KLI, exact-exchange OEP, combinations such as `KLI+PZ` |
| Spin | unpolarised, spin-polarised (Hund's rule) |
| Nucleus | point, uniformly charged sphere (`rel=2`) |

Not ported: scalar-relativistic (`rel=1`), GGA/meta-GGA, SIC, pseudopotential
generation.

## Build and run

```
bazel build //apps:atomx
bazel test //tests:all
./bazel-bin/apps/atomx < input.in
```

Input is an ld1.x `&input` namelist, so ld1.x inputs run unchanged:

```
&input
  zed=80., config='default', rel=2, dft='KLI',
  nucleus='sphere', cau_fact=137.0359895,
  xmin=-7, dx=0.0125, rmax=60, beta=0.3,
/
```

Extra keys:
- `linear_solver = 'lapack' | 'slate' | 'slate-gpu'`
- `output = 'dat' | 'json' | 'dat+json'` (default `'none'`), with ld1.x's
  `prefix` (default `'atomx'`)

### Output files and plots

`output='dat'` writes gnuplot-ready columns with `#` headers:

| File | Columns |
|---|---|
| `<prefix>.wfc` | r, r·R(r) of each orbital (large component); eigenvalues in the header |
| `<prefix>.wfcs` | r, small components (Dirac only) |
| `<prefix>.rho` | r, 4πr²n(r) per spin |
| `<prefix>.pot` | r, V_nucleus, V_Hartree, V_xc per spin, V_KS per spin (Ry) |
| `<prefix>.scf` | iteration, residual |

`output='json'` writes all of this, plus energies and metadata, to
`<prefix>.json` for Python. `tools/plot_atom.py` makes a six-panel figure
from one or more runs:
- orbitals;
- radial density;
- r·V_xc, which compares KLI, OEP and LDA on one axis;
- potential terms;
- level diagram;
- SCF convergence.

```
uv venv .venv && uv pip install --python .venv/bin/python -r tools/requirements.txt
.venv/bin/python tools/plot_atom.py ar_KLI.json ar_OEP.json ar_PZ.json -o ar.png
.venv/bin/python tools/plot_atom.py hg.json --orbitals 1S+,2P-,2P+,6S+ -o hg.png
```

In gnuplot: `plot 'ne.wfc' u 1:2 w l t '1s', '' u 1:4 w l t '2p'` (with
`set logscale x`).

### SLATE backend

```
tools/build_slate.sh              # needs mpicxx, OpenBLAS, CUDA 12.6 for GPUs
bazel build --config=slate //apps:atomx
mpirun -np 2 ./bazel-bin/apps/atomx < input.in   # linear_solver='slate-gpu'
```

## Layout

```
atomx/core      Field/Spinor types, physical constants (identical to QE)
atomx/grid      RadialGrid (log mesh, integration), RadialPoisson
atomx/math      Clebsch-Gordan / Wigner 3j
atomx/linalg    LinearSolver interface; LAPACK and SLATE backends
atomx/atom      Orbital, Configuration, Nucleus, AtomState, element tables
atomx/solvers   RadialSolver (Schroedinger, Dirac), OrbitalResponse (OEP)
atomx/xc        XcModel interface: LDA, KLI, OEP, CompositeXcModel;
                ExchangeOperator (non-relativistic, Dirac)
atomx/scf       ScfDriver, Anderson Mixer, energy terms
atomx/io        namelist reader, CalculationInput
atomx/calc      Calculation: the composition root
atomx/report    ld1.x-style output; .dat/.json data writers
apps/           atomx executable
tests/          unit tests (analytic) and regression tests (published values)
tools/          compare_ld1.py, build_slate.sh, plot_atom.py
```

Every physics component is an interface injected into `ScfDriver`, so new
solvers, functionals, mixers or linear-algebra backends plug in without
changing the loop (see `Composition.KliWithInjectedLinearSolver` in
`tests/regression_test.cpp`).

## Validation

`tools/compare_ld1.py` runs atomx and ld1.x on the same inputs. LDA, LSDA,
KLI, spin-polarised KLI, Dirac LDA, relativistic KLI (point and finite
nucleus) and OEP give the same total energies and eigenvalues as ld1.x to
every printed digit, in the same number of SCF iterations. The exceptions
are runs whose SCF path is sensitive to last-bit differences. In those,
the energies agree to within the SCF threshold:
- Hg relativistic KLI with a finite nucleus;
- Au OEP: Etot agrees to 6e-6 Ry, and the atomx virial discrepancy is
  closer to zero.

The regression tests check published values: KLI and OEP from Ryabinkin
et al., PRL 111, 013001 (2013), and relativistic KLI from Kreibich, Gross,
Engel, PRA 57, 138 (1998).

## License

The numerical kernels are ports of Quantum ESPRESSO (GPL-2.0-or-later), so
this code is distributed under the same license.
