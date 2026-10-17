# Reference papers

The PDFs are © American Physical Society and are not committed. Get local
copies with:

```
docs/papers/fetch_papers.sh
```

## Relativistic exact exchange (used by `rel=2`, `nucleus='sphere'`)

| File | Reference | Used for |
|---|---|---|
| `Kreibich_Gross_Engel_PRA57_138_1998.pdf` | T. Kreibich, E. K. U. Gross, E. Engel, "Approximate relativistic optimized potential method", Phys. Rev. A **57**, 138 (1998) | Relativistic KLI, Eq. (54), implemented in `atomx/xc/kli.cpp` with `DiracExchange`. Table I (18 closed-shell atoms, He–No) is the benchmark in `tests/regression_test.cpp`. |
| `Engel_et_al_PRA58_964_1998.pdf` | E. Engel, A. Facco Bonetti, S. Keller, I. Andrejkovics, R. M. Dreizler, "Relativistic optimized-potential method: Exact transverse exchange and Møller-Plesset-based correlation potential", Phys. Rev. A **58**, 964 (1998) | j-coupled Coulomb exchange between Dirac spinors, Eqs. (A15), (A17), (A20), implemented in `atomx/xc/exchange_operator.cpp`. Nuclear radius R = 1.0793 A^{1/3} + 0.73587 fm, implemented in `atomx/atom/nucleus.cpp`. |
| `Engel_et_al_PRA52_2750_1995.pdf` | E. Engel, S. Keller, A. Facco Bonetti, H. Müller, R. M. Dreizler, "Local and nonlocal relativistic exchange-correlation energy functionals: Comparison to relativistic optimized-potential-model results", Phys. Rev. A **52**, 2750 (1995) | The exact relativistic OEP that the RKLI results above are compared with (scanned PDF). |

## Non-relativistic KLI / OEP (benchmarks, not stored here)

- I. G. Ryabinkin, A. A. Kananenka, V. N. Staroverov, Phys. Rev. Lett.
  **111**, 013001 (2013), Table I. Values used in `tests/regression_test.cpp`.
  [arXiv:1306.3603](https://arxiv.org/abs/1306.3603)
- S. Kümmel, J. P. Perdew, Phys. Rev. Lett. **90**, 043004 (2003).
  [arXiv:cond-mat/0212426](https://arxiv.org/abs/cond-mat/0212426)
- J. B. Krieger, Y. Li, G. J. Iafrate, Phys. Rev. A **45**, 101 (1992) and
  **46**, 5453 (1992).
