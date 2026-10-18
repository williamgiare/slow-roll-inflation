# Canonical slow-roll inflation

[LO](LO/README.md) tests the arbitrary-derivative Bell formula and its
Cortés–Liddle hierarchy, including [large-n checks](LO/large_n_validation/README.md).
[N01LO](N01LO/README.md), [N02LO](N02LO/README.md) and
[N03LO](N03LO/README.md) compare generated spectra with known amplitudes and
independently reconstructed consistency relations.
[N10LO](N10LO/README.md) tests higher-order moments, physical spectra and
observable inversion without requiring known high-order formulas.

The [iterative-versus-closed comparison](iterative_vs_closed/README.md)
checks the independent algebra engines using common generated moment data.

| Check | Command |
|---|---|
| Individual order case (`p=0,1,2,3,10`) | `python -m validation.canonical_SRI.order_cases --p 3` |
| Full slow-roll sweep through weight ten | `python -m validation.canonical_SRI.high_order --p 10` |
| Iterative versus closed algebra | `python -m validation.canonical_SRI.iterative_vs_closed.run` |
| LO hierarchy through `n=1000` | `python -m validation.canonical_SRI.LO.large_n_validation.run --max-n 1000` |

Full high-order regeneration can take hours. The order notebooks read saved
reports by default; set `RERUN=True` to repeat the calculation.

Exact D1 forward-substitution checks cover all 32 published canonical relations:
`n=0,...,3` at `p=0,...,6`, and `n=0` at `p=7,...,10`. The
[coverage guide](../common/consistency_coverage.md) links each case to its saved
zero-residual result and gives commands to inspect or reproduce the checks.

The notebooks display current saved checks or rerun the same validation scripts.
Independent numerical and symbolic tools reside in [../common/](../common/README.md).
