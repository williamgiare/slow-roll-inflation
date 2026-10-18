# Consistency-relation coverage

Exact forward-substitution checks cover all 64 published consistency relations:

| Model | Relative order `p` | Tensor derivative `n` | Relations |
| --- | --- | --- | --- |
| Canonical | 0–6 | 0–3 | 28 |
| Canonical | 7–10 | 0 | 4 |
| K-inflation | 0–6 | 0–3 | 28 |
| K-inflation | 7–10 | 0 | 4 |

For each relation $R_{n,p}$, the D1 check substitutes the forward background
expressions for $r$, the scalar observables and, in k-inflation, the sound
observables. It then subtracts the directly calculated tensor observable.
The residual vanishes identically after truncation at absolute slow-roll
weight $W=n+1+p$.

These checks use numerical mode-generated coefficients represented in exact
rational internal algebra, before public rounding. Exact cancellation tests
the observable inversion; independent moment and analytic-coefficient tests
validate the coefficient values. The table specifies the released range and
does not certify every derivative index at every order or the complete
`--symbolic` pipeline at those orders.

The [coverage report](reports/consistency_coverage.yaml) identifies the saved
zero-residual evidence for every published pair `(n,p)`. The selected canonical
checks at `n=3`, `p=5,6` are in
[canonical_n3_p5_p6.yaml](reports/canonical_n3_p5_p6.yaml).
Reports retain coefficient-bank identifiers, source hashes and measured timings.

From the repository root:

```bash
python -m scripts.validate_consistency_coverage --output /tmp/consistency_coverage.yaml
```

This reads saved evidence and checks that every published relation has an exact
zero-residual result. It performs no relation calculation or moment integration.
To execute only the canonical `n=3`, `p=5,6` checks and refresh the summary:

```bash
python -m scripts.validate_consistency_coverage --run-missing
```

The selected checks reuse the kernel-checked saved coefficient bank and prohibit
new numerical moment integration. The reproduction commands for the other
ranges remain in the [canonical](../canonical_SRI/README.md) and
[k-inflation](../kinflation_SRI/README.md) guides.

The independent iterative-versus-closed comparison has its own tested ranges:
canonical spectra and relations through `p=4`, and generic k-inflation relations
through `p=3`, with `n=0,...,3` and additional LO derivative checks. These ranges
are distinct from the forward-substitution coverage above.
