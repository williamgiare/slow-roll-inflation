# Scientific validation

The package regression suite contains **453 tests**, all passed in the latest
complete run. From the repository root:

```bash
python -m pytest -q tests
```

Independent scientific checks and saved reports are organized as follows:

- **Canonical inflation:** [LO through N3LO](canonical_SRI/README.md), N10LO
  physical-mode tests, high-derivative LO tests and iterative/closed comparisons.
- **K-inflation:** [generic checks through N3LO](kinflation_SRI/README.md),
  N10LO physical-flow tests and the direct literature comparison.
- **Universal coefficients:** [analytic versus mode-equation comparisons](symbolic_coefficients/README.md)
  for all 138 coefficients through weight ten.
- **Background elimination:** [exact forward substitution](common/consistency_coverage.md)
  for all 64 published relations, 32 per model.
- **Shared solvers:** [independent references and mode equations](common/README.md)
  used across the order cases.

Reports record numerical residuals, tolerances and tested orders. Each guide
states its coverage; the published result range and the algebra-engine
comparison range are separate.

## Reproduce selected checks

- Shared symbolic, transport and mode-equation checks:

  ```bash
  python -m validation.run --suite all
  ```

- Canonical N3LO case:

  ```bash
  python -m validation.canonical_SRI.order_cases --p 3
  ```

- All provided k-inflation order cases:

  ```bash
  python -m validation.kinflation_SRI.run --all
  ```

- Summarize saved D1 coverage without rerunning the calculations:

  ```bash
  python -m scripts.validate_consistency_coverage
  ```

## Full coefficient reproduction

Independent analytic regeneration and comparison through weight ten can take
several hours:

```bash
python validation/symbolic_coefficients/run.py --orders 8 9 10 --family-seconds 0 --order-seconds 0 --fresh
```

New reports and checkpoints go to `symbolic_coefficients/reproduced/`;
published results are preserved. The [coefficient guide](symbolic_coefficients/README.md)
explains resuming a run and selecting the output directory.

See the [detailed validation guide](../docs/validation.md) for scientific
assumptions, comparison ranges and commands. Equation documents are in
[results/](../results/README.md).
