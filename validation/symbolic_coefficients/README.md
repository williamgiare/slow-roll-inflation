# Analytic universal coefficients

This suite compares exact analytic coefficients with values from independently
integrated normalized mode equations. A pattern `(j_1,...,j_m)` denotes
$\chi_{m;j_1,\ldots,j_m}$ and has slow-roll weight `m+sum(j_i)`.

The checks cover Hankel/Mellin families and general mixed-insertion families.
They require exact arithmetic, cancellation of negative auxiliary powers and
agreement within the recorded numerical uncertainty bounds.

## Completed comparisons

| Maximum weight | Coefficients | Mixed patterns | Report |
|---|---|---|---|
| 8 | 66 | 21 | [p08.yaml](reports/p08.yaml) |
| 9 | 96 | 39 | [p09.yaml](reports/p09.yaml) |
| 10 | 138 | 67 | [p10.yaml](reports/p10.yaml) |

Each listed set passes the independent comparison. Reports record per-coefficient
errors, uncertainty bounds and pole-cancellation checks. The numerical reference
is the generated mode-equation bank in `validation/common/fixtures/p10_moments.yaml`.
Its integration uncertainty determines the comparison accuracy; evaluating the
analytic expression to more digits does not improve that reference precision.

## Reproduce the complete comparison

From the repository root, using the environment in which SRI is installed:

```bash
python validation/symbolic_coefficients/run.py --orders 8 9 10 --family-seconds 0 --order-seconds 0 --fresh
```

> [!WARNING]
> Full analytic regeneration and comparison can take several hours.

New calculations are written to `reproduced/`; published reports are preserved.
An order report is written only when every required coefficient passes.
`reproduced/reports/sequence.yaml` records the overall status.

- `--orders`: weights to check, completed sequentially.
- `--fresh`: clear the selected reproduction directory's reports and expression
  checkpoints before starting. Omit it to resume a compatible partial run.
- `--family-seconds`, `--order-seconds`: execution budgets; `0` disables the limit.
  A stopped run preserves completed checks and is not marked as passed.
- `--output-root PATH`: choose another reproduction directory.

Family reports and expression checkpoints allow completed work to be reused.
Special-function values are cached for numerical evaluation after analytic
generation. Reading a published report alone does not populate these caches.

## Short regression tests

```bash
python -m pytest -q validation/symbolic_coefficients/test_backend.py
python -m pytest -q tests/test_symbolic.py
```

The first command checks the analytic validation implementation, including
constant-source limits, pole cancellation and a mixed coefficient against fresh
mode integration. The second checks the production `--symbolic` interface,
both models and algebra engines, exact output and representative end-to-end
calculations. These short tests do not regenerate every high-order family.

Complete coefficient coverage is separate from end-to-end spectra and
background-inversion coverage. The [validation guide](../../docs/validation.md#analytic-coefficient-validation)
states the tested ranges. Use the [coefficient archive](../../results/universal_coefficients/README.md)
to read saved expressions and numerical values.
