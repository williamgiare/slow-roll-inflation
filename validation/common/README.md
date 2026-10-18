# Shared validation tools

`symbolic/` contains analytic comparison targets and an independent
observable-flow construction of the consistency hierarchy. `numerical/`
contains separate mode and growing-pump solvers, Hankel/Gamma limit tests,
parameter-Taylor checks and the LO integral/ODE comparison. `fixtures/`
holds reference data used by validation, and `reports/` stores numerical
residuals, execution settings and source identifiers.

These tools are shared by the per-order cases and are independent of the
production calculation. Published coefficients serve only as comparison
targets. These numerical cases generate the coefficients from the mode equations;
the production interface can separately select analytic coefficients with `--symbolic`.

```bash
python -m validation.run --suite all
python -m validation.common.numerical.generic_order --maximum 7
```

The [consistency coverage guide](consistency_coverage.md) documents exact D1
substitution for every published relation.
`python -m scripts.validate_consistency_coverage` builds the coverage summary
from saved reports without executing their scientific checks.

The [validation guide](../../docs/validation.md) describes the tests and ranges.

`python -m validation.common.performance` compares the optimized polynomial
operations with direct exact arithmetic and records timings for both models in
[arithmetic_performance.yaml](reports/arithmetic_performance.yaml).
It also checks mode-solver inputs through order ten. The report records exact
residuals, tested orders and timings.
