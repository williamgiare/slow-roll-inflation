# Validation

Install SRI and the additional validation tools, then run from the repository root:

```bash
python -m pip install -e .
python -m pip install pytest matplotlib
python -m pytest -q tests
```

The package suite contains **453 tests**, all passed in the latest complete run.
Additional analytic-backend and independent physical checks can be run separately:

```bash
python -m pytest -q validation/symbolic_coefficients/test_backend.py
python -m validation.run --suite all
```

The shared suites can also be selected individually with `--suite symbolic`,
`transport`, `modes`, `jost` or `generic`. Assertions enforce the tolerances;
failed checks return a nonzero exit status. Reports retain measured residuals,
integration settings, dependency versions and calculation-source identifiers.

## Symbolic and coefficient checks

The scalar and tensor amplitudes through N3LO are compared monomial by monomial
with independently transcribed reference expressions. Numerical coefficients
are compared at explicit tolerances. Background inverse identities, weight
assignments, flow transformations and zero-`r` limits are checked exactly.

The consistency oracle constructs an observable flow and uses scalar Lagrange
inversion, independently of SRI's direct-`r` inverse. Comparisons include
`n=0,...,3`, `p=0,...,3`, as well as `n=8,p=0` and `n=5,p=1`.
Production coefficient generation is also tested in a process that blocks
validation imports.

Reference amplitudes are documented in
[published_spectra.py](../validation/common/symbolic/published_spectra.py).
The comparisons use the common conformal reference $-k_*\eta_*=1$.
Relevant references are [Stewart & Gong (2001)](https://arxiv.org/abs/astro-ph/0101225),
[Auclair & Ringeval (2022)](https://arxiv.org/abs/2205.12608), and
[Bianchi & Gamonal (2024)](https://arxiv.org/abs/2405.03157).

## Published consistency-relation coverage

Exact D1 forward substitution covers all **64 published relations**, 32 per
physical model: `n=0,...,3` for `p=0,...,6`, and `n=0` for `p=7,...,10`.
Every residual vanishes in rational internal algebra at $W=n+1+p$, using
mode-generated numerical coefficients before public rounding. This checks
observable inversion; coefficient accuracy has separate independent tests.

The [coverage guide](../validation/common/consistency_coverage.md) and
[YAML report](../validation/common/reports/consistency_coverage.yaml) identify
the saved evidence for every published pair. Collect the summary without
running scientific calculations:

```bash
python -m scripts.validate_consistency_coverage --output /tmp/consistency_coverage.yaml
```

Use `--run-missing` to execute only the canonical `n=3`, `p=5,6` checks,
reusing the checked coefficient bank, and refresh the report. These results
do not extend the separate iterative-versus-closed comparison ranges or
establish a complete end-to-end `--symbolic` sweep at all published orders.

## Independent numerical equations

The LO transport suite compares the integral with a tensor ODE on both sides
of the pivot, including logarithmic derivatives and transported reference
normalizations. Query points are shuffled.

The mode suite checks exact constant-`epsilon` Hankel amplitudes and independently
solves canonical scalar and tensor modes for a nonconstant background.
Background rescalings test the expected slow-roll remainder. Initial-time
and freeze-out refinements monitor integration accuracy.

The normalized determinant is compared with separate oscillating-mode and
growing-pump equations for compact sources of either sign. Quadrature refinement
uses 24, 48 and 96 nodes per interval.

## Order cases

| Case | Main check | Command |
| --- | --- | --- |
| [LO](../validation/canonical_SRI/LO/README.md) | Bell formula and Cortés–Liddle hierarchy | `python -m validation.canonical_SRI.order_cases --p 0` |
| [N1LO](../validation/canonical_SRI/N01LO/README.md) | Linear amplitude corrections | `python -m validation.canonical_SRI.order_cases --p 1` |
| [N2LO](../validation/canonical_SRI/N02LO/README.md) | Quadratic amplitude corrections | `python -m validation.canonical_SRI.order_cases --p 2` |
| [N3LO](../validation/canonical_SRI/N03LO/README.md) | Cubic amplitude corrections | `python -m validation.canonical_SRI.order_cases --p 3` |
| [N10LO](../validation/canonical_SRI/N10LO/README.md) | Order-ten coefficients, spectra and first running | `python -m validation.canonical_SRI.order_cases --p 10` |

The LO stress test compares generating-series and recursive evaluations through
`n=1000`, using exact rational arithmetic and three prime fields. Its symbolic
comparison with SRI reaches `n=8`:

```bash
python -m validation.canonical_SRI.LO.large_n_validation.run --max-n 1000
```

## Higher slow-roll orders

```bash
python -m validation.common.numerical.generic_order --maximum 7
python -m validation.canonical_SRI.high_order --p 10
```

Gamma/Mellin response coefficients, the closed two-insertion kernel and exact
Hankel limits supply independent targets beyond published low-order formulas.
The order-ten protocol checks 138 universal moments against a separately
implemented linear mode/pump solver. It tests scalar and tensor spectra at
every `p=0,...,10`, as well as 33 consistency relations for `n=0,1,2`, with
maximum absolute weight `W=13`.

The physical test uses
$a(x)=x^{-1}\exp[-e(x^b-1)/b]$, with $x=-k\eta$ and reference $x=1$.
Independent parameter-Taylor coefficients test the formal slow-roll order;
finite-parameter scalar and tensor modes test the resulting amplitudes.
The latter calculation is repeated at three background rescalings with refined
initial data, integration and freeze-out. The background is a numerical
validation case, rather than a fit to CMB data.

The N10LO notebook displays stored results by default; `RERUN=True` performs
the full calculation. Stored results are accepted only when their recorded
calculation sources match or a recorded compatibility check establishes that
the numerical calculation is unchanged. A converged generated bank may be supplied with
`--bank` to the order case or `--reuse-bank` to the high-order sweep; independent
checks still run and the generation-kernel signature must match.

These finite checks establish the stated tested ranges. The arbitrary-order
construction follows from the finite weighted basis and the derivation.

## Iterative versus closed algebra

`python -m validation.canonical_SRI.iterative_vs_closed.run` compares spectra through `p=4`,
relations through `p=4` for `n=0,...,3`, and LO through `n=8`.
The script checks intermediate maps, closed derivatives, amplitude polynomials,
final relations and forward substitution. All polynomial residuals are exactly
zero before rounding; the numerical moment data are common. This complements,
and does not replace, independent tests of the moment generator.
The notebook in `validation/canonical_SRI/iterative_vs_closed/` runs the same checks.

## Equation documents

The [results collection](../results/README.md) contains standalone LaTeX,
compiled PDF and YAML companions for every `p=0,...,10` in both models. The released expressions include
`n=0,...,3` through N06LO and `n=0` from N07LO to N10LO. These export ranges do
not limit the calculation or change the larger derivative ranges tested above.
Validation commands run the scientific checks; `python -m scripts.export_results`
generates equation documents separately.

## Variable scalar sound speed

The canonical order cases and algebra comparisons are under
`validation/canonical_SRI/`. Variable-sound-speed cases are under
`validation/kinflation_SRI/`:

```bash
python -m validation.kinflation_SRI.run --all
python -m validation.kinflation_SRI.run --p 10
```

Generic exact polynomial tests cover `p=0,...,3`, `n=0,...,3`, and LO through
`n=8`: forward substitution, closed-versus-iterative equality, canonical
recovery and regular zero-`r` limits. The scalar source is checked through a
separately constructed effective canonical tensor background.
N01LO is also compared with an independent analytic sound-speed relation.

The [direct N3LO literature test](../validation/kinflation_SRI/literature/README.md)
compares native scalar and tensor amplitude factors against Bianchi–Gamonal’s
official supplement. Both engines agree symbolically on all 31 scalar and
8 tensor monomials. A separate numerical-backend comparison checks each
coefficient at tolerance $10^{-9}$. This test uses the sector-specific acoustic
or conformal pivot and does not test common-pivot transport or inversion.
Run it with `python -m scripts.validate_kinflation_literature`;
`--output PATH` writes a YAML report.

Order-ten tests compare both engines on two physical flow families and compare
mixed `epsilon_1`/sound-flow Taylor coefficients with exact Hankel solutions.
The evolving profile has `q_j=(j-1)! q_1^j` and
`s_j=(j-1)! q_1^(j-1) s_1`; the constant-flow limit has `q_j=0` and
`s_j=0` for `j>=2`. These tests reach relative amplitude order ten without
using literature tables. They do not replace generic observable-inversion
checks; the [case guide](../validation/kinflation_SRI/README.md) and reports
identify each tested range explicitly.

The additional dense iterative tensor-tilt inverse has been checked by exact
forward substitution through `p=7`. The saved check is in
[generic_inverse.yaml](../validation/kinflation_SRI/generic_inverse.yaml). Reproduce the full sweep with
`python -m validation.kinflation_SRI.run --generic-max-p 7`.
This sweep is optional because its symbolic cost grows rapidly with order.

The public spectra and consistency relations use the same sound scale
parameters `alpha_j^cs`. Independent forward derivatives check both directions
of the polynomial sound-coordinate conversion through relative order ten.
The converted spectra recover the native amplitudes at the same SR cutoff,
and the independent NLO scalar and tensor amplitudes are also checked.

## Exact arithmetic and performance

```bash
python -m validation.common.performance
```

This check compares optimized calculations with direct polynomial arithmetic
using the same generated moment bank. It compares complete background maps,
log amplitudes and final tensor relations before rounding, then verifies exact
forward substitution. The tests also compare sparse products, substitutions,
exponentials and logarithms on rational polynomials with zero-weight constants.
The independent closed engine and the known low-order limits remain part of
the full test suite.

The command also compares mode-solver pump polynomials, initial-data arrays and
product bases through `p=10`, at asymptotic depths 12, 18 and 20. Numerical bank
reuse requires a matching generation kernel; integration settings and measured
convergence differences are recorded in the report.

The [performance report](../validation/common/reports/arithmetic_performance.yaml)
records exact residuals and timings for both physical models. Timings measure
cold background algebra with a common generated moment bank, excluding numerical
mode integration. They depend on order, derivative index and machine.

## Analytic coefficient validation

[The symbolic coefficient suite](../validation/symbolic_coefficients/README.md)
contains the completed independent comparisons through weights eight, nine and
ten: 66, 96 and 138 coefficients respectively. Weight ten includes 67 mixed
patterns. Each expression is exact, all extraction poles cancel identically,
and each numerical comparison satisfies its recorded uncertainty bound.

The public backend integration tests exercise both algebra engines and physical
models, exact low-order spectra and relations, numerical evaluation, mixed-coefficient
definitions, backend-separated caches and all CLI output formats. They forbid
numerical moment generation during symbolic calculations. These tests establish
their listed end-to-end cases; complete coefficient coverage does not imply a
complete end-to-end symbolic inversion sweep through weight ten.

```bash
python -m pytest -q tests/test_symbolic.py
python validation/symbolic_coefficients/run.py --orders 8 9 10 --family-seconds 0 --order-seconds 0 --fresh
```

The second command independently regenerates the coefficient comparisons and
can take several hours. Reproduced checkpoints and reports are written under
`validation/symbolic_coefficients/reproduced/`, preserving the published results.
Use `--output-root PATH` to select another output directory. Numerical reference
coefficients enter only after analytic generation; the reference fixture is from
independently integrated normalized mode equations.

## Archived YAML readers

`tests/test_archive_get.py` checks exact universal coefficient loading, saved
analytic and mode-equation values, model and order selection, and errors for
missing results. Reader calls are tested with calculation entry points disabled.
The tests compare the stored low-order spectra and relations with the production
calculation using a kernel-checked numerical bank, within the precision of the
published decimal coefficients. At relative orders six and ten they also check
that setting `c_s=1` and all sound derivatives to zero recovers the canonical
spectra and relation.

`tests/test_archive_export.py` checks new equation exports, YAML companions,
scientific notation and denominators in the polynomial reader, and rejection of
unsupported syntax. The [coefficient notebook](../results/universal_coefficients/universal_coefficients.ipynb)
and [archive notebook](../notebooks/03_archived_results.ipynb) demonstrate these
getters. The first two [calculation notebooks](../notebooks/README.md) exercise
the public calculation API. Reader verification is recorded
in [the archive reader report](../validation/common/reports/archive_readers.yaml).

```bash
python -m pytest -q tests/test_archive_get.py tests/test_archive_export.py
```
