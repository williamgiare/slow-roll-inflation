# Python API

The two calculation functions return SymPy expressions:

```python
from sri import spectra, consistency

power = spectra(p=2)
relation = consistency(n=2, p=2)
print(power.scalar)
print(power.tensor)
print(relation.expression)
```

`spectra` calculates both primordial spectra. `consistency` calculates one
tensor relation; it accepts any two of `n`, `p` and `W`, with `W=n+1+p`.
`n=0` selects the tensor tilt, `n=1` its running and `n=2` its running of running.
`ConsistencyOrder.resolve(n=..., p=..., W=...)` provides the same order inference
without running a calculation.

## Conventions

- **Reference:** `x=-k*eta=1`. This is the conformal reference used by the formulas;
  `k=aH` is a different convention. Tensor power includes both polarizations.
- **Background:** `epsilon_1=-d ln(H)/dN` and
  `epsilon_(j+1)=d ln|epsilon_j|/dN`, with increasing `N=ln(a)`.
  `H` is the Hubble rate and `Mpl` the reduced Planck mass.
- **Spectral inputs:** `alpha_0^s=n_s-1`, `alpha_0^t=n_t`;
  `alpha_j^s=d^j(n_s-1)/d(ln k)^j` and `alpha_j^t=d^j n_t/d(ln k)^j` for `j>=1`.
- **Sound inputs:** `alpha_0^cs=d ln(c_s)/d ln(k)` and
  `alpha_j^cs=d^j alpha_0^cs/d(ln k)^j` for `j>=1`, at the common tensor reference.
- **Constant:** `C=EulerGamma+log(2)-2`, substituted automatically by `.evaluate(...)`.

## Calculation options

Both functions accept the same independent options:

| Argument | Selection | Default |
|---|---|---|
| `symbolic` | `True`: exact analytic coefficients; `False`: numerical coefficients | `False` |
| `method` | `"iterative"` or `"closed"` background algebra | `"iterative"` |
| `kinflation` | `True`: variable scalar sound speed with `c_t=1`; `False`: canonical inflation | `False` |
| `coefficient_tolerance` | Numerical universal-moment convergence tolerance; ignored with `symbolic=True` | `1e-10` |

```python
power = spectra(p=2, method="closed", symbolic=True)
relation = consistency(n=0, p=2, kinflation=True, symbolic=True)
```

Orders can be arbitrary finite nonnegative integers. Increasing them requires
more computation and memory. All physical parameters remain symbolic in either
coefficient mode.

## ConsistencyResult

- `expression`: the tensor relation in terms of `r`, scalar derivatives and,
  for k-inflation, `c_s` and its derivatives.
- `n`, `p`, `W`, `error_order`: derivative index, relative correction order,
  absolute cutoff and remainder weight (`W+1`).
- `components`: `(absolute weight, expression)` pairs.
- `scalar_symbols`, `r_symbol`, `sound_symbols`, `cs_symbol`: expression symbols.
- `model`, `method`, `coefficient_metadata`: model, algebra method and coefficient information.
- `latex()`: the relation with its formal remainder.
- `evaluate(r=..., scalar=[...], cs=..., sound=[...])`: numerical evaluation;
  `cs` and `sound` are needed only for k-inflation.

The scalar input sequence is `[n_s-1, alpha_s, beta_s, ...]`. Supply every
entry used by the expression. `r` must be finite and nonnegative; the `r=0`
limit is regular.

```python
running = consistency(n=1, p=0)
value = running.evaluate(r=0.02, scalar=[-0.035])
```

## SpectraResult

- `scalar`, `tensor`: full spectra, including leading normalizations.
- `scalar_correction`, `tensor_correction`: relative amplitude factors starting with one.
- `log_scalar_correction`, `log_tensor_correction`: truncated logarithms of those factors.
- `p`, `model`, `method`, `coefficient_metadata`: order and calculation information.
- `sound_symbols`, `cs_symbol`: sound-speed symbols for k-inflation.
- `evaluate(H=..., epsilon=[...], Mpl=1, cs=..., sound=[...])`: numerical scalar
  and tensor powers.

`epsilon` is ordered as `[epsilon_1, epsilon_2, ...]`; at order `p`, the scalar
spectrum can need `p+1` entries. `H`, `Mpl` and `epsilon_1` must be positive
and finite. Higher flow parameters may vanish. Set `Mpl=1` when `H` is given
in reduced Planck units.

```python
values = power.evaluate(H=1e-5, epsilon=[0.005, 0.02, 0.01])
```

## K-inflation API

Use `kinflation=True` for new calculations:

```python
power_cs = spectra(p=2, kinflation=True)
values_cs = power_cs.evaluate(H=1e-5, epsilon=[0.005, 0.02, 0.01],
                              cs=0.7, sound=[0.002, 0.0001])
relation_cs = consistency(n=0, p=2, kinflation=True, method="closed")
value_cs = relation_cs.evaluate(r=0.02, scalar=[-0.035, 0.001],
                                cs=0.7, sound=[0.002, 0.0001])
```

`cs` must be positive and finite. `sound=[alpha_0^cs, alpha_1^cs, ...]` supplies
the sound derivatives actually used by the expression. The result records
`model="kinflation"`; canonical results use `model="canonical"`.
The leading relation is `n_t=-r/(8*c_s)`. Sound-speed data are independent
inputs; the scalar amplitude and tensor spectral parameters are not inputs
to the relation.

## Analytic coefficient selection

`symbolic=True` returns exact coefficients. Short expressions are explicit;
larger mixed expressions use defined `Chi(pattern)` constants. Their weight is
`len(pattern)+sum(pattern)`. Inspect `coefficient_metadata['coefficient_backend']`
for `"numeric"` or `"symbolic"`.

```python
from sri.analytic import moment, expand_coefficients, evaluate

chi = moment((0, 1, 1))
number = evaluate(chi, digits=24)
power_exact = spectra(p=5, symbolic=True)
expanded = expand_coefficients(power_exact.tensor_correction)
```

`moment` returns a complete exact coefficient. `expand_coefficients` replaces
`Chi` constants by their definitions, and `evaluate` evaluates the resulting
special functions at the requested precision. Large families can require hours;
use [archived expressions](#archived-yaml-results) to read existing results.
See the [analytic guide](symbolic.md) for conventions and caching.

## Higher-order coefficient accuracy

Numerical moments enter the algebra as exact rational representations of their
converged values. Public numerical coefficients are rounded to 14 significant
digits. This display precision is separate from the convergence tolerance,
which applies to individual moments rather than every combined output coefficient.
Inspect `coefficient_metadata` for settings and measured convergence differences.
Refined-integration agreement is a numerical estimate, not an interval bound.

For a simplified display copy, use
`sri.presentation.simplify_expression(expression, exact=True)` with exact
coefficients. `exact=False` additionally permits numerical fraction recognition.
The return object contains `expression` and `metadata`. The original calculation
result is preserved; numerical recognition does not certify an analytic identity.

## Archived YAML results

`utils.get` reads the published YAML archives from the repository checkout:

```python
from utils import get

chi = get.coefficient(0, 1, 1)
number = get.value(0, 1, 1)
saved_power = get.spectra(p=3)
saved_relation = get.consistency(n=0, p=3, type="kinflation")
```

- `coefficient`: complete exact SymPy coefficient; accepts separate indices or a tuple/list.
- `value`: stored analytic `Decimal`; use `source="mode_equation"` for the
  independent numerical value.
- `spectra`: stored `scalar`, `tensor`, amplitude corrections, `p`, `type` and `metadata`.
- `consistency`: stored `expression`, `components`, `n`, `p`, `W`, `type` and `metadata`.

`type="canonical"` is the default; `type="kinflation"` selects the sound-speed
archive. Both models have spectra at `p=0,...,10`, relations at `n=0,...,3`
for `p=0,...,6`, and `n=0` for `p=7,...,10`.

Readers load stored expressions without running the generators. Missing results
raise an error. Spectrum and relation coefficients retain the precision of the
published PDFs; universal coefficients are exact. Use `root=PATH` to select
another `results/` directory and `get.clear_cache()` to release loaded data.
See the [archive notebook](../notebooks/03_archived_results.ipynb).

## Notebook equation display

```python
from utils.notebook_display import show_spectra, show_consistency, show_equation

show_spectra(power)
show_consistency(relation)
show_equation(r"\alpha_2^t", relation.expression, digits=None)
```

These helpers display expressions in short aligned rows, grouped by slow-roll
weight. Floating-point coefficients use six significant digits by default;
`digits=None` displays all stored digits. Exact constants and returned
expressions are preserved. The [spectra](../notebooks/01_spectra.ipynb) and
[consistency](../notebooks/02_consistency_relations.ipynb) notebooks demonstrate
calculation, evaluation and simplification.
