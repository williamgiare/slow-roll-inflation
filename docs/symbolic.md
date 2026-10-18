# Exact universal coefficients

Use `--symbolic` in the CLI or `symbolic=True` in Python for exact analytic
universal coefficients. Numerical coefficients remain the default. Coefficient
selection is independent of the physical model and algebra method.

```bash
SRI -p 2 --symbolic
SRI -n 0 -p 2 --kinflation --symbolic --closed
```

## Reading exact output

Short coefficients are printed in terms of rational numbers, `C`, powers of
`pi` and zeta values, with `C=EulerGamma+log(2)-2`. Larger mixed coefficients
appear as `Chi((j_1,...,j_m))`, or $\chi_{m;j_1,\ldots,j_m}$ in LaTeX.
These are exact defined constants, not additional physical parameters.

For the source expansion $g(x)=\sum_j g_j(\ln x)^j$, the pattern
`(j_1,...,j_m)` has weight `m+sum(j_i)`. Patterns are symmetric; summing over
sorted patterns uses multiplicity `m!/product_j count(j)!`. For example,
`(0,1,1)` has weight five and multiplicity three.

Constant-source families follow from Hankel amplitudes, one-varying and
two-insertion families from Mellin generators, and general mixed families
from finite Goncharov multiple polylogarithms. The [module guide](modules.md#analytic-coefficient-backend)
describes their implementation.

## Expand or evaluate a coefficient

```python
from sri import spectra
from sri.analytic import moment, expand_coefficients, evaluate

chi = moment((0, 1, 1))
number = evaluate(chi, digits=24)
power = spectra(p=5, symbolic=True)
expanded = expand_coefficients(power.tensor_correction)
```

`moment` returns the full analytic expression. `expand_coefficients` replaces
compact `Chi` constants by their definitions; `evaluate` evaluates their
special functions. The largest families can require hours and substantial
memory. Compact output avoids that expansion during ordinary calculations.
Generated families and numerical special-function values are cached locally;
`SRI_SYMBOLIC_CACHE` selects the cache directory.

`--symbolic --simplify` collects exact factors. Numerical `.evaluate(...)`
methods on spectra and relations accept the same physical inputs with either
coefficient mode.

## Universal coefficient archive

The [YAML archive](../results/universal_coefficients/README.md) contains all
138 coefficients through weight ten, including complete expressions for
67 mixed patterns. `index.yaml` lists families; `values.yaml` lists saved
analytic and independent mode-equation values.

```python
from utils import get

chi = get.coefficient(0, 1, 1)
number = get.value(0, 1, 1)
independent = get.value(0, 1, 1, source="mode_equation")
```

Run from the repository checkout. These calls read stored expressions and
values without regenerating or evaluating coefficients. See the
[coefficient notebook](../results/universal_coefficients/universal_coefficients.ipynb).

## Scientific checks

All 138 coefficients through weight ten pass the independent mode-equation
comparison, including exact cancellation of the required extraction poles.
The [validation suite](../validation/symbolic_coefficients/README.md) records
per-coefficient errors, uncertainty bounds and reproduction commands.
Full independent reproduction can take several hours.

Weight ten is the completed coefficient validation range; the construction has
no fixed order ceiling. Separate [integration tests](validation.md#analytic-coefficient-validation)
check representative spectra, relations and both algebra engines. Published
spectrum and relation archives use numerical coefficients; the universal
coefficient archive contains exact expressions.
