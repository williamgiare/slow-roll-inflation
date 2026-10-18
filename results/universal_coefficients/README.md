# Universal source coefficients

This collection contains the 138 universal coefficients through slow-roll weight
ten, including 67 mixed patterns. Exact analytic expressions and their numerical
values are available in YAML. The independent mode-equation values are included
for comparison. The coefficients are shared by the canonical and k-inflation models.

- [index.yaml](index.yaml) lists every family and its expression file.
- [values.yaml](values.yaml) lists every coefficient, its weight, permutation
  multiplicity, analytic numerical value and independent mode-equation value.
- [expressions/](expressions/) contains one YAML file per family, grouped by weight.
- [universal_coefficients.ipynb](universal_coefficients.ipynb) demonstrates exact
  coefficient access and saved numerical values.

## Conventions

For `g(x)=sum_j g_j (log x)^j`, a sorted pattern `(j_1,...,j_m)` denotes
`chi_{m;j_1,...,j_m}`. Its slow-roll weight is `w=m+sum(j_i)`. Each coefficient
appears only at its own weight. Lower-order coefficients are included when
assembling a calculation through a larger cutoff. At leading amplitude order
`p=0` no source-insertion coefficient is needed.

The symmetric coefficient has permutation multiplicity
`m!/product_j count(j)!` when summed over sorted patterns. Thus `(0,1,1)` has
weight five and multiplicity three. A family `m03_j07.yaml` has three insertions
and total logarithmic degree seven, hence weight ten.

The constant `C` means `EulerGamma + log(2) - 2`; `pi` is the usual circular
constant. `zeta`, `dirichlet_eta`, `log`, `polylog`, `gamma` and `polygamma` denote
the standard mathematical functions. Rational numbers and integer powers in
`exact_expression` are exact, and `**` denotes exponentiation.

## Reading a family

Each `coefficients` entry contains its pattern, permutation multiplicity, both
numerical values and the complete `exact_expression`. Each subsequent line in
the expression starts with `+` and adds one signed monomial. Factors on a line
are separated by ` * `, and the leading factor is an exact rational number.
The entire block is one finite sum; no derivatives or extraction operations
remain to be performed.

Mixed families also define a local `gpl_words` dictionary. For example,

```yaml
gpl_words:
  G000001: ['0', '-1', '-1/2']
```

means `G000001 = G(0,-1,-1/2;1)`, the Goncharov multiple polylogarithm with
kernels `du/(u-a_i)` and the analytic backend's tangential endpoint convention.
These labels are local to their family file. Every label appearing in a formula
has its complete rational word in that file. They abbreviate standard special
functions, not undetermined source coefficients.

All numerical values are quoted strings to preserve their recorded digits.
`analytic_value` is the saved 24-digit evaluation of the exact expression;
`mode_equation_value` is the saved floating-point integration result. The latter's
accuracy is governed by its integration uncertainty, not the number of analytic
digits. Individual errors and comparison bounds belong to the
[validation reports](../../validation/symbolic_coefficients/README.md).

## Load exact expressions in Python

Run from the repository checkout with SRI installed:

```python
from utils import get

chi = get.coefficient(0, 1, 1)
value = get.value(0, 1, 1)
independent = get.value(0, 1, 1, source="mode_equation")
print(chi)
```

The reader constructs exact SymPy expressions from the documented syntax.
Reading a coefficient does not regenerate its family or evaluate its special
functions. `get.value` returns saved decimal values as `Decimal` objects.
`coefficient` and `value` accept separate integers or a single tuple/list.
Permuting the pattern does not change the coefficient. A missing entry raises
an error without starting a calculation. Use `root=PATH` to select another
`results/` directory. Loaded families and expressions are cached in memory;
`get.clear_cache()` releases them and reloads edited files on the next request.

Large families require more memory and loading time. Use `sri.analytic.evaluate`
for explicit arbitrary-precision evaluation; evaluating the largest families
can take hours.

## Coefficient notebook

Open [universal_coefficients.ipynb](universal_coefficients.ipynb) in this directory
and run its cells. Change `pattern` to select a coefficient. The notebook can be
run from this directory or the repository root; it locates the repository automatically.
The standard `python -m pip install -e .` installation is sufficient when the
selected Python environment already supports notebooks. See
[notebook setup](../../docs/installation.md#notebooks) for optional kernel and
execution dependencies.

## Other archived results

The same utility reads spectra and relations:

```python
power = get.spectra(p=3)
relation = get.consistency(n=0, p=3, type="kinflation")
```

`power.scalar` and `power.tensor` include the normalizations;
`relation.expression` contains the selected relation. These archives keep the
published numerical coefficients. The [API guide](../../docs/api.md#archived-yaml-results)
lists all fields and stored orders; the [archive notebook](../../notebooks/03_archived_results.ipynb)
shows both models.

For new spectra and relations with exact coefficients, add `--symbolic`,
for example `SRI -p 2 --symbolic`.
See the [analytic coefficient guide](../../docs/symbolic.md). The ordinary
[spectrum and relation collections](../README.md) retain their numerical
coefficient presentation. Independent regeneration of the mathematical tests
is described in [validation](../../validation/symbolic_coefficients/README.md).

## Export an independently reproduced bank

After running the coefficient comparison script, the same YAML layout can be
exported from its completed reports and expression checkpoints:

```bash
python -m scripts.export_universal_coefficients \
  --reports validation/symbolic_coefficients/reproduced/reports \
  --checkpoints validation/symbolic_coefficients/reproduced/cache
```

The exporter requires completed comparisons and existing expression checkpoints.
It does not generate missing mixed coefficients or evaluate their constants.
The mathematical regeneration command is documented in the validation guide.
