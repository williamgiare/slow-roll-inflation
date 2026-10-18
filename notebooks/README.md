# Notebook examples

- [Spectra](01_spectra.ipynb): calculate scalar and tensor spectra at low order,
  compare iterative and closed methods, select analytic coefficients or variable
  sound speed, evaluate numerically and format LaTeX/YAML.
- [Consistency relations](02_consistency_relations.ipynb): calculate a relation at
  `n=2`, `p=2`, select the same methods and coefficient options, infer orders,
  simplify a display copy, evaluate numerically and format LaTeX/YAML.
- [Archived results](03_archived_results.ipynb): read and manipulate the saved YAML
  formulas for both physical models.

The [universal coefficient notebook](../results/universal_coefficients/universal_coefficients.ipynb)
and its instructions are kept with the coefficient archive.

Install SRI from the repository root:

```bash
python -m pip install -e .
```

Open the notebooks with this Python environment in Jupyter or a notebook-capable
editor. If kernel or execution tools are missing, see
[notebook setup](../docs/installation.md#notebooks).

The spectra and consistency notebooks perform low-order calculations; the archive
notebook reads saved formulas. See the [API guide](../docs/api.md).

## Reading the formulas

`utils.notebook_display` separates spectral corrections by relative weight and
relations by absolute weight, with short aligned rows and mathematical indices.
Only displayed floating-point coefficients are rounded to six significant
digits. Set `digits=None` on a display helper to keep all stored digits.
The computed SymPy expressions and exported text retain their original precision.
The `simplify` example preserves the exact expression and checks algebraic equality.
