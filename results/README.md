# Published results

- **[Canonical inflation](canonical_SRI/README.md):** scalar and tensor spectra
  from LO through N10LO, with tensor consistency relations.
- **[K-inflation](kinflation_SRI/README.md):** the same orders with variable
  scalar sound speed.
- **[Universal coefficients](universal_coefficients/README.md):** complete exact
  expressions and saved analytic and mode-equation values through weight ten.

Each spectrum/relation order has a **PDF**, **LaTeX source** and **YAML companion**.
Both models provide spectra for `p=0,...,10`, relations for `n=0,...,3` at
`p=0,...,6`, and `n=0` at `p=7,...,10`. The 64 published relations have
[exact forward-substitution checks](../validation/common/consistency_coverage.md).

Spectrum and relation archives contain the numerical coefficients shown in
the PDFs. Universal-coefficient YAML contains exact expressions. Use
`--symbolic` for new calculations with exact coefficients; see the
[analytic guide](../docs/symbolic.md).

## Read results in Python

From the repository checkout, with SRI installed:

```python
from utils import get

chi = get.coefficient(0, 1, 1)
power = get.spectra(p=3)
relation = get.consistency(n=0, p=3, type="kinflation")
```

`type="canonical"` is the default; `type="kinflation"` selects k-inflation.
`power.scalar` and `power.tensor` are the full SymPy spectra;
`relation.expression` is the selected consistency relation. These calls read
stored data without generating new coefficients or spectra. See the
[reader guide](../docs/api.md#archived-yaml-results) and
[archive notebook](../notebooks/03_archived_results.ipynb).

## Generate or export

The model collection guides describe equation generation and PDF compilation.
To create YAML companions from the existing LaTeX documents:

```bash
python -m scripts.archive_results
```

This reads the existing formulas without running a spectrum or coefficient
calculation. `scripts.export_results` also writes YAML when generating new documents.
