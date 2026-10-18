# Canonical spectra and consistency relations

Each order directory contains a standalone LaTeX document, its compiled PDF and a YAML companion.
The provided expressions use the default iterative algebra; the CLI also
supports `--closed`, as described in the [command guide](../../docs/cli.md).
The spectra are evaluated at `x=-k*eta=1`, through relative amplitude order `p`.
Scalar and tensor spectra occupy separate subsections. The consistency relations
are grouped by tensor derivative index and retained through `W=n+1+p`.
Long expressions are split into additive blocks by slow-roll order.

LO through N06LO include `n=0,1,2,3`. N07LO through N10LO include `n=0`.
These are the document export ranges, rather than limits on SRI calculations.
The coefficients are generated numerical values; no analytic identification of
special constants is implied. These equation documents are distinct from the
scientific checks and reports in `validation/`.

## Documents

- LO: [PDF](LO/LO.pdf) · [LaTeX](LO/LO.tex) · [YAML](LO/LO.yaml)
- N01LO: [PDF](N01LO/N01LO.pdf) · [LaTeX](N01LO/N01LO.tex) · [YAML](N01LO/N01LO.yaml)
- N02LO: [PDF](N02LO/N02LO.pdf) · [LaTeX](N02LO/N02LO.tex) · [YAML](N02LO/N02LO.yaml)
- N03LO: [PDF](N03LO/N03LO.pdf) · [LaTeX](N03LO/N03LO.tex) · [YAML](N03LO/N03LO.yaml)
- N04LO: [PDF](N04LO/N04LO.pdf) · [LaTeX](N04LO/N04LO.tex) · [YAML](N04LO/N04LO.yaml)
- N05LO: [PDF](N05LO/N05LO.pdf) · [LaTeX](N05LO/N05LO.tex) · [YAML](N05LO/N05LO.yaml)
- N06LO: [PDF](N06LO/N06LO.pdf) · [LaTeX](N06LO/N06LO.tex) · [YAML](N06LO/N06LO.yaml)
- N07LO: [PDF](N07LO/N07LO.pdf) · [LaTeX](N07LO/N07LO.tex) · [YAML](N07LO/N07LO.yaml)
- N08LO: [PDF](N08LO/N08LO.pdf) · [LaTeX](N08LO/N08LO.tex) · [YAML](N08LO/N08LO.yaml)
- N09LO: [PDF](N09LO/N09LO.pdf) · [LaTeX](N09LO/N09LO.tex) · [YAML](N09LO/N09LO.yaml)
- N10LO: [PDF](N10LO/N10LO.pdf) · [LaTeX](N10LO/N10LO.tex) · [YAML](N10LO/N10LO.yaml)

## Generate

```bash
python -m scripts.export_results --compile
```

PDF export uses an existing `pdflatex` installation. Compiler auxiliaries are
created in a temporary directory; each order directory retains `.tex`, `.pdf` and `.yaml`.
Scientific checks and their reports are documented separately in `validation/`.

A previously generated production moment bank can be reused after its kernel
fingerprint is checked:

```bash
python -m scripts.export_results --compile --bank validation/common/reports/order10_fresh_bank.yaml
```

`--min-p`, `--max-p`, `--max-n`, `--high-order-from` and `--high-order-max-n` control the
export range. The defaults reproduce the documents provided here. Omitting
`--compile` generates the LaTeX sources and YAML companions. Increasing high-order derivative
indices can substantially increase calculation time.

The stored documents use numerical coefficients. For a new calculation with
exact analytic coefficients, use `SRI -p 2 --symbolic`.

## Read the archived expressions

Run from the repository checkout with SRI installed:

```python
from utils import get

power = get.spectra(p=3, type="canonical")
relation = get.consistency(n=0, p=3, type="canonical")
```

`power.scalar` and `power.tensor` are the complete SymPy spectra;
`relation.expression` is the tensor consistency relation. The YAML retains the
published decimal coefficients. See the [reader guide](../../docs/api.md#archived-yaml-results)
and [notebook](../../notebooks/03_archived_results.ipynb).
