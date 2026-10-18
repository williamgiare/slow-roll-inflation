# Installation

Requires Python 3.10 or newer. From the repository root:

```bash
pip install -e .
SRI --version
SRI -p 0
```

The last command calculates LO scalar and tensor spectra. The installation
includes both physical models, both algebra methods and the numerical and
analytic coefficient backends. The package is imported as `sri`;
`python -m sri` is equivalent to `SRI`.

For an isolated environment, create and activate it before installing:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e .
```

## Tests and scientific validation

```bash
python -m pip install pytest matplotlib
python -m pytest -q tests
```

Run commands from the repository root. The [validation guide](validation.md)
separates package tests, saved reports and independent scientific reproduction.
Full coefficient reproduction through weight ten can take several hours.

## Notebooks

Open the [example notebooks](../notebooks/README.md) in Jupyter or a
notebook-capable editor using the environment in which SRI is installed.
If notebook execution tools are missing:

```bash
python -m pip install ipykernel nbformat nbclient matplotlib
```

These install kernel and execution support; use your existing notebook editor.
The [coefficient notebook](../results/universal_coefficients/README.md#coefficient-notebook)
is stored with its archive. Repository utilities such as `utils.get` read
`results/` from the checkout.

## Equation documents

The [published PDFs](../results/README.md) need no TeX installation.
To generate new LaTeX/YAML documents and compile their PDFs:

```bash
python -m scripts.export_results --compile
python -m scripts.export_results --kinflation --compile
```

`--compile` requires `pdflatex`; omit it to generate only LaTeX and YAML.
The full order sweep can be expensive. The collection guides explain how to
select orders or reuse a checked coefficient bank.

## Symbolic coefficients

Exact coefficients are included in the standard installation:

```bash
SRI -p 2 --symbolic
```

Compact results use defined `Chi` constants for large mixed expressions.
Their complete expressions are generated when expanded or evaluated, and
cached under `~/.cache/sri/symbolic`; `SRI_SYMBOLIC_CACHE` selects another
writable location. See the [analytic coefficient guide](symbolic.md).
Use `utils.get.coefficient` to read an existing exact expression from the
[coefficient archive](../results/universal_coefficients/README.md).

## Release version

The release version is defined only in the root `VERSION` file. After changing it,
rerun `pip install -e .` to update package metadata, `sri.__version__` and
`SRI --version` together.
