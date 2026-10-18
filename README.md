# Slow-Roll Inflation at Arbitrary Order

This is a symbolic algorithm that, for both canonical single-field inflation and $k$-inflation with a slowly varying scalar sound speed, computes the following quantities at arbitrary order in the slow-roll expansion:

- **Primordial spectra:** the scalar and tensor primordial power spectra, $\mathcal P_s(k)$ and $\mathcal P_t(k)$, through any specified finite relative order $p$ in the slow-roll expansion.
- **Consistency relations:** the tensor tilt $n_t=\alpha_0^t\equiv d\ln\mathcal P_t/d\ln k$ and its runnings $\alpha_n^t\equiv d^n n_t/d(\ln k)^n$ ($n\geq1$), expressed in terms of the tensor-to-scalar ratio $r=\mathcal P_t/\mathcal P_s$, the scalar tilt $n_s-1=\alpha_0^s\equiv d\ln\mathcal P_s/d\ln k$ and its runnings $\alpha_n^s\equiv d^n(n_s-1)/d(\ln k)^n$ ($n\geq1$). Relations can be calculated for any derivative order $n\geq0$ and slow-roll correction order $p$. For $k$-inflation, they also depend on $c_s$, its logarithmic scale derivative $\alpha_0^{c_s}\equiv d\ln c_s/d\ln k$ and higher derivatives $\alpha_n^{c_s}\equiv d^n\alpha_0^{c_s}/d(\ln k)^n$ ($n\geq1$).

Each relation is derived directly, without requiring lower-order tensor consistency relations as inputs. Calculation time and memory increase with order.

## Installation

Requires Python 3.10 or newer. From the repository root:

```bash
pip install -e .
```

Check the installation and try a leading-order (LO) spectrum calculation:

```bash
SRI --version
SRI -p 0
```

See the [installation guide](docs/installation.md) for optional notebook and validation dependencies.

## Calculation modes

The universal coefficients and the background algebra can be calculated independently in the following ways:

- **Universal coefficients: numerical (default).** The coefficients are determined by integrating the normalized mode and pump equations. Their converged numerical values are represented by rational numbers for the subsequent symbolic algebra.
- **Universal coefficients: symbolic (`--symbolic`).** The same coefficients are obtained from exact analytic expressions. Larger mixed coefficients can appear as defined constants `Chi((j_1,...,j_m))`; their complete expressions are available through `sri.analytic.moment` and the [coefficient archive](results/universal_coefficients/README.md).
- **Background algebra: iterative (default).** Intermediate background quantities are constructed successively, and the background is eliminated from consistency relations using the finite difference inverse.
- **Background algebra: closed (`--closed`).** Background quantities are constructed through explicit finite sums, and consistency relations use the closed perturbative background inverse.

With no mode flags, calculations use numerical coefficients and iterative background algebra. `--symbolic` changes the coefficient mode; `--closed` changes the algebra method. Both algebra methods support spectra and consistency relations for both physical models, with either coefficient mode. Physical parameters remain symbolic in both coefficient modes. See the [analytic coefficient guide](docs/symbolic.md) for exact-output conventions.

## Command Line Usage

| Flag | Meaning | Default |
|---|---|---|
| `-p P` | Relative slow-roll order of the spectra and corrections beyond the leading consistency order: `0` for LO, `1` for NLO, `2` for N2LO, and so on. | Required for spectra; may be inferred for a relation. |
| `-n N` | Tensor derivative index: `0` for the tilt, `1` for its running, `2` for its running of running. | Inferred when `p` and `W` are supplied. |
| `-W W` | Absolute consistency cutoff, $W=n+1+p$; supply any two of `n`, `p`, `W` for a relation. | Inferred when `n` and `p` are supplied. |
| `--kinflation` | Calculate $k$-inflation with a slowly varying scalar sound speed and tensor speed $c_t=1$. | Canonical inflation. |
| `--symbolic` | Use exact analytic universal coefficients. | Numerical coefficients, without a flag. |
| `--closed` | Select closed background algebra and elimination. | Iterative algebra, without a flag. |
| `--coefficient-tolerance T` | Convergence tolerance for numerical universal coefficients; ignored with `--symbolic`. | `1e-10` |
| `--simplify` | Collect display factors; numerical output can also recognize simple fractions. | No optional simplification. |
| `--format FORMAT` | Output as `text`, `latex` or `yaml`. | `text` |
| `--output PATH` | Save the output to a file; its parent directory must exist. | Print to the terminal. |
| `--version` | Show the installed version and exit. | — |
| `-h`, `--help` | Show the available command-line options and exit. | — |

## Command Line Examples

- **Canonical spectra at N5LO:**

  ```bash
  SRI -p 5
  ```

- **The same spectra for k-inflation:**

  ```bash
  SRI -p 5 --kinflation
  ```

- **Tensor running of running with N2LO corrections** ($n=2$, $p=2$, $W=5$), including both spectra:

  ```bash
  SRI -n 2 -p 2
  ```

- **The same relation with exact coefficients, closed algebra, simplification and YAML export:**

  ```bash
  SRI -n 2 -p 2 --symbolic --closed --simplify --format yaml --output relation.yaml
  ```

See the [CLI guide](docs/cli.md) for further details.

## Python Interface

The Python interface returns SymPy expressions:

```python
from sri import spectra, consistency

power = spectra(p=2, symbolic=True)
relation = consistency(n=2, p=2, method="closed", symbolic=True)

print(power.scalar)
print(power.tensor)
print(relation.expression)
```

Use `kinflation=True` for a new $k$-inflation calculation, for example `spectra(p=2, kinflation=True, symbolic=True)`. Use `symbolic=False` (default) for numerical universal coefficients and `method="iterative"` (default) for iterative algebra. See the [Python API guide](docs/api.md) for numerical evaluation and additional methods.

## Notebook Examples

- **[01 — Spectra](notebooks/01_spectra.ipynb):** calculate scalar and tensor spectra at low order, select the algebra and coefficient modes, and evaluate the expressions numerically.
- **[02 — Consistency relations](notebooks/02_consistency_relations.ipynb):** calculate a relation at $n=2$, $p=2$, infer orders, and use `simplify`.
- **[03 — Archived results](notebooks/03_archived_results.ipynb):** read saved spectra, consistency relations and coefficients from YAML without rerunning their calculations.
- **[Universal coefficients](results/universal_coefficients/universal_coefficients.ipynb):** inspect exact analytic coefficients and their saved numerical values.

Open the notebooks in Jupyter or a notebook-capable editor using the environment in which SRI is installed. See the [notebook guide](notebooks/README.md) for setup and formula-display options.

## Output Notation

| Printed symbol | Definition |
|---|---|
| `P_s` | Scalar power spectrum. |
| `P_t` | Tensor power spectrum, both polarizations. |
| `H` | Hubble rate. |
| `Mpl` | Reduced Planck mass. |
| `epsilon_1` | $-d\ln H/dN$, with $N=\ln a$. |
| `epsilon_j` | $d\ln\lvert\epsilon_{j-1}\rvert/dN$, $j\geq2$. |
| `r` | $\mathcal P_t/\mathcal P_s$. |
| `alpha_0^s` | Scalar tilt $n_s-1$. |
| `alpha_n^s` | $d^n(n_s-1)/d(\ln k)^n$, $n\geq1$. |
| `alpha_0^t` | Tensor tilt $n_t$. |
| `alpha_n^t` | $d^n n_t/d(\ln k)^n$, $n\geq1$. |
| `c_s` | Scalar sound speed. |
| `c_t` | Tensor propagation speed; $c_t=1$. |
| `alpha_0^cs` | $d\ln c_s/d\ln k$. |
| `alpha_n^cs` | $d^n\alpha_0^{c_s}/d(\ln k)^n$, $n\geq1$. |
| `C` | $\gamma_E+\ln 2-2$. |
| `EulerGamma` | Euler–Mascheroni constant $\gamma_E$. |
| `Chi((j_1,...,j_m))` | Exact universal coefficient $\chi_{m;j_1,\ldots,j_m}$. |
| `pi` | $\pi$. |
| `log` | Natural logarithm. |
| `zeta(q)` | Riemann zeta function $\zeta(q)$. |
| `n` | Tensor derivative index. |
| `p` | Relative slow-roll correction order. |
| `W` | Absolute cutoff $n+1+p$. |
| `O(SR^q)` | Remainder starting at slow-roll weight $q$. |

The common reference is `x=-k*eta=1`, where `k` is the comoving wavenumber and `eta` is conformal time. For $k$-inflation, `c_s` and its derivatives refer to this same tensor reference; `tau` in the scalar-reference metadata is acoustic conformal time, $d\tau=c_s\,d\eta$.

Text and YAML formulas use `*` for multiplication and `**` for powers. The suffixes `^s`, `^t` and `^cs` label the sector in derivative names.

## Results

- **[Canonical inflation](results/canonical_SRI/README.md):** scalar and tensor spectra and consistency relations through N10LO.
- **[k-inflation](results/kinflation_SRI/README.md):** the corresponding results through N10LO with variable scalar sound speed.
- **[Universal coefficients](results/universal_coefficients/README.md):** all 138 coefficients through slow-roll weight ten, with exact expressions and saved analytic and independent numerical values.

Spectra and consistency relations are available as **PDF, LaTeX and YAML**; universal coefficients are stored in **YAML**. Each collection's README lists the available orders and explains how to read and reuse the results.

## Validation

- **Universal coefficients:** independent analytic and mode-equation determinations agree for all 138 coefficients through weight ten, including 67 mixed patterns. See the [coefficient comparisons](validation/symbolic_coefficients/README.md).
- **Known spectra:** comparisons through N3LO include canonical spectra and a direct $k$-inflation comparison with the Bianchi–Gamonal analytic results. See the [canonical checks](validation/canonical_SRI/README.md) and [direct literature comparison](validation/kinflation_SRI/literature/README.md).
- **Physical mode equations:** exact Hankel limits and independent integrations on nonconstant backgrounds check spectral amplitudes, normalization and the sound-speed reference-scale translation. See the [validation guide](docs/validation.md).
- **Iterative and closed methods:** complete comparisons cover $p=0,\ldots,4$ and $n=0,\ldots,3$ for canonical inflation, and $p=0,\ldots,3$ with the same derivative indices for $k$-inflation.
- **Background elimination:** exact forward-substitution checks cover all 64 published consistency relations, 32 per model. The [coverage guide](validation/common/consistency_coverage.md) lists their orders and saved results.
- **LO derivative hierarchy:** symbolic comparisons through $n=8$ and exact rational and finite-field checks through $n=1000$ test high derivative orders separately from high slow-roll orders. See the [large-$n$ checks](validation/canonical_SRI/LO/large_n_validation/README.md).

Run package tests (`pytest` required):

```bash
python -m pytest -q tests
```

Regenerate and compare all coefficients through weight 10:

```bash
python validation/symbolic_coefficients/run.py --orders 8 9 10 --family-seconds 0 --order-seconds 0 --fresh
```
> [!WARNING]
> Full coefficient reproduction through weight 10 can take several hours.

See the [validation guide](docs/validation.md) for details and saved reports.

## Project Structure

- **`source/`** — Implementation of the `sri` package: CLI, background variables, spectra, consistency relations and polynomial algebra.

- **`utils/`** — Reusable helpers for archived formulas, notebook display and export algebra.

- **`scripts/`** — Commands to export equation documents, create YAML archives and inspect scientific validation coverage.

- **`notebooks/`** — The three numbered usage notebooks for spectra, consistency relations and archived results.

- **`results/`** — Published spectra, consistency relations and universal coefficients.

- **`validation/`** — Scientific checks, independent references, reproducible notebooks and saved reports.

- **`tests/`** — Automated tests of the public interfaces, algebra, physical calculations, output formats and archive readers.

- **`docs/`** — Detailed installation, CLI, Python API, module, analytic-coefficient and validation guides.

## Citation

If you use this code or results in your research, please cite:

> W. Giarè, *Slow-Roll Inflation at Arbitrary Order*, arXiv:2610.xxxxx (placeholder).

## Author & Contributors

- [William Giarè](https://williamgiare.com)
