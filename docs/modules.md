# Module guide

The public entry points are `sri.spectra`, `sri.consistency` and the `SRI` command.
The public interface and shared coefficient selection interface reside in `source/`.
The two algebra engines reside in `source/iterative/` and `source/closed/`.

| Module | Responsibility |
| --- | --- |
| `orders.py` | Resolve compatible nonnegative `n`, `p` and `W` |
| `series.py` | Sparse rational algebra, slow-roll weights and truncation |
| `background.py` | Shared weighted coordinates and physical constants |
| `generic_modes.py` | Normalized mode equations, source basis and numerical convergence |
| `mode_coefficients.py` | Select numerical or analytic universal moments and count insertion multiplicities |
| `analytic/` | Exact Hankel/Mellin and mixed GPL generators, defined Chi constants and local caches |
| `spectra.py`, `consistency.py` | Public calculation functions and result caches |
| `kinflation.py` | Acoustic scalar modes, sound-speed derivatives and k-inflation relations |
| `iterative/` | Intermediate differentiation, amplitude assembly and difference-form inverse |
| `closed/` | Explicit Bell derivatives, partition sums and Lagrange–Good inverse |
| `results.py` | Return objects, numerical evaluation and LaTeX export |
| `cli.py` | Command arguments, combined result selection and output formats |
| `terminal.py` | Aligned equations, term counts and total computation time |
| `latex_output.py` | Standalone spectra and consistency documents grouped by SR order |
| `presentation.py` | Optional fraction recognition and factor collection |
| `serialization.py` | Safe YAML encoding for CLI output, fixtures and reports |

## Spectra

At relative order `p`, the calculation constructs the conformal factor
$h=-\eta aH$ and both pump sources. Source derivatives determine their
logarithmic coefficients. The normalized mode moments multiply finite products
of these coefficients, including the appropriate permutation multiplicities.
The reference normalization contributes $-2\ln h$ to the log amplitude.
Truncated exponentiation gives the scalar and tensor amplitude corrections.

The numerical mode system propagates the logarithmic ratio of the mode and
growing pump on the outgoing contour. Its source basis consists of integer
partitions at the requested weight. Refined integrations test both endpoints,
initial data and integration tolerance. If double precision fails the
convergence criterion, an arbitrary-precision integration is used. Failure
of the refined calculation raises `CoefficientConvergenceError`.

Generated moments enter the weighted algebra as exact rational representations
of their numerical values. Coefficients are rounded for public display.
The formal slow-roll cutoff and numerical convergence tolerance are separate.

## Consistency relations

The requested pair `(n,p)` fixes `W=n+1+p`. Background spectra supply the
scalar observables and the selected tensor derivative. The leading inverse
starts with $\epsilon_1=r/16$ and $\epsilon_2=-r/8-\alpha_0^s$.
The correction to this map raises slow-roll weight, so its finite inverse
terminates at order `W`. Applying that inverse produces a polynomial in
`r` and the required scalar derivatives.

The iterative method builds and differentiates intermediate background
expressions with consistent truncation; the closed method evaluates explicit
finite sums. It does not require lower-order tensor
consistency relations. At fixed relative order `p`, increasing `n` requires
additional derivatives but no higher-order mode moments.

`presentation.py` operates on copies of the output expressions and does not
alter generated moments or background elimination.

## Validation tools

`validation/run.py` runs named shared suites. `validation/canonical_SRI/order_cases.py`
provides the LO, low-order and order-ten cases; `validation/canonical_SRI/high_order.py`
performs the full sweep through a specified order. `scripts/export_results.py`
writes standalone equation documents under `results/`; PDF compilation keeps
auxiliary files in a temporary directory.

Independent symbolic comparisons, mode solvers, fixtures and reports reside in
`validation/common/`. The per-order notebooks call these shared tools and
keep validation results and figures alongside each case. Equation documents
are stored separately under `results/`. Production modules
do not import validation code or reference fixtures.

## Algebra engines

`iterative` is the default. Its conformal factor, background jets, scalar and
selected tensor derivatives are generated successively. Its observable inverse
uses differences that raise slow-roll weight and stop at the requested cutoff.

`closed` evaluates the conformal factor with complete Bell polynomials.
Ordinary background derivatives use the explicit Taylor jets
$E(N+t)=E\exp[\sum_jq_jt^j/j!]$ and
$q_i(N+t)=\sum_jq_{i+j}t^j/j!$. Scale derivatives use finite partial-Bell sums.
The leading observable reconstruction uses the LO Bell solution, and the full
inverse uses the multivariate Lagrange–Good derivative sum. Differentiated
products are kept through `W+sum(m_i*w_i)` before the final truncation to `W`.
The source and amplitude sums contain the same selected universal coefficients
in both engines. The closed engine does not call iterative background derivatives
or the difference-form inverse. Algebraic equivalence is tested exactly before
coefficients are rounded.

`closed` selects the algebraic formulas independently of the coefficient backend.
`symbolic=True` selects the exact analytic provider for either engine. There is no fixed order cap; closed sums
can become expensive as the orders grow.

## Acoustic scalar modes

`kinflation.py` extends the background coordinates with
`s_j=D_N^j ln(c_s)` and `ell=ln(c_s)`. Although `ell` has weight zero,
its derivative is `s_1`; it is not an integration constant.
The scalar clock satisfies `D_N h_s=(1-epsilon_1-s_1)h_s-1`.
The acoustic pump is `a Mpl sqrt(2 epsilon_1/c_s)`. Tensor modes use the
canonical clock and source. Both sectors use the same universal mode moments.

Native scalar amplitudes are transported through a finite Taylor sum from
`k_s=-1/tau` to `k_t=-1/eta` at the same background time. The separation is
`ln(c_s)+ln(h_s)-ln(h)`, with `d tau=c_s d eta`.
The enlarged inverse uses the polynomial coordinate `R=r/c_s` together with
scalar and sound-speed scale derivatives. The output substitutes `R=r/c_s`
and `ell=ln(c_s)`, leaving only observable quantities.

`Operations` selects the successive or finite-sum derivative engine; neither
requires published amplitude coefficients. Shared closed derivatives support
both background jets and the dynamically evolving logarithm of the sound speed.
`validation/kinflation_SRI/run.py` supplies exact algebra, effective-pump and
Hankel checks independent of the production assembly.

Before spectrum export, `sound_scale_replacements` converts the internal
`s_j` into `alpha_(j-1)^cs`. The forward conversion is a triangular linear
map whose diagonal is `h^j`. Reciprocal clock factors are expanded and
truncated, so the exported amplitude factors remain polynomials in the
Hubble-flow and sound scale parameters. The iterative method solves the
triangular system; the closed method sums strictly descending matrix paths.
Both `SpectraResult.evaluate` and `ConsistencyResult.evaluate` accept
`sound=[alpha_0^cs,alpha_1^cs,...]`. The native mode assembly and the
background inverse retain their internal coordinates.

`latex_output.py` groups amplitude terms by their slow-roll weights, treating
`c_s` and `ln(c_s)` as coefficients and `alpha_j^cs` as weight `j+1`.
`scripts/export_results.py` writes the canonical and k-inflation collections under
`results/canonical_SRI/` and `results/kinflation_SRI/`. The `--kinflation` option
selects the latter; `--min-p` and `--max-p` select an order range.

## Polynomial algebra

`series.py` evaluates grouped substitutions, direct multinomial differences,
monomial products, symmetric polynomial squares and weight-aware Horner sums
for exponentials and logarithms. A bounded cache stores monomial weights without
retaining ring instances. All coefficients remain rational, and every operation
uses the same SR cutoff as the direct arithmetic.

The iterative inverse uses `WeightedRing.substitution_difference` to evaluate
`P(x+delta)-P(x)` without constructing the identity contribution. The closed
inverse retains the Lagrange–Good sum. Both physical models and the document
exporter use this shared arithmetic without replacing methods globally.
Every exported relation is checked by exact forward substitution.

`validation/common/performance.py` compares complete observable maps and tensor
relations with the independent reference arithmetic in
`validation/common/symbolic/reference_arithmetic.py`. It also checks all exact and
numerical inputs of the universal mode solver through order ten.

## Analytic coefficient backend

`analytic/fast.py` provides exact constant-source Hankel, one-varying and
two-insertion Mellin families. `analytic/generic_hyperlogs.py` generates the
remaining mixed families through ordered Mellin master systems, exact pole
cancellation and rational directional extraction. `analytic/mixed_hyperlogs.py`
provides Gamma-jet and related exact Mellin utilities.

`SymbolicMomentData` keeps universal-coefficient symbols in the rational
background algebra, resolving short formulas explicitly and mixed constants
as `Chi(pattern)`. Their definitions and arbitrary-precision evaluation belong
to `analytic/`, independently of physical model and algebra engine. Production
never imports the coefficient validation suite or its mode-equation fixture.

## Published archive readers

`utils/get.py` provides `coefficient`, `value`, `spectra` and `consistency` over
the YAML collections. `utils/_expressions.py` decodes the exact coefficient
syntax; result formulas use a restricted arithmetic AST with declared symbols.
Readers do not call the coefficient generators or physical solvers.
`notebooks/03_archived_results.ipynb` demonstrates both physical models.
The coefficient notebook is stored with its archive under
`results/universal_coefficients/universal_coefficients.ipynb`.

`scripts/archive_results.py` exports YAML companions from the finite polynomial
LaTeX syntax used in the published equation documents. Its converter preserves
their displayed numerical precision and rejects unsupported syntax.
`scripts/export_results.py` writes a YAML companion alongside each new LaTeX
document. PDF compilation uses `--compile`.

Reusable export-algebra helpers are in `utils/export_algebra.py`; commands
for generating documents and checking reports are in `scripts/`.

## Notebook display

`utils/notebook_display.py` provides `show_spectra`, `show_consistency` and
`show_equation` for the repository notebooks. Spectra are shown with their
normalizations and relative-weight contributions; relations are grouped by
absolute weight. Short aligned rows and mathematical symbol names improve
readability. Display copies use six significant digits for floating-point
coefficients by default; exact constants and the underlying expressions are
preserved. Pass `digits=None` to retain all stored numerical digits.

The [spectra notebook](../notebooks/01_spectra.ipynb) and
[consistency notebook](../notebooks/02_consistency_relations.ipynb) demonstrate
the calculation API separately. The relation notebook includes exact factor
collection with `simplify_expression`.
