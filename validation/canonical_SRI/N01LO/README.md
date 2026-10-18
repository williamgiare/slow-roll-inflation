# N1LO: recovery of published results

This case tests the first correction to the scalar and tensor amplitudes, and the tensor tilt and first two runnings.

The comparison targets are the pivot amplitudes in
[Auclair & Ringeval (2022), Eqs. (44) and (54)](https://arxiv.org/html/2205.12608v2),
truncated to the order of this case.
The explicitly transcribed targets are in
[`published_spectra.py`](../../common/symbolic/published_spectra.py).
They are imported only by validation. This case regenerates the numerical mode
coefficients independently of those comparison targets.

## What is checked

- Every scalar and tensor monomial coefficient, rather than agreement at one parameter point.
- `n=0,1,2` against the separate observable-flow/Lagrange-inversion oracle.
- Exact forward substitution of the eliminated background, SR weights and the zero-`r` limit.
- Constant-`epsilon` spectra against the exact Hankel solution, including the expected remainder scaling.
- Scalar and tensor parameter-Taylor coefficients for a nonconstant canonical background, obtained from an independent linear mode solver.

The coefficient tolerances are asserted in the shared validation code and recorded
in the report. Published expressions are comparison outputs, never calculation inputs.

## Reproduce

Open [validation.ipynb](validation.ipynb), or run from the repository root:

```bash
python -m validation.canonical_SRI.order_cases --p 1
```

## Spectra

- [Generated equations (PDF)](../../../results/canonical_SRI/N01LO/N01LO.pdf)
- [LaTeX source](../../../results/canonical_SRI/N01LO/N01LO.tex)
- [YAML expressions](../../../results/canonical_SRI/N01LO/N01LO.yaml)
