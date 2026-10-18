# Iterative and closed consistency constructions

The two engines implement the same finite-order slow-roll construction.
The iterative engine differentiates intermediate expressions and applies
the finite inverse in difference form. The closed engine evaluates explicit
Bell/partition sums for the background and scale derivatives and the
multivariate Lagrange–Good derivative sum for background removal.

The comparison checks the conformal factor, both pump sources, log amplitudes,
primordial spectra, the observable map, tensor derivatives, and the final
consistency relations. Forward substitution tests the closed relation as well.
Residuals are compared as exact rational polynomials before output rounding.

Both engines use the same numerically generated universal mode moments.
Equality therefore tests the algebraic construction; independent tests of
the moment generator remain in the per-order and shared validation suites.
No literature formula is used in this comparison.

```bash
python -m validation.canonical_SRI.iterative_vs_closed.run
```

The default run compares scalar and tensor spectra through `p=4`, consistency
relations for `n=0,...,3` through `p=4`, and LO relations through `n=8`.
These are tested ranges, not implementation limits. Increase `--max-p`,
`--max-n`, or `--max-lo-n` to extend them. Closed evaluation can be substantially
more expensive at high orders.

[validation.ipynb](validation.ipynb) runs the comparisons and displays their
zero residuals. Run settings and source fingerprints are recorded in
`results.yaml`. No numerical tolerance is used for the polynomial comparisons.
