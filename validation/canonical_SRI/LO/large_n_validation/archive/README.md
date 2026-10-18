<!-- WG Edit: scientific and software documentation. -->
# Independent projective and recursive implementations

These scripts provide separate LO hierarchy evaluations for the large-n tests.

- `Symbolic_test.py`: symbolic Bell and recursive consistency expressions.
- `High_order_test.py`: exact rational Bell-polynomial and projective-jet checks.
- `Jet_test.py`: projective-jet and recursive-flow coefficient comparison.
- `projective_consistency_hierarchy_tests.py`: rational and finite-field
  generating-series, Bell and Riccati evaluations.
- `projective_consistency_hierarchy 12.29.19.py`: symbolic hierarchy implementation.

The notebooks contain symbolic and high-order rational checks. The common
entry point is `python -m validation.canonical_SRI.LO.large_n_validation.run`, run from the
repository root with SRI installed.
