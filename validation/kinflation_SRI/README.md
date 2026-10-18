# Single-field k-inflation

These cases test a positive, slowly varying scalar sound speed in Einstein
gravity with `c_t=1`. All observables use a common tensor pivot; scalar
amplitudes include acoustic-to-tensor reference transport.

- Run all provided order cases:

  ```bash
  python -m validation.kinflation_SRI.run --all
  ```

- Run N3LO only:

  ```bash
  python -m validation.kinflation_SRI.run --p 3
  ```

The [LO](LO/README.md), [N01LO](N01LO/README.md), [N02LO](N02LO/README.md)
and [N03LO](N03LO/README.md) cases test generic polynomial expressions.
They compare the two algebra engines for `n=0,...,3` (through `n=8` at LO),
substitute each relation into the forward background observables and require
an exactly zero residual. They also test canonical recovery and a regular
zero-`r` limit. N01LO includes an independently derived analytic target.

A [direct external N3LO comparison](literature/README.md) checks the native
scalar and tensor amplitude factors against the official Bianchi–Gamonal
supplement, with exact symbolic equality for both algebra engines. The reference
YAML, pivot conventions and measured numerical residuals are included. Reproduce it with
`python -m scripts.validate_kinflation_literature`.

An independent effective-background construction uses
`a_eff=a sqrt(epsilon_1/c_s)` to map the acoustic scalar pump to a canonical
tensor pump. It checks the scalar source and native amplitude before pivot
transport. This comparison shares the universal moment data; their independent
mode-equation validation belongs to the shared coefficient tests.

[N10LO](N10LO/README.md) tests both spectra on constant-flow and evolving-flow
physical profiles. It compares the two engines, checks the effective-pump map,
and compares mixed scalar/tensor slow-roll Taylor coefficients with exact
Hankel amplitudes for variable `c_s`. These profile checks are distinct from a
full dense order-ten observable inversion. Reports state their precise scope.

Each notebook displays a current saved report by default; `RERUN=True` executes
the case again. Saved reports require matching calculation sources or an explicit
recorded compatibility check establishing an unchanged numerical calculation.
A matching, independently checked mode-generated coefficient bank may be reused
with `--bank validation/common/reports/order10_fresh_bank.yaml`.

Exact D1 forward-substitution checks cover all 32 published k-inflation
relations: `n=0,...,3` at `p=0,...,6`, and `n=0` at `p=7,...,10`.
The [coverage guide](../common/consistency_coverage.md) identifies the saved
evidence for each case and distinguishes it from the algebra-engine comparison.

The additional dense iterative tensor-tilt inverse has been checked by exact
forward substitution through `p=7`. The saved check is in
[generic_inverse.yaml](generic_inverse.yaml). Reproduce the full sweep with
`python -m validation.kinflation_SRI.run --generic-max-p 7`.
This sweep is optional because its symbolic cost grows rapidly with order.

The public spectra and consistency relations use the same sound scale
parameters `alpha_j^cs`. Independent forward derivatives check both directions
of the polynomial sound-coordinate conversion through relative order ten.
The converted spectra recover the native amplitudes at the same SR cutoff,
and the independent NLO scalar and tensor amplitudes are also checked.
