# Direct N3LO literature comparison

The scalar amplitude for single-field k-inflation and the tensor amplitude in
Einstein gravity are compared coefficient by coefficient through relative order
three with [Bianchi & Gamonal, Phys. Rev. D 110, 104032 (2024)](https://arxiv.org/abs/2405.03157).
The reference is `PowerN3LO` at `k=ks` in the authors'
[official Mathematica supplement, v1.0](https://zenodo.org/records/13863577),
released under CC0. Its generic coefficients are stored in
[bianchi_gamonal_n3lo.yaml](bianchi_gamonal_n3lo.yaml).

Run from the repository root after installing SRI:

```bash
python -m scripts.validate_kinflation_literature --output /tmp/kinflation_n3lo.yaml
python -m pytest -q tests/test_kinflation_literature.py
```

The comparison uses the authors' own sector pivot, $-k\tau=1$, with
$d\tau=c\,dt/a$. For scalars this is SRI's native acoustic pivot; for tensors it
is the conformal pivot, with $c_t=1$. The quantities tested are the dimensionless
amplitude factors multiplying the leading normalization. The reference
normalization $\hbar H^2/(4\pi^2 c^3Z)$ gives
$\hbar H^2/(8\pi^2 M_{\rm Pl}^2\epsilon_1c_s)$ for
$Z=2M_{\rm Pl}^2\epsilon_1/c_s^2$. The overall tensor polarization convention
is factored out. This test does not exercise transport to the common tensor
pivot, the public sound-coordinate conversion or consistency-relation inversion.

The reference conversion uses only the definitions $D_N=d/d\ln a$,
$\epsilon_{1H}=\epsilon_1$, $\epsilon_{2H}=-\epsilon_2$,
$\epsilon_{3H}=-\epsilon_3$, $\epsilon_{1Z}=-\epsilon_2+2s_1$ and
$\epsilon_{1c}=-s_1$, where $s_j=D_N^j\ln c_s$.
Higher $Z$ and $c$ flow variables follow from the authors' recursion
$\epsilon_{n+1,X}=-D_N\ln\epsilon_{n,X}$. Apparent denominators cancel before
zero-flow limits are taken. Neither an effective canonical background nor
SRI's generated amplitudes enter the reference expression.

Both iterative and closed algebra give exactly zero symbolic residual for all
31 scalar and 8 tensor monomials, including LO. Numerical-backend coefficients
are checked separately with tolerance $10^{-9}$ on
$|a-b|/(1+|b|)$. The [saved report](report.yaml) records measured errors,
versions and source hashes.
