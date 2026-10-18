# LO: Bell formula and Cortés–Liddle hierarchy

At `p=0`, SRI recovers the closed leading-order relation for tensor scale
derivatives. This case compares that finite Bell expression with the recursive
consistency flow of Cortés and Liddle. Writing
$T_{n+1}=\alpha_n^t$, $S_j=\alpha_{j-1}^s$ and $T_1=-r/8$, the reference flow is

$$T_1'=T_1(T_1-S_1),\qquad S_j'=S_{j+1}.$$

The [notebook](validation.ipynb) compares SRI, the finite Bell formula,
a projective generating-series evaluation and the recursive flow symbolically
for `n=0,...,8`. Every residual vanishes exactly.
The [large-n tests](large_n_validation/README.md) compare exact rational and
finite-field evaluations through `n=1000`.

```bash
python -m validation.canonical_SRI.order_cases --p 0
python -m validation.canonical_SRI.LO.large_n_validation.run --max-n 1000
python -m pytest tests/test_lo_large_n_validation.py tests/test_order_cases.py
```

The large-n execution is recorded in [RESULTS.md](large_n_validation/RESULTS.md).
Leading scalar and tensor normalizations are given in
[the generated expressions](../../../results/canonical_SRI/LO/LO.pdf) ([LaTeX source](../../../results/canonical_SRI/LO/LO.tex)).
