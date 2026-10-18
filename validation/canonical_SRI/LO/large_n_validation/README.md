# LO hierarchy at high derivative order

These tests compare the closed Bell formula, projective generating series and
recursive consistency flow for arbitrary scalar derivative inputs. With
$T_{n+1}=\alpha_n^t$, $S_j=\alpha_{j-1}^s$ and $T_1=-r/8$,

$$T_1'=T_1(T_1-S_1),\qquad S_j'=S_{j+1}.$$

The closed construction follows from

$$K(s)=\exp\!\left[-\sum_{j\geq1}S_j\frac{s^j}{j!}\right],\qquad
T_1(s)=\frac{T_1(0)K(s)}{1-T_1(0)\int_0^sK(u)\,du}.$$

The generating series evaluates its coefficients without expanding large
symbolic polynomials. A separate reference differentiates the Riccati flow.

## Tested ranges

- Symbolic SRI, Bell, projective and recursive expressions: `n=0,...,8`.
- Exact rational generating series and flow: `n=0,...,200`, with direct partial
  Bell sums through `n=25`.
- Finite-field generating series and flow: `n=0,...,1000`, over three primes
  and three scalar input sequences per field.
- Exact rational projective-jet and flow: `n=0,...,1000`.

The finite-field checks are exact at their tested inputs. The symbolic comparison
with SRI establishes polynomial equality through `n=8`. These tests concern
large derivative order `n` at `p=0`, independently of the N10LO tests.

## Execution

```bash
python -m validation.canonical_SRI.LO.large_n_validation.run --max-n 1000
python -m pytest tests/test_lo_large_n_validation.py
```

`run.py` uses the independent implementations in [archive/](archive/README.md)
and writes [RESULTS.md](RESULTS.md). The regular test suite uses shorter stress
ranges while retaining the symbolic comparison through `n=8`.
