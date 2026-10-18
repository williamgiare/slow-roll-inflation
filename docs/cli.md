# Command-line interface

Run `SRI` from the environment in which the package is installed. The defaults
are canonical inflation, numerical universal coefficients and iterative algebra.
Physical parameters remain symbolic in every calculation mode.

## Orders

- **Spectra:** `-p P` returns both scalar and tensor spectra through relative
  order `P`: `0` is LO, `1` is NLO, `2` is N2LO, and so on.
- **Consistency relations:** supply any two of `-n N`, `-p P` and `-W W`.
  The third follows from `W=n+1+p`. The output includes both spectra and the
  selected tensor relation. `n=0` selects the tilt, `n=1` its running and
  `n=2` its running of running.
- **Remainders:** spectra have relative remainder `O(SR^(p+1))`; relations
  have absolute remainder `O(SR^(W+1))`.

All indices must be nonnegative integers. If all three are given, they must
satisfy `W=n+1+p`. Argument order is unrestricted.

- N2LO scalar and tensor spectra:

  ```bash
  SRI -p 2
  ```

- Tensor running of running with N2LO corrections; these commands are equivalent:

  ```bash
  SRI -n 2 -p 2
  SRI -n 2 -W 5
  SRI -p 2 -W 5
  ```

The implementation accepts arbitrary finite orders. Calculation time and memory
increase with order; increasing `n` at fixed `p` requires no higher-weight
universal coefficients.

## Calculation method

- `--symbolic` selects exact analytic coefficients; numerical coefficients are the default.
- `--closed` selects closed background algebra; iterative algebra is the default.
- `--kinflation` selects a slowly varying scalar sound speed in Einstein gravity,
  with tensor speed `c_t=1`; canonical inflation is the default.

These options can be combined for spectra and consistency relations:

```bash
SRI -p 2 --symbolic --closed
SRI -n 0 -p 2 --kinflation --symbolic --closed
```

Both models use the common tensor reference `x=-k*eta=1`. In k-inflation, the
scalar spectrum includes transport from its acoustic reference; `c_s` and
`alpha_j^cs` refer to the common tensor reference. Symbol definitions are in
[Output notation](../README.md#output-notation).

## Export

Use `--format text` (default), `--format latex` or `--format yaml`.
Add `--output PATH` to save the result; otherwise it is printed to the terminal.
The output directory must already exist.

```bash
SRI -p 2 --format latex --output spectra.tex
SRI -n 2 -p 2 --format yaml --output relation.yaml
```

- **Text:** spectra followed by the optional relation, with aligned continuation
  lines. The final report lists term counts, remainders and calculation time.
- **LaTeX:** equation fragments containing both spectra and the optional relation.
  For standalone documents and PDFs, see the [result guide](../results/README.md).
- **YAML:** expression strings, orders, method, model and coefficient metadata.
  A relation result also contains both spectra in a `spectra` object.

CLI YAML describes one calculation. The per-order files in `results/` use the
archive format read by `utils.get`; see [archive access](api.md#archived-yaml-results).
Successful commands exit with status 0; invalid arguments or failed convergence
exit with status 2.

## Optional simplification

`--simplify` collects common factors in both spectra and relations, in every
output format. With numerical coefficients, it can also recognize simple
fractions. With `--symbolic`, it preserves exact coefficients and only collects
factors.

```bash
SRI -n 2 -p 2 --symbolic --simplify --format yaml --output relation.yaml
```

Numerical fraction recognition requires numerator and denominator at most 4096
and relative agreement within `1e-12`. It does not certify an analytic identity
or replace a nonzero coefficient by zero. YAML records the display operations
under `presentation`; the computed expressions and convergence settings are preserved.

## Exact analytic coefficients

Short analytic coefficients are printed explicitly. Larger mixed coefficients
appear as `Chi((j_1,...,j_m))` in text/YAML and as $\chi_{m;j_1,\ldots,j_m}$
in LaTeX. These are defined constants and require no additional physical inputs.
See the [analytic coefficient guide](symbolic.md) to expand or evaluate them.

Compact symbolic output avoids expanding large special-function expressions.
Expanding or numerically evaluating a large mixed family can take hours.
All 138 coefficients through weight ten have passed independent comparison;
this is the validated coefficient range, not an implementation order limit.

## Numerical coefficient accuracy

`--coefficient-tolerance T` sets the numerical moment convergence tolerance
(default `1e-10`, accepted range `1e-13` to `1e-3`). It is independent of the
slow-roll cutoff and is ignored with `--symbolic`.

Refined mode integrations check convergence. The engine uses 40-digit integration
when double precision fails the check, and returns no result if convergence
still fails. The tolerance applies to individual moments; it is not an error
bound for every combined coefficient in a final spectrum or relation.
