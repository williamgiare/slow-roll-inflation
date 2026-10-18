# N03LO: variable scalar sound speed

Cubic scalar and tensor amplitudes and generic relations for `n=0,...,3` pass exact forward-substitution, canonical-limit and iterative-versus-closed checks. The effective-background map checks the acoustic source separately.

```bash
python -m validation.kinflation_SRI.run --p 3
```

[validation.ipynb](validation.ipynb) displays the saved report and can rerun
the case. [report.yaml](report.yaml) records residuals, settings and measured
errors. The [model guide](../README.md) explains the physical assumptions and
the [shared tools](../../common/README.md) check the universal moment generator.
