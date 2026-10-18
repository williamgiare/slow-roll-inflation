# N10LO: variable scalar sound speed

Relative amplitude order ten is tested on constant-flow and evolving-flow physical profiles. Exact iterative-versus-closed polynomial comparisons and mixed-flow Hankel Taylor coefficients test the spectra; the independent effective-background map checks the acoustic source. Generic exact observable-inversion comparisons are covered separately through `p=3`.

```bash
python -m validation.kinflation_SRI.run --p 10
```

[validation.ipynb](validation.ipynb) displays the saved report and can rerun
the case. [report.yaml](report.yaml) records residuals, settings and measured
errors. The [model guide](../README.md) explains the physical assumptions and
the [shared tools](../../common/README.md) check the universal moment generator.
