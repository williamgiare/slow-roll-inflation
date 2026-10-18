# N02LO: variable scalar sound speed

Quadratic scalar and tensor amplitudes, the acoustic pump mapping and generic consistency relations for `n=0,...,3` are checked exactly. The two independent algebra engines agree before coefficient rounding.

```bash
python -m validation.kinflation_SRI.run --p 2
```

[validation.ipynb](validation.ipynb) displays the saved report and can rerun
the case. [report.yaml](report.yaml) records residuals, settings and measured
errors. The [model guide](../README.md) explains the physical assumptions and
the [shared tools](../../common/README.md) check the universal moment generator.
