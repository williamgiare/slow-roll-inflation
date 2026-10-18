# N01LO: variable scalar sound speed

The first corrections include both scalar and sound-speed scale derivatives and the acoustic-to-tensor reference transport. An independently derived analytic tensor-tilt expression tests its coefficients. Generic relations are checked for `n=0,...,3`.

```bash
python -m validation.kinflation_SRI.run --p 1
```

[validation.ipynb](validation.ipynb) displays the saved report and can rerun
the case. [report.yaml](report.yaml) records residuals, settings and measured
errors. The [model guide](../README.md) explains the physical assumptions and
the [shared tools](../../common/README.md) check the universal moment generator.
