#========
# WG Edit
# Require every production moment to come from mode integration, including NLO.
# Analytic formulas are used only after generation as independent comparisons.
#========
import importlib
from fractions import Fraction
import sympy as sp
from validation.common.symbolic.checks import coefficient_comparison
from validation.common.symbolic.reference import reference_consistency


def test_production_regenerates_every_low_order_moment(monkeypatch):
    mc = importlib.import_module('sri.mode_coefficients')
    sm = importlib.import_module('sri.spectra')
    cm = importlib.import_module('sri.consistency')
    gm = importlib.import_module('sri.generic_modes')
    assert not hasattr(mc, 'evaluated_moments')
    calls = []
    original = mc.generate_moments
    def traced(order,tolerance):
        result = original(order,tolerance)
        calls.append(order)
        return result
    mc.coefficient_data.cache_clear()
    sm.mode_calculation.cache_clear()
    sm._spectra.cache_clear()
    cm.calculation.cache_clear()
    gm._generate_moments.cache_clear()
    monkeypatch.setattr(mc,'generate_moments',traced)
    try:
        for p in range(5):
            data = mc.coefficient_data(p,1e-10)
            assert set(data.expressions)==set(gm.source_patterns(p))
            assert data.expressions=={key:sp.Rational(Fraction.from_float(value))
                                     for key,value in data.generated.values.items()}
            assert data.metadata()['analytic_cache_moments']==0
            assert data.metadata()['numerically_generated_moments']==len(data.expressions)
            power = sm.spectra(p=p)
            assert power.coefficient_metadata['analytic_cache_moments']==0
            for n in range(3):
                order,inverse,source,out = cm.calculation(n,p)
                ring = inverse.modes.jet.ring
                assert not ring.add(ring.compose(out,inverse.forward,order.W),ring.scale(source,-1))
                actual = cm.consistency(n=n,p=p)
                if p<=3:
                    coefficient_comparison(actual.expression,reference_consistency(n,p),
                                           (actual.r_symbol,*actual.scalar_symbols))
        assert calls==list(range(5))
    finally:
        mc.coefficient_data.cache_clear()
        sm.mode_calculation.cache_clear()
        sm._spectra.cache_clear()
        cm.calculation.cache_clear()


def test_installed_source_runs_without_validation_or_reference_files():
    import os
    import subprocess
    import sys
    from pathlib import Path
    code = '''
import sys
import math
class BlockReferenceImports:
    def find_spec(self,fullname,path=None,target=None):
        if fullname=='validation' or fullname.startswith('validation.'):
            raise AssertionError('Production attempted to import validation')
sys.meta_path.insert(0,BlockReferenceImports())
from sri import spectra,consistency
for p in range(4):
    assert spectra(p=p).coefficient_metadata['analytic_cache_moments']==0
    for n in (0,1,2):
        result=consistency(n=n,p=p)
        assert result.coefficient_metadata['analytic_cache_moments']==0
        assert math.isfinite(result.evaluate(r=.02,scalar=[-.035]+[0.]*result.W))
print('generated-only production passed')
'''
    # WG Edit: test the installed package independently of checkout paths.
    env = os.environ.copy()
    env.pop('PYTHONPATH', None)
    result = subprocess.run([sys.executable,'-c',code],cwd='/private/tmp',env=env,
                            capture_output=True,text=True,timeout=60)
    assert result.returncode==0,result.stdout+result.stderr
    assert 'generated-only production passed' in result.stdout
