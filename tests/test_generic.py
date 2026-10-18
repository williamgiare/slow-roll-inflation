#========
# WG Edit
# Exercise generated higher-order coefficients, combinatorics, and the complete
# public path beyond N3LO. Analytic oracles do not import the mode generator.
#========
from sri import serialization as structured
from math import factorial, prod
import pytest
import sympy as sp
from sri import spectra, consistency, ConsistencyOrder
from sri.cli import parser, render
from sri.generic_modes import source_patterns, generate_moments, NormalizedModeSystem, ModeRun
from validation.common.numerical.generic_order import analytic_moments
from validation.common.symbolic.checks import coefficient_comparison


def test_partition_basis_is_complete_without_ordered_duplicates():
    from itertools import product
    for p in range(1,9):
        expected=set()
        for rank in range(1,p+1):
            for values in product(range(p-rank+1),repeat=rank):
                if rank+sum(values)<=p:
                    expected.add(tuple(sorted(values)))
        actual=list(source_patterns(p))
        assert len(actual)==len(set(actual))
        assert set(actual)==expected


def test_general_moments_against_independent_gamma_and_hankel():
    generated=generate_moments(7)
    assert len(generated.values)==44
    for key,value in analytic_moments(7).items():
        assert abs(generated.values[key]-value) < 1e-10*(1+abs(value))


def test_growing_endpoint_does_not_change_generated_coefficients():
    system=NormalizedModeSystem(5,12)
    a=system.integrate(ModeRun(lower=1e-6))
    b=system.integrate(ModeRun(lower=1e-9))
    assert max(abs(a[k]-b[k]) for k in a)<2e-8


def test_order_six_spectra_have_new_terms_and_evaluate():
    result=spectra(p=6)
    assert result.p==6 and result.coefficient_metadata['numerically_generated_moments']>0
    assert not any(str(v).startswith('J_') for v in result.scalar.free_symbols)
    values=result.evaluate(H=1e-5,epsilon=[.005,.03,-.02,.01,.01,.01,.01])
    assert values['scalar']>0 and values['tensor']>0
    assert sp.Symbol('epsilon_7') in result.scalar.free_symbols


def test_order_four_consistency_is_not_order_three_with_new_metadata():
    low=consistency(n=0,p=3)
    high=consistency(n=0,p=4)
    assert high.W==5 and high.components[-1][0]==5
    assert high.components[-1][1]!=0
    coefficient_comparison(sum(v for _,v in high.components[:-1]),low.expression,
                           (high.r_symbol,*high.scalar_symbols))
    assert high.evaluate(r=.02,scalar=[-.035,-.0005,.0001,-.00001])<0
    assert high.evaluate(r=0,scalar=[-.035,-.0005,.0001,-.00001])==0


def test_order_four_command_infers_p_and_reports_coefficient_accuracy():
    parsed=structured.loads(render(parser().parse_args(['-W','5','-n','0','--format','yaml'])))
    assert parsed['p']==4
    assert parsed['coefficients']['backend']=='normalized mode equations'
    assert ConsistencyOrder.resolve(n=1,p=7).W==9


@pytest.mark.parametrize('p',[-1,True,1.5])
def test_generator_checks_order_before_cache(p):
    with pytest.raises(ValueError):
        generate_moments(p)


@pytest.mark.parametrize('value',[float('nan'),0,-1,1e-16])
def test_invalid_coefficient_accuracy(value):
    with pytest.raises(ValueError):
        spectra(p=4,coefficient_tolerance=value)
#========
