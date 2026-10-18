#========
# WG Edit
# Fast scientific regressions for the general analytic coefficient generator.
# Full weight-eight to weight-ten comparisons are exercised by run.py.
#========
import pytest
import sympy as sp

from . import run


def test_two_insertion_limit_and_auxiliary_poles():
    jet=run.exact.directional_generator((1,1),1)
    assert jet[-2]==jet[-1]==0
    assert sp.expand(jet[0]-run.fast.constant_moment(2))==0


def test_mixed_coefficient_against_fresh_mode_integration():
    expression=run.exact.family(3,2)[(0,1,1)]
    assert not expression.has(sp.Float,sp.Derivative,sp.Integral)
    assert expression.free_symbols<={run.fast.C}
    reference=run.generate_moments(5,1e-10).values[(0,1,1)]
    value=run.exact.evaluate(expression,24)
    assert abs(float(value)-reference)<1e-9


def test_standard_function_cache_preserves_values(tmp_path,monkeypatch):
    monkeypatch.setattr(run,'ROOT',tmp_path)
    monkeypatch.setattr(run,'_GPL_DATABASE',None)
    expression=run.exact.family(3,2)[(0,1,1)]
    try:
        cold,stats=run.cached_evaluate(expression,24)
        warm,reused=run.cached_evaluate(expression,24)
        uncached=run.exact.evaluate(expression,24)
        assert cold==warm
        assert abs(cold-uncached)<sp.Float('1e-23')
        assert stats['new']>0 and reused['new']==0
    finally:
        if run._GPL_DATABASE is not None: run._GPL_DATABASE.close()


def test_undetermined_coefficients_are_rejected():
    pattern=(0,1,1)
    with pytest.raises(AssertionError,match='Undetermined'):
        run.checked_expression(pattern,sp.Symbol('J_unknown'),{pattern:0.},{pattern:0.})


def test_published_report_without_private_cache_is_recomputed(tmp_path,monkeypatch):
    monkeypatch.setattr(run,'ROOT',tmp_path)
    path=tmp_path/'reports/families/m02_j00.yaml'
    run.atomic_yaml(path,{'signature':run.signature(),'status':'passed',
                         'checks':[{'pattern':[0,0],'analytic_value':'invalid saved value',
                                    'evaluation_seconds':0.,'status':'passed'}]})
    assert run.read_family(2,0,run.signature()) is None
    run.family_worker(2,0)
    report=run.structured.loads(path.read_text())
    assert report['status']=='passed'
    assert report['checks'][0]['analytic_value']!='invalid saved value'
    assert (tmp_path/'cache/m02_j00.pkl').exists()


def test_nonfinite_coefficients_are_rejected():
    with pytest.raises(ArithmeticError,match='Nonfinite'):
        run.checked_expression((0,),sp.nan,{(0,):0.},{(0,):0.})


def test_stale_expression_cache_is_removed_before_generation(tmp_path,monkeypatch):
    monkeypatch.setattr(run,'ROOT',tmp_path)
    report=tmp_path/'reports/families/m02_j00.yaml'
    cache=tmp_path/'cache/m02_j00.pkl'
    cache.parent.mkdir(parents=True)
    cache.write_bytes(run.pickle.dumps({(0,0):sp.Symbol('invalid_old_coefficient')}))
    run.atomic_yaml(report,{'signature':'outdated source','status':'passed','checks':[]})
    original=run.exact.family
    def interrupted(*args): raise RuntimeError('interrupted generation')
    monkeypatch.setattr(run.exact,'family',interrupted)
    with pytest.raises(RuntimeError,match='interrupted'):
        run.family_worker(2,0)
    assert not cache.exists()
    monkeypatch.setattr(run.exact,'family',original)
    run.family_worker(2,0)
    assert run.structured.loads(report.read_text())['status']=='passed'


def test_time_budget_preserves_partial_checks_without_passing_order(tmp_path,monkeypatch):
    monkeypatch.setattr(run,'ROOT',tmp_path)
    partial=tmp_path/'reports/families/m03_j02.yaml'
    run.atomic_yaml(partial,{'signature':run.signature(),'status':'running','checks':[]})
    def timeout(command,**kwargs):
        raise run.subprocess.TimeoutExpired(command,kwargs['timeout'])
    monkeypatch.setattr(run.subprocess,'run',timeout)
    with pytest.raises(TimeoutError,match='exceeded'):
        run.run_order(5,.01,3600)
    assert run.structured.loads(partial.read_text())['status']=='stopped'
    assert not (tmp_path/'reports/p05.yaml').exists()


# WG Edit: verify that uninterrupted runs impose no subprocess timeout.
def test_zero_budgets_disable_execution_timeouts(tmp_path,monkeypatch):
    monkeypatch.setattr(run,'ROOT',tmp_path)
    captured=[]
    def record(command,**kwargs):
        captured.append(kwargs['timeout'])
        raise RuntimeError('stop after inspecting timeout')
    monkeypatch.setattr(run.subprocess,'run',record)
    with pytest.raises(RuntimeError,match='inspecting timeout'):
        run.run_order(5,0,0)
    assert captured==[None]
    assert not (tmp_path/'reports/p05.yaml').exists()


# WG Edit: progress reports must preserve the evaluated coefficient.
def test_evaluation_progress_preserves_values(tmp_path,monkeypatch):
    monkeypatch.setattr(run,'ROOT',tmp_path)
    monkeypatch.setattr(run,'_GPL_DATABASE',None)
    expression=run.exact.family(3,2)[(0,1,1)]
    updates=[]
    try:
        value,stats=run.cached_evaluate(expression,24,
            progress=lambda done,total:updates.append((done,total)))
        expected=run.exact.evaluate(expression,24)
        assert abs(value-expected)<sp.Float('1e-23')
        total=len(expression.atoms(run.exact.GPL))
        assert updates[0]==(0,total)
        assert updates[-1]==(total,total)
        assert all(0<=done<=count==total for done,count in updates)
    finally:
        if run._GPL_DATABASE is not None: run._GPL_DATABASE.close()
