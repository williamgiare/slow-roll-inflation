"""Published equation syntax and new equation exports share one YAML layout."""
from pathlib import Path
import pytest
import sympy as sp

from scripts.archive_results import LatexPolynomial, archive_document
from scripts.export_results import generate
from utils import get


@pytest.mark.parametrize('type,kinflation', [('canonical',False),('kinflation',True)])
def test_new_document_export_includes_readable_yaml(tmp_path,type,kinflation):
    report=generate(root=tmp_path,max_p=1,max_n=1,kinflation=kinflation)
    assert report['status']=='passed'
    for row in report['orders']:
        assert (tmp_path/row['yaml']).is_file()
        power=get.spectra(row['p'],type=type,root=tmp_path/'results')
        relation=get.consistency(1,row['p'],type=type,root=tmp_path/'results')
        assert power.p==row['p'] and relation.W==row['p']+2
        original=(tmp_path/row['tex']).read_bytes()
        archive_document(tmp_path/row['tex'],row['p'],type)
        assert (tmp_path/row['tex']).read_bytes()==original


def test_latex_parser_handles_denominators_powers_and_scientific_notation():
    source=r'- \frac{1.25 \cdot 10^{-12} r^{2}}{c_{s}^{2}} + 0.5 \alpha_{0}^{c_s} \log{\left(c_{s} \right)}^{2}'
    formula=LatexPolynomial(source).formula()
    reader=get._FormulaReader({'r':sp.Symbol('r'),'cs':sp.Symbol('c_s'),'alpha_0_cs':sp.Symbol('alpha_0^cs')})
    expression=reader.formula(formula)
    r,cs,a=sp.symbols('r c_s alpha_0^cs')
    expected=-sp.Float('1.25e-12',17)*r**2/cs**2+sp.Float('.5',17)*a*sp.log(cs)**2
    assert sp.expand(expression-expected)==0


@pytest.mark.parametrize('source',[r'\sin{r}',r'\sqrt{r}',r'\unknown{r}',r'\frac{r}'])
def test_latex_parser_rejects_unsupported_or_incomplete_syntax(source):
    with pytest.raises(ValueError):LatexPolynomial(source).formula()


@pytest.mark.parametrize('kinflation', [False, True])
def test_export_cli_model_selection(tmp_path, monkeypatch, kinflation):
    import sys
    from scripts import export_results
    from sri import serialization as structured

    calls = []
    def fake_generate(**options):
        calls.append(options)
        return {'status': 'passed', 'model': 'kinflation' if options['kinflation'] else 'canonical'}

    reports = tmp_path / 'validation/common/reports'
    reports.mkdir(parents=True)
    monkeypatch.setattr(export_results, 'ROOT', tmp_path)
    monkeypatch.setattr(export_results, 'generate', fake_generate)
    monkeypatch.setattr(sys, 'argv', ['export_results', *(['--kinflation'] if kinflation else [])])
    export_results.main()
    assert calls[0]['kinflation'] is kinflation
    filename = 'results_export_kinflation.yaml' if kinflation else 'results_export.yaml'
    assert structured.loads((reports / filename).read_text())['model'] == ('kinflation' if kinflation else 'canonical')
