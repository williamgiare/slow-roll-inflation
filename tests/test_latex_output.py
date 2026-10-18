#========
# WG Edit
# Protect standalone document sections, algebra-preserving grouping and clean layout.
#========
from pathlib import Path
from types import SimpleNamespace
import pytest
import sympy as sp
from sri import spectra, consistency
from sri.cli import parser, render
from sri.latex_output import results_document, spectrum_groups, order_name

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize('p', [0, 1, 3, 4])
def test_grouped_spectra_reconstruct_the_generated_polynomial(p):
    power = spectra(p=p)
    for field in ('scalar', 'tensor'):
        groups = spectrum_groups(power, field)
        assert sp.expand(sum(expr for _, expr in groups) - getattr(power, field + '_correction')) == 0
        assert [degree for degree, _ in groups] == list(range(p + 1))


def test_document_contains_both_spectrum_subsections_and_requested_derivatives():
    power = spectra(p=2)
    relations = [consistency(n=n, p=2) for n in range(4)]
    text = results_document(power, relations)
    assert text.index(r'\section{Spectra}') < text.index(r'\subsection{Scalar spectrum}')
    assert text.index(r'\subsection{Scalar spectrum}') < text.index(r'\subsection{Tensor spectrum}')
    assert text.index(r'\subsection{Tensor spectrum}') < text.index(r'\section{Consistency relations}')
    for n in range(4):
        assert rf'\subsection{{$n={n}$}}' in text
    assert r'\allowdisplaybreaks' in text and r'\end{document}' in text


@pytest.mark.parametrize('relations', [lambda: [consistency(n=0, p=1)],
                                     lambda: [consistency(n=0, p=0)] * 2])
def test_mixed_orders_or_repeated_indices_are_rejected(relations):
    with pytest.raises(ValueError):
        results_document(spectra(p=0), relations())


def test_total_time_includes_equation_preparation(monkeypatch):
    import sri.cli as cli
    import sri.terminal as terminal
    ticks = iter([100.0, 101.23456])
    clock = SimpleNamespace(perf_counter=lambda: next(ticks))
    monkeypatch.setattr(cli, 'time', clock)
    monkeypatch.setattr(terminal, 'time', clock)
    text = render(parser().parse_args(['-n', '0', '-p', '0']))
    assert text.splitlines()[-1] == 'Total calculation time: 1.235 s'


@pytest.mark.parametrize('model', ['canonical_SRI','kinflation_SRI'])
@pytest.mark.parametrize('p', range(11))
def test_released_documents_have_the_agreed_contents_and_clean_files(p,model):
    name = order_name(p)
    folder = ROOT / 'results' / model / name
    assert {path.name for path in folder.iterdir()} == {name + '.tex', name + '.pdf', name + '.yaml'}
    assert (folder / (name + '.pdf')).read_bytes().startswith(b'%PDF-')
    text = (folder / (name + '.tex')).read_text()
    assert r'\subsection{Scalar spectrum}' in text and r'\subsection{Tensor spectrum}' in text
    count = 4 if p <= 6 else 1
    for n in range(count):
        assert rf'\subsection{{$n={n}$}}' in text
    for n in range(count, 4):
        assert rf'\subsection{{$n={n}$}}' not in text


def test_validation_contains_no_equation_export_directories():
    assert not list((ROOT / 'validation').glob('*/results'))
    assert not (ROOT / 'validation/canonical_SRI/N10LO/results.tex').exists()

#========
# WG Edit
# Verify sound-speed weights, normalization and model-separated document export.
#========
@pytest.mark.parametrize('p', range(4))
def test_kinflation_spectrum_groups_use_scale_derivative_weights(p):
    power=spectra(p=p,kinflation=True)
    for field in ('scalar','tensor'):
        groups=spectrum_groups(power,field)
        assert [degree for degree,_ in groups]==list(range(p+1))
        assert sp.expand(sum(expr for _,expr in groups)-getattr(power,field+'_correction'))==0
        variables=[sp.Symbol(f'epsilon_{j}') for j in range(1,p+2)]+list(power.sound_symbols)
        weights=[1]*(p+1)+list(range(1,len(power.sound_symbols)+1))
        for degree,expr in groups:
            assert all(sum(a*w for a,w in zip(ex,weights))==degree
                       for ex,coefficient in sp.Poly(expr,*variables).terms() if coefficient)


def test_kinflation_document_normalization_symbols_and_model_matching(tmp_path):
    from sri.latex_output import export_document
    power=spectra(p=2,kinflation=True)
    relation=consistency(n=0,p=2,kinflation=True)
    path=export_document(power,[relation],root=tmp_path)
    assert path==tmp_path/'results/kinflation_SRI/N02LO/N02LO.tex'
    text=path.read_text()
    assert r'\epsilon_{1*}c_{s*}' in text
    assert r'\alpha_{0}^{c_s}' in text and r'\alpha_n^{c_s}' in text
    assert 'SRI (k-inflation)' in text and 'acoustic reference' in text
    with pytest.raises(ValueError,match='physical model'):
        results_document(power,[consistency(n=0,p=2)])

# WG Edit: a long equation must keep its terms while using bounded TeX alignments.
def test_long_polynomial_uses_continuation_blocks_without_extra_equalities():
    from sri.latex_output import polynomial_block,polynomial_rows
    x=sp.Symbol('x')
    expression=sum(sp.Rational(j,17)*x**j for j in range(1,401))
    text=polynomial_block('F',expression,symbol_names={})
    assert text.count('\\begin{align*}')>1
    assert text.count(' &= ')==1
    assert r'\phantom{F}' in text
    rows=polynomial_rows(expression,symbol_names={})
    for row in rows:
        assert row in text
