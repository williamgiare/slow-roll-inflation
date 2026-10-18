# ========
# WG Edit
# Check installed package isolation, relocated exports and report invalidation.
# ========
from pathlib import Path
import pytest
import sri
from sri import spectra
from sri.latex_output import export_document
from validation.canonical_SRI import order_cases


def test_package_uses_flat_source_directory():
    root = Path(__file__).resolve().parents[1]
    assert Path(sri.__file__).resolve() == root / 'source' / '__init__.py'


# WG Edit: equation documents belong to results/, independent of validation.
def test_spectrum_exports_follow_results_layout(tmp_path):
    target = export_document(spectra(p=0), root=tmp_path)
    assert target == tmp_path / 'results/canonical_SRI/LO/LO.tex'
    text = target.read_text()
    assert r'\section{Spectra}' in text
    assert r'\subsection{Scalar spectrum}' in text
    assert r'\subsection{Tensor spectrum}' in text
    assert 'P_s(k_*)' in text and 'P_t(k_*)' in text


@pytest.mark.parametrize('p', [0, 1, 2, 3, 10])
def test_stored_order_results_match_calculation_sources(p):
    assert order_cases.saved_report(p)['status'] == 'passed'


def test_stored_results_reject_modified_calculation_sources(monkeypatch):
    monkeypatch.setattr(order_cases, 'fingerprint', lambda: {})
    with pytest.raises(RuntimeError, match='stale'):
        order_cases.saved_report(10)
