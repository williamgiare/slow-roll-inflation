# ========
# WG Edit
# Keep recovered high-n LO checks connected to SRI in the regular test suite.
# ========
import pytest
from validation.canonical_SRI.LO.large_n_validation.run import (verify_archive, symbolic_bridge, rational_bridge,
                              modular_bridge, projective_rational, finite_scale_values)


def test_recovered_files_are_unchanged():
    assert verify_archive() == 11


def test_sri_and_both_archived_constructions_match():
    assert symbolic_bridge(8) == 9


def test_rational_current_formula_matches_archived_flow():
    assert rational_bridge(30) == 31


def test_modular_current_formula_matches_archived_flow():
    assert modular_bridge(100) == 9


def test_recovered_projective_jet():
    assert projective_rational(100) == 101


@pytest.mark.parametrize('prime', [9, 7])
def test_reject_invalid_field(prime):
    with pytest.raises(ValueError):
        finite_scale_values(1, [0]*9, 8, prime)


def test_zero_tensor_normalization():
    assert all(v == 0 for v in finite_scale_values(0, [0, 2, -1, 3], 3))
