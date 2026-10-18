#========
# WG Edit
# Exact differential tests cover arbitrary substitutions, zero-weight constants,
# truncated Taylor sums, coefficient cancellation and bounded cache ownership.
#========
from fractions import Fraction
import gc
import random
import weakref
import pytest
from sri.series import WeightedRing, _weighted_degree
from validation.common.symbolic.reference_arithmetic import (
    reference_arithmetic, reference_add, reference_mul, reference_compose,
    reference_exponential, reference_logarithm_unit)
from sri.consistency import calculation
from sri import spectra


def random_polynomial(rng, size, terms=14):
    result={}
    for _ in range(terms):
        ex=tuple(rng.randrange(4) for _ in range(size))
        result[ex]=Fraction(rng.randrange(-5,6),rng.randrange(1,8))
    return {ex:c for ex,c in result.items() if c}


@pytest.mark.parametrize('weights',[(1,2,0),(0,1,3),(1,1,1),(0,0,2)])
def test_sparse_arithmetic_matches_independent_reference(weights):
    rng=random.Random(173);ring=WeightedRing(weights)
    for _ in range(20):
        a=random_polynomial(rng,3);b=random_polynomial(rng,3)
        maps=[random_polynomial(rng,3,4) for _ in weights]
        maps[rng.randrange(3)]={} if rng.randrange(2) else ring.constant(Fraction(2,3))
        for W in (0,2,5,9):
            assert ring.add(a,b,ring.scale(a,-1))==reference_add(ring,a,b,ring.scale(a,-1))
            for left,right in ((a,b),(a,a),(ring.constant(1),a),(ring.variable(1),a)):
                assert ring.mul(left,right,W)==reference_mul(ring,left,right,W)
            assert ring.compose(a,maps,W)==reference_compose(ring,a,maps,W)


@pytest.mark.parametrize('weights',[(1,2,0),(1,1,3)])
def test_horner_taylor_sums_match_all_retained_coefficients(weights):
    rng=random.Random(71);ring=WeightedRing(weights)
    for _ in range(12):
        poly={ex:c for ex,c in random_polynomial(rng,3,8).items() if ring.degree(ex)>0}
        for W in range(11):
            assert ring.exponential(poly,W)==reference_exponential(ring,poly,W)
            unit=ring.add(ring.constant(1),poly)
            assert ring.logarithm_unit(unit,W)==reference_logarithm_unit(ring,unit,W)


def test_weight_cache_does_not_retain_ring_instances():
    WeightedRing.clear_degree_cache()
    ring=WeightedRing((1,2,0));ref=weakref.ref(ring)
    assert ring.degree([2,3,11])==8
    ring.degree((2,3,11))
    assert _weighted_degree.cache_info().hits>=1
    assert _weighted_degree.cache_info().maxsize==131072
    del ring;gc.collect();assert ref() is None
    WeightedRing.clear_degree_cache()
    assert _weighted_degree.cache_info().currsize==0


def test_substitution_difference_rejects_wrong_size_or_nonraising_shifts():
    ring=WeightedRing((1,));x=ring.variable(0)
    with pytest.raises(ValueError):ring.substitution_difference([],2)
    with pytest.raises(ArithmeticError):ring.substitution_difference([ring.constant(1)],2)(x)


@pytest.mark.parametrize('kinflation',[False,True])
@pytest.mark.parametrize('p',range(5))
def test_public_spectra_and_exact_logs_match_reference(p,kinflation):
    with reference_arithmetic():
        before=spectra(p=p,kinflation=kinflation)
    after=spectra(p=p,kinflation=kinflation)
    for field in ('scalar','tensor','scalar_correction','tensor_correction','scalar_log_correction','tensor_log_correction'):
        if hasattr(before,field):assert getattr(before,field)==getattr(after,field)


@pytest.mark.parametrize('kinflation',[False,True])
@pytest.mark.parametrize('n,p',[(0,1),(1,2),(2,3)])
def test_closed_engine_preserves_its_independent_result(n,p,kinflation):
    with reference_arithmetic():
        before=calculation(n,p,method='closed',kinflation=kinflation)
    after=calculation(n,p,method='closed',kinflation=kinflation)
    assert before[2:]==after[2:]
    assert before[1].forward==after[1].forward
    assert before[1].modes.log_amplitudes==after[1].modes.log_amplitudes
