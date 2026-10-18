#========
# WG Edit
# Document exports use the same optimized arithmetic as the calculation API.
# Keep these entry points free of global method replacement.
#========
from contextlib import contextmanager


def horner_compose(ring, poly, replacements, order):
    return ring.compose(poly, replacements, order)


def difference_operator(ring, replacements, order):
    return ring.substitution_difference(replacements, order)


def fast_eliminate(inverse, poly):
    return inverse.eliminate(poly)


@contextmanager
def export_algebra():
    """Use the native exact arithmetic without altering any calculation method."""
    yield
