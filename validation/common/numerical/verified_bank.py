#========
# WG Edit
# Reuse one converged production bank across all lower-order validation cases.
# The generation step is itself validated independently. This is a scoped
# validation-only cache, never an alternative source of physical coefficients.
#========
from contextlib import contextmanager
from sri.generic_modes import GeneratedMoments


@contextmanager
def reuse_bank(bank):
    import sri.mode_coefficients as mc
    from sri.spectra import mode_calculation, _spectra
    from sri.consistency import calculation
    # WG Edit: all physical-model caches must follow the scoped coefficient bank.
    from sri.kinflation import mode_calculation as sound_modes, calculation as sound_calculation
    original=mc.generate_moments
    def cached(order,tolerance=1e-10):
        if order>bank.order or tolerance<bank.tolerance:return original(order,tolerance)
        keep=lambda k:len(k)+sum(k)<=order
        return GeneratedMoments(order,{k:v for k,v in bank.values.items() if keep(k)},
          {k:v for k,v in bank.differences.items() if keep(k)},bank.tolerance,bank.runs,bank.precision_digits)
    mc.coefficient_data.cache_clear();mode_calculation.cache_clear();_spectra.cache_clear();calculation.cache_clear()
    sound_modes.cache_clear();sound_calculation.cache_clear()
    mc.generate_moments=cached
    try:yield
    finally:
        mc.generate_moments=original
        mc.coefficient_data.cache_clear();mode_calculation.cache_clear();_spectra.cache_clear();calculation.cache_clear()
        sound_modes.cache_clear();sound_calculation.cache_clear()
#========
