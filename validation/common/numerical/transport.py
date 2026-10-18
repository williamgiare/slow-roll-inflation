#========
# WG Edit
# Validate the LO integral against an independent tensor ODE, numerical
# logarithmic derivatives and transported anchors on both sides of the pivot.
#========
import numpy as np
from scipy.integrate import quad, solve_ivp


def run():
    As, r = 2.1e-9, .025
    def scalar(x):
        return As * np.exp(-.035*x + .0003*x*x) * (1 + .035*np.sin(.7*x+.2))
    anchor = r * scalar(0)
    def transported(x, reference=0., normalization=anchor):
        integral, _ = quad(lambda y: 1/scalar(y), reference, x,
                            epsabs=1e-3, epsrel=5e-13, limit=100)
        return normalization / (1 + normalization*integral/8)
    max_ode, max_relation, max_anchor = 0., 0., 0.
    for endpoint in (-8., 8.):
        points = np.linspace(0, endpoint, 101)
        solution = solve_ivp(lambda x,z: [-anchor*z[0]**2/(8*scalar(x))],
            (0,endpoint), [1.], t_eval=points, method='DOP853', rtol=2e-13, atol=1e-14)
        assert solution.success
        prediction = np.array([transported(x) for x in points])
        max_ode = max(max_ode, float(np.max(np.abs(prediction/(anchor*solution.y[0])-1))))
    for x in (-7., -3., 0., 2., 7.):
        step = .003
        values = [np.log(transported(x + j*step)) for j in (-2,-1,1,2)]
        derivative = (values[0] - 8*values[1] + 8*values[2] - values[3])/(12*step)
        max_relation = max(max_relation, abs(-8*derivative/(transported(x)/scalar(x))-1))
    for pivot in (-4., 3.):
        new_normalization = transported(pivot)
        for x in (5., -6., 0., 2., -3.):  # deliberately unordered queries
            max_anchor = max(max_anchor, abs(transported(x,pivot,new_normalization)/transported(x)-1))
    assert max_ode < 2e-11
    assert max_relation < 2e-8
    assert max_anchor < 2e-12
    return {'ode_vs_integral_max_relative_error': max_ode,
            'finite_difference_consistency_max_relative_error': max_relation,
            'transported_anchor_max_relative_error': max_anchor,
            'range_ln_k': [-8,8], 'unordered_queries': 'passed'}
