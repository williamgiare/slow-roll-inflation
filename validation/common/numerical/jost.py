#========
# WG Edit
# Validate the complete normalized determinant against independent mode and pump ODEs.
# Retain the already verified compact-source calculation and refine the quadrature.
#========
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import solve_ivp
from scipy.linalg import lu_factor

def y0(x):
    return (1+1j/x)*np.exp(1j*x)

def dy0(x):
    return np.exp(1j*x)*(1j-1/x-1j/x**2)

def compact_source(x, amplitude):
    lower, upper = .2, 4.
    x = np.asarray(x)
    phase = np.pi*(x-lower)/(upper-lower)
    return np.where((x >= lower) & (x <= upper),
                    amplitude*np.sin(phase)**2*np.exp(-((x-.85)/.65)**2), 0.)

def direct_jost(amplitude):
    lower, upper = .2, 4.
    def mode_rhs(x, state):
        return (state[1], ((2+compact_source(x, amplitude))/x**2-1)*state[0])
    mode = solve_ivp(mode_rhs, (upper, lower), (y0(upper), dy0(upper)),
                     method='DOP853', rtol=2e-13, atol=2e-14)
    assert mode.success
    value, derivative = mode.y[:, -1]
    regular, regular_derivative = y0(lower).real, dy0(lower).real
    jost_bulk = (regular*derivative-regular_derivative*value)/1j
    def pump_rhs(x, state):
        return (state[1], (2+compact_source(x, amplitude))/x**2*state[0])
    pump = solve_ivp(pump_rhs, (lower, 1.), (1/lower, -1/lower**2),
                     method='DOP853', rtol=2e-13, atol=2e-14)
    assert pump.success
    return jost_bulk*pump.y[0, -1]

def fredholm_jost(amplitude, nodes_per_interval):
    t, w = leggauss(nodes_per_interval)
    grids, weights = [], []
    for lo, hi in ((.2, 1.), (1., 4.)):
        grids.append((hi+lo)/2+(hi-lo)/2*t)
        weights.append((hi-lo)/2*w)
    x, weights = np.concatenate(grids), np.concatenate(weights)
    lo, hi = np.minimum.outer(x, x), np.maximum.outer(x, x)
    G_R = -1j*y0(lo).real*y0(hi)
    G_P = np.where(hi < 1., (hi**2-1/hi)/(3*lo), 0.)
    v_weight = compact_source(x, amplitude)/x**2*weights
    A, B = G_R*v_weight[None, :], G_P*v_weight[None, :]
    C = A+B-np.einsum('ik,kj->ij', A, B, optimize=False)
    assert np.all(np.isfinite(C))
    lu, pivots = lu_factor(np.eye(len(x))-C)
    diagonal = np.diag(lu)
    sign = (-1)**int(np.count_nonzero(pivots != np.arange(len(x))))
    determinant = sign*np.prod(diagonal)
    assert np.isfinite(determinant)
    return determinant

def determinant_checks():
    rows = []
    for amplitude in (.04, -.03):
        exact = direct_jost(amplitude)
        errors = []
        for nodes in (24, 48, 96):
            determinant = fredholm_jost(amplitude, nodes)
            error = abs(determinant/exact-1)
            errors.append(error)
            rows.append({'amplitude': amplitude, 'nodes_per_interval': nodes,
                         'mode_J': [float(exact.real), float(exact.imag)],
                         'determinant_J': [float(determinant.real), float(determinant.imag)],
                         'relative_error': float(error)})
        assert errors[-1] < 1e-8, errors
        assert errors[1] < .4*errors[0] and errors[2] < .4*errors[1], errors
    return rows

def run():
    return {"determinant_vs_independent_odes": determinant_checks()}
