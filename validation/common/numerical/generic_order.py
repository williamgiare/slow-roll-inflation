#========
# WG Edit
# Independent analytic moment oracles, exact power-law inflation, and public
# predictions at orders beyond the former table. No SRI moment equation is used
# to construct the Gamma/Hankel oracle coefficients.
#========
from sri import serialization as structured
from math import factorial, comb
import time
import mpmath as mp
import sympy as sp
from sri import spectra, consistency
from sri.generic_modes import generate_moments, high_precision_integrate, NormalizedModeSystem
from ..symbolic.moment_reference import evaluated_moments
from sri.background import CONSTANT_VALUES, C


def analytic_moments(order):
    with mp.workdps(60):
        # Exact Gamma Mellin generator: explicit Bell/Taylor coefficients.
        ell = [mp.mpf(0)]
        for n in range(1,order+1):
            ell.append(-(mp.euler+mp.log(2)-2) if n==1 else
                -(1-2**(1-n))*mp.zeta(n)/n if n%2==0 else (2-mp.zeta(n))/n)
        A = [mp.mpf(1)]
        for n in range(1,order+1):
            A.append(sum(k*ell[k]*A[n-k] for k in range(1,n+1))/n)
        chi = [2*sum(A[k]/mp.mpf(3)**(n+2-k) for k in range(1,n+2))
               for n in range(order)]
        oracle = {(j,):float(factorial(j)*chi[j]) for j in range(order)}
        # Constant source: normalized Hankel amplitude, no source-insertion ODE.
        def hankel(c):
            nu = mp.sqrt(mp.mpf('2.25')+c)
            return 2*(mp.loggamma(nu)-mp.loggamma(mp.mpf('1.5'))) +(2*nu-3)*mp.log(2)
        for n,value in enumerate(mp.taylor(hankel,0,order)):
            if n:
                oracle[(0,)*n] = float(value)
        if order >= 2:
            maximum = order-2
            # Closed two-insertion kernel; removable poles resolved as Taylor jets.
            tangent = mp.taylor(lambda z:mp.tan(mp.pi*z/2)/z if z else mp.pi/2,0,maximum)
            R = [2*sum(A[j]*tangent[k]/mp.mpf(3)**(n-j-k+1)
                    for j in range(n+1) for k in range(n-j+1)) for n in range(maximum+1)]
            denominator = [(-1/mp.mpf(4)**(n+1)-(-1)**n)/5 for n in range(maximum+1)]
            f = [sum(A[j]*denominator[n-j] for j in range(n+1)) for n in range(maximum+1)]
            g = [sum((4*chi[j]-(2 if j==0 else 0))*denominator[n-j]
                    for j in range(n+1)) for n in range(maximum+1)]
            for j in range(maximum+1):
                for k in range(j,maximum-j+1):
                    value = R[j]*R[k]/2
                    value += 2*sum(comb(a+b,a)*f[a+b] for a in range(j+1) for b in range(k+1))
                    value += sum(comb(a+b,a)*g[a+b]/mp.mpf(3)**(j+k-a-b+2)
                                 for a in range(j+1) for b in range(k+1))
                    oracle[(j,k)] = float(value*factorial(j)*factorial(k))
        return oracle


def moment_checks(order=7, tolerance=1e-10):
    generated = generate_moments(order,tolerance)
    expected = analytic_moments(order)
    errors = {str(k):abs(generated.values[k]-value) for k,value in expected.items()}
    for key,value in expected.items():
        assert abs(generated.values[key]-value) <= tolerance*(1+abs(value)),(key,value,generated.values[key])
    for key,value in evaluated_moments().items():
        assert abs(generated.values[key]-float(value.subs(CONSTANT_VALUES))) < tolerance
    return {'order':order,'analytic_oracle_moments':len(errors),
            'maximum_absolute_oracle_error':max(errors.values()),'oracle_errors':errors,
            'generator':generated.metadata()}


def high_precision_check():
    system = NormalizedModeSystem(4,14)
    calculated = high_precision_integrate(system,upper=80,lower=1e-8,step=.25,stages=8)
    errors = {str(k):abs(calculated[k]-v) for k,v in analytic_moments(4).items()}
    assert max(errors.values())<2e-11
    return {'order':4,'precision_digits':40,'maximum_absolute_oracle_error':max(errors.values())}


def power_law_checks(maximum=7):
    rows=[]
    e=sp.Symbol('epsilon_1')
    with mp.workdps(60):
        for p in range(maximum+1):
            result=spectra(p=p)
            values={symbol:0 for symbol in result.tensor_correction.free_symbols | result.scalar_correction.free_symbols
                    if str(symbol).startswith('epsilon_') and symbol!=e}
            tensor=sp.expand(result.tensor_correction.subs(values))
            scalar=sp.expand(result.scalar_correction.subs(values))
            # Identical scalar/tensor pumps in exact power-law inflation.
            difference = sp.Poly(sp.expand(tensor-scalar),e)
            assert all(abs(float(c.subs(CONSTANT_VALUES)))<2e-12 for c in difference.all_coeffs())
            nu=(3-e)/(2*(1-e))
            exact=2**(2*nu-3)*(sp.gamma(nu)/sp.gamma(sp.Rational(3,2)))**2*(1-e)**2
            def exact_function(value):
                index=(3-value)/(2*(1-value))
                return 2**(2*index-3)*(mp.gamma(index)/mp.gamma(mp.mpf('1.5')))**2*(1-value)**2
            exact_coefficient=mp.taylor(exact_function,0,p)[p]
            actual_coefficient=float(tensor.coeff(e,p).subs(CONSTANT_VALUES).evalf(35))
            difference=abs(actual_coefficient-float(exact_coefficient))
            assert difference<3e-8*(1+abs(float(exact_coefficient))),(p,difference)
            errors=[]
            for value in (mp.mpf('.06'),mp.mpf('.03')):
                predicted=mp.mpf(str(tensor.subs(e,sp.Float(str(value),60)).subs(CONSTANT_VALUES).evalf(55)))
                errors.append(float(abs(predicted/exact_function(value)-1)))
            slope=float(mp.log(errors[0]/errors[1],2))
            assert slope>p+.75,(p,errors,slope)
            rows.append({'p':p,'highest_power_coefficient_difference':difference,
                         'relative_errors':errors,'observed_remainder_order':slope})
    return rows


def observable_checks():
    rows=[]
    for n,p in ((0,4),(1,4),(2,4),(0,5)):
        result=consistency(n=n,p=p)
        # Physical power-law inflation: r=16 epsilon, alpha_s0=nt and all runnings vanish.
        errors=[]
        for e in (.008,.004):
            predicted=result.evaluate(r=16*e,scalar=[-2*e/(1-e),*([0]*result.W)])
            exact=-2*e/(1-e) if n==0 else 0.
            errors.append(abs(predicted-exact))
        expected_order=result.W+1
        ratio=errors[0]/max(errors[1],1e-25)
        assert ratio>2**(expected_order-.8),(n,p,errors,ratio)
        assert not any(str(s).startswith(('E','q','J_','epsilon_')) for s in result.expression.free_symbols)
        assert result.W==n+1+p
        rows.append({'n':n,'p':p,'W':result.W,'absolute_errors':errors,
                     'observed_remainder_order':float(mp.log(ratio,2))})
    return rows


def featured_mode_checks():
    # A canonical, nonconstant background is propagated by separate physical
    # scalar/tensor ODEs. Reconstruct higher epsilon_i directly from its profile.
    import numpy as np
    from .physical_modes import integrated_spectra, profile
    u=sp.Symbol('u')
    flow=[sp.Rational(1,10)*(1-u*u)]
    for _ in range(5):
        flow.append(sp.cancel(sp.Rational(2,25)*(1-u*u)*sp.diff(flow[-1],u)/flow[-1]))
    rows=[]
    for lam in (1.,.5):
        coarse,_=integrated_spectra(lam,600.)
        exact,freeze=integrated_spectra(lam,1200.)
        hierarchy=[profile(0,lam)[0],*[lam*float(v.subs(u,np.tanh(.45))) for v in flow]]
        row={'lambda':lam,'independent_mode_amplitudes':exact,'freeze_errors':freeze,'errors':{}}
        for field in ('scalar','tensor'):
            assert abs(coarse[field]/exact[field]-1)<2e-10
            assert freeze[field]<2e-9
        for p in (3,4,5):
            powers=spectra(p=p)
            for field in ('scalar','tensor'):
                expression=getattr(powers,field+'_correction')
                substitutions={sp.Symbol(f'epsilon_{j+1}'):v for j,v in enumerate(hierarchy)}
                predicted=float(expression.evalf(22,subs=substitutions|CONSTANT_VALUES))
                row['errors'][field+str(p)]=abs(predicted/exact[field]-1)
        rows.append(row)
    slopes=[]
    for field in ('scalar','tensor'):
        for p in (4,5):
            slope=float(np.log2(rows[0]['errors'][field+str(p)]/rows[1]['errors'][field+str(p)]))
            assert slope>p+.5,(field,p,slope)
            assert rows[0]['errors'][field+str(p)]<rows[0]['errors'][field+str(p-1)]
            slopes.append({'field':field,'p':p,'observed_remainder_order':slope})
    return {'runs':rows,'remainder_scaling':slopes}


def formal_inverse_checks():
    # Check the complete coordinate substitution, keeping every universal
    # moment symbolic. The residual must be identically zero, not just small.
    from sri.consistency import calculation
    rows=[]
    for n,p in ((0,4),(1,4),(2,4),(0,5)):
        order,inverse,source,observable=calculation(n,p)
        ring=inverse.modes.jet.ring
        residual=ring.add(ring.compose(observable,inverse.forward,order.W),ring.scale(source,-1))
        assert residual=={},(n,p,len(residual))
        rows.append({'n':n,'p':p,'W':order.W,'symbolic_substitution_residual':'0'})
    return rows


def run(maximum=7):
    start=time.monotonic()
    result={'generic_moments':moment_checks(maximum),'arbitrary_precision':high_precision_check(),
            'power_law_amplitudes':power_law_checks(maximum),'observable_relations':observable_checks(),
            'symbolic_inverse':formal_inverse_checks(),'independent_featured_modes':featured_mode_checks()}
    result['seconds']=time.monotonic()-start
    return result


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument('--maximum',type=int,default=7)
    parser.add_argument('--output',type=str)
    args=parser.parse_args()
    report=run(args.maximum)
    rendered=structured.dumps(report,indent=2)+'\n'
    if args.output:
        from pathlib import Path
        Path(args.output).write_text(rendered)
    else:
        print(rendered)
#========
