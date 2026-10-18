#========
# WG Edit
# Compare SRI amplitudes with independently integrated canonical scalar/tensor modes.
# A smooth epsilon_1(N) reconstructs a canonical potential; scaling lambda
# scales all background derivatives by their slow-roll weights.
#========
import numpy as np
import mpmath as mp
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq
from scipy.special import hankel1, h1vp
from sri import spectra
from sri.background import C


def constant_background():
    rows = []
    with mp.workdps(50):
        for p in range(4):
            expression = spectra(p=p).tensor_correction
            errors = []
            for e in (mp.mpf('.04'), mp.mpf('.02'), mp.mpf('.01')):
                nu = (3-e)/(2*(1-e))
                exact = 2**(2*nu-3)*(mp.gamma(nu)/mp.gamma(mp.mpf('1.5')))**2*(1-e)**2
                substitutions = {str(x): 0 for x in expression.free_symbols if str(x).startswith('epsilon_')}
                substitutions['epsilon_1'] = e
                substitutions['C'] = mp.euler+mp.log(2)-2
                # Arbitrary-precision evaluation avoids cancellation in the N3LO remainder.
                import sympy as sp
                evaluated = expression.subs({sp.Symbol(k):sp.Float(str(v),55)
                                             for k,v in substitutions.items()}).evalf(50)
                predicted = mp.mpf(str(evaluated))
                errors.append(abs(predicted/exact-1))
            slopes = [float(mp.log(errors[i]/errors[i+1],2)) for i in range(2)]
            assert min(slopes) > p+.8, (p,slopes)
            rows.append({'p':p, 'relative_errors':list(map(float,errors)),
                         'observed_remainder_orders':slopes})
    return rows


def profile(N, lam):
    e0,b,c,d = .02,.10,.08,.45
    z=c*lam*N+d
    u=np.tanh(z)
    sech=1-u*u
    epsilon=lam*e0*np.exp(b/c*(u-np.tanh(d)))
    q1=lam*b*sech
    q2=-2*lam**2*b*c*sech*u
    q3=-2*lam**3*b*c*c*sech*(1-3*u*u)
    return epsilon,q1,q2,q3


def integrated_spectra(lam, start_x):
    H0=1e-5
    future,past=60.,-13.
    logarithm=-quad(lambda N:profile(N,lam)[0],0,future,
                     epsabs=2e-13,epsrel=2e-13)[0]
    def background_rhs(N,state):
        epsilon=profile(N,lam)[0]
        return (-epsilon,(1-epsilon)*state[1]-1)
    solution=solve_ivp(background_rhs,(future,past),
        (logarithm,1/(1-profile(future,lam)[0])),dense_output=True,
        method='DOP853',rtol=3e-13,atol=2e-14)
    assert solution.success
    def background(N):
        logH,h=solution.sol(N)
        return H0*np.exp(logH),h
    h0=background(0)[1]
    k=H0/h0
    def x(N):
        H,h=background(N)
        return h*k/(np.exp(N)*H)
    start=brentq(lambda N:x(N)-start_x,past,0,xtol=2e-13)
    output={}
    freeze_errors={}
    for field in ('scalar','tensor'):
        def pump_terms(N):
            epsilon,q1,q2,_=profile(N,lam)
            # z_N/z and its derivative, evaluated directly from the physical pump.
            a=1+q1/2 if field=='scalar' else 1.
            da=q2/2 if field=='scalar' else 0.
            return a*a+da+(1-epsilon)*a
        H,h=background(start)
        nu=np.sqrt(.25+h*h*pump_terms(start))
        phase=np.exp(1j*np.pi*(2*nu+1)/4)
        normalization=phase*np.sqrt(np.pi*start_x/2)
        y=normalization*hankel1(nu,start_x)
        yx=normalization*(h1vp(nu,start_x)+hankel1(nu,start_x)/(2*start_x))
        yn=-k/(np.exp(start)*H)*yx
        def rhs(N,state):
            epsilon=profile(N,lam)[0]
            H,_=background(N)
            K=k/(np.exp(N)*H)
            return (state[1],-(1-epsilon)*state[1]-(K*K-pump_terms(N))*state[0])
        mode=solve_ivp(rhs,(start,14.),(y,yn),dense_output=True,
            method='DOP853',rtol=2e-12,atol=2e-13)
        assert mode.success
        def amplitude(N):
            unit=(k/H0)**2*abs(mode.sol(N)[0]/np.exp(N))**2
            if field=='scalar':unit*=profile(0,lam)[0]/profile(N,lam)[0]
            return float(unit)
        output[field]=amplitude(14.)
        freeze_errors[field]=abs(amplitude(12.)/output[field]-1)
    return output,freeze_errors


def featured_background():
    rows=[]
    errors={field:{p:[] for p in range(4)} for field in ('scalar','tensor')}
    for lam in (1.,.5,.25):
        coarse,_=integrated_spectra(lam,600.)
        exact,freeze=integrated_spectra(lam,1200.)
        epsilon,q1,q2,q3=profile(0,lam)
        hierarchy=[epsilon,q1,q2/q1,q3/q2-q2/q1]
        for field in ('scalar','tensor'):
            start_error=abs(coarse[field]/exact[field]-1)
            assert start_error<2e-7,(lam,field,start_error)
            assert freeze[field]<2e-9,(lam,field,freeze[field])
            row={'lambda':lam,'field':field,'mode_amplitude':exact[field],
                 'initial_time_relative_difference':start_error,
                 'freeze_out_relative_difference':freeze[field],'slow_roll_errors':{}}
            for p in range(4):
                import sympy as sp
                expression=getattr(spectra(p=p),field+'_correction')
                substitutions={sp.Symbol(f'epsilon_{j+1}'):value for j,value in enumerate(hierarchy)}
                substitutions[C]=float(mp.euler+mp.log(2)-2)
                predicted=float(expression.evalf(18,subs=substitutions))
                error=abs(predicted/exact[field]-1)
                errors[field][p].append(error)
                row['slow_roll_errors'][str(p)]=error
            rows.append(row)
        print(f'Numerical modes: lambda={lam:g}',flush=True)
    slopes=[]
    for field in errors:
        for p,values in errors[field].items():
            orders=[float(np.log2(values[j]/values[j+1])) for j in range(2)]
            # The smallest N3LO residual may reach the initial-time numerical floor.
            assert orders[0]>p+.4,(field,p,values,orders)
            assert values[2]<values[0],(field,p,values)
            slopes.append({'field':field,'p':p,'observed_remainder_orders':orders})
    return {'amplitudes':rows,'remainder_scaling':slopes,
            'model':'epsilon1=lambda*0.02*exp[(0.10/0.08)(tanh(0.08*lambda*N+0.45)-tanh(0.45))]',
            'canonical_reconstruction':'dphi/dN=Mpl*sqrt(2*epsilon1); V=Mpl^2*H^2*(3-epsilon1)',
            'initial_x_values':[600,1200],'late_N_values':[12,14]}


def run():
    return {'constant_background':constant_background(),
            'independent_featured_modes':featured_background()}
