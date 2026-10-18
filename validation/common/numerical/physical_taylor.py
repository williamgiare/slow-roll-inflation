#========
# WG Edit
# Independent parameter-Taylor modes for an actual canonical background.
# Check coefficients individually instead of inferring an order from potentially
# cancelling finite-parameter residuals. Full physical modes remain a separate test.
#========
from fractions import Fraction as F
from math import factorial
import mpmath as mp
from .linear_modes import add,scale,multiply,derivative,inverse_derivative,LinearModeOracle


def jet_add(*jets):
    n=max(map(len,jets));return [add(*(v[k] for v in jets if k<len(v))) for k in range(n)]
def jet_scale(a,c):return [scale(v,c) for v in a]
def jet_mul(a,b,p):return [add(*(multiply(a[j],b[k-j]) for j in range(k+1) if j<len(a) and k-j<len(b))) for k in range(p+1)]
def unit_inverse(a,p):
    assert a[0]==[F(1)]
    out=[[F(1)]]
    for k in range(1,p+1):out.append(scale(add(*(multiply(a[j],out[k-j]) for j in range(1,k+1))),-1))
    return out

def jet_log(a,p):
    assert a[0]==[F(1)]
    out=[[]]
    for k in range(1,p+1):
        out.append(add(a[k],scale(add(*(scale(multiply(a[j],out[k-j]),k-j) for j in range(1,k))),-F(1,k))))
    return out

def jet_exp(a,p):
    assert not a[0]
    out=[[F(1)]]
    for k in range(1,p+1):out.append(scale(add(*(scale(multiply(a[j],out[k-j]),j) for j in range(1,k+1))),F(1,k)))
    return out


def physical_jets(p,e=F(2,25),b=F(73,1000)):
    v=[[]]+[[F(0)]*(k-1)+[e*b**(k-1)/factorial(k-1)] for k in range(1,p+1)]
    one=[[F(1)]]+[[] for _ in range(p)]
    h=jet_add(one,v);hinv=unit_inverse(h,p)
    bjet=[[],[b]]+[[] for _ in range(p-1)]
    E=jet_mul(jet_mul(v,jet_add(h,bjet),p),jet_mul(hinv,hinv,p),p)
    Ft=[[]]+[[F(0)]*k+[-e*b**(k-1)/factorial(k)] for k in range(1,p+1)]
    logh=jet_log(h,p);log_extra=jet_log(jet_add(h,bjet),p)
    logeps=jet_add([[]]+[[F(0),b]]+[[] for _ in range(p-1)],log_extra,jet_scale(logh,-2))
    # Subtract the reference constants; only the physical growing normalization remains.
    logeps=[[F(0),*a[1:]] if a else [] for a in logeps]
    Fs=jet_add(Ft,jet_scale(logeps,F(1,2)))
    ss=[derivative(a) for a in Fs];st=[derivative(a) for a in Ft]
    q=jet_scale(jet_mul(jet_add(ss,jet_scale(st,-1)),hinv,p),-2)
    qjets=[q]
    for _ in range(p):q=jet_scale(jet_mul([derivative(a) for a in q],hinv,p),-1);qjets.append(q)
    return {'scalar':Fs,'tensor':Ft,'h':h,'E':E,'q':qjets}


def physical_oracle(p,field,initial_terms=20):
    jets=physical_jets(p);Fjet=jets[field];s=[derivative(a) for a in Fjet]
    g=jet_add([derivative(a) for a in s],jet_scale(s,-3),jet_mul(s,s,p))
    Z=jet_exp(Fjet,p);V=unit_inverse(Z,p)
    oracle=object.__new__(LinearModeOracle)
    oracle.order=p;oracle.size=p+1;oracle.basis=tuple((j,) for j in range(p+1));oracle.index={ex:k for k,ex in enumerate(oracle.basis)}
    oracle.weights=list(range(p+1));oracle.pairs=[[(j,k-j) for j in range(1,k)] for k in range(p+1)]
    oracle.multiplicities=[1]*(p+1);oracle.Z=Z;oracle.inverse=V;oracle.s=s
    raw=[]
    for n in range(1,initial_terms+1):
        row=[]
        for k in range(p+1):
            if n==1:rhs=scale(g[k],F(-1,2))
            else:
                prev=raw[-1][k];first=add(derivative(prev),scale(prev,-(n-1)))
                rhs=scale(add(derivative(first),scale(first,-n)),F(1,2))
                rhs=add(rhs,scale(add(*(multiply(g[j],raw[-1][k-j]) for j in range(1,k+1))),F(-1,2)))
                rhs=add(rhs,*(scale(add(derivative(raw[j-1][k]),scale(raw[j-1][k],-j)),-(-1)**(n-2-j)) for j in range(1,n-1)))
            row.append(inverse_derivative(rhs,n))
        raw.append(row)
    oracle.raw=raw
    return oracle,jets


def amplitude_coefficients(log_values,p,e=F(2,25)):
    with mp.workdps(60):
        logs=[mp.mpf(0),*[mp.mpf(log_values[(0,)*k]) for k in range(1,p+1)]]
        units=[mp.mpf(1)]
        for k in range(1,p+1):units.append(sum(j*logs[j]*units[k-j] for j in range(1,k+1))/k)
        e=mp.mpf(e.numerator)/e.denominator
        out=[sum(units[j]*(-1)**(k-j)*(k-j+1)*e**(k-j) for j in range(k+1)) for k in range(p+1)]
        return [mp.nstr(v,50) for v in out]


def sri_parameter_coefficients(modes,field,p):
    from sri.background import CONSTANT_VALUES
    jets=physical_jets(p)
    params=[jets['E'],*jets['q'][:modes.jet.depth]]
    def project(jet):return [a[0] if a else F(0) for a in jet]
    values=[project(v) for v in params]
    for symbol in modes.jet.symbols[len(params):]:
        if symbol in CONSTANT_VALUES:constant=mp.mpf(str(CONSTANT_VALUES[symbol].evalf(65)))
        else:constant=mp.mpf(str(modes.coefficients.values[symbol]))
        values.append([constant]+[mp.mpf(0)]*p)
    with mp.workdps(60):
        values=[[mp.mpf(v.numerator)/v.denominator if isinstance(v,F) else v for v in a] for a in values]
        def mul(a,b):return [mp.fsum(a[j]*b[k-j] for j in range(k+1)) for k in range(p+1)]
        powers={};out=[mp.mpf(0)]*(p+1)
        unit=modes.jet.ring.exponential(modes.log_amplitudes[field],p)
        for ex,coefficient in unit.items():
            term=[mp.mpf(coefficient.numerator)/coefficient.denominator]+[mp.mpf(0)]*p
            for index,count in enumerate(ex):
                if not count:continue
                key=(index,count)
                if key not in powers:
                    v=[mp.mpf(1)]+[mp.mpf(0)]*p
                    for _ in range(count):v=mul(v,values[index])
                    powers[key]=v
                term=mul(term,powers[key])
            out=[a+b for a,b in zip(out,term)]
        return [mp.nstr(v,50) for v in out]
#========
