#========
# WG Edit
# Validate k-inflation by exact inverse identities, canonical recovery, an
# independently constructed effective tensor pump and exact Hankel spectra.
# High-order profiles are physical restrictions, not a dense symbolic scan.
#========
import argparse
from fractions import Fraction
from functools import lru_cache
from math import factorial
from pathlib import Path
import hashlib
from sri import serialization as structured
import time
import mpmath as mp
import sympy as sp
from sri import spectra, consistency
from sri.kinflation import (KinflationJet, Operations, assemble_modes,
                         sound_scale_replacements, scale_log_amplitudes)
from sri.spectra import mode_calculation
from sri.consistency import calculation
from sri.mode_coefficients import coefficient_data, insertion_patterns
from validation.canonical_SRI.high_order import load_bank
from validation.common.numerical.verified_bank import reuse_bank

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent
NAMES={0:'LO',1:'N01LO',2:'N02LO',3:'N03LO',10:'N10LO'}


def fingerprint():
    paths=[*sorted((ROOT/'source').rglob('*.py')),Path(__file__)]
    return {str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


def saved_report(p):
    data=structured.loads((HERE/NAMES[p]/'report.yaml').read_text())
    from validation.common.report_compatibility import compatible_report
    current=fingerprint()
    if data['source_sha256']!=current and not compatible_report(ROOT,data['source_sha256'],current):
        raise RuntimeError('The k-inflation validation is stale; rerun it for the current source.')
    if data['status']!='passed':raise RuntimeError('The k-inflation validation did not pass.')
    return data


def canonical_recovery(p,method='iterative'):
    a=spectra(p=p,kinflation=True,method=method);b=spectra(p=p,method=method)
    zero={a.cs_symbol:1,**dict.fromkeys(a.sound_symbols,0)}
    for field in ('scalar','tensor'):
        assert sp.expand(getattr(a,field).subs(zero)-getattr(b,field))==0
    return {'p':p,'method':method,'scalar_residual':'0','tensor_residual':'0'}


#========
# WG Edit
# Check the public sound-coordinate conversion with independent forward jets.
# Source moments and native amplitudes remain in the original N coordinates.
#========
def sound_coordinate_check(p):
    jet=KinflationJet(p+1);r,v=jet.ring,jet.variables
    op=Operations(jet,'iterative');h=op.clock(v[0],p)
    forward=list(v)
    for j in range(1,p+1):
        forward[jet.sound_start+j-1]=op.scale(v[jet.ell_index],j,h,p)
    inverses=[sound_scale_replacements(jet,h,p,method) for method in ('iterative','closed')]
    assert inverses[0]==inverses[1]
    for j in range(1,p+1):
        index=jet.sound_start+j-1
        assert r.compose(forward[index],inverses[0],p)==v[index]
        assert r.compose(inverses[0][index],forward,p)==v[index]
    return {'p':p,'derivatives_checked':p,'forward_inverse_residual':'0',
            'inverse_forward_residual':'0','iterative_closed_residual':'0',
            'denominators_in_sound_parameters':'none'}


def spectrum_coordinate_check(modes,p):
    jet,r,v=modes.jet,modes.jet.ring,modes.jet.variables
    op=Operations(jet,modes.method);forward=list(v)
    for j in range(1,p+1):
        forward[jet.sound_start+j-1]=op.scale(v[jet.ell_index],j,modes.h,p)
    converted=scale_log_amplitudes(modes,p)
    for field,value in converted.items():
        assert r.compose(value,forward,p)==modes.log_amplitudes[field]
        native_unit=op.exponential(modes.log_amplitudes[field],p)
        public_unit=op.exponential(value,p)
        assert r.compose(public_unit,forward,p)==native_unit
    return {'p':p,'method':modes.method,'log_amplitude_residual':'0','amplitude_residual':'0'}


def profile_coordinate_check(modes,p):
    # WG Edit: apply the generic conversion before restricting to a physical profile.
    jet=KinflationJet(p+1,modes.coefficients.symbols);r=jet.ring;v=jet.variables
    op=Operations(jet,'iterative');h=op.clock(v[0],p)
    inv=sound_scale_replacements(jet,h,p)
    assert inv==sound_scale_replacements(jet,h,p,'closed')
    forward=list(modes.jet.variables);physical=Operations(modes.jet,'iterative')
    for j in range(1,p+1):
        forward[jet.sound_start+j-1]=physical.scale(v[jet.ell_index],j,modes.h,p)
    for field,value in modes.log_amplitudes.items():
        converted=r.compose(value,inv,p)
        assert r.compose(converted,forward,p)==value
    return {'p':p,'log_amplitude_residual':'0','iterative_closed_conversion_residual':'0'}


def relation_check(n,p,compare=True):
    a=calculation(n,p,kinflation=True)
    order,inv,source,out=a;ring=inv.modes.jet.ring
    assert not ring.add(ring.compose(out,inv.forward,order.W),ring.scale(source,-1))
    if compare:
        b=calculation(n,p,kinflation=True,method='closed')
        assert a[1].forward==b[1].forward and a[1].baseline==b[1].baseline
        assert a[1].shifted==b[1].shifted and a[2:]==b[2:]
    value=consistency(n=n,p=p,kinflation=True)
    canon=consistency(n=n,p=p)
    assert sp.expand(value.expression.subs({value.cs_symbol:1,
                        **dict.fromkeys(value.sound_symbols,0)})-canon.expression)==0
    assert value.expression.subs(value.r_symbol,0)==0
    allowed={value.r_symbol,value.cs_symbol,*value.scalar_symbols,*value.sound_symbols}
    assert value.expression.free_symbols <= allowed
    return {'n':n,'p':p,'W':order.W,'terms':len(out),
            'exact_forward_residual':'0','canonical_residual':'0',
            'iterative_closed_residual':'0' if compare else 'not checked',
            'observable_only':'passed','zero_r':'passed'}


def effective_tensor_check(modes,p):
    """Independent scalar-to-tensor mapping using a_eff=a sqrt(E/c_s)."""
    jet,r,v=modes.jet,modes.jet.ring,modes.jet.variables
    Q=r.add(r.constant(1),r.scale(v[1],Fraction(1,2)),r.scale(v[jet.sound_start],Fraction(-1,2)))
    from sri.closed.algebra import unit_power
    invQ=unit_power(r,Q,-1,p)
    D=lambda poly:r.mul(invQ,jet.derivative(poly,p),p)
    effE=r.mul(invQ,r.add(v[0],r.scale(v[1],Fraction(1,2)),
                   r.scale(v[jet.sound_start],Fraction(1,2)),
                   r.scale(r.mul(invQ,jet.derivative(Q,p),p),-1)),p)
    # Canonical conformal clock constructed in the effective background.
    h=term=r.constant(1)
    for _ in range(p):
        term=r.add(r.mul(effE,term,p),D(term));h=r.add(h,term)
    source=r.add(r.mul(r.power(h,2,p),r.add(r.constant(2),r.scale(effE,-1)),p),r.constant(-2))
    assert source==modes.sources['scalar']
    jets=[source]
    for _ in range(1,p):jets.append(r.scale(r.mul(h,D(jets[-1]),p),-1))
    jets=[r.scale(value,Fraction(1,factorial(j))) for j,value in enumerate(jets)]
    log={}
    for pattern in insertion_patterns(p):
        factor=factorial(len(pattern))
        for j in set(pattern):factor//=factorial(pattern.count(j))
        term=r.scale(r.polynomial(modes.coefficients.expressions[pattern],jet.symbols),factor)
        for j in pattern:term=r.mul(term,jets[j],p)
        log=r.add(log,term)
    log=r.add(log,r.scale(r.logarithm_unit(h,p),-2),r.scale(r.logarithm_unit(Q,p),2))
    assert log==modes.native_log_amplitudes['scalar']
    return {'source_residual':'0','native_amplitude_residual':'0',
            'construction':'effective scale factor a sqrt(epsilon_1/c_s); canonical tensor clock and source',
            'moment_backend':'shared; this is an independent source and normalization check'}


def profile_modes(p,data,profile,method):
    jet=KinflationJet(p+1,data.symbols);r=jet.ring;v=jet.variables
    E,q1,s1=v[0],v[1],v[jet.sound_start]
    for j in range(1,p+2):
        v[j]=({} if profile=='constant' else r.scale(r.power(q1,j,p),factorial(j-1)))
        if j>1:
            v[jet.sound_start+j-1]=({} if profile=='constant' else
                r.scale(r.mul(r.power(q1,j-1,p),s1,p),factorial(j-1)))
    return assemble_modes(p,jet,data,method)


def project_value(poly,modes,weight,directions,ell):
    value=mp.mpf('0')
    jet=modes.jet
    for powers,coef in poly.items():
        if jet.ring.degree(powers)!=weight:continue
        term=mp.mpf(coef.numerator)/coef.denominator
        for index,power in enumerate(powers):
            if not power:continue
            if index==jet.ell_index:term*=ell**power
            else:term*=directions.get(index,mp.mpf('0'))**power
        value+=term
    return value


def high_order_checks(bank):
    p=bank.order;data=coefficient_data(p);rows=[]
    for profile in ('constant','evolving'):
        a=profile_modes(p,data,profile,'iterative')
        b=profile_modes(p,data,profile,'closed')
        assert a.h==b.h and a.scalar_h==b.scalar_h
        assert a.sources==b.sources and a.log_amplitudes==b.log_amplitudes
        oracle=effective_tensor_check(a,p)
        coordinate_check=profile_coordinate_check(a,p)
        op=Operations(a.jet,'iterative');r=a.jet.ring;v=a.jet.variables
        zero=list(v)
        for j in range(a.jet.sound_start,a.jet.sound_end):zero[j]={}
        zero[a.jet.ell_index]={}
        canon=mode_calculation(p,p+1)
        # WG Edit: canonical recovery is checked at every weight on both profiles.
        for field in ('scalar','tensor'):
            embedded=r.polynomial(canon.jet.ring.expression(canon.log_amplitudes[field],canon.jet.symbols),a.jet.symbols)
            assert r.compose(a.log_amplitudes[field],zero,p)==r.compose(embedded,zero,p)
        if profile=='constant':
            # Exact Hankel coefficients for mixed E and s1, at one physical pivot.
            with mp.workdps(65):
                for c in ('0.5','0.8','1.2'):
                    cs=mp.mpf(c);ell=mp.log(cs)
                    for sound in ('-.4','.3','.8'):
                        z=mp.mpf(sound)
                        def exact(t,field):
                            E=t;s=t*z;ht=1/(1-E);hs=1/(1-E-s)
                            nu=(mp.mpf('.5')+ht) if field=='tensor' else (mp.mpf('.5')+hs*(1-s/2))
                            h=ht if field=='tensor' else hs
                            M=-2*mp.log(h)+(2*nu-1)*mp.log(2)+2*mp.loggamma(nu)-mp.log(mp.pi)
                            if field=='scalar':M+=-(2*E+s)/(1-E-s)*mp.log(cs*hs/ht)
                            return M
                        errors=[]
                        for field in ('scalar','tensor'):
                            coefficients=mp.taylor(lambda t:exact(t,field),0,p)
                            for weight in range(1,p+1):
                                actual=project_value(a.log_amplitudes[field],a,weight,
                                      {0:mp.mpf('1'),a.jet.sound_start:z},ell)
                                error=abs(actual-coefficients[weight])/(1+abs(coefficients[weight]))
                                assert error<mp.mpf('2e-9'),(c,sound,field,weight,error)
                                errors.append(float(error))
                        rows.append({'profile':profile,'c_s':c,'sound_direction':sound,
                                     'coefficient_weights':[1,p],'maximum_scaled_Hankel_error':max(errors),
                                     'public_sound_conversion':coordinate_check})
        else:
            rows.append({'profile':profile,'background':'E,q1,s1 proportional to 1/(1-q1_star N)',
                         'scalar_log_terms':len(a.log_amplitudes['scalar']),
                         'iterative_closed_residual':'0','effective_tensor_mapping':oracle,
                         'public_sound_conversion':coordinate_check})
    return {'status':'passed','p':p,'profiles':rows,
            'scope':'mixed constant-flow exact Hankel coefficients and nonconstant-flow source/amplitude checks; not a dense order-ten observable-inversion comparison'}


def nlo_check():
    a=consistency(n=0,p=1,kinflation=True);r,cs=a.r_symbol,a.cs_symbol
    C=sp.EulerGamma+sp.log(2)-2
    expected=-r/(8*cs)*(1-r/(16*cs)-(1-sp.log(cs))*a.scalar_symbols[0]
                         -(C+3)*a.sound_symbols[0])
    error=sp.expand(a.expression-expected)
    for params in [(.02,.7,-.035,.002),(.04,.5,-.02,-.003)]:
        sub={r:params[0],cs:params[1],a.scalar_symbols[0]:params[2],a.sound_symbols[0]:params[3]}
        assert abs(float(error.evalf(40,subs=sub)))<1e-12
    return {'status':'passed','equation':'n_t=-r/(8 c_s)[1-r/(16 c_s)-(1-ln c_s)alpha_0^s-(C+3)alpha_0^cs]'}


def nlo_spectrum_check():
    # WG Edit: the independent first-order acoustic amplitude at the tensor pivot.
    power=spectra(p=1,kinflation=True)
    E,q=sp.symbols('epsilon_1 epsilon_2');a=power.sound_symbols[0];cs=power.cs_symbol
    C=sp.EulerGamma+sp.log(2)-2
    scalar=1-2*(C+1)*E-C*q-(C+2)*a-(2*E+q+a)*sp.log(cs)
    tensor=1-2*(C+1)*E
    from validation.common.symbolic.checks import coefficient_comparison
    # Log(c_s) is an independent coefficient coordinate in these polynomial comparisons.
    ell=sp.Symbol('ell_cs')
    checked=[]
    for field,expected in [('scalar',scalar),('tensor',tensor)]:
        actual=getattr(power,field+'_correction')
        checked.append({'field':field,**coefficient_comparison(
            actual.xreplace({sp.log(cs):ell}),expected.xreplace({sp.log(cs):ell}),(E,q,a,ell))})
    return {'status':'passed','comparisons':checked}


def run(p,bank=None):
    if p not in NAMES:raise ValueError('Validation cases are LO, N01LO, N02LO, N03LO and N10LO.')
    start=time.monotonic()
    report={'status':'passed','model':'kinflation','p':p,'checks':{}}
    report['checks']['sound_coordinates']=sound_coordinate_check(p)
    if p<=3:
        report['checks']['canonical_spectra']=[canonical_recovery(p,m) for m in ('iterative','closed')]
        report['checks']['effective_tensor_mapping']=effective_tensor_check(mode_calculation(p,p+1,kinflation=True),p)
        report['checks']['spectrum_coordinates']=[spectrum_coordinate_check(
            mode_calculation(p,p+1,method=method,kinflation=True),p) for method in ('iterative','closed')]
        report['checks']['relations']=[relation_check(n,p) for n in (range(9) if p==0 else range(4))]
        if p:
            report['checks']['NLO_limit']=nlo_check()
            report['checks']['NLO_spectra']=nlo_spectrum_check()
    if p==10:
        if bank is None:
            from sri.generic_modes import generate_moments
            bank=generate_moments(10)
        with reuse_bank(bank):report['checks']['high_order']=high_order_checks(bank)
    report['source_sha256']=fingerprint();report['seconds']=time.monotonic()-start
    target=HERE/NAMES[p]/'report.yaml';target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(structured.dumps(report,indent=2)+'\n')
    print(f'{NAMES[p]} k-inflation: passed ({report["seconds"]:.2f} s).',flush=True)
    return report


# WG Edit: optional generic inversion sweep; high orders are substantially more expensive.
def generic_inverse_checks(max_p, bank=None):
    if max_p < 0:
        raise ValueError('The maximum SR order must be nonnegative.')
    rows=[]
    for p in range(max_p+1):
        start=time.monotonic()
        order,inv,source,out=calculation(0,p,kinflation=True)
        ring=inv.modes.jet.ring
        assert not ring.add(ring.compose(out,inv.forward,order.W),ring.scale(source,-1))
        rows.append({'p':p,'n':0,'W':order.W,'terms':len(out),
                     'exact_forward_residual':'0','seconds':time.monotonic()-start})
        print(f'Generic inverse p={p}: passed.',flush=True)
    report={'status':'passed','model':'kinflation','checks':rows,
            'scope':'dense generic iterative tensor-tilt inverse; no closed-engine comparison',
            'source_sha256':fingerprint()}
    (HERE/'generic_inverse.yaml').write_text(structured.dumps(report,indent=2)+'\n')
    return report


def main():
    parser=argparse.ArgumentParser(description='Independent k-inflation validation cases.')
    parser.add_argument('--p',type=int,choices=NAMES)
    parser.add_argument('--all',action='store_true')
    parser.add_argument('--generic-max-p',type=int,help='Optional dense iterative tensor-tilt inverse sweep.')
    parser.add_argument('--bank',type=Path,help='Reuse a mode-generated bank matching the unchanged coefficient kernel.')
    args=parser.parse_args()
    if not args.all and args.p is None and args.generic_max_p is None:
        parser.error('Supply --p, --all or --generic-max-p.')
    bank=load_bank(args.bank)[0] if args.bank else None
    if args.all or args.p is not None:
        for p in NAMES if args.all else [args.p]:run(p,bank)
    if args.generic_max_p is not None:
        if bank is None:generic_inverse_checks(args.generic_max_p)
        else:
            with reuse_bank(bank):generic_inverse_checks(args.generic_max_p,bank)


if __name__=='__main__':main()
