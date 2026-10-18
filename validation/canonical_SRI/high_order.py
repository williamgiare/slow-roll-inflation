#========
# WG Edit
# Full p-order validation: all universal moments, public amplitudes, n=0,1,2
# inverses, exact Hankel limits, independent physical modes and Taylor coefficients.
#========
import argparse, ast, hashlib, time, platform
from sri import serialization as structured
import importlib.metadata
from pathlib import Path
from math import factorial
import mpmath as mp
import sympy as sp
from sri import spectra,consistency
from sri.consistency import calculation
from sri.spectra import mode_calculation
from sri.generic_modes import generate_moments,GeneratedMoments,NormalizedModeSystem
from sri.background import CONSTANT_VALUES
from ..common.numerical.verified_bank import reuse_bank
from ..common.numerical.linear_modes import LinearModeOracle
from ..common.numerical.generic_order import analytic_moments
from ..common.numerical.physical_high_order import physical_modes,local_flow
from ..common.numerical.physical_taylor import physical_oracle,amplitude_coefficients,sri_parameter_coefficients


def kernel_signature():
    import sri.generic_modes as gm
    import sri.series as sr
    generator=Path(gm.__file__).read_bytes()
    tree=ast.parse(Path(sr.__file__).read_text())
    methods={'constant','variable','degree','cut','add','scale','mul','derivative'}
    arithmetic=[ast.dump(node,include_attributes=False) for cls in tree.body if isinstance(cls,ast.ClassDef)
                for node in cls.body if isinstance(node,ast.FunctionDef) and node.name in methods]
    # WG Edit: include the bounded weight helper in the numerical kernel signature.
    helpers=[ast.dump(node,include_attributes=False) for node in tree.body
             if isinstance(node,ast.FunctionDef) and node.name=='_weighted_degree']
    return hashlib.sha256(generator+'\n'.join(arithmetic+helpers).encode()).hexdigest()


def load_bank(path):
    data=structured.loads(Path(path).read_text())
    if data['producer_kernel_signature']!=kernel_signature():
        raise ValueError('Saved moment bank does not match the current coefficient-generation kernel.')
    meta=data['metadata']
    bank=GeneratedMoments(meta['order'],{tuple(k):v for k,v in data['values']},
          {tuple(k):v for k,v in data['differences']},meta['moment_tolerance'],tuple(meta['runs']),meta['integration_precision_digits'])
    return bank,data


def moment_checks(bank):
    start=time.monotonic()
    expected=bank.values
    oracle=LinearModeOracle(bank.order,18)
    production=NormalizedModeSystem(bank.order,18)
    # Exact comparison of independently generated growing-pump derivatives.
    for pattern,ex in zip(production.patterns,production.basis):
     row=oracle.s[oracle.index[ex]]
     actual=production.exact_s[production.patterns.index(pattern)]
     assert list(actual)==row+[0]*(len(actual)-len(row)),pattern
    print('All growing-pump polynomials agree exactly.',flush=True)
    settings=[dict(upper=160,lower=1e-11,step=.16,stages=10,digits=45),
              dict(upper=200,lower=1e-12,step=.125,stages=12,digits=60)]
    values=[]
    for config in settings:
     v=oracle.integrate(**config);values.append(v);print('Linear integration done',config,flush=True)
    convergence={k:abs(values[0][k]-values[1][k]) for k in expected}
    errors={k:abs(values[1][k]-expected[k]) for k in expected}
    for k,v in expected.items():
     assert convergence[k]<1e-10*(1+abs(v)),('linear convergence',k,convergence[k])
     assert errors[k]<1e-10*(1+abs(v)),('cross method',k,errors[k],v)
    known=analytic_moments(bank.order)
    analytic_errors={k:abs(values[1][k]-v) for k,v in known.items()}
    assert all(analytic_errors[k]<1e-10*(1+abs(v)) for k,v in known.items())
    return {'status':'passed','order':bank.order,'moments_checked':len(expected),
      'exact_pump_polynomials_checked':len(expected),
      'maximum_absolute_cross_method_difference':max(errors.values(),default=0),
      'maximum_scaled_cross_method_difference':max((errors[k]/(1+abs(expected[k])) for k in errors),default=0),
      'maximum_linear_run_difference':max(convergence.values(),default=0),
      'maximum_analytic_oracle_error':max(analytic_errors.values(),default=0),'runs':settings,
      'coefficient_comparisons':[{'log_powers':list(k),'production':expected[k],'linear':values[1][k],
          'cross_method_difference':errors[k],'linear_run_difference':convergence[k]} for k in expected],
      'seconds':time.monotonic()-start}


def algebra_checks(bank):
    start=time.monotonic();rows=[];spec=[]
    for p in range(bank.order+1):
     t=time.monotonic();power=spectra(p=p)
     variables=tuple(sp.Symbol(f'epsilon_{j}') for j in range(1,p+2))
     for field in ('scalar','tensor'):
      expr=getattr(power,field+'_correction')
      assert max((sum(ex) for ex,c in sp.Poly(expr,*variables).terms()),default=0)==p
      assert not any(str(s).startswith('J_') for s in expr.free_symbols)
     spec.append({'p':p,'scalar_terms':len(sp.Poly(power.scalar_correction,*variables).terms()),
                  'tensor_terms':len(sp.Poly(power.tensor_correction,*variables).terms()),
                  'seconds':time.monotonic()-t})
     print('Spectra',p,spec[-1],flush=True)
     for n in range(3):
      t=time.monotonic();order,inv,source,out=calculation(n,p);ring=inv.modes.jet.ring
      residual=ring.add(ring.compose(out,inv.forward,order.W),ring.scale(source,-1))
      assert residual=={},(p,n,len(residual))
      result=consistency(n=n,p=p)
      assert result.W==n+1+p
      assert not any(str(s).startswith(('J_','epsilon_','q')) or str(s)=='E' for s in result.expression.free_symbols)
      for weight,expr in result.components:
       variables=(result.r_symbol,*result.scalar_symbols)
       assert all(ex[0]+sum((i+1)*v for i,v in enumerate(ex[1:]))==weight
                  for ex,c in sp.Poly(expr,*variables).terms() if c)
      assert result.evaluate(r=0,scalar=[0]*result.W)==0
      rows.append({'n':n,'p':p,'W':order.W,'exact_rational_substitution_residual':'0',
                   'terms':len(out),'seconds':time.monotonic()-t})
      print('Consistency',n,p,rows[-1],flush=True)
    return {'status':'passed','p_range':[0,bank.order],'n_indices':[0,1,2],
      'spectra':spec,'relations':rows,
      'coefficient_data':'mode-generated production bank only; exact rational internal algebra without analytic restoration',
      'seconds':time.monotonic()-start}


def gamma_checks(bank):
    rows=[];start=time.monotonic()
    with mp.workdps(65):
        def exact(e):
            nu=(3-e)/(2*(1-e))
            return 2**(2*nu-3)*(mp.gamma(nu)/mp.gamma(mp.mpf('1.5')))**2*(1-e)**2
        coeff=mp.taylor(exact,0,bank.order)
        for p in range(bank.order+1):
            modes=mode_calculation(p,p+1);ring=modes.jet.ring
            projected={field:{ex:c for ex,c in poly.items() if not any(ex[1:modes.jet.depth+1])}
                       for field,poly in modes.log_amplitudes.items()}
            assert projected['scalar']==projected['tensor']
            unit=ring.exponential(projected['tensor'],p)
            expr=ring.expression(unit,modes.jet.symbols).subs(CONSTANT_VALUES)
            E=modes.jet.symbols[0]
            difference=abs(mp.mpf(str(expr.coeff(E,p).evalf(60)))-coeff[p])
            assert difference<mp.mpf('1e-9')*(1+abs(coeff[p]))
            errors=[]
            for e in map(mp.mpf,('.08','.04')):
                pred=mp.mpf(str(expr.subs(E,sp.Float(str(e),70)).evalf(65)))
                errors.append(abs(pred/exact(e)-1))
            observed=float(mp.log(errors[0]/errors[1],2))
            assert observed>p+.6,(p,observed,errors)
            rows.append({'p':p,'coefficient_absolute_error':float(difference),
                         'relative_amplitude_errors':list(map(float,errors)),
                         'observed_remainder_order':observed})
            print('Gamma',p,rows[-1],flush=True)
    return {'status':'passed','p_range':[0,bank.order],'cases':rows,'seconds':time.monotonic()-start}


def taylor_checks(bank):
    start=time.monotonic();rows=[]
    with mp.workdps(65):
     for field in ('scalar','tensor'):
      oracle,jets=physical_oracle(bank.order,field)
      lv=oracle.integrate(upper=150,lower=1e-13,step=.125,stages=12,digits=60)
      physical=amplitude_coefficients(lv,bank.order)
      for p in range(bank.order+1):
       sri=sri_parameter_coefficients(mode_calculation(p,p+1),field,p)
       differences=[float(abs(mp.mpf(sri[j])-mp.mpf(physical[j]))) for j in range(p+1)]
       normalized=[differences[j]/max(abs(float(physical[j])),1e-25) for j in range(p+1)]
       print(field,p,'highest',sri[-1],physical[p],'error',differences[-1],'relative',normalized[-1],flush=True)
       rows.append({'field':field,'p':p,'sri':sri,'linear_physical':physical[:p+1],
                    'absolute_coefficient_differences':differences,'relative_coefficient_differences':normalized})
    # Test coefficients themselves, not a fitted finite-parameter slope.
    for row in rows:
        assert max(row['relative_coefficient_differences'])<1e-8,row
    return {'status':'passed','p_range':[0,bank.order],'comparisons':rows,'seconds':time.monotonic()-start}


def physical_checks(bank):
    start=time.monotonic();rows=[]
    with mp.workdps(70):
     powers={p:spectra(p=p) for p in range(bank.order+1)}
     for lam in ('1','.5','.25'):
      e=mp.mpf('.08')*mp.mpf(lam);b=mp.mpf('.073')*mp.mpf(lam)
      coarse=physical_modes(str(e),str(b),upper=80,lower='1e-12',step='.20',stages=10,digits=50,initial_terms=18)
      fine=physical_modes(str(e),str(b),upper=120,lower='1e-15',step='.10',stages=12,digits=70,initial_terms=24,power_terms=100)
      flow=local_flow(str(e),str(b),bank.order+1,digits=75)
      errors={};predictions={};refinement={}
      for field in ('scalar','tensor'):
       exact=mp.mpf(fine['amplitudes'][field]);refinement[field]=float(abs(mp.mpf(coarse['amplitudes'][field])/exact-1))
       assert refinement[field]<mp.mpf('1e-20'),(lam,field,refinement[field])
      for p,power in powers.items():
       for field in ('scalar','tensor'):
        expr=getattr(power,field+'_correction')
        subs={sp.Symbol(f'epsilon_{j+1}'):sp.Float(v,75) for j,v in enumerate(flow)}
        pred=mp.mpf(str(expr.subs(subs).subs(CONSTANT_VALUES).evalf(70)))
        key=field+str(p);predictions[key]=mp.nstr(pred,50)
        errors[key]=float(abs(pred/mp.mpf(fine['amplitudes'][field])-1))
      row={'lambda':lam,'background_e':mp.nstr(e,20),'background_b':mp.nstr(b,20),
           'epsilon_star':fine['epsilon_star'],'flow_at_pivot':flow,'mode_refinement_errors':refinement,
           'exact_mode_amplitudes':fine['amplitudes'],'relative_errors':errors,'predictions':predictions,
           'mode_runs':[coarse['settings'],fine['settings']]}
      rows.append(row);print('Physical',lam,'refinement',refinement,'final order',errors['scalar'+str(bank.order)],errors['tensor'+str(bank.order)],flush=True)
    # Finite-parameter errors can cancel. The formal order is tested directly
    # by taylor_checks; here we require resolved modes and final-order accuracy.
    for row in rows:
        for field in ('scalar','tensor'):
            error=row['relative_errors'][field+str(bank.order)]
            if bank.order==10:
                assert error < 5e-9,row
                assert error<row['relative_errors'][field+'9'],row
    orders=[]
    for field in ('scalar','tensor'):
        for p in range(bank.order+1):
            errors=[row['relative_errors'][field+str(p)] for row in rows]
            slopes=[float(mp.log(errors[j]/errors[j+1],2)) for j in range(2)]
            orders.append({'field':field,'p':p,'observed_finite_parameter_slopes':slopes})
    return {'status':'passed','p_range':[0,bank.order],
      'canonical_background':'a(x)=x^(-1) exp[-e (x^b-1)/b]',
      'physical_reconstruction':'epsilon=v(1+v+b)/(1+v)^2, v=e*x^b; dphi/dN=sqrt(2 epsilon), V=H^2(3-epsilon)',
      'criterion':'mode refinement and finite-parameter residuals; formal SR order checked by independent Taylor coefficients',
      'runs':rows,'finite_parameter_error_slopes':orders,'seconds':time.monotonic()-start}


def main(argv=None):
    parser=argparse.ArgumentParser(description='Independent validation through a requested slow-roll correction order.')
    parser.add_argument('--p',type=int,default=10)
    parser.add_argument('--reuse-bank',type=Path,help='Reuse a production bank matching the coefficient kernel; all independent checks still run.')
    parser.add_argument('--output',type=Path)
    args=parser.parse_args(argv)
    if not 0 <= args.p <= 10:
        parser.error('This validation protocol is calibrated for p=0..10; the production engine has no fixed order ceiling.')
    start=time.monotonic()
    if args.reuse_bank:
        bank,generation=load_bank(args.reuse_bank)
        if bank.order!=args.p:raise ValueError('Bank order must equal requested p.')
    else:
        bank=generate_moments(args.p)
        generation={'metadata':bank.metadata(),'producer_kernel_signature':kernel_signature()}
    report={'status':'passed','python':platform.python_version(),
      'dependencies':{name:importlib.metadata.version(name) for name in ('sympy','numpy','scipy','mpmath')},
      'p_range':[0,args.p],'n_indices':[0,1,2],
      'maximum_W':args.p+3,'generation':generation,'checks':{}}
    report['checks']['universal_moments']=moment_checks(bank)
    with reuse_bank(bank):
        for name,check in [('algebra',algebra_checks),('gamma',gamma_checks),
                ('physical_taylor',taylor_checks),('physical_modes',physical_checks)]:
            print('Running',name,flush=True)
            report['checks'][name]=check(bank)
    root=Path(__file__).resolve().parents[2]
    files=[*sorted((root/'source').rglob('*.py')),*sorted((root/'validation').rglob('*.py'))]
    report['source_sha256']={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    report['seconds']=time.monotonic()-start
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(structured.dumps(report,indent=2)+'\n')
    print('All high-order checks passed.',flush=True)
    return report


if __name__=='__main__':main()
#========
