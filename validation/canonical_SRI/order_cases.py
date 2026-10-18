# ========
# WG Edit
# Thin per-order entry points over independent shared validation machinery.
# Literature targets live exclusively in validation, never in source.
# ========
from pathlib import Path
import argparse, hashlib, time, platform
from sri import serialization as structured
import importlib.metadata
import numpy as np
import sympy as sp
from sri import spectra, consistency
from sri.consistency import calculation
from sri.spectra import mode_calculation
from sri.generic_modes import generate_moments
from sri.background import CONSTANT_VALUES
from .high_order import load_bank, kernel_signature, moment_checks, gamma_checks, taylor_checks, physical_checks
from ..common.numerical.verified_bank import reuse_bank
from ..common.symbolic.checks import coefficient_comparison
from ..common.symbolic.reference import reference_consistency
from ..common.symbolic.published_spectra import amplitudes, REFERENCE

ROOT=Path(__file__).resolve().parents[2]
ORDER_NAMES={0:'LO',1:'NLO',2:'N2LO',3:'N3LO',10:'N10LO'}

def fingerprint():
    # WG Edit: include the recovered LO validation sources in report fingerprints.
    paths=[*sorted((ROOT/'source').rglob('*.py')),*sorted((ROOT/'validation').rglob('*.py')),*sorted((ROOT/'validation/common/fixtures').glob('*.yaml'))]
    return {str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}

def save(report):
    report['source_sha256']=fingerprint()
    report['python']=platform.python_version()
    report['dependencies']={name:importlib.metadata.version(name) for name in ('numpy','scipy','sympy','mpmath','matplotlib')}
    target=ROOT/'validation/common/reports'/f"{ORDER_NAMES[report['p']].lower()}_case.yaml"
    target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(structured.dumps(report,indent=2)+'\n')
    return report

def saved_report(p):
    data=structured.loads((ROOT/'validation/common/reports'/f'{ORDER_NAMES[p].lower()}_case.yaml').read_text())
    from validation.common.report_compatibility import compatible_report
    current=fingerprint()
    if data['source_sha256']!=current and not compatible_report(ROOT,data['source_sha256'],current):
        raise RuntimeError('The saved report is stale. Rerun validation for the current source.')
    if data['status']!='passed': raise RuntimeError('The saved validation did not pass.')
    return data

def lo_bell(n):
    """Direct finite Bell sum from -log[1+(r/8) integral exp(-log Ps)]."""
    result=consistency(n=n,p=0)
    a=result.scalar_symbols
    inner=[sp.Integer(1)]+[sum(sp.bell(j,m,tuple(-a[i] for i in range(j-m+1))) for m in range(1,j+1)) for j in range(1,n+1)]
    expression=sum((-result.r_symbol/8)**m*sp.factorial(m-1)*sp.bell(n+1,m,tuple(inner[:n+2-m])) for m in range(1,n+2))
    return sp.expand(expression)

def relation_checks(n,p):
    order,inv,source,out=calculation(n,p)
    ring=inv.modes.jet.ring
    residual=ring.add(ring.compose(out,inv.forward,order.W),ring.scale(source,-1))
    assert residual=={},(n,p,'forward substitution failed')
    result=consistency(n=n,p=p)
    symbols=(result.r_symbol,*result.scalar_symbols)
    assert result.W==n+1+p
    assert not any(str(s).startswith(('epsilon_','q','J_')) or str(s)=='E' for s in result.expression.free_symbols)
    for weight,expr in result.components:
        assert all(ex[0]+sum((j+1)*v for j,v in enumerate(ex[1:]))==weight for ex,c in sp.Poly(expr,*symbols).terms() if c)
    assert result.expression.subs(result.r_symbol,0)==0
    # Exact power-law relation, expanded at the requested absolute weight.
    r=result.r_symbol
    nt=-r/(8*(1-r/16))
    value=result.expression.subs({result.scalar_symbols[0]:nt,**dict.fromkeys(result.scalar_symbols[1:],0)})
    remainder=sp.series(value-(nt if n==0 else 0),r,0,result.W+1).removeO().expand()
    limit=coefficient_comparison(remainder,sp.Integer(0),(r,))
    return {'n':n,'p':p,'W':order.W,'terms':len(out),'exact_forward_substitution_residual':'0',
            'homogeneous_weights':'passed','zero_r':'passed','constant_epsilon':limit}

# WG Edit: LO now validates Bell/Cortes--Liddle equivalence independently of native solvers.
def run_lo():
    from validation.canonical_SRI.LO.large_n_validation.run import symbolic_bridge, verify_archive
    start=time.monotonic()
    verify_archive()
    symbolic_bridge(8)
    rows=[]
    for n in range(9):
        actual=consistency(n=n,p=0)
        assert sp.expand(actual.expression-lo_bell(n))==0,n
        rows.append({'n':n,'W':n+1,'bell_vs_engine_residual':'0',
                     'archived_bell_vs_riccati_residual':'0'})
    return save({'status':'passed','p':0,'checks':{'bell_n0_to_n8':rows,
        'archived_symbolic_bridge':'passed','archive_integrity':'passed'},
        'seconds':time.monotonic()-start})

def run_low(p):
    start=time.monotonic()
    if p not in (1,2,3):raise ValueError('Use run_low for p=1,2,3.')
    bank=generate_moments(p)
    report={'status':'passed','p':p,'generation':bank.metadata(),'checks':{}}
    with reuse_bank(bank):
        expected=amplitudes(p);power=spectra(p=p);e=sp.symbols(f'epsilon_1:{p+3}')
        comparisons=[]
        for field in ('scalar','tensor'):
            comparisons.append({'field':field,'reference_equation':54 if field=='scalar' else 44,
                **coefficient_comparison(getattr(power,field+'_correction'),expected[field],e)})
        report['checks']['published_spectra']={'source':REFERENCE,'pivot':'-k_* eta_*=1','comparisons':comparisons}
        relations=[]
        for n in (0,1,2):
            row=relation_checks(n,p);actual=consistency(n=n,p=p)
            row['independent_observable_flow']=coefficient_comparison(actual.expression,reference_consistency(n,p),(actual.r_symbol,*actual.scalar_symbols))
            relations.append(row)
        report['checks']['consistency']=relations
        report['checks']['hankel_limit']=gamma_checks(bank)
        report['checks']['physical_taylor']=taylor_checks(bank)
    report['seconds']=time.monotonic()-start
    return save(report)

def run_high(*,bank_path=None):
    start=time.monotonic()
    if bank_path is None:
        bank=generate_moments(10);generation={'metadata':bank.metadata(),'producer_kernel_signature':kernel_signature()}
    else:
        bank,generation=load_bank(bank_path)
        if bank.order!=10:raise ValueError('The order-ten case requires a p=10 bank.')
    report={'status':'passed','p':10,'generation':generation,'checks':{}}
    report['checks']['independent_linear_modes']=moment_checks(bank)
    with reuse_bank(bank):
        report['checks']['hankel_limit']=gamma_checks(bank)
        report['checks']['physical_taylor']=taylor_checks(bank)
        report['checks']['physical_modes']=physical_checks(bank)
        report['checks']['consistency']=relation_checks(1,10)
        power=spectra(p=10);relation=consistency(n=1,p=10)
        # WG Edit: equation documents are exported separately under results/.

    report['seconds']=time.monotonic()-start
    plot_high(report);return save(report)

def figure_style():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family':'serif','mathtext.fontset':'stix','font.size':15,
      'axes.linewidth':1.15,'xtick.direction':'in','ytick.direction':'in',
      'xtick.top':True,'ytick.right':True,'xtick.major.size':6,'ytick.major.size':6,
      'xtick.minor.size':3,'ytick.minor.size':3,'savefig.dpi':200})
    return plt

def decorate(ax):
    ax.minorticks_on();ax.grid(True,which='major',ls='--',lw=.65,alpha=.32,color='0.45')
    ax.set_axisbelow(True)

def plot_high(report):
    plt=figure_style();rows=report['checks']['physical_modes']['runs']
    case=rows[0]
    fig,ax=plt.subplots(figsize=(6.4,4.6),layout='constrained')
    for field,col,ls in [('scalar','darkgreen','-'),('tensor','darkred','--')]:
        ax.semilogy(range(11),[case['relative_errors'][field+str(p)] for p in range(11)],color=col,ls=ls,lw=2,marker='o',ms=4.5)
    ax.set_xlabel(r'$p$');ax.set_ylabel(r'$|F_A^{(p)}/F_A^{\mathrm{num}}-1|$')
    ax.set_xticks(range(0,11,2));decorate(ax)
    directory=ROOT/'validation/canonical_SRI/N10LO/figures';directory.mkdir(exist_ok=True)
    fig.savefig(directory/'mode_convergence.png');plt.close(fig)

def main():
    parser=argparse.ArgumentParser(description='Run a documented per-order SRI validation case.')
    parser.add_argument('--p',type=int,choices=ORDER_NAMES,required=True)
    parser.add_argument('--bank',type=Path,help='Reuse a kernel-matched generated bank; independent checks still run.')
    args=parser.parse_args()
    if args.p==0:
        report=run_lo()
    elif args.p==10:report=run_high(bank_path=args.bank)
    else:report=run_low(args.p)
    print(ORDER_NAMES[args.p]+': all requested checks passed.',flush=True)
    return report

if __name__=='__main__':main()
# ========
