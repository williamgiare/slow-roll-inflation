#========
# WG Edit
# Sequential, checkpointed comparison of exact and numerical universal moments.
# Each mixed family is generated independently before the ODE fixture is read.
#========
from __future__ import annotations

import argparse
from hashlib import sha256
from sri import serialization as structured
from pathlib import Path
import pickle
import resource
import sqlite3
import shutil
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parent
BACKEND=ROOT/'backend'
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
sys.path.insert(0,str(BACKEND))

import sympy as sp
import mpmath as mp
import sri
from sri.generic_modes import source_patterns,generate_moments

import coefficients as fast
import generic_hyperlogs as exact

_GPL_DATABASE=None


def cached_evaluate(expression,digits=24,progress=None):
    """Cache standard-function evaluations, never physical source coefficients."""
    global _GPL_DATABASE
    atoms=sorted(expression.atoms(exact.GPL),key=lambda a:(len(a.args[0]),str(a)))
    # WG Edit: expose evaluation progress without changing coefficient arithmetic.
    if progress is not None: progress(0,len(atoms))
    if not atoms: return exact.evaluate(expression,digits),{'new':0,'reused':0}
    if _GPL_DATABASE is None:
        path=ROOT/'cache/gpl_values.sqlite'
        path.parent.mkdir(parents=True,exist_ok=True)
        _GPL_DATABASE=sqlite3.connect(path)
        _GPL_DATABASE.execute('CREATE TABLE IF NOT EXISTS values_cache (signature TEXT, digits INTEGER, word TEXT, value TEXT, PRIMARY KEY(signature,digits,word))')
    fingerprint=signature()
    working=digits+15
    replacements={};statistics={'new':0,'reused':0}
    for index,atom in enumerate(atoms,1):
        word=tuple(exact.Q(int(a.p),int(a.q)) for a in atom.args[0])
        key=','.join(map(str,word))
        record=_GPL_DATABASE.execute('SELECT value FROM values_cache WHERE signature=? AND digits=? AND word=?',(fingerprint,working,key)).fetchone()
        if record is None:
            number=mp.nstr(exact.gpl_value(word,working),working)
            _GPL_DATABASE.execute('INSERT INTO values_cache VALUES (?,?,?,?)',(fingerprint,working,key,number))
            statistics['new']+=1
            if statistics['new']%64==0: _GPL_DATABASE.commit()
        else:
            number=record[0];statistics['reused']+=1
        replacements[atom]=sp.Float(number,working)
        if progress is not None and (index%128==0 or index==len(atoms)):
            progress(index,len(atoms))
    _GPL_DATABASE.commit()
    return expression.xreplace(replacements).subs(fast.C,fast.C_VALUE).evalf(digits),statistics


def atomic_yaml(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(structured.dumps(value,indent=2)+'\n')
    tmp.replace(path)


def fixture_path():
    checkout=Path(__file__).resolve().parents[1]/'common/fixtures/p10_moments.yaml'
    if checkout.exists(): return checkout
    return Path(sri.__file__).resolve().parent.parent/'validation/common/fixtures/p10_moments.yaml'


def signature():
    digest=sha256()
    for file in sorted(BACKEND.glob('*.py')):
        digest.update(file.name.encode());digest.update(file.read_bytes())
    digest.update(fixture_path().read_bytes())
    from validation.common.report_compatibility import preserved_coefficient_signature
    return preserved_coefficient_signature(Path(__file__).resolve().parents[2],digest.hexdigest())


def numerical_reference():
    data=structured.loads(fixture_path().read_text())
    return ({tuple(k):v for k,v in data['values']},
            {tuple(k):v for k,v in data['differences']},data['metadata'])


def checked_expression(pattern,expression,numbers,uncertainties,progress=None):
    if expression.has(sp.Float,sp.Derivative,sp.Integral):
        raise AssertionError(f'Nonexact or unevaluated expression for {pattern}.')
    if expression.free_symbols-{fast.C}:
        raise AssertionError(f'Undetermined symbols for {pattern}: {expression.free_symbols}')
    start=time.perf_counter()
    value,cache_statistics=cached_evaluate(expression,24,progress=progress)
    if value.is_finite is not True or value.is_real is not True:
        raise ArithmeticError(f'Nonfinite or complex coefficient for {pattern}: {value}')
    error=abs(value-sp.Float(numbers[pattern],24))
    difference=abs(uncertainties[pattern])
    bound=max(5e-12,10*difference)
    if error>sp.Float(bound,24):
        raise AssertionError(f'{pattern}: error {error} exceeds numerical uncertainty bound {bound}.')
    return {'pattern':list(pattern),'analytic_value':str(value),
            'absolute_error':float(error),'numerical_run_difference':difference,
            'comparison_bound':bound,'GPL_constants':len(expression.atoms(exact.GPL)),
            'GPL_evaluations':cache_statistics,
            'evaluation_digits':24,'evaluation_seconds':time.perf_counter()-start,
            'status':'passed'}


def family_worker(rank,degree):
    label=f'm{rank:02d}_j{degree:02d}'
    report_path=ROOT/'reports/families'/f'{label}.yaml'
    cache_path=ROOT/'cache'/f'{label}.pkl'
    if cache_path.exists() and not report_path.exists():
        raise ValueError('Checkpoint has no readable validation metadata; choose a separate output directory.')
    fingerprint=signature()
    previous=structured.loads(report_path.read_text()) if report_path.exists() else {}
    if previous.get('signature')!=fingerprint:
        previous={}
        cache_path.unlink(missing_ok=True)
    resuming=bool(previous) and cache_path.exists()
    report={'rank':rank,'log_degree':degree,'SR_weight':rank+degree,
            'signature':fingerprint,'status':'running','checks':previous.get('checks',[]) if resuming else [],
            'generation_seconds':previous.get('generation_seconds',0),
            'directions':previous.get('directions',[])}
    atomic_yaml(report_path,report)
    start=time.perf_counter()
    if previous and cache_path.exists():
        # Read only checkpoints produced by this test suite for this fingerprint.
        with cache_path.open('rb') as handle: expressions=pickle.load(handle)
    else:
        original=exact.directional_generator
        def record_direction(slopes,maximum=2):
            begin=time.perf_counter()
            # WG Edit: report the active direction before a potentially long calculation.
            report['phase']='generation'
            report['active_direction']=list(slopes)
            report['direction_started_unix']=time.time()
            atomic_yaml(report_path,report)
            direction_path=cache_path.with_suffix('')/('-'.join(map(str,slopes))+'.pkl')
            direction_meta=direction_path.with_suffix('.yaml')
            saved=structured.loads(direction_meta.read_text()) if direction_meta.exists() else {}
            if saved.get('signature')==fingerprint and direction_path.exists():
                with direction_path.open('rb') as handle: jet=pickle.load(handle)
                elapsed=saved['seconds']
            else:
                jet=original(slopes,maximum)
                elapsed=time.perf_counter()-begin
                direction_path.parent.mkdir(parents=True,exist_ok=True)
                temp=direction_path.with_suffix('.tmp')
                with temp.open('wb') as handle: pickle.dump(jet,handle)
                temp.replace(direction_path)
                atomic_yaml(direction_meta,{'signature':fingerprint,'seconds':elapsed})
            poles={n:sp.expand(jet[n]) for n in range(-rank,0)}
            if any(v!=0 for v in poles.values()):
                raise AssertionError(f'Uncancelled pole in direction {slopes}.')
            report['directions']=[d for d in report['directions'] if d['slopes']!=list(slopes)]
            report['directions'].append({'slopes':list(slopes),'negative_powers':list(poles),
                                         'pole_residuals':['0']*len(poles),'seconds':elapsed})
            report['phase']='generation'
            report['active_direction']=None
            report['direction_started_unix']=None
            rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
            report['peak_RSS_bytes']=rss if sys.platform=='darwin' else rss*1024
            atomic_yaml(report_path,report)
            return jet
        exact.directional_generator=record_direction
        expressions=exact.family(rank,degree)
        report['generation_seconds']=time.perf_counter()-start
        report['phase']='generated'
        cache_path.parent.mkdir(parents=True,exist_ok=True)
        temp=cache_path.with_suffix('.tmp')
        with temp.open('wb') as handle: pickle.dump(expressions,handle)
        temp.replace(cache_path)
        atomic_yaml(report_path,report)
    # The numerical data do not enter symbolic generation.
    numbers,uncertainties,metadata=numerical_reference()
    completed={tuple(c['pattern']) for c in report['checks']}
    report['phase']='evaluation'
    for pattern,expression in expressions.items():
        if pattern in completed: continue
        # WG Edit: checkpoint progress through standard-function evaluation.
        report['active_pattern']=list(pattern)
        report['evaluation_started_unix']=time.time()
        report['evaluation_progress']={}
        atomic_yaml(report_path,report)
        def record_evaluation(done,total):
            report['evaluation_progress']={'evaluated_GPL_constants':done,
                                           'GPL_constants_total':total}
            atomic_yaml(report_path,report)
        result=checked_expression(pattern,expression,numbers,uncertainties,
                                  progress=record_evaluation)
        report['checks'].append(result)
        rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        report['peak_RSS_bytes']=rss if sys.platform=='darwin' else rss*1024
        atomic_yaml(report_path,report)
        print(f'Checked {pattern}: {result["analytic_value"]}',flush=True)
    report['status']='passed';report['phase']='complete'
    report['active_pattern']=None
    report['evaluation_started_unix']=None
    report['patterns']=len(expressions)
    report['evaluation_seconds']=sum(c['evaluation_seconds'] for c in report['checks'])
    report['numerical_backend']=metadata['backend']
    report['numerical_integration_digits']=metadata.get('integration_precision_digits')
    atomic_yaml(report_path,report)
    print(f'Family ({rank},{degree}) passed.',flush=True)


def read_family(rank,degree,fingerprint):
    path=ROOT/'reports/families'/f'm{rank:02d}_j{degree:02d}.yaml'
    cache=ROOT/'cache'/f'm{rank:02d}_j{degree:02d}.pkl'
    if not path.exists() or not cache.exists(): return None
    report=structured.loads(path.read_text())
    return report if report.get('signature')==fingerprint and report.get('status')=='passed' else None


def run_order(order,family_limit,order_limit):
    start=time.perf_counter()
    fingerprint=signature()
    keys=list(source_patterns(order))
    groups=sorted({(len(k),sum(k)) for k in keys if not fast.supported(k)})
    reports=[]
    for rank,degree in groups:
        result=read_family(rank,degree,fingerprint)
        if result is None:
            # WG Edit: zero disables a budget for uninterrupted higher-order runs.
            remaining=order_limit-(time.perf_counter()-start) if order_limit else None
            if remaining is not None and remaining<=0:
                raise TimeoutError(f'Order {order} reached its time budget.')
            timeout=family_limit or None
            if remaining is not None:
                timeout=remaining if timeout is None else min(timeout,remaining)
            print(f'p={order}: generating family ({rank},{degree}).',flush=True)
            try:
                subprocess.run([sys.executable,'-u',str(Path(__file__).resolve()),
                                '--worker',str(rank),str(degree),'--output-root',str(ROOT)],
                               check=True,timeout=timeout)
            except subprocess.TimeoutExpired as error:
                progress=ROOT/'reports/families'/f'm{rank:02d}_j{degree:02d}.yaml'
                if progress.exists():
                    partial=structured.loads(progress.read_text())
                    partial['status']='stopped'
                    partial['stop_reason']='execution time budget exceeded'
                    atomic_yaml(progress,partial)
                raise TimeoutError(f'Family ({rank},{degree}) exceeded {timeout:.0f} seconds.') from error
            result=read_family(rank,degree,fingerprint)
            if result is None: raise AssertionError(f'Missing completed family ({rank},{degree}).')
        reports.append(result)
    # Generate every fast-family expression before loading the numerical data.
    begin=time.perf_counter()
    fast_expressions={k:fast.moment(k) for k in keys if fast.supported(k)}
    fast_seconds=time.perf_counter()-begin
    numbers,uncertainties,metadata=numerical_reference()
    checks=[checked_expression(k,value,numbers,uncertainties) for k,value in fast_expressions.items()]
    mixed={tuple(c['pattern']):c for f in reports for c in f['checks']}
    checks.extend(mixed[k] for k in keys if not fast.supported(k))
    if {tuple(c['pattern']) for c in checks}!=set(keys):
        raise AssertionError(f'Incomplete coefficient coverage at p={order}.')
    report={'status':'passed','p':order,'signature':fingerprint,'patterns':len(keys),
            'mixed_patterns':sum(not fast.supported(k) for k in keys),
            'family_coefficient_patterns':len(mixed),'checks':checks,'family_reports':[
                {'rank':f['rank'],'degree':f['log_degree'],'generation_seconds':f['generation_seconds'],
                 'evaluation_seconds':f['evaluation_seconds']} for f in reports],
            'fast_generation_seconds':fast_seconds,'this_run_seconds':time.perf_counter()-start,
            'cumulative_family_seconds':sum(f['generation_seconds']+f['evaluation_seconds'] for f in reports),
            'numerical_backend':metadata,'environment':{
                'python':sys.version.split()[0],'sympy':sp.__version__,'mpmath':mp.__version__}}
    atomic_yaml(ROOT/'reports'/f'p{order:02d}.yaml',report)
    print(f'p={order}: PASSED, {len(keys)} patterns; saved p{order:02d}.yaml.',flush=True)
    return report


def main():
    global ROOT
    parser=argparse.ArgumentParser()
    parser.add_argument('--output-root',type=Path,default=ROOT/'reproduced',
                        help='Directory for reproduced reports and checkpoints; published reports are preserved.')
    parser.add_argument('--orders',nargs='+',type=int,default=[8,9,10])
    # WG Edit: allow the operator to continue through checkpoints without a timeout.
    parser.add_argument('--family-seconds',type=float,default=900,
                        help='Per-family execution budget; 0 disables the limit.')
    parser.add_argument('--order-seconds',type=float,default=3600,
                        help='Per-order execution budget; 0 disables the limit.')
    parser.add_argument('--worker',nargs=2,type=int,help=argparse.SUPPRESS)
    parser.add_argument('--fresh',action='store_true',help='Regenerate all local checkpoints and reports.')
    args=parser.parse_args()
    ROOT=args.output_root.resolve()
    if ROOT==Path(__file__).resolve().parent:
        parser.error('Choose a separate reproduction directory; published reports are preserved.')
    if args.family_seconds<0 or args.order_seconds<0:
        parser.error('Execution budgets must be nonnegative.')
    if args.worker:
        family_worker(*args.worker);return
    if args.fresh:
        shutil.rmtree(ROOT/'cache',ignore_errors=True)
        shutil.rmtree(ROOT/'reports',ignore_errors=True)
    if any(not 1<=p<=10 for p in args.orders): parser.error('The independent fixture covers p=1,...,10.')
    summary={'status':'running','requested_orders':args.orders,'completed_orders':[]}
    atomic_yaml(ROOT/'reports/sequence.yaml',summary)
    try:
        for order in args.orders:
            result=run_order(order,args.family_seconds,args.order_seconds)
            summary['completed_orders'].append({'p':order,'patterns':result['patterns'],
                                                'this_run_seconds':result['this_run_seconds']})
            atomic_yaml(ROOT/'reports/sequence.yaml',summary)
    except (TimeoutError,subprocess.CalledProcessError) as error:
        summary['status']='stopped';summary['reason']=str(error)
        atomic_yaml(ROOT/'reports/sequence.yaml',summary)
        raise SystemExit(str(error))
    summary['status']='passed';atomic_yaml(ROOT/'reports/sequence.yaml',summary)


if __name__=='__main__': main()
