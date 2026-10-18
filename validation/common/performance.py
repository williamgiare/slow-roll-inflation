#========
# WG Edit
# Exact arithmetic comparisons and reproducible performance measurements.
# Reference calculations use the unoptimized arithmetic, not literature inputs.
#========
import argparse
import ast
import copy
import hashlib
from sri import serialization as structured
from pathlib import Path
import platform
import time
import numpy as np
from sri.consistency import calculation
from sri.spectra import spectra
from sri.generic_modes import NormalizedModeSystem
from .symbolic import reference_arithmetic as reference
from .numerical.verified_bank import reuse_bank
from ..canonical_SRI.high_order import kernel_signature, load_bank

ROOT=Path(__file__).resolve().parents[2]


def reference_kernel_signature():
    """Recover the numerical kernel signature with the reference arithmetic."""
    source=ROOT/'source/series.py'
    tree=ast.parse(source.read_text())
    prior=ast.parse(Path(reference.__file__).read_text())
    methods={'constant','variable','degree','cut','add','scale','mul','derivative'}
    nodes={node.name.removeprefix('reference_'):node for node in prior.body
           if isinstance(node,ast.FunctionDef) and node.name.startswith('reference_')}
    arithmetic=[]
    for cls in tree.body:
        if not isinstance(cls,ast.ClassDef):continue
        for node in cls.body:
            if isinstance(node,ast.FunctionDef) and node.name in methods:
                value=copy.deepcopy(nodes.get(node.name,node));value.name=node.name
                arithmetic.append(ast.dump(value,include_attributes=False))
    generator=(ROOT/'source/generic_modes.py').read_bytes()
    return hashlib.sha256(generator+'\n'.join(arithmetic).encode()).hexdigest()


def mode_input_checks():
    """Compare every exact and numerical ODE input at all orders through ten."""
    rows=[]
    for p in range(1,11):
        for depth in (12,18,20):
            with reference.reference_arithmetic():before=NormalizedModeSystem(p,depth)
            after=NormalizedModeSystem(p,depth)
            for key in ('exact_s','exact_pump','s','pump','out','left','right','multiplicities'):
                assert np.array_equal(getattr(before,key),getattr(after,key)),(p,depth,key)
            for key in ('exact_w','exact_u','w_initial','u_initial'):
                a,b=getattr(before,key),getattr(after,key)
                assert len(a)==len(b)
                assert all(np.array_equal(x,y) for x,y in zip(a,b)),(p,depth,key)
            for key in ('patterns','basis','exact_multiplicities','size'):
                assert getattr(before,key)==getattr(after,key),(p,depth,key)
            rows.append({'p':p,'asymptotic_terms':depth,'exact_and_numerical_inputs':'identical'})
        print(f'Mode input comparison p={p}: passed.',flush=True)
    return rows


def certify_bank(path,checks):
    """Retain a generated bank only after proving its mode inputs unchanged."""
    data=structured.loads(path.read_text());previous=data['producer_kernel_signature']
    if previous==kernel_signature():return
    if previous!=reference_kernel_signature():
        raise ValueError('The bank does not match the independent reference kernel.')
    assert len(checks)==30 and all(row['exact_and_numerical_inputs']=='identical' for row in checks)
    before=data['values'],data['differences'],data['metadata']
    data['producer_kernel_signature']=kernel_signature()
    data.setdefault('exact_arithmetic_verifications',[]).append({
        'previous_kernel_signature':previous,'kernel_signature':kernel_signature(),
        'verification':'exact and numerical pump, initial data, product basis and multiplicities equal through p=10 at asymptotic depths 12,18,20',
        'report':'validation/common/reports/arithmetic_performance.yaml'})
    assert before==(data['values'],data['differences'],data['metadata'])
    path.write_text(structured.dumps(data,indent=2)+'\n')


def relation_benchmark(n,p,kinflation):
    with reference.reference_arithmetic():
        start=time.perf_counter();a=calculation(n,p,kinflation=kinflation);old=time.perf_counter()-start
    reference.clear_calculation_caches()
    start=time.perf_counter();b=calculation(n,p,kinflation=kinflation);new=time.perf_counter()-start
    assert a[2:]==b[2:]
    for field in ('forward','baseline','shifted'):assert getattr(a[1],field)==getattr(b[1],field)
    for field in ('h','sources','log_amplitudes'):assert getattr(a[1].modes,field)==getattr(b[1].modes,field)
    ring=b[1].modes.jet.ring
    assert not ring.add(ring.compose(b[3],b[1].forward,b[0].W),ring.scale(b[2],-1))
    result={'model':'kinflation' if kinflation else 'canonical','n':n,'p':p,'W':b[0].W,
            'reference_seconds':old,'optimized_seconds':new,'speedup':old/new,
            'terms':len(b[3]),'exact_result_residual':'0','exact_forward_residual':'0'}
    print(f"{result['model']} n={n}, p={p}: {old:.3f} -> {new:.3f} s ({old/new:.2f}x).",flush=True)
    reference.clear_calculation_caches()
    return result


def run(bank_path,*,benchmarks=True):
    start=time.perf_counter();checks=mode_input_checks()
    certify_bank(bank_path,checks)
    # This fixture contains generated mode moments; its numbers are never seeded in production.
    fixture=ROOT/'validation/common/fixtures/p10_moments.yaml'
    if fixture.exists():certify_bank(fixture,checks)
    rows=[]
    if benchmarks:
        bank,_=load_bank(bank_path)
        with reuse_bank(bank):
            for n,p,sound in ((1,7,False),(0,6,True),(1,5,True)):
                rows.append(relation_benchmark(n,p,sound))
    files=[*sorted((ROOT/'source').rglob('*.py')),Path(__file__),Path(reference.__file__)]
    return {'status':'passed','comparison':'exact rational polynomials before rounding',
            'mode_inputs':checks,'benchmarks':rows,'python':platform.python_version(),
            'timing_scope':'cold background algebra using one common generated moment bank; numerical integration excluded',
            'seconds':time.perf_counter()-start,
            'source_sha256':{str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in files}}


def main():
    parser=argparse.ArgumentParser(description='Verify optimized arithmetic and measure both inflation models.')
    parser.add_argument('--bank',type=Path,default=ROOT/'validation/common/reports/order10_fresh_bank.yaml')
    parser.add_argument('--certify-only',action='store_true')
    parser.add_argument('--output',type=Path,default=ROOT/'validation/common/reports/arithmetic_performance.yaml')
    args=parser.parse_args();report=run(args.bank,benchmarks=not args.certify_only)
    args.output.write_text(structured.dumps(report,indent=2)+'\n')
    print('Exact arithmetic validation passed.',flush=True)


if __name__=='__main__':main()
