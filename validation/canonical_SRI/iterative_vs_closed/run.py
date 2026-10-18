# ========
# WG Edit
# Compare both engines before rounding on one common generated moment bank.
# Check each construction step and forward-substitute the closed relation.
# ========
from pathlib import Path
import argparse
import hashlib
from sri import serialization as structured
import time
import sympy as sp
from sri import spectra
from sri.spectra import mode_calculation
from sri.consistency import calculation
from sri.closed.algebra import ClosedBackground
from sri.closed.spectra import exponential


ROOT = Path(__file__).resolve().parents[3]


def spectra_check(p):
    iterative = mode_calculation(p,p+1,method='iterative')
    closed = mode_calculation(p,p+1,method='closed')
    ring = iterative.jet.ring
    assert iterative.h == closed.h
    assert iterative.sources == closed.sources
    assert iterative.log_amplitudes == closed.log_amplitudes
    bg = ClosedBackground(closed.jet)
    for field in ('scalar','tensor'):
        left = ring.exponential(iterative.log_amplitudes[field],p)
        right = exponential(ring,closed.log_amplitudes[field],p)
        assert left == right
        assert sp.expand(iterative.jet.to_hubble(left)-bg.to_hubble(right)) == 0
    a = spectra(p=p,method='iterative')
    b = spectra(p=p,method='closed')
    assert a.scalar == b.scalar and a.tensor == b.tensor
    return {'p':p,'h_residual':'0','source_residual':'0','log_amplitude_residual':'0',
            'amplitude_residual':'0','Hubble_flow_residual':'0','public_spectra_residual':'0'}


def relation_check(n,p):
    a = calculation(n,p,method='iterative')
    b = calculation(n,p,method='closed')
    assert a[1].forward == b[1].forward
    assert a[1].baseline == b[1].baseline
    assert a[1].shifted == b[1].shifted
    assert a[2] == b[2]
    assert a[3] == b[3]
    ring = b[1].modes.jet.ring
    residual = ring.add(ring.compose(b[3],b[1].forward,b[0].W),ring.scale(b[2],-1))
    assert not residual
    return {'n':n,'p':p,'W':n+1+p,'forward_map_residual':'0',
            'LO_inverse_residual':'0','shift_residual':'0','tensor_derivative_residual':'0',
            'observable_relation_residual':'0','forward_substitution_residual':'0'}


def run(*,max_p=4,max_n=3,max_lo_n=8):
    if any(not isinstance(v,int) or isinstance(v,bool) or v<0 for v in (max_p,max_n,max_lo_n)):
        raise ValueError('Validation limits must be nonnegative integers.')
    started = time.monotonic()
    powers = []; relations = []
    for p in range(max_p+1):
        powers.append(spectra_check(p))
        for n in range((max(max_n,max_lo_n) if p==0 else max_n)+1):
            relations.append(relation_check(n,p))
        print(f'p={p}: spectra and consistency comparisons passed.',flush=True)
    files = sorted((ROOT/'source').rglob('*.py'))
    return {'status':'passed','max_p':max_p,'max_n':max_n,'max_lo_n':max_lo_n,
            'spectra':powers,'consistency':relations,'seconds':time.monotonic()-started,
            'comparison':'exact rational polynomial equality before public rounding',
            'moment_backend':'shared normalized mode equations; not an independent test of moment generation',
            'source_sha256':{str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest() for f in files}}


def main():
    parser = argparse.ArgumentParser(description='Exact iterative/closed algebra comparisons.')
    parser.add_argument('--max-p',type=int,default=4)
    parser.add_argument('--max-n',type=int,default=3)
    parser.add_argument('--max-lo-n',type=int,default=8)
    parser.add_argument('--output',type=Path,default=Path(__file__).with_name('results.yaml'))
    args = parser.parse_args()
    result = run(max_p=args.max_p,max_n=args.max_n,max_lo_n=args.max_lo_n)
    args.output.write_text(structured.dumps(result,indent=2)+'\n')
    print(f"Passed: {len(result['spectra'])} spectrum cases and {len(result['consistency'])} relations.")


if __name__=='__main__':
    main()
