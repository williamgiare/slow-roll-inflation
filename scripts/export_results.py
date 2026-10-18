#========
# WG Edit
# Generate standalone equation documents independently of scientific validation.
# Reusing a bank is optional and checks the production-kernel fingerprint.
#========
import argparse
from contextlib import nullcontext
from sri import serialization as structured
from pathlib import Path
import shutil
import subprocess
import tempfile
import time
from sri import spectra, consistency
from sri.latex_output import export_document, order_name
from sri.terminal import term_count
from scripts.archive_results import archive_document

ROOT = Path(__file__).resolve().parents[1]


def compile_pdf(path):
    compiler = shutil.which('pdflatex')
    if compiler is None:
        raise RuntimeError('PDF export requires an existing pdflatex installation.')
    with tempfile.TemporaryDirectory(prefix='sri-latex-') as temporary:
        run = subprocess.run([compiler, '-interaction=nonstopmode', '-halt-on-error',
            '-output-directory', temporary, str(path.resolve())],
            cwd=temporary, capture_output=True, text=True, timeout=120)
        if run.returncode:
            raise RuntimeError(run.stdout[-4000:] + run.stderr[-1000:])
        log = (Path(temporary) / (path.stem + '.log')).read_text(errors='replace')
        if 'Overfull \\hbox' in log or 'Overfull \\vbox' in log:
            raise RuntimeError(f'{path.name}: an equation exceeds the page area.')
        target = path.with_suffix('.pdf')
        shutil.copyfile(Path(temporary) / target.name, target)
    return target


def generate(*, root=ROOT, max_p=10, max_n=3, high_order_from=7,
             high_order_max_n=0, bank_path=None, compile_documents=False, kinflation=False, min_p=0):
    if any(not isinstance(value, int) or isinstance(value, bool) or value < 0
           for value in (max_p, max_n, high_order_from, high_order_max_n, min_p)):
        raise ValueError('Orders must be nonnegative integers.')
    if min_p > max_p:
        raise ValueError('The minimum order cannot exceed the maximum order.')
    if bank_path is None:
        scope = nullcontext()
    else:
        # WG Edit: use the canonical validation bank loader.
        from validation.canonical_SRI.high_order import load_bank
        from validation.common.numerical.verified_bank import reuse_bank
        bank, _ = load_bank(bank_path)
        if bank.order < max_p:
            raise ValueError('The bank does not cover the requested spectral orders.')
        scope = reuse_bank(bank)
    rows = []
    # WG Edit: production and export share the optimized exact arithmetic.
    with scope:
        # WG Edit: generate either model and allow a selected order range.
        for p in range(min_p, max_p + 1):
            started = time.perf_counter()
            power = spectra(p=p, kinflation=kinflation)
            limit = max_n if p < high_order_from else high_order_max_n
            relations = []
            relation_checks = {}
            for n in range(limit + 1):
                t = time.perf_counter()
                relations.append(consistency(n=n, p=p, kinflation=kinflation))
                # WG Edit: verify each newly generated relation by exact forward substitution.
                from sri.consistency import calculation
                order,inverse,source,out=calculation(n,p,kinflation=kinflation)
                ring=inverse.modes.jet.ring
                residual=ring.add(ring.compose(out,inverse.forward,order.W),ring.scale(source,-1))
                if residual:raise ArithmeticError(f'Forward substitution failed at n={n}, p={p}.')
                relation_checks[str(n)]='0'
                print(f'{order_name(p)}: n={n} generated in {time.perf_counter()-t:.2f} s', flush=True)
            path = export_document(power, relations, root=root)
            yaml_path = archive_document(path, p, 'kinflation' if kinflation else 'canonical')
            pdf = compile_pdf(path) if compile_documents else None
            rows.append({'p': p, 'n': list(range(limit + 1)),
                'spectral_terms': {'scalar': term_count(power.scalar_correction),
                                  'tensor': term_count(power.tensor_correction)},
                'relation_terms': {str(value.n): term_count(value.expression) for value in relations},
                'exact_forward_residuals': relation_checks,
                'tex': str(path.relative_to(root)),
                'yaml': str(yaml_path.relative_to(root)),
                'pdf': str(pdf.relative_to(root)) if pdf else None,
                'seconds': time.perf_counter() - started})
            print(f'{order_name(p)}: document generated in {rows[-1]["seconds"]:.2f} s', flush=True)
            # Release background jets after each complete document.
            from sri.spectra import mode_calculation, _spectra
            from sri.consistency import calculation
            mode_calculation.cache_clear(); _spectra.cache_clear(); calculation.cache_clear()
            if kinflation:
                from sri.kinflation import mode_calculation as sound_modes, calculation as sound_calculation
                sound_modes.cache_clear(); sound_calculation.cache_clear()
            from sri.series import WeightedRing
            WeightedRing.clear_degree_cache()
    return {'status': 'passed', 'method': 'iterative',
            'model': 'kinflation' if kinflation else 'canonical', 'orders': rows,
            'moment_bank': 'kernel-verified generated bank' if bank_path else 'fresh mode generation'}


def main():
    parser = argparse.ArgumentParser(description='Export spectra and tensor consistency relations by SR order.')
    # WG Edit: select the physical model without changing canonical defaults.
    parser.add_argument('--kinflation', dest='kinflation', action='store_true', help='Export variable-sound-speed spectra and relations.')
    parser.add_argument('--min-p', type=int, default=0)
    parser.add_argument('--max-p', type=int, default=10)
    parser.add_argument('--max-n', type=int, default=3)
    parser.add_argument('--high-order-from', type=int, default=7)
    parser.add_argument('--high-order-max-n', type=int, default=0)
    parser.add_argument('--bank', type=Path, help='Reuse a kernel-matched production moment bank.')
    parser.add_argument('--compile', action='store_true', help='Also export PDF using existing pdflatex.')
    args = parser.parse_args()
    result = generate(max_p=args.max_p, max_n=args.max_n,
        high_order_from=args.high_order_from, high_order_max_n=args.high_order_max_n,
        bank_path=args.bank, compile_documents=args.compile,
        kinflation=args.kinflation, min_p=args.min_p)
    filename = 'results_export_kinflation.yaml' if args.kinflation else 'results_export.yaml'
    report = ROOT / 'validation/common/reports' / filename
    report.write_text(structured.dumps(result, indent=2) + '\n')


if __name__ == '__main__':
    main()
