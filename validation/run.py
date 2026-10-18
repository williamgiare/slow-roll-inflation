#========
# WG Edit
# Run named validation suites with thresholds, reproducible metadata and failing exit codes.
#========
import argparse
import hashlib
import importlib.metadata
from sri import serialization as structured
import platform
from pathlib import Path
import time


def main(argv=None):
    parser=argparse.ArgumentParser(description='Run independent SRI validation suites.')
    parser.add_argument('--suite',choices=('all','symbolic','transport','modes','jost','generic'),default='all')
    parser.add_argument('--output',type=Path)
    args=parser.parse_args(argv)
    from .common.symbolic import checks
    from .common.numerical import transport,physical_modes,jost,generic_order
    suites={'symbolic':checks.run,'transport':transport.run,
            'modes':physical_modes.run,'jost':jost.run,'generic':generic_order.run}
    start=time.monotonic()
    root=Path(__file__).resolve().parent.parent
    report={'status':'passed','python':platform.python_version(),
        'dependencies':{name:importlib.metadata.version(name)
            for name in ('sympy','numpy','scipy','mpmath')},'suites':{}}
    for name,run in suites.items():
        if args.suite not in ('all',name):continue
        print(f'Running {name}',flush=True)
        report['suites'][name]=run()
    files=[*sorted((root/'source').rglob('*.py')),
           *sorted((root/'validation').rglob('*.py')),
           *sorted((root/'validation'/'common'/'fixtures').glob('*.yaml'))]
    report['source_sha256']={str(path.relative_to(root)):hashlib.sha256(path.read_bytes()).hexdigest()
                             for path in files}
    report['seconds']=time.monotonic()-start
    if args.output is not None:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(structured.dumps(report,indent=2)+'\n')
    print(f'Passed: {", ".join(report["suites"])} ({report["seconds"]:.1f} s)',flush=True)
    return 0


if __name__=='__main__':
    raise SystemExit(main())
