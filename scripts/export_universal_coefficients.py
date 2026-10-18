"""Export verified universal coefficients as exact, independently readable YAML."""
import argparse
from collections import defaultdict
from functools import lru_cache
from hashlib import sha256
from sri import serialization as structured
from math import factorial
from pathlib import Path
import pickle

import sympy as sp
import yaml
from sri.analytic import C, fast
from sri.analytic.generic_hyperlogs import GPL
from utils._expressions import ExactReader, read_family
from sri.generic_modes import source_patterns


class Formula(str):
    pass


class Dumper(getattr(yaml, 'CSafeDumper', yaml.SafeDumper)):
    def ignore_aliases(self, value):
        return True


Dumper.add_representer(Formula, lambda dumper, value: dumper.represent_scalar('tag:yaml.org,2002:str', str(value), style='|'))
Loader = getattr(yaml, 'CSafeLoader', yaml.SafeLoader)


def write_yaml(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w') as handle:
        yaml.dump(data, handle, Dumper=Dumper, sort_keys=False, allow_unicode=True, width=100)


class ValidatedExpressions(pickle.Unpickler):
    def find_class(self, module, name):
        if module == 'generic_hyperlogs' and name == 'GPL':
            return GPL
        return super().find_class(module, name)


def printer(names):
    @lru_cache(maxsize=16384)
    def render(value):
        if value in names:
            return names[value]
        if value.is_Rational:
            return str(value.p) if value.q == 1 else f'({value.p}/{value.q})'
        if value.is_Symbol:
            return str(value)
        if value == sp.pi:
            return 'pi'
        if value == sp.EulerGamma:
            return 'EulerGamma'
        if value.is_Add:
            return '('+'+'.join(render(a) for a in value.args)+')'
        if value.is_Mul:
            return '('+'*'.join(render(a) for a in value.args)+')'
        if value.is_Pow:
            return f'({render(value.base)})**({render(value.exp)})'
        if value.is_Function:
            return value.func.__name__+'('+','.join(render(a) for a in value.args)+')'
        raise TypeError(f'Unsupported exact constant: {value!r}')
    return render


def expression_text(expression, render):
    lines = []
    for term in sp.Add.make_args(expression):
        q, rest = term.as_coeff_Mul()
        assert q.is_Rational
        scalar = str(q.p) if q.q == 1 else f'{q.p}/{q.q}'
        factors = [] if rest == 1 else [render(f) for f in sp.Mul.make_args(rest)]
        lines.append(('' if not lines else '+ ')+' * '.join([scalar, *factors]))
    return Formula('\n'.join(lines)+'\n')




def export(*, reports, checkpoints, reference, output):
    report = structured.loads((reports/'p10.yaml').read_text())
    assert report['status'] == 'passed' and report['patterns'] == 138
    checks = {tuple(row['pattern']): row for row in report['checks']}
    assert set(checks) == set(source_patterns(10))
    reference_data = structured.loads(reference.read_text())
    numerical = {tuple(key): value for key, value in reference_data['values']}
    groups = defaultdict(list)
    for key in checks:
        groups[len(key), sum(key)].append(key)
    index = {'schema_version': 1, 'maximum_weight': 10, 'coefficient_count': 138,
             'mixed_coefficient_count': 67, 'numerical_summary': 'values.yaml', 'families': []}
    summary = {'schema_version': 1, 'coefficients': []}
    audited = []
    for (rank, degree), keys in sorted(groups.items(), key=lambda item: (sum(item[0]), item[0][0])):
        keys.sort()
        weight = rank+degree
        if any(not fast.supported(key) for key in keys):
            name = f'm{rank:02d}_j{degree:02d}'
            saved = structured.loads((reports/'families'/f'{name}.yaml').read_text())
            assert saved['status'] == 'passed' and saved['signature'] == report['signature']
            path = checkpoints/f'{name}.pkl'
            if not path.exists(): raise FileNotFoundError(f'Missing verified checkpoint {path}; no regeneration was requested.')
            with path.open('rb') as handle:
                expressions = ValidatedExpressions(handle).load()
        else:
            expressions = {key: fast.moment(key) for key in keys}
        assert set(expressions) == set(keys)
        atoms = set().union(*(value.atoms(GPL) for value in expressions.values()))
        ordered = sorted(atoms, key=lambda value: (len(value.args[0]), tuple((int(a.p),int(a.q)) for a in value.args[0])))
        names = {atom: f'G{i:06d}' for i, atom in enumerate(ordered,1)}
        words = {names[atom]: [str(a) for a in atom.args[0]] for atom in ordered}
        render = printer(names)
        data = {'schema_version': 1, 'weight': weight, 'rank': rank, 'log_degree': degree,
                'gpl_words': words, 'coefficients': []}
        relative = f'expressions/weight_{weight:02d}/m{rank:02d}_j{degree:02d}.yaml'
        for key in keys:
            check = checks[key]
            assert check['status'] == 'passed' and check['absolute_error'] <= check['comparison_bound']
            expression = expressions[key]
            assert not expression.has(sp.Float,sp.Integral,sp.Derivative)
            assert expression.free_symbols <= {C}
            row = {'pattern': list(key), 'permutation_multiplicity': factorial(rank)//sp.prod(factorial(key.count(j)) for j in set(key)),
                   'analytic_value': check['analytic_value'], 'mode_equation_value': repr(numerical[key])}
            row['permutation_multiplicity'] = int(row['permutation_multiplicity'])
            data['coefficients'].append({**row,'exact_expression':expression_text(expression,render)})
            summary['coefficients'].append({'weight':weight,**row,'expression_file':relative})
        path = output/relative
        write_yaml(path,data)
        # Read the actual YAML and reconstruct all exact expressions before accepting a family.
        decoded = read_family(path)
        assert decoded == expressions, relative
        index['families'].append({'weight':weight,'rank':rank,'log_degree':degree,'coefficient_count':len(keys),'file':relative})
        audited.append({'file':relative,'coefficients':len(keys),'sha256':sha256(path.read_bytes()).hexdigest(),'exact_roundtrip':'passed'})
        print(f'Archived weight {weight}, family ({rank},{degree}): {len(keys)} coefficients.',flush=True)
        del expressions, decoded, ordered, names, words, atoms, data, render
    write_yaml(output/'values.yaml',summary)
    write_yaml(output/'index.yaml',index)
    assert sum(row['coefficients'] for row in audited) == 138
    return {'status':'passed','coefficients':138,'families':audited,'validation_signature':report['signature']}


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reports',type=Path,default=root/'validation/symbolic_coefficients/reports')
    parser.add_argument('--checkpoints',type=Path,required=True,help='Existing verified family checkpoints; mixed coefficients are never regenerated.')
    parser.add_argument('--reference',type=Path,default=root/'validation/common/fixtures/p10_moments.yaml')
    parser.add_argument('--output',type=Path,default=root/'results/universal_coefficients')
    parser.add_argument('--audit',type=Path,help='Write serialization verification outside the scientific results directory.')
    args = parser.parse_args()
    result = export(reports=args.reports,checkpoints=args.checkpoints,reference=args.reference,output=args.output)
    if args.audit:
        args.audit.parent.mkdir(parents=True,exist_ok=True)
        write_yaml(args.audit,result)


if __name__ == '__main__':
    main()
