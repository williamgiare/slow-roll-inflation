"""Read published YAML results without generating or numerically integrating them."""
import ast
from copy import deepcopy
from dataclasses import dataclass
from decimal import Decimal
from functools import lru_cache
from pathlib import Path
import sympy as sp
import yaml

from ._expressions import ExactReader

RESULTS = Path(__file__).resolve().parents[1] / 'results'
_MODELS = {'canonical': 'canonical_SRI', 'kinflation': 'kinflation_SRI'}
_LOADER = getattr(yaml, 'CSafeLoader', yaml.SafeLoader)


def _integer(value, name):
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError(f'{name} must be a nonnegative integer.')
    return value


def _pattern(powers):
    if len(powers) == 1 and isinstance(powers[0], (tuple, list)):
        powers = tuple(powers[0])
    if not powers:
        raise ValueError('Provide a nonempty coefficient pattern.')
    return tuple(sorted(_integer(j, 'Logarithmic power') for j in powers))


def _root(root):
    return Path(root).expanduser().resolve() if root is not None else RESULTS


@lru_cache(maxsize=4)
def _load(path):
    if not path.is_file():
        raise FileNotFoundError(f'No archived YAML result at {path}. No calculation was started.')
    data = yaml.load(path.read_text(), Loader=_LOADER)
    if not isinstance(data, dict) or data.get('schema_version') != 1:
        raise ValueError(f'Unsupported result archive schema in {path}.')
    return data


def _coefficient_row(pattern, root):
    directory = root / 'universal_coefficients'
    for row in _load(directory / 'values.yaml')['coefficients']:
        if tuple(row['pattern']) == pattern:
            return row, directory
    raise LookupError(f'Coefficient {pattern} is not in the published archive. No generation was started.')


@lru_cache(maxsize=8)
def _coefficient(pattern, root):
    row, directory = _coefficient_row(pattern, root)
    path = (directory / row['expression_file']).resolve()
    if not path.is_relative_to(directory.resolve()):
        raise ValueError('Coefficient file must belong to its archive.')
    family = _load(path)
    for item in family['coefficients']:
        if tuple(item['pattern']) == pattern:
            return ExactReader(family.get('gpl_words', {})).expression(item['exact_expression'])
    raise ValueError(f'Coefficient {pattern} is missing from its archived family.')


def coefficient(*powers, root=None):
    """Return the exact SymPy coefficient; accepts (0,1,1) or a tuple/list."""
    return _coefficient(_pattern(powers), _root(root))


def value(*powers, source='analytic', root=None):
    """Return a stored Decimal value, without reevaluating special functions.

    source='analytic' returns the saved analytic evaluation;
    source='mode_equation' returns the independent numerical integration.
    """
    if source not in ('analytic', 'mode_equation'):
        raise ValueError("source must be 'analytic' or 'mode_equation'.")
    row, _ = _coefficient_row(_pattern(powers), _root(root))
    return Decimal(row[source + '_value'])


class _FormulaReader(ExactReader):
    def __init__(self, symbols):
        super().__init__({})
        self.names.update(symbols)

    def node(self, node):
        if isinstance(node, ast.Constant) and isinstance(node.value, float):
            return sp.Float(str(node.value), 17)
        return super().node(node)

    def factor(self, text):
        return self.node(ast.parse(text, mode='eval').body)

    def formula(self, text):
        # Each line is one signed monomial. Parse separately to bound AST depth.
        return sp.Add(*(self.factor(line) for line in text.splitlines() if line.strip()))


def _result_file(p, type, root):
    _integer(p, 'p')
    if not isinstance(type, str) or type not in _MODELS:
        raise ValueError("type must be 'canonical' or 'kinflation'.")
    name = 'LO' if p == 0 else f'N{p:02d}LO'
    path = root / _MODELS[type] / name / f'{name}.yaml'
    data = _load(path)
    if data.get('p') != p or data.get('type') != type:
        raise ValueError(f'Order or physical model does not match {path}.')
    return path, data


def _reader(data):
    symbols = {}
    for alias, name in data['symbols'].items():
        if not isinstance(alias, str) or not alias.isidentifier() or not isinstance(name, str):
            raise ValueError('Invalid symbol declaration in result archive.')
        symbols[alias] = sp.Symbol(name)
    return _FormulaReader(symbols)


@dataclass(frozen=True)
class ArchivedSpectra:
    """Stored scalar/tensor expressions with their physical order and metadata."""
    scalar: sp.Expr
    tensor: sp.Expr
    scalar_correction: sp.Expr
    tensor_correction: sp.Expr
    p: int
    type: str
    metadata: dict


@dataclass(frozen=True)
class ArchivedConsistency:
    """Stored tensor consistency relation and its absolute slow-roll cutoff."""
    expression: sp.Expr
    n: int
    p: int
    W: int
    type: str
    components: tuple
    metadata: dict


@lru_cache(maxsize=8)
def _spectra(p, type, root):
    _, data = _result_file(p, type, root)
    reader = _reader(data)
    power = data['spectra']
    corrections = {field: reader.formula(power[field]['correction']) for field in ('scalar', 'tensor')}
    amplitudes = {field: reader.formula(power[field]['normalization']) * corrections[field]
                  for field in ('scalar', 'tensor')}
    return amplitudes, corrections


def spectra(p, type='canonical', *, root=None):
    """Read scalar and tensor spectra at stored relative order p.

    type is 'canonical' (default) or 'kinflation'. Missing results raise
    FileNotFoundError; no solver or coefficient generator is ever called.
    """
    root = _root(root)
    _, data = _result_file(p, type, root)
    amplitudes, corrections = _spectra(p, type, root)
    return ArchivedSpectra(amplitudes['scalar'], amplitudes['tensor'],
        corrections['scalar'], corrections['tensor'], p, type, deepcopy(data['metadata']))


@lru_cache(maxsize=16)
def _consistency(n, p, type, root):
    _, data = _result_file(p, type, root)
    rows = data['consistency']
    if str(n) not in rows:
        raise LookupError(f'n={n}, p={p}, type={type!r} is not archived; available n: {", ".join(rows)}. No calculation was started.')
    row = rows[str(n)]
    if row['W'] != n + 1 + p:
        raise ValueError('Invalid absolute cutoff in archived consistency relation.')
    reader = _reader(data)
    components = tuple((part['weight'], reader.formula(part['expression'])) for part in row['components'])
    return components, sp.Add(*(expression for _, expression in components))


def consistency(n, p, type='canonical', *, root=None):
    """Read the archived relation alpha_n^t at W=n+1+p; never calculate it."""
    _integer(n, 'n')
    root = _root(root)
    _, data = _result_file(p, type, root)
    components, expression = _consistency(n, p, type, root)
    return ArchivedConsistency(expression, n, p, n + 1 + p, type,
                               components, deepcopy(data['metadata']))


def clear_cache():
    """Release loaded YAML families and expressions, or reload edited archives."""
    for function in (_load, _coefficient, _spectra, _consistency):
        function.cache_clear()
    ExactReader.factor.cache_clear()
