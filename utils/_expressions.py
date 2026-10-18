"""Safe decoder for exact universal-coefficient YAML expressions."""
import ast
from functools import lru_cache
from pathlib import Path
import sympy as sp
import yaml
from sri.analytic import C
from sri.analytic.generic_hyperlogs import GPL

Loader = getattr(yaml, "CSafeLoader", yaml.SafeLoader)

class ExactReader:
    """Decode the documented exact syntax without executing file contents."""
    functions = {'zeta': sp.zeta, 'dirichlet_eta': sp.dirichlet_eta,
                 'polylog': sp.polylog, 'log': sp.log, 'gamma': sp.gamma,
                 'polygamma': sp.polygamma}

    def __init__(self, words):
        self.words = words
        self.names = {'C': C, 'pi': sp.pi, 'EulerGamma': sp.EulerGamma}

    @lru_cache(maxsize=131072)
    def factor(self, text):
        return self.node(ast.parse(text, mode='eval').body)

    def node(self, node):
        if isinstance(node, ast.Constant) and type(node.value) is int:
            return sp.Integer(node.value)
        if isinstance(node, ast.Name):
            if node.id in self.words:
                return GPL(sp.Tuple(*(sp.Rational(v) for v in self.words[node.id])))
            if node.id in self.names:
                return self.names[node.id]
            raise ValueError(f'Unknown exact constant {node.id}.')
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.USub, ast.UAdd)):
            value = self.node(node.operand)
            return -value if isinstance(node.op, ast.USub) else value
        if isinstance(node, ast.BinOp):
            left, right = self.node(node.left), self.node(node.right)
            if isinstance(node.op, ast.Add): return left+right
            if isinstance(node.op, ast.Sub): return left-right
            if isinstance(node.op, ast.Mult): return left*right
            if isinstance(node.op, ast.Div): return left/right
            if isinstance(node.op, ast.Pow): return left**right
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and not node.keywords:
            function = self.functions.get(node.func.id)
            if function is not None: return function(*(self.node(a) for a in node.args))
        raise ValueError('Unsupported node in an exact coefficient expression.')

    def expression(self, text):
        terms = []
        for line in text.splitlines():
            parts = line.removeprefix('+ ').split(' * ')
            coefficient = sp.Rational(parts[0])
            terms.append(sp.Mul(coefficient, *(self.factor(part) for part in parts[1:])))
        return sp.Add(*terms)


def read_family(path):
    """Load a public YAML family; definitions remain exact SymPy expressions."""
    with Path(path).open() as handle:
        data = yaml.load(handle, Loader=Loader)
    if data['schema_version'] != 1:
        raise ValueError('Unsupported coefficient archive schema.')
    reader = ExactReader(data.get('gpl_words', {}))
    return {tuple(row['pattern']): reader.expression(row['exact_expression']) for row in data['coefficients']}

