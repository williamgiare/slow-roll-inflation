"""Export the published equation collections to readable YAML without recalculation.

Only the polynomial LaTeX syntax written by sri.latex_output is accepted.
PDF and LaTeX files are read-only inputs and retain their published precision.
"""
import argparse
from hashlib import sha256
from pathlib import Path
import re
import yaml

ROOT = Path(__file__).resolve().parents[1]


class Formula(str):
    pass


class Dumper(getattr(yaml, 'CSafeDumper', yaml.SafeDumper)):
    def ignore_aliases(self, value):
        return True


Dumper.add_representer(Formula, lambda dumper, value:
    dumper.represent_scalar('tag:yaml.org,2002:str', str(value), style='|'))


class LatexPolynomial:
    """Translate the exporter's finite polynomial grammar; never execute input."""
    token = re.compile(r'\s*(\\frac|(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?|[A-Za-z_][A-Za-z_0-9]*|[{}()+*/^\-])')

    def __init__(self, text):
        text = text.replace(r'\log{\left(c_{s} \right)}', 'log(cs)')
        text = re.sub(r'\\epsilon_\{(\d+)\}', r'epsilon_\1', text)
        text = re.sub(r'\\alpha_\{(\d+)\}\^\{c_s\}', r'alpha_\1_cs', text)
        text = re.sub(r'\\alpha_\{(\d+)\}\^s', r'alpha_\1_s', text)
        text = text.replace('c_{s}', 'cs').replace(r'\left', '').replace(r'\right', '').replace(r'\cdot', '*')
        self.tokens = []
        position = 0
        while position < len(text.rstrip()):
            match = self.token.match(text, position)
            if match is None:
                raise ValueError(f'Unsupported polynomial LaTeX near {text[position:position+60]!r}.')
            self.tokens.append(match.group(1)); position = match.end()
        self.position = 0

    def peek(self):
        return self.tokens[self.position] if self.position < len(self.tokens) else None

    def take(self, expected=None):
        value = self.peek()
        if value is None or expected is not None and value != expected:
            raise ValueError(f'Expected {expected!r}, found {value!r}.')
        self.position += 1
        return value

    def atom(self):
        value = self.take()
        if value == r'\frac':
            self.take('{'); numerator = self.sum(); self.take('}')
            self.take('{'); denominator = self.sum(); self.take('}')
            return f'({numerator}) / ({denominator})'
        if value in ('{', '('):
            result = self.sum(); self.take('}' if value == '{' else ')')
            return f'({result})'
        if value == 'log':
            self.take('('); argument = self.sum(); self.take(')')
            return f'log({argument})'
        if not re.fullmatch(r'(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?|r|cs|epsilon_\d+|alpha_\d+_(?:s|cs)', value):
            raise ValueError(f'Unsupported polynomial atom {value!r}.')
        return value

    def power(self):
        value = self.atom()
        if self.peek() == '^':
            self.take('^'); exponent = self.signed_power()
            value = f'({value})**({exponent})'
        return value

    def signed_power(self):
        if self.peek() in ('+', '-'):
            sign = self.take(); return sign + self.power()
        return self.power()

    def product(self):
        result = self.signed_power()
        while self.peek() not in (None, '+', '-', '}', ')'):
            operation = self.take() if self.peek() in ('*', '/') else '*'
            result += ' ' + operation + ' ' + self.signed_power()
        return result

    def sum(self):
        return ' '.join(self.terms())

    def terms(self):
        parts = [self.product()]
        while self.peek() in ('+', '-'):
            sign = self.take(); parts.append(sign + ' ' + self.product())
        return parts

    def formula(self):
        result = self.terms()
        if self.peek() is not None:
            raise ValueError(f'Unexpected trailing LaTeX token {self.peek()!r}.')
        return Formula('\n'.join(result) + '\n')


def read_document(path, p, type):
    text = path.read_text()
    pattern = r'\\begin\{align\*\}\n(.*?)\n\\end\{align\*\}'
    groups = {}
    current = None
    for block in re.findall(pattern, text, re.S):
        pieces = []
        for line in block.splitlines():
            if '&=' in line:
                label, expression = line.split('&=', 1)
                match = re.fullmatch(r'([FR])_\{([st]|\d+),(\d+)\}\s*', label)
                if match is None:
                    raise ValueError(f'Unsupported equation label {label!r}.')
                current = match.groups()
            elif line.startswith(r'\phantom{'):
                expression = line.split(r'&\phantom{={}}{}', 1)[1]
            elif line.startswith('&{}'):
                expression = line[len('&{}'):]
            else:
                raise ValueError(f'Unsupported alignment line {line!r}.')
            pieces.append(expression.removesuffix(' \\\\'))
        if current is None:
            raise ValueError('Equation continuation without an initial label.')
        formula = LatexPolynomial(' '.join(pieces)).formula()
        groups.setdefault(current, []).append(formula)
    if not groups:
        raise ValueError(f'No exported equations found in {path}.')
    symbols = {'H': 'H', 'Mpl': 'Mpl', 'r': 'r', 'cs': 'c_s'}
    for formulae in groups.values():
        for formula in formulae:
            for alias in re.findall(r'epsilon_\d+|alpha_\d+_(?:cs|s)', formula):
                symbols[alias] = re.sub(r'alpha_(\d+)_(cs|s)', r'alpha_\1^\2', alias)
    spectra = {}
    for field, letter in (('scalar', 's'), ('tensor', 't')):
        parts = [(int(degree), Formula(''.join(values))) for (kind, index, degree), values in groups.items()
                 if kind == 'F' and index == letter]
        if sorted(degree for degree, _ in parts) != list(range(p + 1)):
            raise ValueError(f'Incomplete spectral weights in {path}.')
        correction = Formula(''.join(('' if j == 0 or expression.startswith(('+', '-')) else '+ ') + expression
                                    for j, (_, expression) in enumerate(sorted(parts))))
        normalization = ('H**2 / (8*pi**2*Mpl**2*epsilon_1' + ('*cs' if type == 'kinflation' else '') + ')'
                         if field == 'scalar' else '2*H**2 / (pi**2*Mpl**2)')
        symbols.setdefault('epsilon_1', 'epsilon_1')
        spectra[field] = {'normalization': normalization, 'correction': correction}
    relations = {}
    for n in sorted({int(index) for kind, index, degree in groups if kind == 'R'}):
        parts = sorted((int(degree), Formula(''.join(values))) for (kind, index, degree), values in groups.items()
                       if kind == 'R' and int(index) == n)
        if [degree for degree, _ in parts] != list(range(n + 1, n + 2 + p)):
            raise ValueError(f'Incomplete consistency weights for n={n} in {path}.')
        relations[str(n)] = {'W': n + 1 + p, 'components': [
            {'weight': degree, 'expression': expression} for degree, expression in parts]}
    return {'schema_version': 1, 'type': type, 'p': p,
            'metadata': {'method': 'iterative', 'coefficient_backend': 'numeric',
                         'coefficient_precision': 'published decimal coefficients',
                         'pivot': '-k_* eta_*=1', 'C': 'EulerGamma + log(2) - 2',
                         'spectrum_remainder_weight': p + 1},
            'symbols': symbols, 'spectra': spectra, 'consistency': relations}


def archive_document(source, p, type, *, target=None):
    """Write a YAML companion to one published LaTeX document."""
    source = Path(source)
    target = Path(target) if target is not None else source.with_suffix('.yaml')
    data = read_document(source, p, type)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(yaml.dump(data, Dumper=Dumper, sort_keys=False, allow_unicode=True, width=100))
    return target


def export(*, root=ROOT, output=None):
    root = Path(root)
    output = Path(output) if output is not None else root
    rows = []
    for type, directory in (('canonical', 'canonical_SRI'), ('kinflation', 'kinflation_SRI')):
        for source in sorted((root/'results'/directory).glob('*/*.tex')):
            p = 0 if source.stem == 'LO' else int(source.stem[1:-2])
            target = output/source.relative_to(root).with_suffix('.yaml')
            archive_document(source, p, type, target=target)
            data = yaml.load(target.read_text(), Loader=getattr(yaml, 'CSafeLoader', yaml.SafeLoader))
            rows.append({'type': type, 'p': p, 'path': str(target.relative_to(output)),
                         'tex_sha256': sha256(source.read_bytes()).hexdigest(),
                         'yaml_sha256': sha256(target.read_bytes()).hexdigest(),
                         'relations': list(data['consistency'])})
    if not rows:
        raise FileNotFoundError('No published equation documents found.')
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=ROOT)
    parser.add_argument('--output', type=Path, help='Alternate repository output root.')
    args = parser.parse_args()
    rows = export(root=args.root, output=args.output)
    print(f'Archived {len(rows)} order documents as YAML; no spectra or coefficients regenerated.')


if __name__ == '__main__':
    main()
