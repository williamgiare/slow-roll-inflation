#========
# WG Edit
# One safe YAML encoding for command output, scientific fixtures and reports.
#========
import yaml

_LOADER = getattr(yaml, 'CSafeLoader', yaml.SafeLoader)
class _DUMPER(getattr(yaml, 'CSafeDumper', yaml.SafeDumper)):
    def ignore_aliases(self, value):
        return True


# NumPy float64 and other float subclasses have the same numerical payload.
_DUMPER.add_multi_representer(
    float, lambda dumper, value: dumper.represent_float(float(value)))


def loads(text):
    """Read standard YAML data without constructing arbitrary Python objects."""
    return yaml.load(text, Loader=_LOADER)


def dumps(value, *, indent=2):
    """Write readable YAML while preserving string types and floating-point values."""
    return yaml.dump(value, Dumper=_DUMPER, sort_keys=False,
                     allow_unicode=True, indent=indent)
