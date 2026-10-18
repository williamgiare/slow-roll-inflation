"""Scientific values and command output use standard, lossless YAML."""
import struct
import numpy as np
import yaml

from sri import serialization
from sri.cli import parser, render


def test_scientific_values_and_string_types_survive_yaml():
    data = {'analytic_value': '0.00424504144990432475268140',
            'words': ['0', '-1', '-1/2'],
            'labels': ['yes', 'no', '2026-10-14', '1e-12'],
            'values': [1e-300, 1e300, -0.0, 1.2345678901234567],
            'status': 'passed', 'formula': 'C = γ_E + log(2) - 2'}
    text = serialization.dumps(data)
    decoded = yaml.safe_load(text)
    assert decoded == data
    assert all(type(v) is str for v in decoded['labels'])
    assert type(decoded['analytic_value']) is str
    assert [struct.pack('!d',v) for v in decoded['values']] == [struct.pack('!d',v) for v in data['values']]


def test_command_yaml_is_readable_by_an_independent_yaml_loader():
    text = render(parser().parse_args(['-n','0','-p','2','--symbolic','--format','yaml']))
    data = yaml.safe_load(text)
    assert data['mode'] == 'consistency'
    assert data['coefficients']['coefficient_backend'] == 'symbolic'
    assert data['spectra']['p'] == 2


def test_numpy_float64_retains_its_binary_value():
    value = np.float64(1.2345678901234567e-14)
    decoded = yaml.safe_load(serialization.dumps({'difference': value}))
    assert type(decoded['difference']) is float
    assert struct.pack('!d', decoded['difference']) == struct.pack('!d', value)


def test_yaml_does_not_construct_arbitrary_objects():
    import pytest
    with pytest.raises(yaml.constructor.ConstructorError):
        serialization.loads('!!python/object/apply:os.system [echo invalid]')
