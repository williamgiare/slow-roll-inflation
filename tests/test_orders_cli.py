#========
# WG Edit
# Protect order inference and the exact public commands requested for SRI.
#========
from sri import serialization as structured
import pytest
from sri import ConsistencyOrder, consistency, spectra
from sri.cli import main, parser, render


@pytest.mark.parametrize('args', [dict(n=2, p=3), dict(n=2, W=6),
                                dict(p=3, W=6), dict(n=2, p=3, W=6)])
def test_equivalent_orders(args):
    assert ConsistencyOrder.resolve(**args) == ConsistencyOrder(2, 3, 6)


@pytest.mark.parametrize('args', [dict(n=2, p=3, W=5), dict(n=2, W=2),
    dict(p=3, W=3), dict(n=-1, p=0), dict(n=2), dict(n=True, p=1), dict(n=0.0, p=1)])
def test_invalid_orders(args):
    with pytest.raises(ValueError):
        ConsistencyOrder.resolve(**args)


def test_no_fixed_slow_roll_cap():
    assert ConsistencyOrder.resolve(n=2,p=20).W == 23


@pytest.mark.parametrize('value', [True, 1.0, -1])
def test_spectra_do_not_bypass_validation_with_a_cache_hit(value):
    spectra(p=1)
    with pytest.raises(ValueError):
        spectra(p=value)


@pytest.mark.parametrize('flags', [('-n', '2', '-p', '3'), ('-W', '6', '-n', '2'),
                                 ('-p', '3', '-W', '6')])
def test_cli_consistency_metadata(flags):
    value = structured.loads(render(parser().parse_args([*flags, '--format', 'yaml'])))
    assert (value['n'], value['p'], value['W'], value['error_order']) == (2, 3, 6, 7)


def test_cli_spectra_output_file(tmp_path):
    target = tmp_path / 'spectra.yaml'
    # WG Edit: spectra are selected by -p alone.
    assert main(['-p', '3', '--format', 'yaml', '--output', str(target)]) == 0
    value = structured.loads(target.read_text())
    assert value['mode'] == 'spectra' and value['relative_error_order'] == 4


def test_bad_cli_exits_with_failure(capsys):
    with pytest.raises(SystemExit) as result:
        main(['-n', '2', '-p', '3', '-W', '5'])
    assert result.value.code == 2
    assert 'Incompatible orders' in capsys.readouterr().err


#========
# WG Edit
# Protect mode inference, including n=0, simplify and removal of both spectra flags.
#========
@pytest.mark.parametrize('p',[0,3,4])
@pytest.mark.parametrize('simplify',[False,True])
@pytest.mark.parametrize('fmt',['text','latex','yaml'])
def test_p_alone_returns_spectra(p,simplify,fmt):
    flags=['-p',str(p),'--format',fmt,*(['--simplify'] if simplify else [])]
    output=render(parser().parse_args(flags))
    if fmt=='yaml':
        value=structured.loads(output)
        assert value['mode']=='spectra' and value['p']==p
        assert ('presentation' in value)==simplify
        assert value['coefficients']['analytic_cache_moments']==0
        if simplify:
            assert set(value['presentation'])=={'scalar','tensor','scalar_correction','tensor_correction'}
    elif fmt=='latex':
        assert r'\mathcal P_s = ' in output and r'\mathcal P_t = ' in output
    else:
        # WG Edit: the recap precedes aligned scalar and tensor equations.
        assert output.startswith('SRI | Slow-roll inflation')
        assert f'Orders: p = {p}' in output
        assert output.index('P_s =') < output.index('P_t =')
        assert 'Consistency relation' not in output


@pytest.mark.parametrize('flags',[
    ['-n','0','-p','0'],['-n','2','-p','4'],['-n','2','-W','7'],
    ['-p','4','-W','7'],['-n','2','-p','4','-W','7'],
])
@pytest.mark.parametrize('simplify',[False,True])
def test_consistency_selection_is_unchanged(flags,simplify):
    value=structured.loads(render(parser().parse_args([*flags,'--format','yaml',
                                               *(['--simplify'] if simplify else [])])))
    assert value['mode']=='consistency'
    assert value['W']==value['n']+1+value['p']
    assert ('presentation' in value)==simplify


@pytest.mark.parametrize('removed_flag',['-spectra','--spectra'])
def test_removed_spectra_flags_are_not_accepted(removed_flag,capsys):
    with pytest.raises(SystemExit) as result:
        main(['-p','3',removed_flag])
    assert result.value.code==2
    assert 'unrecognized arguments' in capsys.readouterr().err


def test_negative_spectra_order_fails(capsys):
    with pytest.raises(SystemExit) as result:
        main(['-p','-1'])
    assert result.value.code==2


@pytest.mark.parametrize('closed', [False, True])
@pytest.mark.parametrize('symbolic', [False, True])
@pytest.mark.parametrize('kinflation', [False, True])
@pytest.mark.parametrize('selection', [['-p', '2'], ['-n', '0', '-p', '2']])
def test_cli_default_and_alternative_modes(closed, symbolic, kinflation, selection):
    flags = [*selection, '--format', 'yaml']
    if closed:
        flags.append('--closed')
    if symbolic:
        flags.append('--symbolic')
    if kinflation:
        flags.append('--kinflation')
    actual = structured.loads(render(parser().parse_args(flags)))
    method = 'closed' if closed else 'iterative'
    assert actual['method'] == method
    assert actual['coefficients']['coefficient_backend'] == ('symbolic' if symbolic else 'numeric')
    expected = spectra(p=2, method=method, symbolic=symbolic, kinflation=kinflation)
    power = actual if actual['mode'] == 'spectra' else actual['spectra']
    assert power['scalar'] == str(expected.scalar)
    assert power['tensor'] == str(expected.tensor)
    if actual['mode'] == 'consistency':
        relation = consistency(n=0, p=2, method=method, symbolic=symbolic, kinflation=kinflation)
        assert actual['expression'] == str(relation.expression)
        assert power['coefficients']['coefficient_backend'] == actual['coefficients']['coefficient_backend']


@pytest.mark.parametrize('flag', ['--numerical', '--iterative'])
def test_default_modes_have_no_cli_flags(flag, capsys):
    with pytest.raises(SystemExit) as result:
        main(['-p', '2', flag])
    assert result.value.code == 2
    assert 'unrecognized arguments' in capsys.readouterr().err



def test_kinflation_cli_help_and_model_selection(capsys):
    arguments = parser()
    assert arguments.parse_args(['-p', '2']).kinflation is False
    assert arguments.parse_args(['-p', '2', '--kinflation']).kinflation is True
    help_text = arguments.format_help()
    assert '--kinflation' in help_text
    assert '--kessence' not in help_text
    with pytest.raises(SystemExit) as result:
        main(['-p', '2', '--kessence'])
    assert result.value.code == 2
    assert 'unrecognized arguments' in capsys.readouterr().err
