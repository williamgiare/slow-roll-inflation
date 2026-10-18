#========
# WG Edit
# Preserve historical source fingerprints after a checked additive backend extension.
# Any later source change invalidates the recorded compatibility automatically.
#========
from pathlib import Path
from sri import serialization as structured


def compatible_report(root, previous, current):
    entries = []
    for name in ('symbolic_backend_compatibility.yaml', 'yaml_format_compatibility.yaml',
                 'release_metadata_compatibility.yaml'):
        path = Path(root)/'validation/common/reports'/name
        if not path.exists():
            continue
        report = structured.loads(path.read_text())
        if report.get('status') == 'passed':
            entries.extend(report.get('accepted_fingerprints', []))
    pending = [previous]
    seen = []
    while pending:
        candidate = pending.pop()
        if candidate == current:
            return True
        if candidate in seen:
            continue
        seen.append(candidate)
        pending.extend(entry['current_source_sha256'] for entry in entries
                       if entry.get('previous_source_sha256') == candidate)
    return False


def preserved_coefficient_signature(root, current):
    """Retain a verified mathematical signature across an exact format conversion."""
    path = Path(root)/'validation/common/reports/yaml_format_compatibility.yaml'
    if path.exists():
        report = structured.loads(path.read_text())
        entry = report.get('coefficient_signature', {})
        if report.get('status') == 'passed' and entry.get('current_raw_signature') == current:
            return entry['preserved_signature']
    return current
