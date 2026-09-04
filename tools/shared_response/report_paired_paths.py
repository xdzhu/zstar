"""Compare two retained response ensembles without changing either source."""

import argparse
import json
from pathlib import Path

from report_benchmark import result, compare


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--unified', type=Path, required=True)
    parser.add_argument('--cartesian', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    a, b = result(args.unified), result(args.cartesian)
    comparison = compare(a, b)
    comparison.pop('measured_shared_over_cartesian_core_hours', None)
    if a.get('dimension') == 0:
        comparison['static_comparison'] = 'molecular_rigid_rotation_subspace_requires_separate_treatment'
    report = {'unified': a, 'cartesian': b, 'comparison': comparison}
    args.output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(comparison, indent=2))
