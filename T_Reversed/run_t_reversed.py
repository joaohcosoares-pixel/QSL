#!/usr/bin/env python3
"""Compile, execute, audit and package fresh T-REVERSED results.

Run with python3 -W error run_t_reversed.py. Every child uses warnings as
errors; scientific files must remain unchanged throughout the run.
"""
from __future__ import annotations

import json
import os
from pathlib import Path
import py_compile
import subprocess
import sys
import zipfile

os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('MPLCONFIGDIR', '/tmp/qsl_t_mpl')

ROOT = Path(__file__).resolve().parent


def main():
    import quadrupolar_spin_liquid_t_reversed as model
    initial = model.metadata()
    scripts = sorted(ROOT.glob('*.py'))
    for script in scripts:
        py_compile.compile(str(script), doraise=True)
    print(f'Compilation PASS: {len(scripts)} scripts', flush=True)
    (ROOT/'logs').mkdir(exist_ok=True)

    def run(script, *args, optimized=False, log_name=None):
        command = [sys.executable, '-W', 'error']
        if optimized:
            command.append('-O')
        command.extend([str(ROOT/script), *map(str, args)])
        print('Running: '+script+(' (-O)' if optimized else ''), flush=True)
        result = subprocess.run(command, cwd=ROOT, text=True, stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, check=False)
        (ROOT/'logs'/(log_name or script.replace('.py', '.log'))).write_text(result.stdout)
        model.require(result.returncode == 0, f'{script} failed; see logs.\n{result.stdout[-5000:]}')
        print(result.stdout[-3500:], flush=True)

    run('quadrupolar_spin_liquid_t_reversed.py')
    run('audit_explicit_majoranas.py', '--output', ROOT/'explicit_majorana_validation.json')
    run('draw_t_lattice.py')
    run('compare_old_vs_t_reversed.py')
    # Full suite under optimization: critical conditions must survive -O.
    run('quadrupolar_spin_liquid_t_reversed.py', optimized=True, log_name='optimized_core.log')

    current = model.metadata()
    for key in ('script_sha256', 'local_source_sha256', 'source_sha256'):
        model.require(initial[key] == current[key], f'Source changed during run: {key}')

    # Structural PNG checks complement the human visual inspection recorded
    # separately. They do not claim to assess graphic meaning or typography.
    from PIL import Image
    pngs = {}
    for filename in ('qsl_lattice_t_reversed.png', 'dispersion_old_vs_t_reversed.png',
                     'Fq_old_vs_t_reversed.png', 'all_gauge_sectors_comparison.png'):
        path = ROOT/filename
        with Image.open(path) as im:
            im.load()
            dpi = im.info.get('dpi', (0., 0.))
            model.require(min(dpi) > 299., f'PNG resolution below 300 dpi: {filename}')
            model.require(min(im.size) >= 1000, f'PNG unexpectedly small: {filename}')
            pngs[filename] = {'size': list(im.size), 'dpi': list(dpi), 'sha256': model.sha256(path)}
    for script in scripts:
        model.require(chr(956) not in script.read_text(), f'Forbidden gauge notation in {script.name}')
    checks = {'status': 'PASS', 'metadata': current, 'compiled_scripts': [p.name for p in scripts],
              'warnings_as_errors': True, 'optimized_core_full_suite': 'PASS',
              'png_structural_checks': pngs, 'visual_review': 'Requires inspection of rendered PNGs',
              'current_sector': 'T-REVERSED'}
    (ROOT/'execution_checks.json').write_text(json.dumps(checks, indent=2)+'\n')

    # Run in-process so its final manifest does not hash a still-changing log.
    import build_report_t_reversed as builder
    manifest = builder.build_report(ROOT)
    for relative, expected in manifest['artifacts_sha256'].items():
        model.require(model.sha256(ROOT/relative) == expected, f'Final artifact mismatch: {relative}')
    target = ROOT.parent/'T_Reversed_entregaveis.zip'
    with zipfile.ZipFile(target, 'w', compression=zipfile.ZIP_DEFLATED) as archive:
        for path in sorted(ROOT.rglob('*')):
            if path.is_file() and '__pycache__' not in path.parts:
                archive.write(path, arcname=str(Path(ROOT.name)/path.relative_to(ROOT)))
    print(f'ALL CHECKS PASS\nPackage: {target}', flush=True)


if __name__ == '__main__':
    main()
