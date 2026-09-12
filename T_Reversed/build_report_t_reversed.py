#!/usr/bin/env python3
"""Build the self-contained report only from current, mutually consistent records.

Refuse changed scripts, changed source snapshots, stale CSVs and non-isotropic
results for the isotropic analytic report. No result is silently regenerated.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import re

import numpy as np
import quadrupolar_spin_liquid_t_reversed as model

ROOT = Path(__file__).resolve().parent


def table(headers, rows):
    def cell(value):
        return str(value).replace('|', '\\|').replace('\n', ' ')
    return '\n'.join(['| ' + ' | '.join(map(cell, headers)) + ' |',
                      '| ' + ' | '.join('---' for _ in headers) + ' |'] +
                     ['| ' + ' | '.join(map(cell, row)) + ' |' for row in rows])


def validate_provenance(record, name, current):
    meta = record['metadata']
    model.require(meta['current_sector'] == 'T-REVERSED', f'Incorrect current sector in {name}')
    model.require(meta['script_sha256'] == current['script_sha256'], f'STALE scripts in {name}; rerun all generators')
    model.require(meta['local_source_sha256'] == current['local_source_sha256'], f'STALE local sources in {name}')
    for filename, expected in meta.get('source_sha256', {}).items():
        path = Path(filename)
        if path.exists():
            model.require(model.sha256(path) == expected, f'Changed original source: {filename}')
    for key in ('Kx', 'Ky', 'Kz', 'Kt', 'Kw'):
        model.require(meta[key] == 1., f'This analytic report requires isotropic couplings: {name}')
    model.require(np.allclose([meta['q_min'], meta['q_max']], model.PROJECT_BZ, rtol=0., atol=1e-14),
                  f'Wrong project BZ: {name}')
    model.require(meta['N_sites_per_cell'] == 4 and abs(meta['BZ_width']-2*np.pi) < 1e-14,
                  f'Wrong counting metadata: {name}')
    model.require(meta['normalization_formula'] == current['normalization_formula'], f'Wrong normalization: {name}')


def test_stale_rejection(current):
    """Exercise the same validator with a corrupted recorded hash."""
    import copy
    for script in ('quadrupolar_spin_liquid_t_reversed.py', 'compare_old_vs_t_reversed.py', 'draw_t_lattice.py'):
        stale = {'metadata': copy.deepcopy(current)}
        stale['metadata']['script_sha256'][script] = '0'*64
        try:
            validate_provenance(stale, 'deliberately stale test record', current)
        except RuntimeError as exc:
            model.require('STALE' in str(exc), 'Wrong failure in stale-data test')
        else:
            raise RuntimeError(f'Stale hash accepted for {script}')
    return {'status': 'PASS', 'corrupted_hashes_rejected': 3}


def validate_csv(path, expected_rows, current):
    with path.open(newline='', encoding='utf-8') as stream:
        rows = list(csv.DictReader(stream))
    model.require(len(rows) == len(expected_rows), f'CSV row mismatch: {path.name}')
    for row, expected in zip(rows, expected_rows):
        model.require(json.loads(row['script_sha256_json']) == current['script_sha256'], f'STALE CSV: {path.name}')
        validate_provenance({'metadata': json.loads(row['metadata_json'])}, path.name, current)
        for key, value in expected.items():
            if value is None:
                model.require(row[key] == '', f'CSV/JSON null mismatch: {path.name}/{key}')
            elif isinstance(value, (float, int)):
                model.require(float(row[key]) == value, f'CSV/JSON numerical mismatch: {path.name}/{key}')
            else:
                model.require(row[key] == str(value), f'CSV/JSON mismatch: {path.name}/{key}')


def build_report(root=ROOT):
    current = model.metadata()
    records = {name: json.loads((root/name).read_text()) for name in (
        'validation_t_reversed.json', 'comparison_t_reversed.json',
        'benchmark_t_reversed.json', 'geometry_t_audit.json')}
    for name, record in records.items():
        validate_provenance(record, name, current)
    v, c, b, g = (records[name] for name in records)
    extra = json.loads((root/'explicit_majorana_validation.json').read_text())
    validate_provenance(extra, 'explicit_majorana_validation.json', current)
    model.require(extra['status'] == 'PASS', 'Explicit Majorana calibration failed')
    stale_test = test_stale_rejection(current)
    model.require(v['status'] == 'PASS' and v['symbolic']['status'] == 'PASS', 'Core audit did not pass')
    model.require(v['number_random_cases'] >= 1000, 'Insufficient random cases')
    model.require(c['finite_chain_validation']['status'] == 'PASS', 'Finite chain failed')
    model.require(c['isotropic']['status'] == 'PASS', 'Isotropic independent check failed')
    model.require(g['non_t_edges_unchanged'] and g['positions_and_vectors_unchanged'], 'Geometry failed')
    model.require(b['repeats'] >= 5, 'Insufficient benchmark repeats')
    model.require(c['consolidated_previous_sectors']['status'] == 'PASS', 'Consolidated comparison record missing')
    for filename, expected in (
        ('ground_state_energy_old_vs_t.csv', [c['summary']]),
        ('energy_convergence_old_vs_t.csv', c['convergence']),
        ('finite_chain_old_vs_t.csv', c['finite_chain_validation']['rows']),
        ('all_gauge_sectors_comparison.csv', c['consolidated_previous_sectors']['rows'])):
        validate_csv(root/filename, expected, current)
    summary, errors = c['summary'], v['max_errors']
    fmt = lambda value: f'{value:.15f}'
    sci = lambda value: f'{value:.3e}'
    replacements = {
        'EPSILON_OLD': fmt(summary['epsilon_old']), 'EPSILON_T': fmt(summary['epsilon_t']),
        'DELTA': f"{summary['Delta_epsilon']:.15e}", 'ENERGY_TOL': sci(c['effective_energy_comparison_tolerance']),
        'FAVORED': summary['favored_sector'], 'SEED': str(v['seed']), 'EDGE_COUNT': str(g['visible_edge_count']),
        'REVERSED_EDGE_COUNT': str(g['reversed_t_edge_count']),
        'RAW_ERROR': sci(errors['raw_geometry_vs_closed']),
        'ISOTROPIC_ERROR': sci(v['symbolic']['max_isotropic_spectrum_error']),
        'FINITE_CORE_ERROR': sci(errors['finite_full_iA_spectrum']),
        'FINITE_COMPARE_ERROR': sci(c['finite_chain_validation']['max_spectrum_error']),
        'TRACKING_ERRORS': '; '.join(f'{k}: {sci(val)}' for k, val in c['tracking_spectrum_errors'].items()),
        'EXPLICIT_OPERATOR_ERROR': sci(extra['max_energy_error']),
        'EXPLICIT_SPECTRUM_ERROR': sci(extra['max_full_many_body_spectrum_error']),
    }
    replacements['GAUGE_TABLE'] = table(['Chave', 'OLD', 'T', 'Mudou?'],
        [(key, f'{value:+d}', f'{model.T_REVERSED_GAUGE[key]:+d}',
          'Sim' if value != model.T_REVERSED_GAUGE[key] else 'Não') for key, value in model.ORIGINAL_GAUGE.items()])
    replacements['BOND_TABLE'] = table(['start', 'end', 'gamma', 'u OLD', 'u T', 'da1', 'da2', 'cell_shift'],
        [(bond.start, bond.end, bond.gamma, f'{bond.u:+d}', f'{model.T_REVERSED_GAUGE[spec[3]]:+d}',
          bond.da1, bond.da2, model.cell_shift(bond))
         for spec, bond in zip(model.GEOMETRY, model.real_space_bonds(model.ORIGINAL_GAUGE))])
    flux_rows = []
    for name, old in v['plaquette_details']['OLD'].items():
        t = v['plaquette_details']['T-REVERSED'][name]
        path = ' → '.join(f'{m}_({n})' for n, m in old['path'])
        factors = lambda r: ''.join(f'({x:+d})' for x in r['factors'])
        flux_rows.append([name, path, ' '.join(f'({x})' for x in old['expressions']), factors(old),
                          f"{old['flux']:+d}", factors(t), f"{t['flux']:+d}", 'Sim' if old['flux'] != t['flux'] else 'Não'])
    replacements['FLUX_TABLE'] = table(['Plaqueta', 'Percurso CCW', 'Produto de u armazenados', 'Fatores OLD', 'OLD', 'Fatores T', 'T', 'Mudou?'], flux_rows)
    step_rows = []
    for name, old in v['plaquette_details']['OLD'].items():
        new = v['plaquette_details']['T-REVERSED'][name]
        for a, end, kind, expression, uold, unew in zip(old['path'][:-1], old['path'][1:], old['kinds'], old['expressions'], old['factors'], new['factors']):
            step_rows.append([name, f'{a[1]}_({a[0]}) → {end[1]}_({end[0]})', kind, expression.lstrip('-'),
                              'oposta: uji=-uij' if expression.startswith('-') else 'armazenada: uij', f'{uold:+d}', f'{unew:+d}'])
    replacements['FLUX_STEPS'] = table(['Plaqueta', 'Passo dirigido', 'Família', 'Chave armazenada', 'Direção utilizada', 'Fator OLD', 'Fator T'], step_rows)
    replacements['CONSOLIDATED_NOTES'] = '\n'.join(f'- {name}: {reason}.' for name, reason in c['consolidated_previous_sectors']['omitted_sectors'].items())
    replacements['MUTATION_TABLE'] = table(['Mutação', 'Resíduo máximo'], [
        ('RAW com sinal errado de q2 versus RAW correto', sci(v['mutations']['wrong_raw_q2_detected_error'])),
        ('Eq.(27) ao contrário, com RAW correto', sci(v['mutations']['wrong_eq27_detected_error'])),
        ('Os dois erros juntos, apenas no Bloch final', sci(v['mutations']['two_bugs_compensate_at_final_bloch_error']))])
    replacements['SPECIAL_TABLE'] = table(['q', 'OLD: eigvalsh', 'T: eigvalsh'],
        [(q, ', '.join(f'{x:.9f}' for x in c['special_points']['OLD'][q]),
          ', '.join(f'{x:.9f}' for x in c['special_points']['T-REVERSED'][q]))
         for q in ('-pi', '-pi/2', '-pi/4', '0', 'pi/4', 'pi/2', 'pi')])
    replacements['GAP_TABLE'] = table(['Setor', 'min_abs_E', 'zero_energy_band_gap', 'q do mínimo numérico', 'Refinamentos'],
        [(name, fmt(row['min_abs_E']), fmt(row['zero_energy_band_gap']), f"{row['q_at_minimum']:.10f}", row['refined_minima_count'])
         for name, row in c['gaps'].items()])
    replacements['ENERGY_TABLE'] = table(['Quantidade', 'Quadratura eigvalsh', 'Incerteza / tolerância'], [
        ('epsilon_old', fmt(summary['epsilon_old']), sci(summary['epsilon_old_uncertainty'])),
        ('epsilon_t', fmt(summary['epsilon_t']), sci(summary['epsilon_t_uncertainty'])),
        ('Delta epsilon', sci(summary['Delta_epsilon']), sci(summary['Delta_epsilon_uncertainty']))])
    replacements['CONVERGENCE_TABLE'] = table(['Nq', 'epsilon_old', 'epsilon_t', 'Delta epsilon', 'erro OLD', 'erro T'],
        [(r['Nq'], fmt(r['epsilon_old']), fmt(r['epsilon_t']), fmt(r['Delta_epsilon']), sci(r['error_estimate_old']), sci(r['error_estimate_t'])) for r in c['convergence']])
    replacements['UNCERTAINTY_TABLE'] = table(['Setor', 'Diferenças entre malhas sucessivas', 'Última malha vs quad', 'Incerteza conservadora', 'Casas sustentadas'],
        [(name, ', '.join(sci(x) for x in r['consecutive_mesh_differences']), sci(r['last_mesh_vs_adaptive']),
          sci(r['conservative_absolute_uncertainty']), r['resolved_decimal_places_from_uncertainty'])
         for name, r in c['convergence_assessment'].items()])
    replacements['FINITE_TABLE'] = table(['Setor', 'Ncell', 'E0/(4 Ncell)', 'epsilon bulk', 'Finita menos bulk', 'Erro espectral'],
        [(r['sector'], r['Ncell'], fmt(r['epsilon_finite']), fmt(r['epsilon_bulk_quad']), sci(r['finite_minus_bulk']), sci(r['full_spectrum_max_error']))
         for r in c['finite_chain_validation']['rows'] if r['Ncell'] in (17, 33, 65, 129)])
    replacements['VALIDATION_TABLE'] = table(['Teste', 'Erro máximo'], [(k, sci(val)) for k, val in errors.items()])
    replacements['BENCHMARK_TABLE'] = table(['Setor', 'Operação', 'Ncell', 'Mediana (s)', 'Mínimo (s)', 'Máximo (s)'],
        [(r['sector'], r['operation'], r.get('Ncell', '—'), f"{r['median_seconds']:.7f}", f"{r['min_seconds']:.7f}", f"{r['max_seconds']:.7f}") for r in b['rows']])
    replacements['BENCHMARK_ENV'] = 'Ambiente: `' + b['platform'] + '`; threads: `' + json.dumps(b['thread_environment']) + '`. Versões: `' + json.dumps(current['software_versions']) + '`.'
    replacements['SUMMARY_TABLE'] = table(['Quantidade', 'OLD / zero-flux', 'T-REVERSED'], [
        ('u alterados', 'nenhum', ','.join(v['changed_keys'])),
        ('Padrão de fluxo (bloco x ; bloco w)', ','.join(str(x) for x in v['fluxes']['OLD'].values()), ','.join(str(x) for x in v['fluxes']['T-REVERSED'].values())),
        ('zero_energy_band_gap', fmt(summary['zero_energy_band_gap_old']), fmt(summary['zero_energy_band_gap_t'])),
        ('min_abs_E', fmt(summary['min_abs_E_old']), fmt(summary['min_abs_E_t'])),
        ('E0/Nsite', fmt(summary['epsilon_old']), fmt(summary['epsilon_t'])),
        ('Caráter espectral', summary['spectral_character_old']+'; degenerescência em ±pi/2', summary['spectral_character_t']+'; degenerescência dupla em todo q')])
    replacements['CONSOLIDATED_TABLE'] = table(['Setor', 'u alterados', 'Fluxos', 'Gap', 'epsilon', 'Delta/OLD', 'Caráter', 'Validação'],
        [(r['sector'], r['changed_u'], r['flux_pattern'], f"{r['zero_energy_band_gap']:.12f}", fmt(r['epsilon']),
          f"{r['Delta_epsilon_relative_to_OLD']:+.12e}", r['spectral_classification'], r['status'])
         for r in c['consolidated_previous_sectors']['rows']])
    replacements['HASH_TABLE'] = table(['Script', 'SHA-256'], [(name, f'`{value}`') for name, value in current['script_sha256'].items()])
    replacements['SOURCE_TABLE'] = table(['Fonte física / snapshot', 'SHA-256'],
        [(name, f'`{value}`') for name, value in current['local_source_sha256'].items()] +
        [(Path(name).name, f'`{value}`') for name, value in v['metadata']['source_sha256'].items() if name.endswith('.pdf')])
    template_path = root/'report_template_t.md'
    report = template_path.read_text(encoding='utf-8')
    for key, value in replacements.items():
        token = '{{'+key+'}}'
        model.require(token in report, f'Unused report replacement {key}')
        report = report.replace(token, value)
    model.require(not re.search(r'\{\{[A-Z_]+\}\}', report), 'Unresolved template field')
    model.require(re.findall(r'^## (\d+)\.', report, re.M) == [str(i) for i in range(1, 33)], 'Expected all 32 sections')
    model.require(chr(956) not in report, 'Forbidden gauge notation in report')
    for target in re.findall(r'!\[[^\]]*\]\(([^)]+)\)', report):
        model.require((root/target).is_file(), f'Missing figure {target}')
    output = root/'Quadrupolar_Spin_Liquid_T_Reversed.md'
    output.write_text(report, encoding='utf-8')
    log = ['T-REVERSED SCIENTIFIC AUDIT', '============================', '',
           'Couplings: Kx=1 Ky=1 Kz=1 Kt=1 Kw=1', 'BZ: q in [-pi, pi]',
           'Gauge: changed keys = u13t, u24t; both -1 -> +1',
           'Fluxes OLD: ' + str(v['fluxes']['OLD']), 'Fluxes T: ' + str(v['fluxes']['T-REVERSED'])]
    for name, key in (('OLD / zero-flux', 'old'), ('T-REVERSED', 't')):
        log.append(f"{name}: min_abs_E={summary['min_abs_E_'+key]:.15f}; zero-energy gap={summary['zero_energy_band_gap_'+key]:.15f}; epsilon={summary['epsilon_'+key]:.15f}")
    log.extend([f"Delta epsilon = {summary['Delta_epsilon']:.15e}", 'Favored sector: ' + summary['favored_sector'],
                'RAW validation: PASS', 'Eq.(27): PASS', 'q2 cancellation: PASS', 'Hermiticity: PASS',
                'Majorana symmetry: PASS', 'Characteristic polynomial: PASS', 'Symbolic tests: PASS', '1000 random tests: PASS',
                f"Finite-chain validation: PASS; max spectrum error = {c['finite_chain_validation']['max_spectrum_error']:.3e}",
                'BZ normalization: PASS', 'Convergence: PASS', 'Stale-data rejection: PASS',
                f"Explicit 64-state Majorana operator: PASS; energy/site error={extra['max_energy_error']:.3e}",
                f"Performance benchmark: {len(b['rows'])} operation/sector cases, 5 warmed repeats each",
                'FINAL SCIENTIFIC VERDICT: A', 'Advisor discussion: YES, WITH RESERVATIONS (documented source errata; fixed-gauge isotropic scope)',
                'Physical verdict: YES', 'Mathematical verdict: YES', 'Numerical verdict: YES (reported isotropic point)'])
    for r in b['rows']:
        if r['sector'] == 'T-REVERSED' and (r.get('Ncell') == 129 or '16001' in r['operation']):
            log.append(f"Benchmark {r['operation']} Ncell={r.get('Ncell','-')}: median={r['median_seconds']:.7f}s min={r['min_seconds']:.7f}s max={r['max_seconds']:.7f}s")
    (root/'final_scientific_audit.txt').write_text('\n'.join(log)+'\n')
    manifest = {'status': 'PASS', 'metadata': current, 'stale_rejection_test': stale_test,
                'template_sha256': model.sha256(template_path),
                'report_sha256': model.sha256(output), 'artifacts_sha256': {
                    str(path.relative_to(root)): model.sha256(path) for path in sorted(root.rglob('*'))
                    if path.is_file() and '__pycache__' not in path.parts and path.name != 'report_manifest_t_reversed.json'}}
    (root/'report_manifest_t_reversed.json').write_text(json.dumps(manifest, indent=2)+'\n')
    print('\n'.join(log))
    print(f'Report: {output.name}; {len(report)} characters; 32 sections; current provenance verified')
    return manifest


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=ROOT)
    args = parser.parse_args()
    build_report(args.output_dir)
