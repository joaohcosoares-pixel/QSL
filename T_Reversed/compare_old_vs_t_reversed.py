#!/usr/bin/env python3
"""OLD/T physical numerics on PROJECT_BZ=(-pi,pi), independent finite checks.

All physical quantities use eigvalsh; eigh and overlap tracking are exclusively
for dispersion figures. Run with OPENBLAS_NUM_THREADS=1 python3 -W error ... .
The core derives the physics from the real bonds. This module adds independent
quadrature, finite-size comparisons, convergence and reproducible timings.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
from itertools import permutations
import json
import os
from pathlib import Path
import platform
import time

os.environ.setdefault("MPLCONFIGDIR", "/tmp/qsl_t_mpl")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq, minimize_scalar
from scipy.special import ellipe

import quadrupolar_spin_liquid_t_reversed as model

CONVERGENCE_NQ = (1001, 2001, 4001, 8001, 16001, 32001)
FINITE_SIZES = (5, 8, 17, 33, 65, 129)
FINITE_CONVERGENCE_SIZES = (17, 33, 65, 129)
PERMUTATIONS = np.asarray(tuple(permutations(range(4))), dtype=int)
GAP_ZERO_TOL = model.GAP_TOLERANCE
ENERGY_FLOOR = 1e-11
COUPLING_NAMES = ("Kx", "Ky", "Kz", "Kt", "Kw")


def require(condition, message):
    """Critical checks stay active under python -O."""
    if not condition:
        raise RuntimeError(message)


def script_hashes():
    root = Path(__file__).resolve().parent
    names = ("quadrupolar_spin_liquid_t_reversed.py", "compare_old_vs_t_reversed.py", "draw_t_lattice.py")
    require(all((root/name).is_file() for name in names), "Required scientific script missing")
    return {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(root.glob("*.py"))}


def physical_bands(q, couplings, gauge):
    values = np.asarray(model.bands_sorted(q, couplings, gauge), dtype=float)
    require(values.shape == (len(q), 4), "Expected four physical bands at every q")
    require(np.isfinite(values).all(), "Non-finite physical eigenvalues")
    return values


def tracked_bands(q_values, couplings, gauge):
    """Choose the largest total squared overlap among all 24 permutations.

    Tracking is only a visual assignment. Within exactly degenerate eigenspaces,
    a unitary Procrustes transport avoids arbitrary basis rotations. Degenerate
    labels are a plotting convention and have no independent physical meaning.
    """
    result = np.empty((len(q_values), 4))
    values, previous = np.linalg.eigh(model.h_bloch(q_values[0], couplings, gauge))
    result[0] = values
    for k, q in enumerate(q_values[1:], 1):
        values, current = np.linalg.eigh(model.h_bloch(q, couplings, gauge))
        overlaps = abs(previous.conj().T @ current) ** 2
        scores = overlaps[np.arange(4)[None, :], PERMUTATIONS].sum(axis=1)
        chosen = PERMUTATIONS[int(np.argmax(scores))]
        values, current = values[chosen], current[:, chosen]
        remaining = set(range(4))
        tolerance = 2e-13 * max(1., float(np.max(abs(values))))
        while remaining:
            first = min(remaining)
            group = sorted(i for i in remaining if abs(values[i] - values[first]) < tolerance)
            remaining.difference_update(group)
            if len(group) > 1:
                left, _, right = np.linalg.svd(current[:, group].conj().T @ previous[:, group])
                current[:, group] = current[:, group] @ (left @ right)
        result[k] = values
        previous = current
    return result


def gap_analysis(couplings, gauge, nq=16001):
    """Dense grid, all resolved local minima, shifted-coordinate refinement.

    Global Majorana symmetry on the symmetric project interval makes the global
    gap 2*min|E|, even without same-q pairing. Generic globality is numerical;
    the core's isotropic analytic spectrum supplies the independent proof.
    """
    q = np.linspace(*model.PROJECT_BZ, nq)
    values = physical_bands(q, couplings, gauge)
    distances = np.min(abs(values), axis=1)
    scale = max(1., float(np.max(abs(values))))
    noise = 64 * np.finfo(float).eps * scale
    local = np.where((distances[1:-1] <= distances[:-2] + noise)
                     & (distances[1:-1] <= distances[2:] + noise)
                     & (np.maximum(distances[:-2], distances[2:]) - distances[1:-1] > noise))[0] + 1
    candidates = [{"q": float(q[j]), "min_abs_E": float(distances[j]), "method": "endpoint"}
                  for j in (0, -1)]
    j_global = int(np.argmin(distances))
    candidates.append({"q": float(q[j_global]), "min_abs_E": float(distances[j_global]), "method": "dense-grid"})

    def distance(x):
        return float(np.min(abs(np.linalg.eigvalsh(model.h_bloch(x, couplings, gauge)))))

    for j in local:
        center = float(q[j])
        candidates.append({"q": center, "min_abs_E": float(distances[j]), "method": "grid-local-minimum"})
        # Centering the optimization variable avoids the sqrt(eps)*abs(q)
        # stopping floor of bounded minimization near nonzero cusp positions.
        optimum = minimize_scalar(lambda shift: distance(center + shift),
                                  bounds=(float(q[j - 1] - center), float(q[j + 1] - center)),
                                  method="bounded", options={"xatol": 2e-14, "maxiter": 200})
        require(optimum.success, "Continuous gap minimization failed")
        candidates.append({"q": center + float(optimum.x), "min_abs_E": float(optimum.fun),
                           "method": "bounded-continuous-refinement"})
    crossings = []
    for band in range(4):
        changes = np.where(np.signbit(values[:-1, band]) != np.signbit(values[1:, band]))[0]
        for j in changes:
            # Sorted eigenvalues are continuous, including at degeneracies.
            root = brentq(lambda x: float(np.linalg.eigvalsh(model.h_bloch(x, couplings, gauge))[band]),
                          float(q[j]), float(q[j+1]), xtol=1e-14)
            entry = {"q": float(root), "min_abs_E": distance(root), "method": "bracketed-band-sign-crossing"}
            candidates.append(entry)
            crossings.append({"band_sorted_index": band, "left": float(q[j]), "right": float(q[j+1]), **entry})
    near_zero_grid = np.where(distances < GAP_ZERO_TOL*scale)[0]
    # Keep representative touching points without manufacturing sign changes:
    # a paired ±E cone touches zero without a sorted band's sign switching.
    touching = [entry for entry in candidates if entry["min_abs_E"] < GAP_ZERO_TOL*scale]
    best = min(candidates, key=lambda row: row["min_abs_E"])
    symmetry_samples = q[::max(1, (nq - 1) // 80)]
    majorana_error = max(float(np.max(abs(model.h_bloch(-x, couplings, gauge)
                                         + model.h_bloch(x, couplings, gauge).conj())))
                         for x in symmetry_samples)
    require(majorana_error < 2e-11 * scale, "Global particle-hole gap relation not validated")
    pairing_error = float(np.max(abs(values + values[:, ::-1])))
    raw_gap = 2 * best["min_abs_E"]
    classified_zero = best["min_abs_E"] < GAP_ZERO_TOL * scale
    # A plateau may be fully gapless. Interior dense-grid zero candidates are
    # represented compactly; physical values are never rounded before analysis.
    roots = sorted({row["q"] for row in candidates
                    if row["min_abs_E"] < 100 * GAP_ZERO_TOL * scale
                    and model.PROJECT_BZ[0] < row["q"] < model.PROJECT_BZ[1]})
    roots = deduplicate_points(roots, 2e-7)
    return {"min_abs_E": best["min_abs_E"], "zero_energy_band_gap": raw_gap,
            "raw_twice_min_abs_E": raw_gap, "q_at_minimum": best["q"],
            "classification": "gapless" if classified_zero else "gapped",
            "zero_tolerance": GAP_ZERO_TOL * scale, "dense_grid_size": nq,
            "refined_minima_count": int(len(local)), "candidates": candidates,
            "band_sign_crossings": crossings, "band_sign_crossing_count": len(crossings),
            "zero_touching_candidates": touching, "near_zero_grid_samples": int(len(near_zero_grid)),
            "crossing_method": "brentq on sorted-band sign brackets plus absolute-energy minima for tangential/paired zero touchings",
            "possible_cusp_points": roots, "majorana_max_error": majorana_error,
            "same_q_pairing_max_error": pairing_error,
            "definition": "Global zero-energy band gap = 2 min_{q,n}|E_n(q)| by Majorana symmetry on [-pi,pi]",
            "generic_global_minimum_status": "dense-grid plus continuous refinement; not an interval-arithmetic proof"}


def deduplicate_points(points, tolerance=1e-9):
    result = []
    for point in sorted(points):
        if not result or point - result[-1] > tolerance:
            result.append(float(point))
    return result


def adaptive_energy(couplings, gauge, gap):
    """Fresh eigvalsh evaluations, independently of the trapezoidal grid.

    Split at symmetry points and continuously refined potential zero-energy
    cusps. Quadrature is repeated at tighter tolerances. Its reported error is
    supplemented by mesh and quadrature-difference checks in the final bound.
    """
    points = deduplicate_points([model.PROJECT_BZ[0], model.PROJECT_BZ[1],
                                -np.pi / 2, 0., np.pi / 2] + gap["possible_cusp_points"])
    require(points[0] == model.PROJECT_BZ[0] and points[-1] == model.PROJECT_BZ[1],
            "Quadrature endpoints differ from project BZ")

    def integrand(q):
        energies = np.linalg.eigvalsh(model.h_bloch(q, couplings, gauge))
        # Deliberately independent of model.occupied_sum and energy_per_site.
        return float(energies[energies < 0].sum())

    repetitions = []
    for tolerance in (3e-10, 2e-12):
        segments = []
        for left, right in zip(points[:-1], points[1:]):
            integral, error = quad(integrand, left, right, epsabs=tolerance,
                                   epsrel=tolerance, limit=300)
            segments.append({"left": left, "right": right, "integral_F": integral,
                             "estimated_error_integral": error})
        repetitions.append({"epsabs_epsrel": tolerance,
                            "epsilon": sum(s["integral_F"] for s in segments) / (8 * np.pi),
                            "estimated_error_per_site": sum(s["estimated_error_integral"] for s in segments) / (8 * np.pi),
                            "segments": segments})
    result = dict(repetitions[-1])
    result.update({"normalization": "integral_{-pi}^{pi} F(q) dq / (8*pi)",
                   "split_points": points, "tolerance_repeats": repetitions,
                   "repeat_difference": abs(repetitions[-1]["epsilon"] - repetitions[-2]["epsilon"])})
    return result


def finite_chain_comparison(couplings, sectors, bulk):
    """Build A from real bonds and compare complete eigvalsh(iA/2) spectra.

    Ncell independent blocks contain 4*Ncell eigenvalues. Sampling 2*Ncell
    momenta over [-pi,pi) repeats each block once; each such sample receives
    half the original block weight. Both counts reproduce the real-chain E0.
    """
    rows = []
    for name, gauge in sectors.items():
        for ncell in FINITE_SIZES:
            A = np.asarray(model.finite_majorana_A(ncell, couplings, gauge))
            require(A.shape == (4*ncell, 4*ncell), "Finite A shape mismatch")
            require(np.max(abs(A.imag)) < 1e-13, "Finite A is not real")
            antisymmetric_error = float(np.max(abs(A+A.T)))
            require(antisymmetric_error < 1e-13, "Finite A is not antisymmetric")
            real_values = np.linalg.eigvalsh(.5j*A)
            allowed_q = model.finite_momenta(ncell)
            require(len(allowed_q) == ncell, "Independent momentum count mismatch")
            block_values = np.sort(physical_bands(allowed_q, couplings, gauge).ravel())
            spectrum_error = float(np.max(abs(real_values-block_values)))
            scale = max(1., float(np.max(abs(real_values))))
            require(spectrum_error < 2e-11*scale, f"Finite/Bloch spectrum mismatch at {name}, {ncell}")
            block_error = completeness_error = 0.
            if ncell in (5, 8, 17):
                basis = []
                for point in allowed_q:
                    V = model.finite_fourier_matrix(ncell, point)
                    basis.append(V)
                    block_error = max(block_error, float(np.max(abs(V.conj().T@(.5j*A)@V
                                                                   -model.h_bloch(point, couplings, gauge)))))
                Vfull = np.concatenate(basis, axis=1)
                completeness_error = float(np.max(abs(Vfull.conj().T@Vfull-np.eye(4*ncell))))
                require(max(block_error, completeness_error) < 2e-11*scale, "Finite Fourier basis/block mismatch")
            E0 = -.5*float(np.sum(abs(real_values)))
            epsilon = E0/(4*ncell)
            epsilon_blocks = float(np.sum(np.minimum(block_values, 0.)))/(4*ncell)
            doubled_q = np.linspace(*model.PROJECT_BZ, 2*ncell, endpoint=False)
            doubled_F = model.occupied_sum(physical_bands(doubled_q, couplings, gauge))
            epsilon_double = float(doubled_F.sum())/(8*ncell)
            discrete_error = max(abs(epsilon-epsilon_blocks), abs(epsilon-epsilon_double))
            require(discrete_error < 2e-12*scale, "Discrete Majorana/double-cover energy factor mismatch")
            rows.append({"sector": name, "Ncell": ncell, "Nsite": 4*ncell,
                         "independent_momentum_blocks": ncell, "independent_eigenvalues": 4*ncell,
                         "project_interval_samples": 2*ncell, "weight_per_sample_in_Ecell": 1./(2*ncell),
                         "weight_per_independent_block_in_Ecell": 1./ncell,
                         "antisymmetric_max_error": antisymmetric_error,
                         "full_spectrum_max_error": spectrum_error,
                         "fourier_block_max_error": block_error if ncell in (5,8,17) else None,
                         "fourier_completeness_max_error": completeness_error if ncell in (5,8,17) else None,
                         "E0_finite": E0, "epsilon_finite": epsilon,
                         "epsilon_discrete_bloch": epsilon_blocks, "epsilon_discrete_double_cover": epsilon_double,
                         "finite_discrete_bloch_error": discrete_error,
                         "epsilon_bulk_quad": bulk[name], "finite_minus_bulk": epsilon-bulk[name]})
    conclusions = {}
    reduced_checks = {}
    for name, gauge in sectors.items():
        selected = [r for r in rows if r["sector"] == name and r["Ncell"] in FINITE_CONVERGENCE_SIZES]
        residuals = [abs(r["finite_minus_bulk"]) for r in selected]
        noise = 2e-11*max(1., abs(bulk[name]))
        decreasing = residuals[-1] <= max(noise, .4*residuals[0])
        require(decreasing, f"Finite energies do not approach bulk for {name}")
        tolerance = 5e-4*max(1., abs(bulk[name]))
        require(residuals[-1] < tolerance, f"Finite-chain normalization unresolved for {name}")
        conclusions[name] = {"sizes": list(FINITE_CONVERGENCE_SIZES), "absolute_bulk_residuals": residuals,
                             "largest_chain_bulk_residual": residuals[-1], "finite_size_acceptance_tolerance": tolerance,
                             "converges_toward_bulk": decreasing, "normalization_factor_2_or_4_excluded": True}
        q_check = np.linspace(*model.PROJECT_BZ, 257)
        period_error = float(np.max(abs(physical_bands(q_check+np.pi, couplings, gauge)
                                       -physical_bands(q_check, couplings, gauge))))
        require(period_error < 2e-11, "Reduced spectral interval has not been justified")
        def occupied(point):
            E = np.linalg.eigvalsh(model.h_bloch(point, couplings, gauge))
            return float(E[E < 0].sum())
        reduced_integral, reduced_error = quad(occupied, -np.pi/2., np.pi/2., epsabs=2e-11,
                                               epsrel=2e-11, points=[0.], limit=400)
        epsilon_reduced = reduced_integral/(4*np.pi)
        require(abs(epsilon_reduced-bulk[name]) < 2e-9, "Reduced/project interval normalization mismatch")
        reduced_checks[name] = {"reduced_interval": [-np.pi/2., np.pi/2.],
                                "formula": "epsilon=integral_reduced F/(4*pi)",
                                "epsilon_reduced": epsilon_reduced, "estimated_error": reduced_error/(4*np.pi),
                                "project_vs_reduced_error": abs(epsilon_reduced-bulk[name]),
                                "spectral_period_pi_max_error": period_error}
    return {"status": "PASS", "rows": rows, "convergence": conclusions,
            "max_spectrum_error": max(r["full_spectrum_max_error"] for r in rows),
            "finite_energy_definition": "h_real=iA/2; E0=-sum(abs(eigvalsh(h_real)))/2; epsilon=E0/(4*Ncell)",
            "allowed_q": "pi*ell/Ncell, Ncell independent representatives; translation 2*a1",
            "double_coverage": {"status": "PASS", "reduced_interval_checks": reduced_checks,
                                 "counting": "4*Ncell sites; Ncell independent 4x4 blocks; 2*Ncell samples on [-pi,pi) with weight 1/(2*Ncell) in Ecell"}}


def benchmark(couplings, sectors, repeats=5):
    require(repeats >= 5, "Benchmark requires at least five repeats")
    rows = []

    def measure(label, sector, action, ncell=None):
        action()  # Each operation gets its own warm-up.
        measurements = []
        for _ in range(repeats):
            start = time.perf_counter()
            action()
            measurements.append(time.perf_counter() - start)
        row = {"operation": label, "sector": sector, "repeats": repeats,
               "median_seconds": float(np.median(measurements)),
               "min_seconds": float(min(measurements)), "max_seconds": float(max(measurements)),
               "samples_seconds": measurements}
        if ncell is not None:
            row["Ncell"] = ncell
        rows.append(row)

    q1000 = np.linspace(*model.PROJECT_BZ, 1000)
    for name, gauge in sectors.items():
        matrices = [model.h_bloch(q, couplings, gauge) for q in q1000]
        measure("1000 H(q) evaluations", name,
                lambda: [model.h_bloch(q, couplings, gauge) for q in q1000])
        measure("1000 eigvalsh 4x4 (H prebuilt)", name,
                lambda: [np.linalg.eigvalsh(H) for H in matrices])
        for nq in (2001, 16001):
            q = np.linspace(*model.PROJECT_BZ, nq)
            measure(f"physical bands {nq} (H build + eigvalsh)", name,
                    lambda q=q: model.bands_sorted(q, couplings, gauge))
        q_integral = np.linspace(*model.PROJECT_BZ, 16001)
        E_integral = physical_bands(q_integral, couplings, gauge)
        measure("F(q) 16001 (H build + eigvalsh + negative sum)", name,
                lambda: model.occupied_sum(model.bands_sorted(q_integral, couplings, gauge)))
        measure("integral 16001 (physical bands precomputed)", name,
                lambda: model.energy_per_site(q_integral, E_integral))
        measure("energy 16001 (H build + eigvalsh + integral)", name,
                lambda: model.energy_per_site(q_integral, model.bands_sorted(q_integral, couplings, gauge)))
        for ncell in FINITE_SIZES:
            measure("finite A build + full eigvalsh(iA)", name,
                    lambda ncell=ncell: np.linalg.eigvalsh(1j * model.finite_majorana_A(ncell, couplings, gauge)),
                    ncell=ncell)
    measure("adversarial suite 1000 cases (both sectors, 500 mixed-sign cases)", "OLD + T-REVERSED",
            lambda: model.randomized_suite(random_cases=1000, seed=model.SEED))
    return {"clock": "time.perf_counter", "warmup": "one untimed call per operation and sector",
            "repeats": repeats, "rows": rows, "thread_environment":
            {key: os.environ.get(key) for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS")},
            "platform": platform.platform(), "processor": platform.processor(),
            "interpretation": "Wall-clock timings on this machine; no precision reduction; finite timings include matrix build and diagonalization"}


def configure_axis(ax):
    ax.set_xlim(*model.PROJECT_BZ)
    ax.set_xticks([-np.pi, -np.pi / 2, 0, np.pi / 2, np.pi])
    ax.set_xticklabels(["−π", "−π/2", "0", "π/2", "π"])
    ax.set_xlabel("q = k·a₁")
    ax.grid(alpha=.18, linewidth=.6)
    ax.spines[["top", "right"]].set_visible(False)


def draw_comparison(out, q, sorted_values, tracked, couplings):
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11, "axes.titlesize": 13})
    colors = ["#21498d", "#d5752d", "#21866b", "#854898"]
    fig, axes = plt.subplots(1, 2, figsize=(12., 4.4), sharey=True, layout="constrained")
    scale = max(float(np.max(abs(v))) for v in sorted_values.values())
    for ax, (name, values) in zip(axes, tracked.items()):
        for band in range(4):
            ax.plot(q, values[:, band], color=colors[band], linewidth=1.5,
                    linestyle="--" if band in (1, 3) else "-")
        ax.axhline(0, color="#555555", linewidth=.7)
        ax.set_ylim(-1.12 * max(scale, 1.), 1.12 * max(scale, 1.))
        configure_axis(ax)
        ax.set_title("OLD / zero-flux" if name == "OLD" else "T bonds reversed")
    axes[0].set_ylabel("Band energy E")
    fig.suptitle("Kx = Ky = Kz = Kt = Kw = 1" if couplings == (1.,) * 5
                 else "Couplings (x, y, z, t, w): " + str(couplings), fontsize=11)
    fig.savefig(out / "dispersion_old_vs_t_reversed.png", dpi=300)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(8.7, 4.7), layout="constrained")
    for (name, values), color in zip(sorted_values.items(), ("#21498d", "#bc4831")):
        ax.plot(q, model.occupied_sum(values), color=color, linewidth=2., label=name)
    configure_axis(ax)
    ax.set_ylabel("F(q) = sum of negative band energies")
    ax.set_title("Occupied-band sum on the project Brillouin interval")
    ax.legend(frameon=False)
    fig.savefig(out / "Fq_old_vs_t_reversed.png", dpi=300)
    plt.close(fig)


def write_csv(path, rows, hashes, metadata):
    encoded_hashes = json.dumps(hashes, sort_keys=True, separators=(",", ":"))
    enriched = [dict(row, script_sha256_json=encoded_hashes,
                     metadata_json=json.dumps(metadata, sort_keys=True, separators=(",", ":"))) for row in rows]
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(enriched[0]))
        writer.writeheader()
        writer.writerows(enriched)


def write_json(path, data):
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")


def benchmark_markdown(data):
    lines = ["# T-REVERSED numerical benchmark", "", "One warm-up and at least five timed repeats per operation. Units: seconds.", "",
             "| Sector | Operation | Ncell | Median | Min | Max |", "|---|---|---:|---:|---:|---:|"]
    for row in data["rows"]:
        lines.append(f"| {row['sector']} | {row['operation']} | {row.get('Ncell', '—')} | "
                     f"{row['median_seconds']:.7f} | {row['min_seconds']:.7f} | {row['max_seconds']:.7f} |")
    lines.extend(["", "Finite-chain timings include building A and diagonalizing iA. The dedicated F(q) timing includes matrix construction and eigvalsh. The 1000-case randomized suite is repeated in full.", "",
                  "Thread environment: `" + json.dumps(data["thread_environment"], sort_keys=True) + "`.", "",
                  "Script SHA-256: `" + json.dumps(data["script_sha256"], sort_keys=True) + "`.", ""])
    return "\n".join(lines)


def old_regression(reference, couplings, epsilon, gap, uncertainty):
    if reference is None or not reference.is_file():
        return {"status": "NOT_RUN", "reason": "No legacy numerical reference supplied"}
    previous = json.loads(reference.read_text(encoding="utf-8"))
    summary = previous.get("summary", {})
    previous_K = tuple(summary.get(name) for name in COUPLING_NAMES)
    if previous_K != couplings:
        return {"status": "NOT_APPLICABLE", "reason": "Reference couplings differ from requested couplings"}
    reference_epsilon, reference_gap = float(summary["epsilon_old"]), float(summary["gap_old"])
    errors = {"epsilon_old": abs(epsilon - reference_epsilon),
              "zero_energy_band_gap_old": abs(gap - reference_gap)}
    require(errors["epsilon_old"] <= max(1e-9, 5 * uncertainty), "OLD energy regression changed")
    require(errors["zero_energy_band_gap_old"] < 2e-7, "OLD gap regression changed")
    return {"status": "PASS", "source": str(reference.resolve()),
            "source_sha256": hashlib.sha256(reference.read_bytes()).hexdigest(),
            "previous_epsilon_old": reference_epsilon, "previous_zero_energy_band_gap_old": reference_gap,
            "absolute_errors": errors,
            "scope": "Numerical regression only; previous CSV/JSON are not reused as T physical input"}


def isotropic_independent_checks(q, values, energies, couplings):
    """Validate the independently derived isotropic polynomial and elliptic E.

    OLD follows S=5,R=4 cos²(q); T follows S=1+8 sin²(q),R=0.
    Constants in these identities are algebraic derivations from the core's
    geometry-based polynomial, never stored numerical energy benchmarks.
    """
    if couplings != (1.,) * 5:
        return {"status": "NOT_APPLICABLE", "reason": "Couplings are not isotropic unit couplings"}
    invariants = {name: np.asarray([model.reduced_invariants(x, couplings, gauge) for x in q])
                  for name, gauge in (("OLD", model.ORIGINAL_GAUGE), ("T-REVERSED", model.T_REVERSED_GAUGE))}
    errors = {
        "old_S": float(np.max(abs(invariants["OLD"][:, 0] - 5.))),
        "old_R": float(np.max(abs(invariants["OLD"][:, 1] - 4. * np.cos(q)**2))),
        "t_S": float(np.max(abs(invariants["T-REVERSED"][:, 0] - (1. + 8.*np.sin(q)**2)))),
        "t_R": float(np.max(abs(invariants["T-REVERSED"][:, 1]))),
    }
    require(max(errors.values()) < 3e-12, "Isotropic polynomial invariant identities failed")
    t_radius = np.sqrt(1. + 8.*np.sin(q)**2)
    t_exact = np.column_stack((-t_radius, -t_radius, t_radius, t_radius))
    t_band_error = float(np.max(abs(values["T-REVERSED"] - t_exact)))
    require(t_band_error < 2e-11, "Derived isotropic T roots fail eigvalsh comparison")
    old_exact_F = -np.sqrt(5.+4.*np.cos(q))-np.sqrt(5.-4.*np.cos(q))
    t_exact_F = -2.*t_radius
    F_errors = {"OLD": float(np.max(abs(model.occupied_sum(values["OLD"])-old_exact_F))),
                "T-REVERSED": float(np.max(abs(model.occupied_sum(values["T-REVERSED"])-t_exact_F)))}
    require(max(F_errors.values()) < 2e-11, "Isotropic occupied sums fail exact formula")
    # OLD's half-angle reduction and T's sin² integral independently give the
    # same elliptic parameter, proving equality without assuming equal bands.
    old_parameter = 8./9.
    t_parameter = 8./9.
    elliptic_old = -3./np.pi * float(ellipe(old_parameter))
    elliptic_t = -3./np.pi * float(ellipe(t_parameter))
    energy_errors = {"OLD": abs(energies["OLD"]-elliptic_old),
                     "T-REVERSED": abs(energies["T-REVERSED"]-elliptic_t)}
    require(max(energy_errors.values()) < 2e-11, "Independent elliptic energy check failed")
    # Explicit change-of-variable identity supplies an independent energy
    # equality argument. This is not same-q spectral equality.
    F_old_2q = model.occupied_sum(physical_bands(2.*q, couplings, model.ORIGINAL_GAUGE))
    F_t_shift = model.occupied_sum(physical_bands(q+np.pi/2., couplings, model.T_REVERSED_GAUGE))
    momentum_identity_error = float(np.max(abs(F_old_2q - .5*(t_exact_F+F_t_shift))))
    require(momentum_identity_error < 2e-11, "Isotropic occupied-sum momentum identity failed")
    return {"status": "PASS", "elliptic_energy_old": elliptic_old, "elliptic_energy_t": elliptic_t,
            "elliptic_parameter_old": old_parameter, "elliptic_parameter_t": t_parameter,
            "invariant_max_errors": errors, "t_exact_spectrum_max_error": t_band_error,
            "occupied_sum_max_errors": F_errors, "elliptic_vs_adaptive_errors": energy_errors,
            "momentum_identity": "F_OLD(2q)=[F_T(q)+F_T(q+pi/2)]/2",
            "momentum_identity_max_error": momentum_identity_error,
            "exact_bulk_energy_equality": True,
            "exact_t_min_abs_E": 1., "exact_t_zero_energy_band_gap": 2.,
            "exact_old_min_abs_E": 1., "exact_old_zero_energy_band_gap": 2.,
            "energy_formula": "epsilon_OLD=epsilon_T=-3 E(m=8/9)/pi; E is the complete elliptic integral with parameter m",
            "derivation": "OLD F=-sqrt(5+4cos q)-sqrt(5-4cos q), half-angle substitution; T F=-2sqrt(1+8sin²q); both reduce to E(8/9)"}


def consolidate_previous(couplings, comparison_rows, legacy_reference):
    """Recompute each supplied prior sector from a fresh copy of ORIGINAL.

    Previous JSON numbers are regression targets only. They never determine
    the new T bonds, the Hamiltonian or its integration normalization.
    """
    root = Path(__file__).resolve().parent
    candidates = (("Z-REVERSED", ("u42z", "u31z"), legacy_reference, "epsilon_z", "gap_z"),
                  ("W-REVERSED", ("u43w", "u12w"), root/"sources/legacy_comparison_w_reversed.json",
                   "epsilon_w", "zero_energy_band_gap_w"))
    rows = [dict(row) for row in comparison_rows]
    gauges = {"OLD": model.ORIGINAL_GAUGE, "T-REVERSED": model.T_REVERSED_GAUGE}
    audits, omitted = {}, {"upper-horizontal-reversed": "No compatible audited numerical artifact supplied"}
    for name, changed_keys, reference, energy_key, gap_key in candidates:
        if reference is None or not reference.is_file():
            omitted[name] = "Prior audited comparison artifact absent"
            continue
        previous = json.loads(reference.read_text(encoding="utf-8"))
        summary = previous.get("summary", {})
        if tuple(summary.get(key) for key in COUPLING_NAMES) != couplings or energy_key not in summary:
            omitted[name] = "Prior artifact couplings differ or energy result is absent"
            continue
        gauge = dict(model.ORIGINAL_GAUGE)
        for key in changed_keys:
            gauge[key] *= -1
        require({key for key in gauge if gauge[key] != model.ORIGINAL_GAUGE[key]} == set(changed_keys),
                "Previous-sector gauge reconstruction contaminated")
        gauges[name] = gauge
        q = np.linspace(*model.PROJECT_BZ, CONVERGENCE_NQ[-1])
        mesh_energy = float(model.energy_per_site(q, physical_bands(q, couplings, gauge)))
        gap = gap_analysis(couplings, gauge)
        adaptive = adaptive_energy(couplings, gauge, gap)
        finite = finite_chain_comparison(couplings, {name: gauge}, {name: adaptive["epsilon"]})
        energy_error = abs(adaptive["epsilon"]-float(summary[energy_key]))
        gap_error = abs(gap["zero_energy_band_gap"]-float(summary[gap_key]))
        # Legacy values only serve regression when the stored BZ and energy
        # conventions explicitly match. Fresh values remain comparable.
        meta = previous.get("metadata", {})
        prior_interval = previous.get("PROJECT_BZ", [meta.get("q_min"), meta.get("q_max")])
        interval_matches = len(prior_interval) == 2 and all(
            point is not None and abs(point-target) < 1e-12
            for point, target in zip(prior_interval, model.PROJECT_BZ))
        normalization_matches = (meta.get("N_sites_per_cell") == 4
                                 and meta.get("normalization_formula") == "epsilon = integral[-pi,pi] F(q) dq / (8*pi)")
        conventions_match = interval_matches and normalization_matches
        if conventions_match:
            require(energy_error < 1e-9 and gap_error < 2e-7, f"Recomputed {name} differs from supplied compatible regression")
        rows.append({"sector": name, "epsilon": adaptive["epsilon"],
                     "absolute_energy_uncertainty": max(ENERGY_FLOOR, 2*abs(mesh_energy-adaptive["epsilon"]),
                                                          10*adaptive["estimated_error_per_site"]),
                     "min_abs_E": gap["min_abs_E"], "zero_energy_band_gap": gap["zero_energy_band_gap"],
                     "q_min": float(model.PROJECT_BZ[0]), "q_max": float(model.PROJECT_BZ[1]),
                     "status": "recomputed and finite-chain validated"})
        audits[name] = {"gauge": gauge, "adaptive_quadrature": adaptive, "gap": gap,
                       "finite_chain_validation": finite, "legacy_reference": str(reference.resolve()),
                       "legacy_reference_sha256": hashlib.sha256(reference.read_bytes()).hexdigest(),
                       "legacy_conventions_match": conventions_match,
                       "legacy_energy_difference": energy_error, "legacy_gap_difference": gap_error}
    for row in rows:
        gauge = gauges[row["sector"]]
        row.update({"changed_u": ",".join(key for key in model.ORIGINAL_GAUGE
                                          if gauge[key] != model.ORIGINAL_GAUGE[key]) or "none",
                    "flux_pattern": ",".join(f"{v:+d}" for v in model.fluxes(gauge).values()),
                    "Delta_epsilon_relative_to_OLD": row["epsilon"]-comparison_rows[0]["epsilon"],
                    "spectral_classification": "gapped" if row["min_abs_E"] > model.GAP_TOLERANCE else "gapless"})
    return {"status": "PASS", "rows": rows, "recomputed_sectors": audits,
            "compatibility": "All displayed sectors freshly recomputed on [-pi,pi], epsilon=integral F/(8pi), full finite iA/2 spectra checked",
            "omitted_sectors": omitted}


def draw_consolidated(out, rows):
    fig, axes = plt.subplots(1, 2, figsize=(12., 4.4), layout="constrained")
    names = [row["sector"] for row in rows]
    colors = ["#21498d", "#bc4831", "#318268", "#854898"]
    for ax, key, title, ylabel in ((axes[0], "epsilon", "Ground-state energy per site", "ε"),
                                    (axes[1], "zero_energy_band_gap", "Zero-energy band gap", "Full gap")):
        heights = [row[key] for row in rows]
        bars = ax.bar(names, heights, color=colors, width=.58)
        ax.bar_label(bars, labels=[f"{value:.6f}" for value in heights], padding=5, fontsize=10)
        ax.set_title(title)
        ax.set_ylabel(ylabel)
        ax.axhline(0., color="#555555", linewidth=.7)
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="y", alpha=.15)
        ax.margins(y=.20)
    fig.suptitle("Compatible sectors recomputed on q ∈ [−π, π]", fontsize=12)
    fig.savefig(out / "all_gauge_sectors_comparison.png", dpi=300)
    plt.close(fig)


def run_analysis(out, couplings=(1.,) * 5, benchmark_repeats=5, regression_reference=None):
    couplings = tuple(float(x) for x in couplings)
    require(len(couplings) == 5 and np.isfinite(couplings).all(), "Five finite couplings required")
    require(tuple(model.PROJECT_BZ) == (-np.pi, np.pi), "Principal BZ must remain [-pi,pi]")
    initial_hashes = script_hashes()
    initial_sources = model.metadata(couplings=couplings)
    out.mkdir(parents=True, exist_ok=True)
    sectors = {"OLD": model.ORIGINAL_GAUGE, "T-REVERSED": model.T_REVERSED_GAUGE}
    changed = {key for key in model.ORIGINAL_GAUGE if model.ORIGINAL_GAUGE[key] != model.T_REVERSED_GAUGE[key]}
    require(changed == {"u13t", "u24t"}, "Gauge contamination: changes other than t")
    for key in ("u12x", "u34x", "u23y", "u42z", "u31z", "u43w", "u13t", "u24t"):
        require(model.T_REVERSED_GAUGE[key] == 1, f"Forbidden change of {key}")
    print("T-REVERSED: convergence and continuous gap search on [-pi,pi]", flush=True)
    convergence = []
    for nq in CONVERGENCE_NQ:
        q = np.linspace(*model.PROJECT_BZ, nq)
        energies = {name: float(model.energy_per_site(q, physical_bands(q, couplings, gauge)))
                    for name, gauge in sectors.items()}
        convergence.append({"Nq": nq, "epsilon_old": energies["OLD"], "epsilon_t": energies["T-REVERSED"],
                            "Delta_epsilon": energies["T-REVERSED"] - energies["OLD"]})
        print(f"  Nq={nq:5d}: OLD={energies['OLD']:.13f}; T={energies['T-REVERSED']:.13f}", flush=True)
    # Each reported mesh error derives from an adjacent refinement. The first
    # row uses the next grid because it has no previous refinement.
    for index, record in enumerate(convergence):
        neighbor = convergence[index-1 if index else 1]
        for suffix in ("old", "t"):
            record[f"error_estimate_{suffix}"] = max(ENERGY_FLOOR,
                2.*abs(record[f"epsilon_{suffix}"]-neighbor[f"epsilon_{suffix}"]))
    gaps = {name: gap_analysis(couplings, gauge) for name, gauge in sectors.items()}
    adaptive = {name: adaptive_energy(couplings, gauge, gaps[name]) for name, gauge in sectors.items()}
    estimates = {}
    for name, key in (("OLD", "epsilon_old"), ("T-REVERSED", "epsilon_t")):
        adjacent = [abs(b[key] - a[key]) for a, b in zip(convergence[:-1], convergence[1:])]
        quad_mesh = abs(convergence[-1][key] - adaptive[name]["epsilon"])
        bound = max(ENERGY_FLOOR, 2 * max(adjacent[-2:]),
                    2 * quad_mesh, 10 * adaptive[name]["estimated_error_per_site"],
                    2 * adaptive[name]["repeat_difference"])
        estimates[name] = {"consecutive_mesh_differences": adjacent,
                           "last_mesh_vs_adaptive": quad_mesh, "conservative_absolute_uncertainty": bound,
                           "uncertainty_kind": "empirical conservative estimate, not a rigorous mathematical bound",
                           "resolved_decimal_places_from_uncertainty": max(0, int(np.floor(-np.log10(2 * bound))))}
    old_energy, t_energy = (adaptive[name]["epsilon"] for name in sectors)
    delta = t_energy - old_energy
    comparison_uncertainty = max(model.ENERGY_TOLERANCE,
                                 sum(v["conservative_absolute_uncertainty"] for v in estimates.values()))
    favored = "OLD" if delta > comparison_uncertainty else "T-REVERSED" if delta < -comparison_uncertainty else "DEGENERATE / UNRESOLVED NUMERICALLY"
    row = dict(zip(COUPLING_NAMES, couplings))
    row.update({"epsilon_old": old_energy, "epsilon_t": t_energy, "Delta_epsilon": delta,
                "min_abs_E_old": gaps["OLD"]["min_abs_E"], "min_abs_E_t": gaps["T-REVERSED"]["min_abs_E"],
                "zero_energy_band_gap_old": gaps["OLD"]["zero_energy_band_gap"],
                "zero_energy_band_gap_t": gaps["T-REVERSED"]["zero_energy_band_gap"],
                "gap_old": gaps["OLD"]["zero_energy_band_gap"], "gap_t": gaps["T-REVERSED"]["zero_energy_band_gap"],
                "spectral_character_old": gaps["OLD"]["classification"],
                "spectral_character_t": gaps["T-REVERSED"]["classification"],
                "epsilon_old_uncertainty": estimates["OLD"]["conservative_absolute_uncertainty"],
                "epsilon_t_uncertainty": estimates["T-REVERSED"]["conservative_absolute_uncertainty"],
                "Delta_epsilon_uncertainty": comparison_uncertainty,
                "favored_sector": favored, "energy_estimator": "independent adaptive quadrature"})
    regression = old_regression(regression_reference, couplings, old_energy,
                                gaps["OLD"]["zero_energy_band_gap"], comparison_uncertainty)
    print("T-REVERSED: independent real-space finite-chain validation", flush=True)
    finite = finite_chain_comparison(couplings, sectors, {name: adaptive[name]["epsilon"] for name in sectors})
    qplot = np.linspace(*model.PROJECT_BZ, 2001)
    values = {name: physical_bands(qplot, couplings, gauge) for name, gauge in sectors.items()}
    tracked = {name: tracked_bands(qplot, couplings, gauge) for name, gauge in sectors.items()}
    tracking_errors = {name: float(np.max(abs(np.sort(tracked[name], axis=1) - values[name]))) for name in sectors}
    require(max(tracking_errors.values()) < 2e-11 * max(1., max(np.abs(couplings))), "Tracking changed physical eigenvalues")
    draw_comparison(out, qplot, values, tracked, couplings)
    special_q = (("-pi", -np.pi), ("-pi/2", -np.pi / 2), ("-pi/4", -np.pi / 4),
                 ("0", 0.), ("pi/4", np.pi / 4), ("pi/2", np.pi / 2), ("pi", np.pi))
    special_points = {name: {label: np.linalg.eigvalsh(model.h_bloch(q, couplings, gauge)).tolist()
                             for label, q in special_q} for name, gauge in sectors.items()}
    analytic_points = {name: {label: model.analytic_spectrum(q, couplings, gauge).tolist()
                              for label, q in special_q} for name, gauge in sectors.items()}
    analytic_errors = {name: max(float(np.max(abs(np.asarray(special_points[name][label])
                                                 - analytic_points[name][label]))) for label, _ in special_q)
                       for name in sectors}
    require(max(analytic_errors.values()) < 2e-11 * max(1., max(np.abs(couplings))), "Closed-form special-point spectrum mismatch")
    isotropic = isotropic_independent_checks(qplot, values, {"OLD": old_energy, "T-REVERSED": t_energy}, couplings)
    base_rows = [{"sector": name, "epsilon": adaptive[name]["epsilon"],
                  "absolute_energy_uncertainty": estimates[name]["conservative_absolute_uncertainty"],
                  "min_abs_E": gaps[name]["min_abs_E"], "zero_energy_band_gap": gaps[name]["zero_energy_band_gap"],
                  "q_min": float(model.PROJECT_BZ[0]), "q_max": float(model.PROJECT_BZ[1]),
                  "status": "recomputed and finite-chain validated"} for name in sectors]
    consolidated = consolidate_previous(couplings, base_rows, regression_reference)
    if consolidated["status"] == "PASS":
        draw_consolidated(out, consolidated["rows"])
    print("T-REVERSED: warmed performance benchmark (five or more repeats)", flush=True)
    timings = benchmark(couplings, sectors, benchmark_repeats)
    # Provenance is collected only after all artifacts' computations. Refuse to
    # stamp fresh hashes on a run whose source files changed while it executed.
    hashes = script_hashes()
    require(hashes == initial_hashes, "Scientific scripts changed during this run; rerun before publishing results")
    metadata = model.metadata(couplings=couplings)
    metadata.update({"script_sha256": hashes, "random_cases": 1000,
                     "timestamp": metadata["generated_at_utc"],
                     "Majorana_convention": metadata["majorana_convention"]})
    for kind in ("source_sha256", "local_source_sha256"):
        require(initial_sources[kind] == metadata[kind], f"Sources changed during comparison: {kind}")
    timings.update({"metadata": metadata, "script_sha256": hashes, "couplings": dict(zip(COUPLING_NAMES, couplings)),
                    "PROJECT_BZ": list(model.PROJECT_BZ)})
    results = {"metadata": metadata, "script_sha256": hashes,
               "couplings": dict(zip(COUPLING_NAMES, couplings)), "PROJECT_BZ": list(model.PROJECT_BZ),
               "summary": row, "convergence": convergence, "convergence_assessment": estimates,
               "convergence_row_error_definition": "max(1e-11,2*absolute difference to preceding refinement); first row uses next refinement; empirical estimates",
               "adaptive_quadrature": adaptive, "gaps": gaps,
               "effective_energy_comparison_tolerance": comparison_uncertainty,
               "finite_chain_validation": finite, "OLD_regression": regression,
               "tracking_spectrum_errors": tracking_errors, "tracking_used_for_physics": False,
               "tracking_permutations_per_step": len(PERMUTATIONS), "special_points": special_points,
               "special_points_analytic": analytic_points, "special_points_analytic_max_errors": analytic_errors,
               "isotropic": isotropic, "consolidated_previous_sectors": consolidated,
               "negative_band_counts_on_plot_grid": {name: sorted(set(np.sum(E < -1e-10, axis=1).tolist())) for name, E in values.items()},
               "figures_dpi": 300,
               "normalization": "Nsite=4*Ncell; iA=2h; E0=-sum|eig(iA)|/4; epsilon=integral F/(8pi) on [-pi,pi]"}
    write_csv(out / "ground_state_energy_old_vs_t.csv", [row], hashes, metadata)
    write_csv(out / "energy_convergence_old_vs_t.csv", convergence, hashes, metadata)
    write_csv(out / "finite_chain_old_vs_t.csv", finite["rows"], hashes, metadata)
    if consolidated["status"] == "PASS":
        write_csv(out / "all_gauge_sectors_comparison.csv", consolidated["rows"], hashes, metadata)
    write_json(out / "comparison_t_reversed.json", results)
    write_json(out / "benchmark_t_reversed.json", timings)
    (out / "benchmark_t_reversed.md").write_text(benchmark_markdown(timings), encoding="utf-8")
    print("\nT-REVERSED NUMERICAL COMPARISON\n================================")
    print("Couplings: " + ", ".join(f"{name}={value:g}" for name, value in zip(COUPLING_NAMES, couplings)))
    print("BZ: q in [-pi, pi]\nChanged keys: " + ", ".join(sorted(changed)))
    for name, key in (("OLD", "epsilon_old"), ("T-REVERSED", "epsilon_t")):
        print(f"{name}: min_abs_E={gaps[name]['min_abs_E']:.12g}; "
              f"zero_energy_band_gap={gaps[name]['zero_energy_band_gap']:.12g}; epsilon={row[key]:.15f}")
    print(f"Delta epsilon = {delta:.15g} +/- {comparison_uncertainty:.3g}\nFavored sector: {favored}")
    print(f"Finite-chain validation: PASS; max spectrum error={finite['max_spectrum_error']:.3g}")
    print("Energy normalization: PASS (independent real-space finite-chain convergence)")
    print(f"Benchmark: {len(timings['rows'])} operations; {benchmark_repeats} repeats after warm-up")
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--couplings", type=float, nargs=5, default=[1.] * 5, metavar="K")
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--benchmark-repeats", type=int, default=5)
    parser.add_argument("--regression-reference", type=Path,
                        default=Path(__file__).resolve().parent / "sources" / "legacy_comparison_z_reversed.json")
    args = parser.parse_args()
    run_analysis(args.output_dir, args.couplings, args.benchmark_repeats, args.regression_reference)


if __name__ == "__main__":
    main()
