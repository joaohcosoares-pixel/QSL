#!/usr/bin/env python3
"""Independent 64-state matter-operator calibration of the T/OLD energy scale.

Twelve Majoranas are represented by explicit Pauli tensor products. The
many-body Hamiltonian uses the literal real-space PDF Eq.(11) transcription
below, independently of the production geometry and momentum formulas.

This checks the unprojected matter Hilbert space, not the physical spin
projection. Run with OPENBLAS_NUM_THREADS=1 python3 -W error this_file.py.
"""
from __future__ import annotations

import argparse
from itertools import product
import json
from pathlib import Path

import numpy as np

import quadrupolar_spin_liquid_t_reversed as model


# Independent literal source evidence, PDF pages 3-4 (Eq.11) and 5 (Eq.15).
# start, end, coupling index in (Kx,Ky,Kz,Kt,Kw), OLD u, end-minus-start cell.
# This intentional redundancy belongs to the auditor, never to production.
PDF_REAL_BONDS = (
    (1, 2, 0, +1, 0),
    (3, 4, 0, +1, 0),
    (2, 3, 1, +1, 0),
    (4, 1, 1, -1, 0),
    (4, 2, 2, +1, 0),
    (3, 1, 2, +1, +1),
    (4, 3, 4, +1, -1),
    (1, 2, 4, -1, -1),
    (1, 3, 3, -1, 0),
    (2, 4, 3, -1, +1),
)
COUPLING_CASES = (
    (1., 1., 1., 1., 1.),
    (1.13, .79, 1.31, .92, .67),
    (-1.21, .71, 1.41, -.93, .53),
)
N_CELL = 3
N_SITE = 4 * N_CELL
N_FERMION = N_SITE // 2
DIMENSION = 2 ** N_FERMION
TOLERANCE = 2e-12


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def tensor_product(factors):
    result = np.ones((1, 1), dtype=complex)
    for factor in factors:
        result = np.kron(result, factor)
    return result


def explicit_majoranas():
    identity = np.eye(2, dtype=complex)
    pauli_x = np.array([[0, 1], [1, 0]], dtype=complex)
    pauli_y = np.array([[0, -1j], [1j, 0]], dtype=complex)
    pauli_z = np.diag([1., -1.]).astype(complex)
    result = []
    for fermion in range(N_FERMION):
        for local in (pauli_x, pauli_y):
            result.append(tensor_product(
                [pauli_z] * fermion + [local] +
                [identity] * (N_FERMION - fermion - 1)))
    return tuple(result)


def source_many_body_hamiltonian(majoranas, couplings, sector):
    """Literal once-counted real bond sum, without importing model.GEOMETRY."""
    require(sector in ("OLD", "T-REVERSED"), "Unknown audit sector")
    matrix = np.zeros((DIMENSION, DIMENSION), dtype=complex)
    for cell in range(N_CELL):
        for start, end, coupling_index, original_u, shift in PDF_REAL_BONDS:
            used_u = (-original_u if sector == "T-REVERSED" and
                      coupling_index == 3 else original_u)
            i = 4 * cell + start - 1
            j = 4 * ((cell + shift) % N_CELL) + end - 1
            matrix += (1j * couplings[coupling_index] * used_u *
                       (majoranas[i] @ majoranas[j]))
    return matrix


def run_explicit_majorana_audit():
    initial_provenance = model.metadata()
    majoranas = explicit_majoranas()
    identity = np.eye(DIMENSION)
    clifford_error = 0.
    for i, first in enumerate(majoranas):
        for j, second in enumerate(majoranas):
            expected = 2 * identity if i == j else np.zeros_like(identity)
            clifford_error = max(clifford_error, float(np.max(abs(
                first @ second + second @ first - expected))))
    require(clifford_error < TOLERANCE, "Explicit Majoranas fail the Clifford algebra")

    cases = []
    for couplings in COUPLING_CASES:
        for sector, gauge in (("OLD", model.ORIGINAL_GAUGE),
                              ("T-REVERSED", model.T_REVERSED_GAUGE)):
            matrix = source_many_body_hamiltonian(majoranas, couplings, sector)
            hermiticity_error = float(np.max(abs(matrix - matrix.conj().T)))
            require(hermiticity_error < TOLERANCE, "Non-Hermitian explicit matter operator")
            exact_levels = np.linalg.eigvalsh(matrix)
            exact_epsilon = float(exact_levels[0] / N_SITE)

            # Production A is used only as the comparison target. No bond,
            # phase, Fourier coefficient or energy helper comes from it above.
            finite_A = model.finite_majorana_A(N_CELL, couplings, gauge)
            finite_frequencies = np.linalg.eigvalsh(1j * finite_A)
            finite_epsilon = float(-np.sum(abs(finite_frequencies)) / (4 * N_SITE))
            energy_error = abs(exact_epsilon - finite_epsilon)
            require(energy_error < TOLERANCE, "Explicit operator disagrees with vacuum normalization")

            # Each canonical complex fermion contributes +/- e, where e is
            # a positive eigenvalue of h=iA/2. Compare all 64 many-body levels.
            positive_h = finite_frequencies[N_FERMION:] / 2
            require(float(np.min(positive_h)) > 0., "Calibration case unexpectedly has zero modes")
            quadratic_levels = np.sort([
                np.dot(signs, positive_h)
                for signs in product((-1., 1.), repeat=N_FERMION)
            ])
            spectrum_error = float(np.max(abs(exact_levels - quadratic_levels)))
            require(spectrum_error < TOLERANCE,
                    "Explicit and canonical many-body spectra disagree")
            cases.append({
                "sector": sector,
                "couplings_Kx_Ky_Kz_Kt_Kw": couplings,
                "explicit_E0_per_site": exact_epsilon,
                "finite_iA_E0_per_site": finite_epsilon,
                "energy_error": energy_error,
                "hermiticity_error": hermiticity_error,
                "full_64_level_spectrum_error": spectrum_error,
            })

    final_provenance = model.metadata()
    for category in ("script_sha256", "source_sha256", "local_source_sha256"):
        require(initial_provenance[category] == final_provenance[category],
                f"Source changed during explicit Majorana audit: {category}")
    return {
        "status": "PASS",
        "metadata": final_provenance,
        "method": "Pauli tensor products; PDF-literal real bonds; explicit eigvalsh of 64x64 H",
        "source": "Original PDF pp.3-4 Eq.(11), p.5 Eq.(15)",
        "independent_of_production_geometry_and_fourier": True,
        "Ncell": N_CELL,
        "Nsite": N_SITE,
        "N_complex_fermions": N_FERMION,
        "Hilbert_dimension": DIMENSION,
        "Majorana_convention": "{theta_i,theta_j}=2 delta_ij",
        "clifford_max_error": clifford_error,
        "tolerance": TOLERANCE,
        "max_energy_error": max(row["energy_error"] for row in cases),
        "max_full_many_body_spectrum_error": max(row["full_64_level_spectrum_error"] for row in cases),
        "cases": cases,
        "scope": "Unprojected matter Hilbert space; physical spin projection and parity selection not implemented",
        "finite_size_note": "Equal thermodynamic isotropic OLD/T densities do not imply equality at each finite PBC size",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        help="Optional JSON destination; includes model.metadata provenance")
    args = parser.parse_args()
    result = run_explicit_majorana_audit()
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print("T-REVERSED explicit 12-Majorana / 64-state audit: PASS")
    print(f"Clifford algebra maximum error: {result['clifford_max_error']:.3e}")
    print(f"Vacuum energy per site maximum error: {result['max_energy_error']:.3e}")
    print(f"Complete many-body spectrum maximum error: {result['max_full_many_body_spectrum_error']:.3e}")
    for row in result["cases"]:
        print(f"{row['sector']} K={row['couplings_Kx_Ky_Kz_Kt_Kw']} "
              f"E0/N={row['explicit_E0_per_site']:.15f}")


if __name__ == "__main__":
    main()
