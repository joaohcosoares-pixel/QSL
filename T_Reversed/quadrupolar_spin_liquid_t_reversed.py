#!/usr/bin/env python3
"""OLD -> T only, from the original PDF real bonds (pp. 3-7, Eqs. 11,15,17).

Physics path: RealSpaceBond -> RAW -> U.dagger RAW U -> eigvalsh.
Closed coefficients and reduced polynomial are separate validation targets.
Eq.(26) is preserved in the source and documented as inconsistent with Eq.(11).
All project integrals and plots use PROJECT_BZ=(-pi,pi), including both
equivalent representatives of the cell reciprocal interval, with correct weight.
Run: OPENBLAS_NUM_THREADS=1 python3 -W error quadrupolar_spin_liquid_t_reversed.py
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import platform
from datetime import datetime, timezone
from itertools import product
from math import prod
from pathlib import Path
from typing import Mapping, NamedTuple

import numpy as np

ROOT = Path(__file__).resolve().parent
PROJECT_BZ = (-np.pi, np.pi)
N_SITES_PER_CELL = 4
SEED = 20260911
GAP_TOLERANCE = 1e-8
ENERGY_TOLERANCE = 1e-10
COUPLING_ORDER = ("x", "y", "z", "t", "w")
ORIGINAL_GAUGE = {
    "u12x": +1, "u34x": +1, "u23y": +1, "u41y": -1,
    "u42z": +1, "u31z": +1, "u43w": +1, "u12w": -1,
    "u13t": -1, "u24t": -1,
}
T_REVERSED_GAUGE = dict(ORIGINAL_GAUGE)
for _key in ("u13t", "u24t"):
    T_REVERSED_GAUGE[_key] = -ORIGINAL_GAUGE[_key]
GAUGES = {"OLD": ORIGINAL_GAUGE, "T-REVERSED": T_REVERSED_GAUGE}
INTERNAL_A1 = (0, 1, 1, 0)
INTERNAL_A2 = (0, 0, 1, 1)
G_PI = np.diag([1., -1., -1., 1.])
J_PAIR = np.array([[0., 0., -1., 0.], [0., 0., 0., -1.],
                   [1., 0., 0., 0.], [0., 1., 0., 0.]])

# Stored indices and end-minus-start displacements, transcribed from Eq.(11).
# Translation is 2*a1. These are not inferred from a Fourier coefficient.
GEOMETRY = (
    (1, 2, "x", "u12x", +1, 0),
    (3, 4, "x", "u34x", -1, 0),
    (2, 3, "y", "u23y", 0, +1),
    (4, 1, "y", "u41y", 0, -1),
    (4, 2, "z", "u42z", +1, -1),
    (3, 1, "z", "u31z", +1, -1),
    (4, 3, "w", "u43w", -1, 0),
    (1, 2, "w", "u12w", -1, 0),
    (1, 3, "t", "u13t", +1, +1),
    (2, 4, "t", "u24t", +1, +1),
)
PLAQUETTES = {
    "W124": (((0, 1), (0, 2), (0, 4), (0, 1)), "xzy"),
    "W234": (((0, 2), (0, 3), (0, 4), (0, 2)), "yxz"),
    "W123": (((0, 1), (0, 2), (0, 3), (0, 1)), "xyt"),
    "W134": (((0, 1), (0, 3), (0, 4), (0, 1)), "txy"),
    "W213": (((-1, 2), (0, 1), (-1, 3), (-1, 2)), "wzy"),
    "W143": (((0, 1), (0, 4), (-1, 3), (0, 1)), "ywz"),
    "W214": (((-1, 2), (0, 1), (0, 4), (-1, 2)), "wyt"),
    "W243": (((-1, 2), (0, 4), (-1, 3), (-1, 2)), "twy"),
}


def require(condition, message):
    """Critical checks remain active with python -O."""
    if not condition:
        raise RuntimeError(message)


def validate_gauge(gauge: Mapping[str, int]):
    if set(gauge) != set(ORIGINAL_GAUGE):
        raise ValueError("Gauge requires exactly the ten original keys.")
    for value in gauge.values():
        if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)) or value not in (-1, 1):
            raise ValueError("Every stored u must be integer +1 or -1.")


def assert_only_t_reversed(gauge=None):
    gauge = T_REVERSED_GAUGE if gauge is None else gauge
    validate_gauge(ORIGINAL_GAUGE)
    validate_gauge(gauge)
    changed_keys = {key for key in ORIGINAL_GAUGE if ORIGINAL_GAUGE[key] != gauge[key]}
    require(changed_keys == {"u13t", "u24t"}, f"Gauge contamination: {changed_keys}")
    for key in ("u34x", "u43w", "u42z", "u31z", "u13t", "u24t"):
        require(gauge[key] == +1, f"Forbidden change of {key}")
    return sorted(changed_keys)


class RealSpaceBond(NamedTuple):
    start: int
    end: int
    gamma: str
    u: int
    da1: int
    da2: int


def real_space_bonds(gauge=T_REVERSED_GAUGE):
    validate_gauge(gauge)
    return tuple(RealSpaceBond(i, j, kind, gauge[key], d1, d2)
                 for i, j, kind, key, d1, d2 in GEOMETRY)


def cell_shift(bond):
    delta = bond.da1 - (INTERNAL_A1[bond.end-1] - INTERNAL_A1[bond.start-1])
    require(delta % 2 == 0, "Incompatible a1 displacement")
    require(bond.da2 == INTERNAL_A2[bond.end-1] - INTERNAL_A2[bond.start-1],
            "Incompatible a2 displacement")
    return delta // 2


def directed_edge(start, end, gamma, gauge):
    """Resolve actual translated edge; nodes are (cell,sublattice)."""
    for entry, bond in zip(GEOMETRY, real_space_bonds(gauge)):
        if gamma != bond.gamma:
            continue
        shift = cell_shift(bond)
        if start[1] == bond.start and end == (start[0]+shift, bond.end):
            return bond.u, entry[3]
        if end[1] == bond.start and start == (end[0]+shift, bond.end):
            return -bond.u, "-"+entry[3]
    raise ValueError(f"No real {gamma} edge from {start} to {end}")


def directed_u(start, end, gamma, gauge):
    return directed_edge(start, end, gamma, gauge)[0]


def plaquette_details(gauge):
    result = {}
    for name, (nodes, kinds) in PLAQUETTES.items():
        require(nodes[0] == nodes[-1], f"Open plaquette {name}")
        xy = [(2*n+INTERNAL_A1[m-1]+0.5*INTERNAL_A2[m-1], INTERNAL_A2[m-1])
              for n, m in nodes]
        area2 = sum(a[0]*b[1]-b[0]*a[1] for a, b in zip(xy[:-1], xy[1:]))
        require(area2 > 0, f"Non-anticlockwise traversal: {name}")
        factors = [directed_edge(a, b, g, gauge)
                   for a, b, g in zip(nodes[:-1], nodes[1:], kinds)]
        result[name] = {"path": nodes, "kinds": kinds, "factors": [f[0] for f in factors],
                        "expressions": [f[1] for f in factors], "signed_area": area2/2,
                        "flux": prod(f[0] for f in factors)}
    return result


def plaquette_factors(gauge):
    return {p: d["factors"] for p, d in plaquette_details(gauge).items()}


def fluxes(gauge):
    return {p: d["flux"] for p, d in plaquette_details(gauge).items()}


def gauge_equivalences(gauge_a, gauge_b):
    validate_gauge(gauge_a)
    validate_gauge(gauge_b)
    return [signs for signs in product((-1, +1), repeat=4)
            if all(signs[i-1]*gauge_a[key]*signs[j-1] == gauge_b[key]
                   for i, j, _, key, _, _ in GEOMETRY)]


def coupling_map(couplings):
    values = np.asarray(couplings, dtype=float)
    if values.shape != (5,) or not np.all(np.isfinite(values)):
        raise ValueError("Expected finite Kx, Ky, Kz, Kt, Kw.")
    return dict(zip(COUPLING_ORDER, values))


def raw_from_bonds(q1, q2, couplings=(1.,)*5, gauge=T_REVERSED_GAUGE):
    K = coupling_map(couplings)
    if not np.all(np.isfinite([q1, q2])):
        raise ValueError("Momenta must be finite.")
    H = np.zeros((4, 4), dtype=complex)
    for b in real_space_bonds(gauge):
        value = 1j*K[b.gamma]*b.u*np.exp(1j*(q1*b.da1+q2*b.da2))
        H[b.start-1, b.end-1] += value
        H[b.end-1, b.start-1] += value.conjugate()
    return H


def coefficients(q, couplings=(1.,)*5, gauge=T_REVERSED_GAUGE, q2=0.):
    """Independent closed targets, collected by hand from each oriented bond."""
    validate_gauge(gauge)
    Kx, Ky, Kz, Kt, Kw = coupling_map(couplings).values()
    ep, em, cross = np.exp(1j*q), np.exp(-1j*q), np.exp(1j*q2)
    return {
        "f12": 1j*(gauge["u12x"]*Kx*ep+gauge["u12w"]*Kw*em),
        "f13": cross*1j*(gauge["u13t"]*Kt*ep-gauge["u31z"]*Kz*em),
        "f14": cross*(-1j*gauge["u41y"]*Ky),
        "f23": cross*(1j*gauge["u23y"]*Ky),
        "f24": cross*1j*(gauge["u24t"]*Kt*ep-gauge["u42z"]*Kz*em),
        "f34": 1j*(gauge["u34x"]*Kx*em-gauge["u43w"]*Kw*ep),
    }


def matrix_from_coefficients(c):
    H = np.zeros((4, 4), dtype=complex)
    for key, value in c.items():
        i, j = int(key[1])-1, int(key[2])-1
        H[i, j], H[j, i] = value, np.conjugate(value)
    return H


def basis_matrix(q2):
    return np.diag([1., 1., np.exp(-1j*q2), np.exp(-1j*q2)])


def h_bloch(q, couplings=(1.,)*5, gauge=T_REVERSED_GAUGE, q2=0.):
    U = basis_matrix(q2)
    return U.conj().T @ raw_from_bonds(q, q2, couplings, gauge) @ U


def bloch_grid(q_values, couplings=(1.,)*5, gauge=T_REVERSED_GAUGE):
    """Vectorized bond summation, not a closed-form physical shortcut."""
    q = np.asarray(q_values, dtype=float)
    if q.ndim != 1 or q.size == 0 or not np.all(np.isfinite(q)):
        raise ValueError("Expected a nonempty finite one-dimensional q grid.")
    K = coupling_map(couplings)
    H = np.zeros((len(q), 4, 4), dtype=complex)
    # Eq.(27) was independently checked: q2=0 is permitted for final Bloch.
    for b in real_space_bonds(gauge):
        v = 1j*K[b.gamma]*b.u*np.exp(1j*q*b.da1)
        H[:, b.start-1, b.end-1] += v
        H[:, b.end-1, b.start-1] += v.conjugate()
    return H


def bands_sorted(q_values, couplings=(1.,)*5, gauge=T_REVERSED_GAUGE):
    return np.linalg.eigvalsh(bloch_grid(q_values, couplings, gauge))


def occupied_sum(eigenvalues):
    E = np.asarray(eigenvalues, dtype=float)
    if E.ndim < 1 or E.shape[-1] != 4 or not np.all(np.isfinite(E)):
        raise ValueError("Expected finite four-band eigenvalues.")
    return np.where(E < 0., E, 0.).sum(axis=-1)


def energy_per_site(q_values, eigenvalues):
    """epsilon=(1/(8*pi))*integral[-pi,pi] sum(E<0) E dq.

    The width-2pi interval repeats cell momenta twice; the normalization is
    1/(N_sublattices*BZ_width). See finite iA spectrum for an independent audit.
    """
    q, E = np.asarray(q_values, dtype=float), np.asarray(eigenvalues, dtype=float)
    if q.ndim != 1 or len(q) < 2 or not np.all(np.isfinite(q)) or not np.all(np.diff(q) > 0):
        raise ValueError("Integration grid must be finite and strictly increasing.")
    if E.shape != (len(q), 4) or not np.allclose(q[[0, -1]], PROJECT_BZ, rtol=0, atol=1e-12):
        raise ValueError("Principal integration must span PROJECT_BZ=[-pi,pi].")
    integrate = np.trapezoid if hasattr(np, "trapezoid") else np.trapz
    return float(integrate(occupied_sum(E), q)/(4*(PROJECT_BZ[1]-PROJECT_BZ[0])))


def reduced_components(q, couplings=(1.,)*5, gauge=T_REVERSED_GAUGE):
    validate_gauge(gauge)
    require(gauge["u34x"] == gauge["u12x"] and gauge["u43w"] == -gauge["u12w"]
            and gauge["u24t"] == gauge["u13t"] and gauge["u42z"] == gauge["u31z"]
            and gauge["u23y"] == -gauge["u41y"], "Gauge lacks required reduced structure")
    Kx, Ky, Kz, Kt, Kw = coupling_map(couplings).values()
    ux, uw, ut, uz = (gauge[k] for k in ("u12x", "u12w", "u13t", "u31z"))
    return ((uw*Kw-ux*Kx)*np.sin(q), (ux*Kx+uw*Kw)*np.cos(q),
            -(ut*Kt+uz*Kz)*np.sin(q), (ut*Kt-uz*Kz)*np.cos(q), gauge["u23y"]*Ky)


def reduced_invariants(q, couplings=(1.,)*5, gauge=T_REVERSED_GAUGE):
    ar, ai, br, bi, y = reduced_components(q, couplings, gauge)
    S = ar*ar+ai*ai+br*br+bi*bi+y*y
    R = ai*ai*(br*br+bi*bi)+bi*bi*y*y
    return S, R


def polynomial_coefficients(q, couplings=(1.,)*5, gauge=T_REVERSED_GAUGE):
    S, R = reduced_invariants(q, couplings, gauge)
    return np.array([1., 0., -2*S, 0., S*S-4*R])


def analytic_spectrum(q, couplings=(1.,)*5, gauge=T_REVERSED_GAUGE):
    """Validation target from the independently proven even characteristic polynomial."""
    S, R = reduced_invariants(q, couplings, gauge)
    high2, low2 = S+2*np.sqrt(R), S-2*np.sqrt(R)
    require(low2 >= -1e-12*max(1., S), "Negative squared eigenvalue")
    high, low = np.sqrt(max(0., high2)), np.sqrt(max(0., low2))
    return np.array([-high, -low, low, high])


def finite_majorana_A(ncell, couplings=(1.,)*5, gauge=T_REVERSED_GAUGE):
    if isinstance(ncell, bool) or not isinstance(ncell, (int, np.integer)) or ncell < 3:
        raise ValueError("Ncell must be an integer >=3.")
    K = coupling_map(couplings)
    A = np.zeros((4*ncell, 4*ncell))
    for cell in range(ncell):
        for b in real_space_bonds(gauge):
            i = 4*cell+b.start-1
            j = 4*((cell+cell_shift(b)) % ncell)+b.end-1
            A[i, j] += 2*K[b.gamma]*b.u
            A[j, i] -= 2*K[b.gamma]*b.u
    require(np.array_equal(A.T, -A), "Finite A is not antisymmetric")
    return A


def finite_momenta(ncell):
    # One representative for each physical cell momentum, always inside project BZ.
    return np.pi*np.arange(-(ncell//2), ncell-(ncell//2))/ncell


def finite_fourier_matrix(ncell, q, q2=0.):
    require(abs(np.exp(2j*ncell*q)-1) < 1e-10, "Momentum incompatible with PBC")
    V = np.zeros((4*ncell, 4), dtype=complex)
    for cell in range(ncell):
        for m in range(4):
            V[4*cell+m, m] = np.exp(1j*(q*(2*cell+INTERNAL_A1[m])+q2*INTERNAL_A2[m]))/np.sqrt(ncell)
    return V


def finite_fourier_block(A, q, q2=0.):
    ncell = len(A)//4
    require(A.shape == (4*ncell, 4*ncell), "Invalid A shape")
    V = finite_fourier_matrix(ncell, q, q2)
    return V.conj().T @ (0.5j*A) @ V


def delta_from_t_bonds(q, couplings=(1.,)*5):
    """Independent delta using only the two real t bonds, never f13/f24."""
    K = coupling_map(couplings)
    H = np.zeros((4, 4), dtype=complex)
    for b in real_space_bonds(ORIGINAL_GAUGE):
        if b.gamma == "t":
            value = -2j*K[b.gamma]*b.u*np.exp(1j*q*b.da1)
            H[b.start-1, b.end-1] += value
            H[b.end-1, b.start-1] += value.conjugate()
    return H


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def metadata(couplings=(1.,)*5, random_cases=1000, seed=SEED):
    timestamp = datetime.now(timezone.utc).isoformat()
    majorana_convention = "{gamma_i,gamma_j}=2 delta_ij; H=(i/4)gamma^T A gamma; h=iA/2"
    sources = [
        Path('/home/joaohenrique/Pesquisa/Vanuildo (Mat. Condensada)/Quadrupolar Spin Liquid in 1D (cópia).pdf'),
        Path('/home/joaohenrique/Pesquisa/Pesquisa/QSL/Quadrupolar_Spin_Liquid_OLD_reference_obsidian.md'),
        Path('/home/joaohenrique/.codex/attachments/bab48a6d-bb75-4528-9c03-7d1ae00f3f7a/pasted-text.txt'),
    ]
    sources += [ROOT.parent/'comparison_z_reversed.json', ROOT.parent/'validation_z_reversed.json']
    sources += sorted(Path('/home/joaohenrique/Computing/QSL/Z_Reversed_entregaveis').glob('*.py'))
    return {
        "generated_at_utc": timestamp, "timestamp": timestamp,
        "current_sector": "T-REVERSED", "reference_sector": "OLD",
        **{f"K{k}": float(v) for k, v in coupling_map(couplings).items()},
        "q_min": float(PROJECT_BZ[0]), "q_max": float(PROJECT_BZ[1]),
        "bz_width": float(2*np.pi), "random_cases": random_cases,
        "BZ_width": float(2*np.pi), "N_sites_per_cell": N_SITES_PER_CELL,
        "normalization_formula": "epsilon = integral[-pi,pi] F(q) dq / (8*pi)",
        "integration_method": "trapezoid on nested grids; independent adaptive scipy.quad",
        "gap_tolerance": GAP_TOLERANCE, "energy_tolerance": ENERGY_TOLERANCE,
        "software_versions": {"python": platform.python_version(), **{
            name: importlib.metadata.version(name) for name in ('numpy', 'scipy', 'sympy', 'matplotlib')}},
        "seed": seed, "number_random_cases": random_cases,
        "script_sha256": {p.name: sha256(p) for p in sorted(ROOT.glob('*.py'))},
        "source_sha256": {str(p): sha256(p) for p in sources if p.exists()},
        "local_source_sha256": {str(p.relative_to(ROOT)): sha256(p) for p in sorted((ROOT/'sources').glob('*')) if p.is_file()},
        "cell_translation": "2*a1", "q_definition": "k dot a1",
        "finite_momentum_rule": "q=pi*ell/Ncell; Ncell independent representatives, each inside PROJECT_BZ",
        "majorana_convention": majorana_convention, "Majorana_convention": majorana_convention,
        "finite_energy_formula": "E0=-sum(abs(eigvalsh(i*A)))/4; epsilon=E0/(4*Ncell)",
        "warnings": "Run with -W error; none suppressed",
    }


def randomized_suite(random_cases=1000, seed=SEED):
    """Only the full per-case suite, independently repeatable for benchmarking."""
    if isinstance(random_cases, bool) or not isinstance(random_cases, int) or random_cases < 1000:
        raise ValueError("Audit requires an integer >=1000 random cases.")
    maxima = {}
    def check(name, left, right, tol=3e-11):
        error = float(np.max(np.abs(np.asarray(left)-np.asarray(right))))
        require(np.isfinite(error), f"Nonfinite error: {name}")
        maxima[name] = max(maxima.get(name, 0.), error)
        require(error <= tol, f"{name}: error {error} > {tol}")

    rng = np.random.default_rng(seed)
    signed_count = 0
    min_spacing = float('inf')
    for case in range(random_cases):
        while True:
            K = rng.uniform(.15, 2., 5)
            spacing = float(np.min(np.diff(np.sort(K))))
            if spacing > .025:
                break
        min_spacing = min(min_spacing, spacing)
        if case % 2:
            # Every such case contains both signs; no seed selection by performance.
            signs = rng.choice([-1, 1], 5)
            signs[0], signs[1] = -1, +1
            K *= signs
            signed_count += 1
        q, q2 = rng.uniform(*PROJECT_BZ, size=2)
        U = basis_matrix(q2)
        check('unitarity', U.conj().T@U, np.eye(4))
        assert_only_t_reversed()
        for gauge in GAUGES.values():
            raw = raw_from_bonds(q, q2, K, gauge)
            H = h_bloch(q, K, gauge, q2)
            E = np.linalg.eigvalsh(H)
            check('raw_geometry_vs_closed', raw, matrix_from_coefficients(coefficients(q, K, gauge, q2)))
            # A direct phase assertion independent of Eq.(27) and the coefficient formulas.
            raw0 = raw_from_bonds(q, 0., K, gauge)
            for i, j in ((0, 2), (0, 3), (1, 2), (1, 3)):
                check('raw_cross_phase_positive_q2', raw[i, j], np.exp(1j*q2)*raw0[i, j])
            check('eq27_vs_independent_closed_bloch', U.conj().T@raw@U,
                  matrix_from_coefficients(coefficients(q, K, gauge)))
            check('q2_cancellation', H, h_bloch(q, K, gauge, 0.))
            check('vectorized_bonds', H, bloch_grid([q], K, gauge)[0])
            check('hermiticity', H, H.conj().T)
            check('majorana', h_bloch(-q, K, gauge), -H.conjugate())
            check('period_2pi', h_bloch(q+2*np.pi, K, gauge), H)
            check('covariance_pi', h_bloch(q+np.pi, K, gauge), G_PI@H@G_PI)
            check('spectrum_period_pi', np.linalg.eigvalsh(h_bloch(q+np.pi, K, gauge)), E)
            check('raw_bloch_spectrum', np.linalg.eigvalsh(raw), E)
            check('same_q_antiunitary', J_PAIR@H.conjugate()@J_PAIR.T, -H)
            check('same_q_pairing', E, -E[::-1])
            check('closed_spectrum', E, analytic_spectrum(q, K, gauge), tol=2e-9)
            poly = polynomial_coefficients(q, K, gauge)
            scale = max(1., float(np.max(np.abs(poly))))
            check('characteristic_polynomial_relative', np.poly(H)/scale, poly/scale)
            root_scale = np.maximum(1., np.sum(abs(poly)[None, :]*abs(E[:, None])**np.arange(4, -1, -1), axis=1))
            check('polynomial_at_eigenvalues_relative', np.polyval(poly, E)/root_scale, np.zeros(4))
            # Direct determinant at arbitrary lambda avoids using eigenvalues on both sides.
            lam = rng.uniform(-4, 4)
            pval = np.polyval(poly, lam)
            check('determinant_polynomial_relative', np.linalg.det(lam*np.eye(4)-H)/max(1., abs(pval)),
                  pval/max(1., abs(pval)), tol=1e-10)
        check('only_t_matrix_difference', h_bloch(q, K, T_REVERSED_GAUGE)-h_bloch(q, K, ORIGINAL_GAUGE),
              delta_from_t_bonds(q, K))

    return {"status": "PASS", "number_random_cases": random_cases, "sectors_per_case": 2,
            "signed_cases": signed_count, "minimum_absolute_coupling_spacing": min_spacing,
            "seed": seed, "max_errors": maxima}


def self_test(random_cases=1000, seed=SEED):
    if random_cases < 1000:
        raise ValueError("Audit requires at least 1000 random cases.")
    initial_provenance = metadata(random_cases=random_cases, seed=seed)
    changed = assert_only_t_reversed()
    # Independent symbolic file derives its matrices before comparing this module.
    from symbolic_t_audit import run_symbolic_audit
    symbolic = run_symbolic_audit()
    f_old, f_t = fluxes(ORIGINAL_GAUGE), fluxes(T_REVERSED_GAUGE)
    require(all(v == 1 for v in f_old.values()), "PDF zero-flux reference failed")
    for name, (_, kinds) in PLAQUETTES.items():
        require(f_t[name] == f_old[name]*(-1)**kinds.count('t'), "Flux edge parity failed")
    require(gauge_equivalences(ORIGINAL_GAUGE, T_REVERSED_GAUGE) == [], "Unexpected gauge equivalence")
    require(len(gauge_equivalences(ORIGINAL_GAUGE, ORIGINAL_GAUGE)) == 2, "Gauge search regression")
    for gauge in GAUGES.values():
        for signs in product((-1, 1), repeat=4):
            transformed = {key: signs[i-1]*gauge[key]*signs[j-1]
                           for i, j, _, key, _, _ in GEOMETRY}
            require(fluxes(transformed) == fluxes(gauge), "Flux is not gauge invariant")

    maxima = {}
    def check(name, left, right, tol=3e-11):
        error = float(np.max(np.abs(np.asarray(left)-np.asarray(right))))
        require(np.isfinite(error), f"Nonfinite error: {name}")
        maxima[name] = max(maxima.get(name, 0.), error)
        require(error <= tol, f"{name}: error {error} > {tol}")

    random_result = randomized_suite(random_cases, seed)
    maxima.update(random_result["max_errors"])
    signed_count = random_result["signed_cases"]
    min_spacing = random_result["minimum_absolute_coupling_spacing"]

    finite_details = []
    for ncell in (5, 8, 17):
        for K in ((1.13, .79, 1.31, .92, .67), (-1.21, .71, 1.41, -.93, .53)):
            for name, gauge in GAUGES.items():
                A = finite_majorana_A(ncell, K, gauge)
                momenta = finite_momenta(ncell)
                blocks, basis = [], []
                for q in momenta:
                    V = finite_fourier_matrix(ncell, q)
                    basis.append(V)
                    block = V.conj().T @ (0.5j*A) @ V
                    check('finite_bloch_block', block, h_bloch(q, K, gauge))
                    check('finite_raw_block', finite_fourier_block(A, q, .371), raw_from_bonds(q, .371, K, gauge))
                    blocks.extend(np.linalg.eigvalsh(block))
                Vfull = np.concatenate(basis, axis=1)
                check('finite_fourier_completeness', Vfull.conj().T@Vfull, np.eye(4*ncell))
                Er = np.linalg.eigvalsh(1j*A)
                check('finite_full_iA_spectrum', np.sort(2*np.array(blocks)), Er)
                e_real = -np.sum(np.abs(Er))/(16*ncell)
                Eb = np.array(blocks)
                e_bloch = Eb[Eb < 0].sum()/(4*ncell)
                check('finite_energy_normalization', e_real, e_bloch)
                finite_details.append({'sector': name, 'Ncell': ncell, 'Nsite': 4*ncell,
                                       'couplings': K, 'epsilon_real': e_real, 'epsilon_bloch': e_bloch})

    # Exactly solvable real-space dimers anchor the absolute energy convention.
    # Each isolated family is a perfect matching, hence two dimers per cell.
    # E0/Nsite = -(2*Ncell*abs(K))/(4*Ncell), independent of every Fourier formula.
    dimer_checks = []
    for family_index, family in enumerate(COUPLING_ORDER):
        K = np.zeros(5)
        K[family_index] = -1.7
        expected_energy = -abs(K[family_index])/2
        q_dimer = np.linspace(*PROJECT_BZ, 33)
        for name, gauge in GAUGES.items():
            A_dimer = finite_majorana_A(5, K, gauge)
            finite_energy = -np.sum(abs(np.linalg.eigvalsh(1j*A_dimer)))/(16*5)
            bulk_energy = energy_per_site(q_dimer, bands_sorted(q_dimer, K, gauge))
            check('isolated_dimers_finite_absolute_normalization', finite_energy, expected_energy)
            check('isolated_dimers_bulk_absolute_normalization', bulk_energy, expected_energy)
            dimer_checks.append({'family': family, 'sector': name, 'exact_epsilon': expected_energy,
                                 'finite_epsilon': finite_energy, 'bulk_epsilon': bulk_energy})

    # Adversarial mutation tests explicitly show why final Bloch alone is insufficient.
    K, q, q2 = (1.17, .83, 1.49, .62, .39), .417, .693
    U = basis_matrix(q2)
    raw = raw_from_bonds(q, q2, K)
    wrong_raw = raw_from_bonds(q, -q2, K)
    H = h_bloch(q, K)
    phase_error = float(np.max(np.abs(wrong_raw-raw)))
    wrong_eq27_error = float(np.max(np.abs(U@raw@U.conj().T-H)))
    compensated_error = float(np.max(np.abs(U@wrong_raw@U.conj().T-H)))
    require(phase_error > .1 and wrong_eq27_error > .1, 'Mutation did not expose q2 errors')
    require(compensated_error < 1e-12, 'Expected compensating double-bug demonstration failed')
    accepted_sign_assignments = []
    for signs in product((-1, +1), repeat=10):
        candidate = dict(zip(ORIGINAL_GAUGE, signs))
        try:
            assert_only_t_reversed(candidate)
        except RuntimeError:
            require(candidate != T_REVERSED_GAUGE, 'Correct T gauge was rejected')
        else:
            require(candidate == T_REVERSED_GAUGE, 'Gauge mutation escaped detection')
            accepted_sign_assignments.append(candidate)
    require(len(accepted_sign_assignments) == 1, 'Only the exact T sector may pass')

    # Numerical regression values belong to prior audited data, not the literal Eq.(26).
    old_data = json.loads((ROOT/'sources/legacy_comparison_z_reversed.json').read_text())
    q_grid = np.linspace(*PROJECT_BZ, 8001)
    old_energy = energy_per_site(q_grid, bands_sorted(q_grid, gauge=ORIGINAL_GAUGE))
    old_expected = old_data['summary']['epsilon_old']
    check('OLD_regression_energy', old_energy, old_expected, tol=1e-11)
    for label, old_spectrum in old_data['special_points']['OLD'].items():
        point = {'-pi/2': -np.pi/2, '-pi/4': -np.pi/4, '0': 0.,
                 'pi/4': np.pi/4, 'pi/2': np.pi/2}[label]
        check('OLD_regression_special_spectra', np.linalg.eigvalsh(h_bloch(point, gauge=ORIGINAL_GAUGE)), old_spectrum)
    final_provenance = metadata(random_cases=random_cases, seed=seed)
    for kind in ('script_sha256', 'source_sha256', 'local_source_sha256'):
        require(initial_provenance[kind] == final_provenance[kind],
                f'Sources changed during validation: {kind}; rerun the audit')
    return {'status': 'PASS', 'metadata': final_provenance,
            'symbolic': symbolic, 'changed_keys': changed, 'number_random_cases': random_cases,
            'sectors_per_case': 2, 'signed_cases': signed_count, 'minimum_absolute_coupling_spacing': min_spacing,
            'seed': seed, 'max_errors': maxima, 'finite_chains': finite_details,
            'isolated_dimers_absolute_normalization': dimer_checks,
            'fluxes': {name: fluxes(gauge) for name, gauge in GAUGES.items()},
            'plaquette_details': {name: plaquette_details(gauge) for name, gauge in GAUGES.items()},
            'diagonal_gauges_tested': 16, 'gauge_equivalence_solutions': [],
            'mutations': {'wrong_raw_q2_detected_error': phase_error,
                          'wrong_eq27_detected_error': wrong_eq27_error,
                          'two_bugs_compensate_at_final_bloch_error': compensated_error,
                          'raw_test_rejects_compensating_pair': True, 'non_t_mutations_rejected': True,
                          'partial_t_reversal_rejected': True, 'stored_sign_assignments_tested': 1024,
                          'accepted_sign_assignments': accepted_sign_assignments},
            'OLD_regression': {'epsilon_old': old_energy, 'prior_audited_epsilon': old_expected,
                               'new_OLD_recomputed_on_project_BZ_with_equivalent_normalization': True}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--random-cases', type=int, default=1000)
    parser.add_argument('--seed', type=int, default=SEED)
    parser.add_argument('--output-dir', type=Path, default=ROOT)
    args = parser.parse_args()
    result = self_test(args.random_cases, args.seed)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir/'validation_t_reversed.json').write_text(json.dumps(result, indent=2)+'\n')
    print(f"T-REVERSED core audit: PASS; {args.random_cases} random cases; finite chains 5,8,17")
    print(json.dumps(result['max_errors'], indent=2))


if __name__ == '__main__':
    main()
