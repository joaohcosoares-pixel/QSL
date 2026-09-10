#!/usr/bin/env python3
"""Geometry-derived OLD and Z-REVERSED Majorana Hamiltonians.

Required OLD reference: Quadrupolar_Spin_Liquid_OLD_reference.md.
Its literal Eq.(26) is documented as an erratum: the geometry / Eq.(11) /
Eq.(17) give f34=i*(Kx*exp(-iq)-Kw*exp(iq)). No source file is modified.

Pipeline: original directed bonds -> RAW Fourier -> Eq.(27) -> Bloch.
Eq.(27) is the active substitution H_raw[U Theta], giving U^dagger H_raw U.
Equivalently, for coordinates of the same operator, Theta_raw=U Theta_Bloch.
Finite-chain audit: H=(i/4) gamma^T A gamma, h=iA/2, Aij=2*K*u.

Run this file for symbolic tests, 1000 asymmetric random cases and finite chains.
Only NumPy and SymPy are needed by this core module.
"""
from __future__ import annotations

import argparse
import json
from itertools import product
from math import prod
from pathlib import Path
from typing import Mapping, NamedTuple

import numpy as np
import sympy as sp

ORIGINAL_GAUGE = {
    "u12x": +1, "u34x": +1, "u23y": +1, "u41y": -1,
    "u42z": +1, "u31z": +1, "u43w": +1, "u12w": -1,
    "u13t": -1, "u24t": -1,
}
Z_REVERSED_GAUGE = {
    **ORIGINAL_GAUGE,
    "u42z": -ORIGINAL_GAUGE["u42z"],
    "u31z": -ORIGINAL_GAUGE["u31z"],
}
GAUGES = {"OLD": ORIGINAL_GAUGE, "Z-REVERSED": Z_REVERSED_GAUGE}
COUPLING_ORDER = ("x", "y", "z", "t", "w")
INTERNAL_A1 = (0, 1, 1, 0)
INTERNAL_A2 = (0, 0, 1, 1)
G_PI = np.diag([1., -1., -1., 1.])

# Direct transcription of stored orientations and end-minus-start displacements
# from Eq.(11), independent of the closed Fourier formulas below.
GEOMETRY = (
    (1, 2, "x", "u12x", +1, 0),
    (2, 3, "y", "u23y", 0, +1),
    (3, 4, "x", "u34x", -1, 0),
    (4, 1, "y", "u41y", 0, -1),
    (4, 2, "z", "u42z", +1, -1),
    (4, 3, "w", "u43w", -1, 0),
    (3, 1, "z", "u31z", +1, -1),
    (1, 2, "w", "u12w", -1, 0),
    (1, 3, "t", "u13t", +1, +1),
    (2, 4, "t", "u24t", +1, +1),
)

# Each node is (cell, sublattice). All paths are genuinely closed and CCW.
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


def validate_gauge(gauge: Mapping[str, int]) -> None:
    if set(gauge) != set(ORIGINAL_GAUGE):
        raise ValueError("Gauge must contain exactly the ten original bond keys.")
    for value in gauge.values():
        if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)) or value not in (-1, 1):
            raise ValueError("All u must be integer +1 or -1.")


def assert_only_z_reversed() -> None:
    validate_gauge(ORIGINAL_GAUGE)
    validate_gauge(Z_REVERSED_GAUGE)
    assert ORIGINAL_GAUGE["u34x"] == +1
    assert ORIGINAL_GAUGE["u43w"] == +1
    assert Z_REVERSED_GAUGE["u34x"] == +1
    assert Z_REVERSED_GAUGE["u43w"] == +1
    assert Z_REVERSED_GAUGE["u42z"] == -ORIGINAL_GAUGE["u42z"]
    assert Z_REVERSED_GAUGE["u31z"] == -ORIGINAL_GAUGE["u31z"]
    for key in ORIGINAL_GAUGE:
        if key not in ("u42z", "u31z"):
            if Z_REVERSED_GAUGE[key] != ORIGINAL_GAUGE[key]:
                raise AssertionError(f"Non-z bond changed: {key}")


class RealSpaceBond(NamedTuple):
    start: int
    end: int
    gamma: str
    u: int
    da1: int
    da2: int


def real_space_bonds(gauge=Z_REVERSED_GAUGE):
    validate_gauge(gauge)
    return tuple(RealSpaceBond(i, j, g, gauge[key], d1, d2)
                 for i, j, g, key, d1, d2 in GEOMETRY)


def cell_shift(bond):
    numerator = bond.da1 - (INTERNAL_A1[bond.end-1]-INTERNAL_A1[bond.start-1])
    if numerator % 2 or bond.da2 != INTERNAL_A2[bond.end-1]-INTERNAL_A2[bond.start-1]:
        raise ValueError("Bond displacement incompatible with original unit cell.")
    return numerator // 2


def directed_u(start, end, gamma, gauge):
    """Find an actual periodic edge, reversing u only for a reversed traversal."""
    for bond in real_space_bonds(gauge):
        shift = cell_shift(bond)
        if gamma != bond.gamma:
            continue
        if start[1] == bond.start and end == (start[0]+shift, bond.end):
            return bond.u
        if end[1] == bond.start and start == (end[0]+shift, bond.end):
            return -bond.u
    raise ValueError(f"No {gamma} bond between {start} and {end}")


def plaquette_factors(gauge):
    result = {}
    for name, (nodes, kinds) in PLAQUETTES.items():
        assert nodes[0] == nodes[-1]
        result[name] = [directed_u(a, b, g, gauge)
                        for a, b, g in zip(nodes[:-1], nodes[1:], kinds)]
        xy = [(2*n+INTERNAL_A1[m-1]+.5*INTERNAL_A2[m-1], INTERNAL_A2[m-1])
              for n, m in nodes]
        signed_area_twice = sum(a[0]*b[1]-b[0]*a[1] for a, b in zip(xy[:-1], xy[1:]))
        assert signed_area_twice > 0, f"Non-CCW path {name}"
    return result


def fluxes(gauge):
    return {name: prod(factors) for name, factors in plaquette_factors(gauge).items()}


def gauge_equivalences(gauge_a, gauge_b):
    validate_gauge(gauge_a)
    validate_gauge(gauge_b)
    solutions = []
    for signs in product((-1, +1), repeat=4):
        if all(signs[i-1]*gauge_a[key]*signs[j-1] == gauge_b[key]
               for i, j, g, key, d1, d2 in GEOMETRY):
            solutions.append(signs)
    return solutions


def coupling_map(couplings):
    values = np.asarray(couplings, dtype=float)
    if values.shape != (5,) or not np.all(np.isfinite(values)):
        raise ValueError("Expected five finite real couplings: Kx,Ky,Kz,Kt,Kw.")
    return dict(zip(COUPLING_ORDER, values))


def raw_from_bonds(q1, q2, couplings, gauge=Z_REVERSED_GAUGE):
    K = coupling_map(couplings)
    H = np.zeros((4, 4), dtype=complex)
    for b in real_space_bonds(gauge):
        value = 1j*K[b.gamma]*b.u*np.exp(1j*(q1*b.da1+q2*b.da2))
        H[b.start-1, b.end-1] += value
        H[b.end-1, b.start-1] += value.conjugate()
    return H


def coefficients(q, couplings, gauge=Z_REVERSED_GAUGE, q2=0.):
    """Closed formulas as an audit target, not the geometry source."""
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


def h_bloch(q, couplings=(1.,)*5, gauge=Z_REVERSED_GAUGE, q2=0.):
    """Physics path always starts with the real bond list."""
    U = basis_matrix(q2)
    return U.conj().T @ raw_from_bonds(q, q2, couplings, gauge) @ U


def bands_sorted(q_values, couplings=(1.,)*5, gauge=Z_REVERSED_GAUGE):
    q = np.asarray(q_values, dtype=float)
    if q.ndim != 1 or q.size == 0 or not np.all(np.isfinite(q)):
        raise ValueError("q_values must be a nonempty finite one-dimensional grid.")
    return np.linalg.eigvalsh(np.stack([h_bloch(v, couplings, gauge) for v in q]))


def occupied_sum(eigenvalues):
    E = np.asarray(eigenvalues, dtype=float)
    if not np.all(np.isfinite(E)) or E.shape[-1] != 4:
        raise ValueError("Expected finite four-band eigenvalues.")
    return np.where(E < 0., E, 0.).sum(axis=-1)


def energy_per_site(q_values, eigenvalues, extended=False):
    """Physical BZ width pi, or equivalent extended BZ width 2*pi."""
    q = np.asarray(q_values, dtype=float)
    E = np.asarray(eigenvalues, dtype=float)
    width = 2*np.pi if extended else np.pi
    if q.ndim != 1 or len(q) < 2 or not np.all(np.isfinite(q)) or not np.all(np.diff(q) > 0):
        raise ValueError("Integration grid must be finite and strictly increasing.")
    if E.shape != (len(q), 4) or not np.isclose(q[-1]-q[0], width, rtol=0., atol=1e-12):
        raise ValueError("Wrong band shape or Brillouin-zone integration width.")
    integrate = np.trapezoid if hasattr(np, "trapezoid") else np.trapz
    return float(integrate(occupied_sum(E), q)/(4*width))


def reduced_components(q, couplings, gauge=Z_REVERSED_GAUGE):
    if gauge not in (ORIGINAL_GAUGE, Z_REVERSED_GAUGE):
        raise ValueError("Reduced polynomial is defined only for the two audited gauges.")
    Kx, Ky, Kz, Kt, Kw = coupling_map(couplings).values()
    sigma = gauge["u31z"]
    return (-(Kx+Kw)*np.sin(q), (Kx-Kw)*np.cos(q),
            (Kt-sigma*Kz)*np.sin(q), -(Kt+sigma*Kz)*np.cos(q), Ky)


def polynomial_coefficients(q, couplings, gauge=Z_REVERSED_GAUGE):
    ar, ai, br, bi, y = reduced_components(q, couplings, gauge)
    S = ar*ar+ai*ai+br*br+bi*bi+y*y
    R = ai*ai*(br*br+bi*bi)+bi*bi*y*y
    return np.array([1., 0., -2*S, 0., S*S-4*R])


def finite_majorana_A(ncell, couplings, gauge=Z_REVERSED_GAUGE):
    if isinstance(ncell, bool) or not isinstance(ncell, int) or ncell < 3:
        raise ValueError("ncell must be an integer >=3.")
    K = coupling_map(couplings)
    A = np.zeros((4*ncell, 4*ncell))
    for cell in range(ncell):
        for b in real_space_bonds(gauge):
            i = 4*cell+b.start-1
            j = 4*((cell+cell_shift(b)) % ncell)+b.end-1
            A[i, j] += 2*K[b.gamma]*b.u
            A[j, i] -= 2*K[b.gamma]*b.u
    assert np.array_equal(A.T, -A)
    return A


def finite_fourier_block(A, q, q2=0.):
    ncell = len(A)//4
    if A.shape != (4*ncell, 4*ncell) or ncell < 3:
        raise ValueError("Invalid finite Majorana matrix shape.")
    if abs(np.exp(2j*ncell*q)-1) > 1e-10:
        raise ValueError("Momentum is incompatible with periodic boundary conditions.")
    V = np.zeros((4*ncell, 4), dtype=complex)
    for cell in range(ncell):
        for m in range(4):
            V[4*cell+m, m] = np.exp(1j*(q*(2*cell+INTERNAL_A1[m])+q2*INTERNAL_A2[m]))/np.sqrt(ncell)
    return V.conj().T @ (0.5j*A) @ V


def symbolic_tests():
    q, q2 = sp.symbols("q q2", real=True)
    Kx, Ky, Kz, Kt, Kw = sp.symbols("Kx Ky Kz Kt Kw", real=True)
    K = dict(zip(COUPLING_ORDER, (Kx, Ky, Kz, Kt, Kw)))
    U = sp.diag(1, 1, sp.exp(-sp.I*q2), sp.exp(-sp.I*q2))
    simp = lambda M: M.applyfunc(sp.simplify)
    assert simp(U.conjugate().T*U) == sp.eye(4)
    ar, ai, br, bi, y, lam = sp.symbols("a_R a_I b_R b_I y lambda", real=True)
    a, b = ar+sp.I*ai, br+sp.I*bi
    reduced = sp.Matrix([[0, a, b, sp.I*y], [sp.conjugate(a), 0, sp.I*y, b],
                         [sp.conjugate(b), -sp.I*y, 0, -sp.conjugate(a)],
                         [-sp.I*y, sp.conjugate(b), -a, 0]])
    S = ar**2+ai**2+br**2+bi**2+y**2
    R = ai**2*(br**2+bi**2)+bi**2*y**2
    characteristic = reduced.charpoly()
    polynomial = characteristic.as_expr().subs(characteristic.gen, lam)
    compact = lam**4-2*S*lam**2+S**2-4*R
    assert sp.expand(polynomial-compact) == 0
    J = sp.Matrix([[0,0,-1,0],[0,0,0,-1],[1,0,0,0],[0,1,0,0]])
    assert simp(J*sp.conjugate(reduced)*J.T+reduced) == sp.zeros(4)
    assert J*J == -sp.eye(4)
    for gauge in GAUGES.values():
        raw = sp.zeros(4)
        for bond in real_space_bonds(gauge):
            value = sp.I*K[bond.gamma]*bond.u*sp.exp(sp.I*(bond.da1*q+bond.da2*q2))
            raw[bond.start-1, bond.end-1] += value
            raw[bond.end-1, bond.start-1] += sp.conjugate(value)
        ep, em, phase = sp.exp(sp.I*q), sp.exp(-sp.I*q), sp.exp(sp.I*q2)
        vals = [sp.I*(gauge["u12x"]*Kx*ep+gauge["u12w"]*Kw*em),
                phase*sp.I*(gauge["u13t"]*Kt*ep-gauge["u31z"]*Kz*em),
                -phase*sp.I*gauge["u41y"]*Ky, phase*sp.I*gauge["u23y"]*Ky,
                phase*sp.I*(gauge["u24t"]*Kt*ep-gauge["u42z"]*Kz*em),
                sp.I*(gauge["u34x"]*Kx*em-gauge["u43w"]*Kw*ep)]
        formula = sp.zeros(4)
        for (i, j), value in zip(((0,1),(0,2),(0,3),(1,2),(1,3),(2,3)), vals):
            formula[i,j], formula[j,i] = value, sp.conjugate(value)
        assert simp(raw-formula) == sp.zeros(4)
        H = simp(U.conjugate().T*raw*U)
        assert not H.has(q2)
        assert simp(H-H.conjugate().T) == sp.zeros(4)
        assert simp(H.subs(q, -q)+sp.conjugate(H)) == sp.zeros(4)
        assert simp(H.subs(q, q+2*sp.pi)-H) == sp.zeros(4)
        D = sp.diag(1,-1,-1,1)
        assert simp(H.subs(q, q+sp.pi)-D*H*D) == sp.zeros(4)
        assert simp(sp.Matrix([H[1,2]-H[0,3], H[1,3]-H[0,2], H[2,3]+sp.conjugate(H[0,1])])) == sp.zeros(3,1)
        sigma = gauge["u31z"]
        substitution = {ar:-(Kx+Kw)*sp.sin(q), ai:(Kx-Kw)*sp.cos(q),
                        br:(Kt-sigma*Kz)*sp.sin(q), bi:-(Kt+sigma*Kz)*sp.cos(q), y:Ky}
        assert simp((H-reduced.subs(substitution)).expand(complex=True)) == sp.zeros(4)
        Hiso = simp(H.subs({Kx:1,Ky:1,Kz:1,Kt:1,Kw:1}))
        if gauge == Z_REVERSED_GAUGE:
            assert simp((Hiso*Hiso-(1+8*sp.sin(q)**2)*sp.eye(4)).expand(complex=True)) == sp.zeros(4)
        else:
            isotropic_polynomial = compact.subs(substitution).subs({Kx:1,Ky:1,Kz:1,Kt:1,Kw:1})
            assert sp.trigsimp(isotropic_polynomial-(lam**4-10*lam**2+25-16*sp.cos(q)**2)) == 0
    return {"symbolic_assertions": "passed for OLD and Z-REVERSED",
            "characteristic_polynomial_expanded": str(sp.expand(polynomial)),
            "characteristic_polynomial_compact": str(compact),
            "linear_lambda_coefficient": 0}


def self_test(random_cases=1000, seed=20260908):
    if random_cases < 1000:
        raise ValueError("Scientific audit requires at least 1000 random cases.")
    assert_only_z_reversed()
    symbolic = symbolic_tests()
    assert tuple(fluxes(ORIGINAL_GAUGE).values()) == (1,)*8
    assert tuple(fluxes(Z_REVERSED_GAUGE).values()) == (-1,-1,1,1,-1,-1,1,1)
    assert gauge_equivalences(ORIGINAL_GAUGE, Z_REVERSED_GAUGE) == []
    assert len(gauge_equivalences(ORIGINAL_GAUGE, ORIGINAL_GAUGE)) == 2
    # Gauge flux invariance independently checked under all 16 local choices.
    for gauge in GAUGES.values():
        for signs in product((-1,1), repeat=4):
            transformed = {key:signs[i-1]*gauge[key]*signs[j-1] for i,j,g,key,d1,d2 in GEOMETRY}
            assert fluxes(transformed) == fluxes(gauge)
    maxima = {}
    def check(name, left, right, tol=2e-11):
        error = float(np.max(np.abs(np.asarray(left)-np.asarray(right))))
        maxima[name] = max(maxima.get(name, 0.), error)
        if error > tol:
            raise AssertionError(f"{name}: {error} exceeds {tol}")
    rng = np.random.default_rng(seed)
    for case in range(random_cases):
        while True:
            K = rng.uniform(.15, 2., 5)
            if np.min(np.diff(np.sort(K))) > .025:
                break
        if case % 2:
            K *= rng.choice([-1,1], 5)
        q, q2 = rng.uniform(-3*np.pi, 3*np.pi, 2)
        U = basis_matrix(q2)
        check("unitarity", U.conj().T@U, np.eye(4))
        for gauge in GAUGES.values():
            raw = raw_from_bonds(q, q2, K, gauge)
            H = h_bloch(q, K, gauge, q2)
            E = np.linalg.eigvalsh(H)
            check("raw_formula_vs_geometry", raw, matrix_from_coefficients(coefficients(q,K,gauge,q2)))
            check("eq27", H, U.conj().T@raw@U)
            check("bloch_vs_closed", H, matrix_from_coefficients(coefficients(q,K,gauge)))
            check("q2_cancellation", H, h_bloch(q,K,gauge,0.))
            check("hermiticity", H, H.conj().T)
            check("majorana", h_bloch(-q,K,gauge), -H.conjugate())
            check("matrix_period_2pi", h_bloch(q+2*np.pi,K,gauge), H)
            check("matrix_covariance_pi", h_bloch(q+np.pi,K,gauge), G_PI@H@G_PI)
            check("spectrum_period_pi", np.linalg.eigvalsh(h_bloch(q+np.pi,K,gauge)), E)
            check("raw_bloch_spectrum", np.linalg.eigvalsh(raw), E)
            check("same_q_pairing", E, -E[::-1])
            poly = polynomial_coefficients(q,K,gauge)
            scale = max(1., float(np.max(np.abs(poly))))
            check("polynomial_relative", np.poly(H)/scale, poly/scale)
        difference = h_bloch(q,K,Z_REVERSED_GAUGE)-h_bloch(q,K,ORIGINAL_GAUGE)
        expected = np.zeros((4,4), dtype=complex)
        for i,j in ((0,2),(1,3)):
            expected[i,j] = 2j*K[2]*np.exp(-1j*q)
            expected[j,i] = expected[i,j].conjugate()
        check("only_z_matrix_difference", difference, expected)
    finite_details = []
    for ncell in (5, 8, 17):
        K = (1.13,.79,1.31,.92,.67)
        for name, gauge in GAUGES.items():
            A = finite_majorana_A(ncell,K,gauge)
            blocks = []
            for ell in range(-(ncell//2), ncell-(ncell//2)):
                q = np.pi*ell/ncell
                block = finite_fourier_block(A,q)
                check("finite_chain_bloch", block, h_bloch(q,K,gauge))
                check("finite_chain_raw", finite_fourier_block(A,q,.371), raw_from_bonds(q,.371,K,gauge))
                blocks.extend(np.linalg.eigvalsh(block))
            Ereal = np.linalg.eigvalsh(.5j*A)
            check("finite_chain_full_spectrum", np.sort(blocks), Ereal)
            Ecell = np.array(blocks)
            check("finite_chain_energy_per_site", Ereal[Ereal<0].sum()/(4*ncell), Ecell[Ecell<0].sum()/(4*ncell))
            finite_details.append({"sector": name, "cells": ncell, "sites": 4*ncell, "momenta": ncell})
    return {**symbolic, "random_cases": random_cases, "sectors_per_case": 2,
            "seed": seed, "asymmetric_couplings": True, "signed_cases": random_cases//2,
            "max_errors": maxima, "finite_chains": finite_details,
            "fluxes": {name:fluxes(gauge) for name,gauge in GAUGES.items()},
            "gauge_equivalence_solutions": [], "diagonal_gauges_tested": 16,
            "anti_contamination": "passed", "numpy_version": np.__version__, "sympy_version": sp.__version__}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--random-cases", type=int, default=1000)
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    result = self_test(args.random_cases)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir/"validation_z_reversed.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
