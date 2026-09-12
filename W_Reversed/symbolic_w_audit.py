#!/usr/bin/env python3
"""Independent exact audit, transcribed from PDF Eq.(11), gauge Eq.(15).

This audit deliberately owns a small source transcription and does not import
the production geometry or closed coefficients. PDF pages 3--5 were inspected
in auditoria_fontes/original-{3,4,5}.png. No previous sector defines W.
The primary integration and plotting interval is always [-pi, pi].
Run: python3 -W error symbolic_w_audit.py
"""
from __future__ import annotations

import json
from itertools import product

import numpy as np
import sympy as sp


# (start, end, family, original u, end-minus-start in a1, a2).
# This independent transcription is a test fixture, not production constants.
PDF_EQ11_BONDS = (
    (1, 2, "x", +1, +1, 0),
    (2, 3, "y", +1, 0, +1),
    (3, 4, "x", +1, -1, 0),
    (4, 1, "y", -1, 0, -1),
    (4, 2, "z", +1, +1, -1),
    (4, 3, "w", +1, -1, 0),
    (3, 1, "z", +1, +1, -1),
    (1, 2, "w", -1, -1, 0),
    (1, 3, "t", -1, +1, +1),
    (2, 4, "t", -1, +1, +1),
)


def require_zero(value, label):
    values = list(value) if isinstance(value, sp.MatrixBase) else [value]
    residuals = [sp.simplify(sp.expand(v)) for v in values]
    if any(v != 0 for v in residuals):
        raise AssertionError(f"{label}: {residuals}")


def symbolic_raw(q, q2, K, reverse_w=False, only_w=False):
    result = sp.zeros(4)
    for i, j, family, u, da1, da2 in PDF_EQ11_BONDS:
        if only_w and family != "w":
            continue
        if reverse_w and family == "w":
            u = -u
        term = sp.I * K[family] * u * sp.exp(sp.I * (q * da1 + q2 * da2))
        result[i-1, j-1] += term
        result[j-1, i-1] += sp.conjugate(term)
    return result


def run_symbolic_audit():
    q, q2 = sp.symbols("q q2", real=True)
    lam = sp.Symbol("lambda")
    Kx, Ky, Kz, Kt, Kw = sp.symbols("Kx Ky Kz Kt Kw", real=True)
    K = dict(zip(("x", "y", "z", "t", "w"), (Kx, Ky, Kz, Kt, Kw)))
    U = sp.diag(1, 1, sp.exp(-sp.I*q2), sp.exp(-sp.I*q2))
    require_zero(U.H * U - sp.eye(4), "Eq27 unitarity")
    G = sp.diag(1, -1, -1, 1)
    J = sp.Matrix([[0, 0, -1, 0], [0, 0, 0, -1],
                   [1, 0, 0, 0], [0, 1, 0, 0]])
    require_zero(J * J + sp.eye(4), "antiunitary square")
    derived = {}
    matrices = {}
    for label, reverse in (("OLD", False), ("W-REVERSED", True)):
        raw = symbolic_raw(q, q2, K, reverse)
        H = (U.H * raw * U).applyfunc(sp.simplify)
        matrices[label] = H
        require_zero(H.diff(q2), f"{label} q2 cancellation")
        require_zero(H.H - H, f"{label} Hermiticity")
        require_zero(H.subs(q, -q) + sp.conjugate(H), f"{label} Majorana")
        require_zero(H.subs(q, q + 2*sp.pi) - H, f"{label} 2pi period")
        require_zero(H.subs(q, q + sp.pi) - G * H * G.H, f"{label} pi covariance")
        require_zero(J * sp.conjugate(H) * J.H + H, f"{label} same-q pairing")
        require_zero(H[1, 2] - H[0, 3], f"{label} f23=f14")
        require_zero(H[1, 3] - H[0, 2], f"{label} f24=f13")
        require_zero(H[2, 3] + sp.conjugate(H[0, 1]), f"{label} f34=-conj(f12)")
        # RAW is checked before the basis change, including the sign of q2.
        require_zero(raw[0, 2] - sp.exp(sp.I*q2) * H[0, 2], f"{label} RAW f13 phase")
        require_zero(raw[1, 3] - sp.exp(sp.I*q2) * H[1, 3], f"{label} RAW f24 phase")
        require_zero(raw[0, 3] - sp.I*Ky*sp.exp(sp.I*q2), f"{label} RAW f14 phase")
        require_zero(raw[1, 2] - sp.I*Ky*sp.exp(sp.I*q2), f"{label} RAW f23 phase")
        wrong_basis = (U * raw * U.H).applyfunc(sp.simplify)
        if sp.simplify(wrong_basis[0, 3].diff(q2)) == 0:
            raise AssertionError("Wrong Eq27 must fail q2 independence")
        derived[label] = {f"f{i+1}{j+1}": str(H[i, j])
                          for i in range(4) for j in range(i+1, 4)}

    # Determinant is derived in independent real/imaginary coordinates.
    A, B, C, D, y = sp.symbols("A B C D y", real=True)
    a, b = A + sp.I*B, C + sp.I*D
    reduced = sp.Matrix([[0, a, b, sp.I*y],
                         [sp.conjugate(a), 0, sp.I*y, b],
                         [sp.conjugate(b), -sp.I*y, 0, -sp.conjugate(a)],
                         [-sp.I*y, sp.conjugate(b), -a, 0]])
    S = A*A+B*B+C*C+D*D+y*y
    R = B*B*(C*C+D*D)+D*D*y*y
    charpoly = reduced.charpoly()
    determinant = sp.Poly(charpoly.as_expr().subs(charpoly.gen, lam), lam)
    expected = lam**4 - 2*S*lam**2 + S**2 - 4*R
    require_zero(determinant.as_expr()-expected, "characteristic determinant")
    require_zero(J * sp.conjugate(reduced) * J.H + reduced, "general antiunitary")
    reductions = {}
    for label, w_sign in (("OLD", -1), ("W-REVERSED", +1)):
        substitutions = {
            A: -(Kx-w_sign*Kw)*sp.sin(q),
            B: (Kx+w_sign*Kw)*sp.cos(q),
            C: (Kt-Kz)*sp.sin(q),
            D: -(Kt+Kz)*sp.cos(q), y: Ky,
        }
        require_zero((matrices[label]-reduced.subs(substitutions)).applyfunc(sp.expand_complex),
                     f"{label} reduced matrix from bonds")
        actual_S = S.subs(substitutions)
        target_S = sum(value**2 for value in K.values()) + 2*(Kt*Kz+w_sign*Kx*Kw)*sp.cos(2*q)
        require_zero(sp.trigsimp(actual_S-target_S), f"{label} S")
        reductions[label] = {"S": str(target_S), "R": str(R.subs(substitutions))}

    # This difference uses only the two w bonds, never final f12/f34 formulas.
    delta_from_two_bonds = U.H * (symbolic_raw(q, q2, K, True, True)
                                 - symbolic_raw(q, q2, K, False, True)) * U
    require_zero(matrices["W-REVERSED"] - matrices["OLD"] - delta_from_two_bonds,
                 "DeltaH directly from the two w bonds")
    support = {(i+1, j+1) for i in range(4) for j in range(4)
               if sp.simplify(delta_from_two_bonds[i, j]) != 0}
    if support != {(1, 2), (2, 1), (3, 4), (4, 3)}:
        raise AssertionError(f"Unexpected DeltaH support {support}")
    literal_eq26 = sp.I*(Kx*sp.exp(sp.I*q)-Kw*sp.exp(-sp.I*q))
    source_discrepancy = sp.trigsimp(sp.expand_complex(matrices["OLD"][2, 3] - literal_eq26))
    require_zero(source_discrepancy - 2*(Kx+Kw)*sp.sin(q), "documented literal Eq26 discrepancy")

    isotropic = matrices["W-REVERSED"].subs({value: 1 for value in K.values()}).applyfunc(sp.simplify)
    isotropic = isotropic.applyfunc(sp.expand_complex).applyfunc(sp.simplify)
    c = sp.symbols("c", real=True)
    isotropic_c = isotropic.subs(sp.cos(q), c)
    iso_charpoly = isotropic_c.charpoly()
    iso_P = sp.Poly(iso_charpoly.as_expr().subs(iso_charpoly.gen, lam), lam)
    require_zero(iso_P.as_expr()-(lam**4-2*(1+8*c*c)*lam**2+1), "isotropic polynomial")
    require_zero(isotropic_c.det()-1, "isotropic determinant")
    # Each proposed eigenvalue is independently substituted into the determinant.
    for sign_a, sign_b in product((-1, 1), repeat=2):
        root = sign_a*(sp.sqrt(1+4*c*c) + sign_b*2*sp.Abs(c))
        require_zero(iso_P.as_expr().subs(lam, root), "closed isotropic root")
    H_func = sp.lambdify(q, isotropic, "numpy")
    special = []
    max_spectrum_error = 0.0
    momenta = (-sp.pi, -sp.pi/2, -sp.pi/4, sp.Integer(0), sp.pi/4, sp.pi/2, sp.pi)
    for momentum in momenta:
        c_exact = abs(sp.cos(momentum))
        big, small = sp.sqrt(1+4*c_exact*c_exact)+2*c_exact, sp.sqrt(1+4*c_exact*c_exact)-2*c_exact
        exact = [-big, -small, small, big]
        numeric = np.linalg.eigvalsh(np.asarray(H_func(float(momentum)), dtype=complex))
        error = float(np.max(np.abs(numeric - np.array([float(value) for value in exact]))))
        max_spectrum_error = max(max_spectrum_error, error)
        special.append({"q": str(momentum), "exact": [str(x) for x in exact],
                        "eigvalsh": numeric.tolist(), "max_abs_error": error})
    qgrid = np.linspace(-np.pi, np.pi, 1001)
    actual = np.linalg.eigvalsh(np.stack([np.asarray(H_func(v), dtype=complex) for v in qgrid]))
    cgrid = np.abs(np.cos(qgrid))
    large, small = np.sqrt(1+4*cgrid*cgrid)+2*cgrid, np.sqrt(1+4*cgrid*cgrid)-2*cgrid
    closed = np.column_stack((-large, -small, small, large))
    max_spectrum_error = max(max_spectrum_error, float(np.max(np.abs(actual-closed))))
    if max_spectrum_error > 2e-13:
        raise AssertionError(f"Isotropic spectrum error {max_spectrum_error}")
    # Only this comparison stage imports production; all reference quantities
    # above have already been reconstructed from the independent PDF fixture.
    import quadrupolar_spin_liquid_w_reversed as production
    independent_old = {f"u{i}{j}{family}": u for i, j, family, u, _, _ in PDF_EQ11_BONDS}
    independent_w = {f"u{i}{j}{family}": -u if family == "w" else u
                     for i, j, family, u, _, _ in PDF_EQ11_BONDS}
    if production.ORIGINAL_GAUGE != independent_old or production.W_REVERSED_GAUGE != independent_w:
        raise AssertionError("Production gauge differs from independent PDF transcription")
    rng = np.random.default_rng(20260910)
    production_errors = {"RAW": 0.0, "Bloch": 0.0, "DeltaH": 0.0}
    symbols = (q, q2, Kx, Ky, Kz, Kt, Kw)
    delta_func = sp.lambdify(symbols, delta_from_two_bonds, "numpy")
    for label, reverse in (("OLD", False), ("W-REVERSED", True)):
        raw_func = sp.lambdify(symbols, symbolic_raw(q, q2, K, reverse), "numpy")
        bloch_func = sp.lambdify(symbols, matrices[label], "numpy")
        gauge = production.W_REVERSED_GAUGE if reverse else production.ORIGINAL_GAUGE
        for _ in range(64):
            qval, q2val = rng.uniform(-np.pi, np.pi, 2)
            couplings = rng.uniform(-3.0, 3.0, 5)
            args = (qval, q2val, *couplings)
            ref_raw = np.asarray(raw_func(*args), dtype=complex)
            ref_bloch = np.asarray(bloch_func(*args), dtype=complex)
            actual_raw = production.raw_from_bonds(qval, q2val, couplings, gauge)
            actual_bloch = production.h_bloch(qval, couplings, gauge, q2val)
            ref_delta = np.asarray(delta_func(*args), dtype=complex)
            actual_delta = production.delta_from_w_bonds(qval, couplings)
            for key, residual in (("RAW", actual_raw-ref_raw), ("Bloch", actual_bloch-ref_bloch),
                                  ("DeltaH", actual_delta-ref_delta)):
                production_errors[key] = max(production_errors[key], float(np.max(np.abs(residual))))
    if max(production_errors.values()) > 2e-13:
        raise AssertionError(f"Independent PDF fixture vs production: {production_errors}")
    return {
        "status": "PASS", "primary_BZ": [-float(sp.pi), float(sp.pi)],
        "source": "Independent Eq.(11) bonds and Eq.(15) signs, PDF pp.3-5",
        "coefficients": derived, "reduced_forms": reductions,
        "characteristic_coefficients": [str(sp.factor(v)) for v in determinant.all_coeffs()],
        "isotropic_polynomial": str(iso_P.as_expr()),
        "isotropic_H": str(isotropic),
        "isotropic_min_abs_E": str(sp.sqrt(5)-2),
        "isotropic_zero_energy_band_gap": str(2*sp.sqrt(5)-4),
        "isotropic_degeneracies": "q = -pi/2, pi/2; eigenvalues -1,-1,1,1",
        "special_momenta": special,
        "max_isotropic_spectrum_error": max_spectrum_error,
        "independent_production_comparison": {"cases_per_sector": 64, "max_errors": production_errors},
        "delta_from_w_bonds": str(delta_from_two_bonds.applyfunc(sp.simplify)),
        "literal_eq26_difference": str(source_discrepancy),
        "normalization_formula": "E0/Nsite = integral[-pi,pi] F(q) dq / (8*pi)",
        "pairing": "J conjugate(H(q)) J^dagger = -H(q); J^2=-I",
        "pi_covariance": "H(q+pi) = diag(1,-1,-1,1) H(q) diag(1,-1,-1,1)",
    }


if __name__ == "__main__":
    print(json.dumps(run_symbolic_audit(), indent=2, ensure_ascii=False))
