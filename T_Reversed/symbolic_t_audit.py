#!/usr/bin/env python3
"""Independent exact audit, transcribed from PDF Eq.(11), gauge Eq.(15).

This audit deliberately owns a small source transcription and does not import
the production geometry or closed coefficients. PDF pages 3--5 were inspected
in auditoria_fontes/original-{3,4,5}.png. No previous modified sector defines T.
The primary integration and plotting interval is always [-pi, pi].
Run: python3 -W error symbolic_t_audit.py
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
        raise RuntimeError(f"{label}: {residuals}")


def symbolic_raw(q, q2, K, reverse_t=False, only_t=False):
    result = sp.zeros(4)
    for i, j, family, u, da1, da2 in PDF_EQ11_BONDS:
        if only_t and family != "t":
            continue
        if reverse_t and family == "t":
            u = -u
        term = sp.I * K[family] * u * sp.exp(sp.I * (q * da1 + q2 * da2))
        result[i-1, j-1] += term
        result[j-1, i-1] += sp.conjugate(term)
    return result


def run_symbolic_audit():
    """Derive references before importing the production module.

    A source transcription is intentionally duplicated only in this independent
    audit. Scientific production and all figures use the single core geometry.
    """
    q, q2 = sp.symbols("q q2", real=True)
    lam = sp.Symbol("lambda")
    Kx, Ky, Kz, Kt, Kw = sp.symbols("Kx Ky Kz Kt Kw", real=True)
    K = dict(zip(("x", "y", "z", "t", "w"), (Kx, Ky, Kz, Kt, Kw)))
    U = sp.diag(1, 1, sp.exp(-sp.I*q2), sp.exp(-sp.I*q2))
    require_zero(U.H*U-sp.eye(4), "Eq27 unitarity")
    G = sp.diag(1, -1, -1, 1)
    J = sp.Matrix([[0, 0, -1, 0], [0, 0, 0, -1],
                   [1, 0, 0, 0], [0, 1, 0, 0]])
    require_zero(J*J+sp.eye(4), "antiunitary square")
    derived, raw_derived, matrices = {}, {}, {}
    for label, reverse in (("OLD", False), ("T-REVERSED", True)):
        raw = symbolic_raw(q, q2, K, reverse)
        H = (U.H*raw*U).applyfunc(sp.simplify)
        matrices[label] = H
        require_zero(H.diff(q2), f"{label} q2 cancellation")
        require_zero(H.H-H, f"{label} Hermiticity")
        require_zero(H.subs(q, -q)+sp.conjugate(H), f"{label} Majorana")
        require_zero(H.subs(q, q+2*sp.pi)-H, f"{label} 2pi period")
        require_zero(H.subs(q, q+sp.pi)-G*H*G.H, f"{label} pi covariance")
        require_zero(J*sp.conjugate(H)*J.H+H, f"{label} same-q pairing")
        require_zero(H[1, 2]-H[0, 3], f"{label} f23=f14")
        require_zero(H[1, 3]-H[0, 2], f"{label} f24=f13")
        require_zero(H[2, 3]+sp.conjugate(H[0, 1]), f"{label} f34=-conj(f12)")
        for i, j in ((0, 2), (0, 3), (1, 2), (1, 3)):
            require_zero(raw[i, j]-sp.exp(sp.I*q2)*H[i, j],
                         f"{label} RAW independently tested positive q2 phase {i}{j}")
        # The y amplitudes anchor the phase without consulting final H.
        require_zero(raw[0, 3]-sp.I*Ky*sp.exp(sp.I*q2), f"{label} RAW y14")
        require_zero(raw[1, 2]-sp.I*Ky*sp.exp(sp.I*q2), f"{label} RAW y23")
        wrong_basis = U*raw*U.H
        if sp.simplify(wrong_basis[0, 3].diff(q2)) == 0:
            raise RuntimeError("Wrong Eq27 must fail q2 independence")
        derived[label] = {f"f{i+1}{j+1}": str(H[i, j])
                          for i in range(4) for j in range(i+1, 4)}
        raw_derived[label] = {f"f{i+1}{j+1}": str(raw[i, j])
                              for i in range(4) for j in range(i+1, 4)}

    # Compute the determinant in independent real coordinates. The exact
    # reduction back to the bond-derived matrices is checked below.
    A, B, C, D, y = sp.symbols("A B C D y", real=True)
    a, b = A+sp.I*B, C+sp.I*D
    reduced = sp.Matrix([[0, a, b, sp.I*y],
                         [sp.conjugate(a), 0, sp.I*y, b],
                         [sp.conjugate(b), -sp.I*y, 0, -sp.conjugate(a)],
                         [-sp.I*y, sp.conjugate(b), -a, 0]])
    S = A*A+B*B+C*C+D*D+y*y
    R = B*B*(C*C+D*D)+D*D*y*y
    determinant = sp.Poly(reduced.charpoly(lam).as_expr(), lam)
    expected = lam**4-2*S*lam**2+S**2-4*R
    require_zero(determinant.as_expr()-expected, "independent exact characteristic determinant")
    require_zero(determinant.as_expr().subs(lam, -lam)-determinant.as_expr(), "even polynomial")
    require_zero(J*sp.conjugate(reduced)*J.H+reduced, "general antiunitary pairing")
    # Polynomial in squared energies verifies the general nested-root formula.
    ss, rr = sp.symbols("S R", nonnegative=True)
    for sign_a, sign_b in product((-1, 1), repeat=2):
        root = sign_a*sp.sqrt(ss+sign_b*2*sp.sqrt(rr))
        require_zero(root**4-2*ss*root**2+ss**2-4*rr, "general exact spectral root")

    reductions, sector_polynomials = {}, {}
    for label, t_sign in (("OLD", -1), ("T-REVERSED", +1)):
        substitutions = {
            A: -(Kx+Kw)*sp.sin(q), B: (Kx-Kw)*sp.cos(q),
            C: -(t_sign*Kt+Kz)*sp.sin(q), D: (t_sign*Kt-Kz)*sp.cos(q), y: Ky,
        }
        require_zero((matrices[label]-reduced.subs(substitutions)).applyfunc(sp.expand_complex),
                     f"{label} reduced matrix independently recovered from bonds")
        actual_S, actual_R = S.subs(substitutions), R.subs(substitutions)
        compact_S = sum(value**2 for value in K.values())-2*(Kx*Kw+t_sign*Kt*Kz)*sp.cos(2*q)
        require_zero(sp.trigsimp(actual_S-compact_S), f"{label} compact S")
        sector_P = lam**4-2*compact_S*lam**2+compact_S**2-4*actual_R
        # Avoid hiding any odd coefficient by reporting every coefficient.
        sector_polynomials[label] = [str(v) for v in sp.Poly(sector_P, lam).all_coeffs()]
        reductions[label] = {"A": str(substitutions[A]), "B": str(substitutions[B]),
                             "C": str(substitutions[C]), "D": str(substitutions[D]),
                             "S": str(compact_S), "R": str(actual_R)}

    delta = (U.H*(symbolic_raw(q, q2, K, True, True)
                  -symbolic_raw(q, q2, K, False, True))*U).applyfunc(sp.simplify)
    require_zero(matrices["T-REVERSED"]-matrices["OLD"]-delta,
                 "DeltaH derived directly from only two real t bonds")
    support = {(i+1, j+1) for i in range(4) for j in range(4) if sp.simplify(delta[i, j]) != 0}
    if support != {(1, 3), (3, 1), (2, 4), (4, 2)}:
        raise RuntimeError(f"Unexpected DeltaH support: {support}")
    for i, j in ((0, 1), (0, 3), (1, 2), (2, 3)):
        require_zero(matrices["T-REVERSED"][i, j]-matrices["OLD"][i, j], "Non-t coefficient unchanged")
    literal_eq26 = sp.I*(Kx*sp.exp(sp.I*q)-Kw*sp.exp(-sp.I*q))
    discrepancy = sp.trigsimp(sp.expand_complex(matrices["OLD"][2, 3]-literal_eq26))
    require_zero(discrepancy-2*(Kx+Kw)*sp.sin(q), "documented literal Eq26 discrepancy")

    isotropic = {name: H.subs({value: 1 for value in K.values()}).applyfunc(sp.expand_complex).applyfunc(sp.simplify)
                 for name, H in matrices.items()}
    s, c = sp.symbols("s c", real=True)
    iso_t = isotropic["T-REVERSED"].subs(sp.sin(q), s)
    iso_poly = sp.Poly(iso_t.charpoly(lam).as_expr(), lam)
    require_zero(iso_poly.as_expr()-(lam**2-1-8*s*s)**2, "isotropic T determinant")
    require_zero(iso_t*iso_t-(1+8*s*s)*sp.eye(4), "all-q isotropic T double degeneracy")
    for sign in (-1, 1):
        require_zero(iso_poly.as_expr().subs(lam, sign*sp.sqrt(1+8*s*s)), "isotropic T exact root")
    iso_old_P = lam**4-10*lam**2+25-16*c*c
    old_reduction = {A: -2*s, B: 0, C: 0, D: -2*c, y: 1}
    old_actual_P = sp.Poly(determinant.as_expr().subs(old_reduction), lam).as_expr()
    # s^2+c^2=1; reduce in s rather than assume them independent.
    require_zero(sp.rem(sp.Poly(old_actual_P-iso_old_P, s), sp.Poly(s*s+c*c-1, s)).as_expr(),
                 "isotropic OLD polynomial from independent determinant")
    for sign_a, sign_b in product((-1, 1), repeat=2):
        require_zero(iso_old_P.subs(lam, sign_a*sp.sqrt(5+sign_b*4*c)), "OLD isotropic exact root")

    # Exact equality of integrals does not imply equality of spectra at a fixed q.
    # Each radicand is >=1 for real q, so no sign/absolute-value branch is lost.
    require_zero(sp.trigsimp(1+8*sp.sin(q)**2-(5-4*sp.cos(2*q))), "energy identity first radicand")
    require_zero(sp.trigsimp(1+8*sp.sin(q+sp.pi/2)**2-(5+4*sp.cos(2*q))),
                 "energy identity second radicand")
    energy_identity = "F_OLD(2q) = (F_T(q)+F_T(q+pi/2))/2; pi periodicity gives equal bulk integrals"
    exact_energy = -3*sp.elliptic_e(sp.Rational(8, 9))/sp.pi

    H_funcs = {name: sp.lambdify(q, H, "numpy") for name, H in isotropic.items()}
    special = {name: [] for name in matrices}
    max_iso_error = 0.0
    momenta = (-sp.pi, -sp.pi/2, -sp.pi/4, sp.Integer(0), sp.pi/4, sp.pi/2, sp.pi)
    for momentum in momenta:
        t_energy = sp.sqrt(1+8*sp.sin(momentum)**2)
        old_high = sp.sqrt(5+4*sp.Abs(sp.cos(momentum)))
        old_low = sp.sqrt(5-4*sp.Abs(sp.cos(momentum)))
        exact_by_sector = {"OLD": [-old_high, -old_low, old_low, old_high],
                           "T-REVERSED": [-t_energy, -t_energy, t_energy, t_energy]}
        for name, exact in exact_by_sector.items():
            numeric = np.linalg.eigvalsh(np.asarray(H_funcs[name](float(momentum)), dtype=complex))
            error = float(np.max(abs(numeric-np.array([float(v) for v in exact]))))
            max_iso_error = max(max_iso_error, error)
            special[name].append({"q": str(momentum), "exact": [str(v) for v in exact],
                                  "eigvalsh": numeric.tolist(), "max_abs_error": error})
    grid = np.linspace(-np.pi, np.pi, 1001)
    for name in matrices:
        actual = np.linalg.eigvalsh(np.stack([np.asarray(H_funcs[name](v), dtype=complex) for v in grid]))
        if name == "OLD":
            large, small = np.sqrt(5+4*abs(np.cos(grid))), np.sqrt(5-4*abs(np.cos(grid)))
        else:
            large = small = np.sqrt(1+8*np.sin(grid)**2)
        closed = np.column_stack((-large, -small, small, large))
        max_iso_error = max(max_iso_error, float(np.max(abs(actual-closed))))
    if max_iso_error > 2e-13:
        raise RuntimeError(f"Exact isotropic spectra disagree with eigvalsh: {max_iso_error}")

    # Production is first imported here; all reference expressions above were
    # derived from the independent documentary bond transcription.
    import quadrupolar_spin_liquid_t_reversed as production
    independent_old = {f"u{i}{j}{family}": u for i, j, family, u, _, _ in PDF_EQ11_BONDS}
    independent_t = {f"u{i}{j}{family}": -u if family == "t" else u
                     for i, j, family, u, _, _ in PDF_EQ11_BONDS}
    if production.ORIGINAL_GAUGE != independent_old or production.T_REVERSED_GAUGE != independent_t:
        raise RuntimeError("Production gauge differs from independent PDF transcription")
    production_geometry = {(i, j, g, da1, da2) for i, j, g, _, da1, da2 in production.GEOMETRY}
    reference_geometry = {(i, j, g, da1, da2) for i, j, g, _, da1, da2 in PDF_EQ11_BONDS}
    if production_geometry != reference_geometry or len(production.GEOMETRY) != len(PDF_EQ11_BONDS):
        raise RuntimeError("Production geometry differs from independent PDF transcription")
    symbols = (q, q2, Kx, Ky, Kz, Kt, Kw)
    delta_func = sp.lambdify(symbols, delta, "numpy")
    production_errors = {"RAW": 0.0, "Bloch": 0.0, "DeltaH": 0.0}
    rng = np.random.default_rng(20260911)
    for name, reverse in (("OLD", False), ("T-REVERSED", True)):
        raw_func = sp.lambdify(symbols, symbolic_raw(q, q2, K, reverse), "numpy")
        bloch_func = sp.lambdify(symbols, matrices[name], "numpy")
        gauge = independent_t if reverse else independent_old
        for _ in range(64):
            qval, q2val = rng.uniform(-np.pi, np.pi, 2)
            couplings = rng.uniform(-3., 3., 5)
            args = (qval, q2val, *couplings)
            residuals = {
                "RAW": production.raw_from_bonds(qval, q2val, couplings, gauge)-np.asarray(raw_func(*args), complex),
                "Bloch": production.h_bloch(qval, couplings, gauge, q2val)-np.asarray(bloch_func(*args), complex),
                "DeltaH": production.delta_from_t_bonds(qval, couplings)-np.asarray(delta_func(*args), complex),
            }
            for key, residual in residuals.items():
                production_errors[key] = max(production_errors[key], float(np.max(abs(residual))))
    if max(production_errors.values()) > 2e-13:
        raise RuntimeError(f"Independent PDF fixture vs production: {production_errors}")
    return {
        "status": "PASS", "primary_BZ": [-float(sp.pi), float(sp.pi)],
        "metadata": production.metadata(),
        "source": "Independent Eq.(11) bonds and Eq.(15) signs, PDF pp.3-5",
        "coefficients": derived, "raw_coefficients": raw_derived, "reduced_forms": reductions,
        "characteristic_coefficients": [str(sp.factor(v)) for v in determinant.all_coeffs()],
        "sector_characteristic_coefficients": sector_polynomials,
        "general_eigenvalues": ["-sqrt(S+2sqrt(R))", "-sqrt(S-2sqrt(R))", "sqrt(S-2sqrt(R))", "sqrt(S+2sqrt(R))"],
        "isotropic_polynomial": str(iso_poly.as_expr()), "isotropic_polynomial_factored": str(sp.factor(iso_poly.as_expr())),
        "isotropic_polynomial_coordinate": "s = sin(q)",
        "isotropic_H": str(isotropic["T-REVERSED"]), "isotropic_min_abs_E": "1",
        "isotropic_zero_energy_band_gap": "2", "isotropic_minima": "q=-pi,0,pi",
        "isotropic_degeneracies": "Twofold positive and twofold negative eigenvalues for every q",
        "isotropic_eigenvalues": ["-sqrt(1+8*sin(q)^2)"]*2+["sqrt(1+8*sin(q)^2)"]*2,
        "isotropic_OLD_polynomial": str(iso_old_P),
        "isotropic_exact_bulk_energy": str(exact_energy), "isotropic_exact_energy_equality": energy_identity,
        "energy_equality_proof": "Both radicand identities were symbolically proved. F has period pi. Integrate the identity over [-pi,pi]; the doubled argument and shifted term each preserve the integral. The common integral reduces to the complete elliptic integral E(8/9).",
        "same_q_isospectral": False, "special_momenta": special["T-REVERSED"],
        "special_momenta_by_sector": special, "max_isotropic_spectrum_error": max_iso_error,
        "independent_production_comparison": {"cases_per_sector": 64, "max_errors": production_errors},
        "delta_from_t_bonds": str(delta), "delta_support": sorted(support),
        "literal_eq26_difference": str(discrepancy),
        "normalization_formula": "E0/Nsite = integral[-pi,pi] F(q) dq / (8*pi)",
        "pairing": "J conjugate(H(q)) J^dagger = -H(q); J^2=-I",
        "pi_covariance": "H(q+pi) = diag(1,-1,-1,1) H(q) diag(1,-1,-1,1)",
    }


if __name__ == "__main__":
    print(json.dumps(run_symbolic_audit(), indent=2, ensure_ascii=False))
