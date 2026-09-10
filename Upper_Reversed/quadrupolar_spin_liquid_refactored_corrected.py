"""
Quadrupolar spin liquid in 1D — corrected right-oriented upper-chain gauge.

This module implements the full derivation pipeline without the accidental
cancellation of two intermediate sign/convention errors that was present in an
earlier version.

Physical gauge
--------------
Arrow i -> j  <=>  mu_ij = +1, with mu_ji = -mu_ij.

    mu_12,x = +1      mu_34,x = -1
    mu_23,y = +1      mu_41,y = -1
    mu_42,z = +1      mu_31,z = +1
    mu_43,w = -1      mu_12,w = -1
    mu_13,t = -1      mu_24,t = -1

The upper horizontal x and w links are reversed relative to the original
zero-flux representative.  This gauge is NOT changed anywhere in this file.

Correct derivation pipeline
---------------------------
    gauge / real-space bonds
        -> Fourier RAW matrix H_raw(q1,q2), with cross-chain phase exp(+i q2)
        -> Eq.(27): U(q2) = diag(1,1,exp(-i q2),exp(-i q2))
        -> active spinor convention Theta' = U Theta
        -> H_Bloch = U^dagger H_raw U
        -> q2-independent 1D Bloch Hamiltonian H(q), q = q1 = k.a1

For the physical gauge above the final independent coefficients are

    f12 =  i [Kx exp(+iq) - Kw exp(-iq)]
    f13 = -i [Kt exp(+iq) + Kz exp(-iq)]
    f14 =  i Ky
    f23 =  i Ky
    f24 = -i [Kt exp(+iq) + Kz exp(-iq)]
    f34 = -i [Kx exp(-iq) - Kw exp(+iq)]

and therefore

    f23 = f14,
    f24 = f13,
    f34 = conjugate(f12).

The physical primitive cell has length 2 a1.  In the Bloch gauge used here,
internal a1-coordinates are r = (0,1,1,0), which implies

    H(q + pi) = G_pi H(q) G_pi^dagger,
    G_pi = diag(1,-1,-1,1).

Only Python, NumPy and SymPy are required.
"""

from __future__ import annotations

from itertools import product
from typing import Mapping, NamedTuple

import numpy as np
import sympy as sp


# =============================================================================
# 1. Symbolic variables and constants
# =============================================================================

f12, f13, f14, f23, f24, f34 = sp.symbols(
    "f12 f13 f14 f23 f24 f34", complex=True
)
lambda_ = sp.symbols("lambda", real=True)
q1, q2 = sp.symbols("q1 q2", real=True)
Kx, Ky, Kz, Kt, Kw = sp.symbols("Kx Ky Kz Kt Kw", real=True)

COUPLING_SYMBOLS: dict[str, sp.Expr] = {
    "x": Kx,
    "y": Ky,
    "z": Kz,
    "t": Kt,
    "w": Kw,
}

G_PI_NUMERIC = np.diag([1.0, -1.0, -1.0, 1.0]).astype(np.complex128)
G_PI_SYMBOLIC = sp.diag(1, -1, -1, 1)
INTERNAL_A1_POSITIONS = np.array([0, 1, 1, 0], dtype=int)


# =============================================================================
# 2. Gauge configuration and antisymmetry helpers
# =============================================================================

RIGHT_GAUGE: dict[str, int] = {
    "mu12x": +1,
    "mu34x": -1,
    "mu23y": +1,
    "mu41y": -1,
    "mu42z": +1,
    "mu31z": +1,
    "mu43w": -1,
    "mu12w": -1,
    "mu13t": -1,
    "mu24t": -1,
}

OLD_GAUGE: dict[str, int] = {
    **RIGHT_GAUGE,
    "mu34x": +1,
    "mu43w": +1,
}


def validate_gauge(gauge: Mapping[str, int]) -> None:
    """Require every stored Z2 gauge variable to be exactly +1 or -1."""
    required = set(RIGHT_GAUGE)
    missing = required.difference(gauge)
    if missing:
        raise ValueError(f"Gauge is missing bond variables: {sorted(missing)}")

    invalid = {key: value for key, value in gauge.items() if value not in (-1, +1)}
    if invalid:
        raise ValueError(f"Gauge variables must be +/-1. Invalid entries: {invalid}")


def mu_from_stored(
    gauge: Mapping[str, int],
    stored_key: str,
    requested_i: int,
    requested_j: int,
) -> int:
    """Return mu_ij from a stored key muABgamma, using mu_ji = -mu_ij."""
    validate_gauge(gauge)
    if not stored_key.startswith("mu") or len(stored_key) < 5:
        raise ValueError(f"Invalid stored gauge key: {stored_key!r}")
    stored_i = int(stored_key[2])
    stored_j = int(stored_key[3])
    if (requested_i, requested_j) == (stored_i, stored_j):
        return int(gauge[stored_key])
    if (requested_i, requested_j) == (stored_j, stored_i):
        return -int(gauge[stored_key])
    raise ValueError(
        f"Requested link ({requested_i},{requested_j}) does not match {stored_key}."
    )


# =============================================================================
# 3. Gauge-invariant plaquette fluxes
# =============================================================================


def x_block_fluxes(gauge: Mapping[str, int] = RIGHT_GAUGE) -> dict[str, int]:
    """
    Counterclockwise triangle fluxes in the x block.

    Paths:
        1 -> 2 -> 4 -> 1
        2 -> 3 -> 4 -> 2
        1 -> 2 -> 3 -> 1
        1 -> 3 -> 4 -> 1
    """
    validate_gauge(gauge)
    return {
        "W124": (
            mu_from_stored(gauge, "mu12x", 1, 2)
            * mu_from_stored(gauge, "mu42z", 2, 4)
            * mu_from_stored(gauge, "mu41y", 4, 1)
        ),
        "W234": (
            mu_from_stored(gauge, "mu23y", 2, 3)
            * mu_from_stored(gauge, "mu34x", 3, 4)
            * mu_from_stored(gauge, "mu42z", 4, 2)
        ),
        "W123": (
            mu_from_stored(gauge, "mu12x", 1, 2)
            * mu_from_stored(gauge, "mu23y", 2, 3)
            * mu_from_stored(gauge, "mu13t", 3, 1)
        ),
        "W134": (
            mu_from_stored(gauge, "mu13t", 1, 3)
            * mu_from_stored(gauge, "mu34x", 3, 4)
            * mu_from_stored(gauge, "mu41y", 4, 1)
        ),
    }


def w_block_fluxes(gauge: Mapping[str, int] = RIGHT_GAUGE) -> dict[str, int]:
    """
    Counterclockwise triangle fluxes in the w block.

    Paths:
        2 -> 1 -> 3 -> 2
        1 -> 4 -> 3 -> 1
        2 -> 1 -> 4 -> 2
        2 -> 4 -> 3 -> 2

    Every reversed traversal is evaluated from mu_ji = -mu_ij.  No flux value
    is hard-coded.
    """
    validate_gauge(gauge)
    return {
        "W213": (
            mu_from_stored(gauge, "mu12w", 2, 1)
            * mu_from_stored(gauge, "mu31z", 1, 3)
            * mu_from_stored(gauge, "mu23y", 3, 2)
        ),
        "W143": (
            mu_from_stored(gauge, "mu41y", 1, 4)
            * mu_from_stored(gauge, "mu43w", 4, 3)
            * mu_from_stored(gauge, "mu31z", 3, 1)
        ),
        "W214": (
            mu_from_stored(gauge, "mu12w", 2, 1)
            * mu_from_stored(gauge, "mu41y", 1, 4)
            * mu_from_stored(gauge, "mu24t", 4, 2)
        ),
        "W243": (
            mu_from_stored(gauge, "mu24t", 2, 4)
            * mu_from_stored(gauge, "mu43w", 4, 3)
            * mu_from_stored(gauge, "mu23y", 3, 2)
        ),
    }


# =============================================================================
# 4. Generic 4x4 Hermitian Hamiltonians
# =============================================================================


def ham_qsl_symbolic(
    f12_: sp.Expr,
    f13_: sp.Expr,
    f14_: sp.Expr,
    f23_: sp.Expr,
    f24_: sp.Expr,
    f34_: sp.Expr,
) -> sp.Matrix:
    """Generic Hermitian 4x4 matrix in the (1,2,3,4) sublattice basis."""
    return sp.Matrix(
        [
            [0, f12_, f13_, f14_],
            [sp.conjugate(f12_), 0, f23_, f24_],
            [sp.conjugate(f13_), sp.conjugate(f23_), 0, f34_],
            [sp.conjugate(f14_), sp.conjugate(f24_), sp.conjugate(f34_), 0],
        ]
    )


def ham_qsl_numeric(
    f12_: complex,
    f13_: complex,
    f14_: complex,
    f23_: complex,
    f24_: complex,
    f34_: complex,
) -> np.ndarray:
    """Numerical Hermitian 4x4 matrix in the same basis and convention."""
    return np.array(
        [
            [0.0, f12_, f13_, f14_],
            [np.conj(f12_), 0.0, f23_, f24_],
            [np.conj(f13_), np.conj(f23_), 0.0, f34_],
            [np.conj(f14_), np.conj(f24_), np.conj(f34_), 0.0],
        ],
        dtype=np.complex128,
    )


H_GENERIC = ham_qsl_symbolic(f12, f13, f14, f23, f24, f34)


def generic_characteristic_polynomial() -> sp.Expr:
    """Return det(lambda I - H) for the fully generic Hermitian matrix."""
    return sp.factor((lambda_ * sp.eye(4) - H_GENERIC).det())


# =============================================================================
# 5. Explicit real-space bond geometry
# =============================================================================


class RealSpaceBond(NamedTuple):
    """One oriented hopping used to derive H_raw from geometry."""

    start: int
    end: int
    gamma: str
    mu: int
    da1: int
    da2: int


def real_space_bonds(
    gauge: Mapping[str, int] = RIGHT_GAUGE,
) -> tuple[RealSpaceBond, ...]:
    """
    Explicit bond list in a canonical upper-triangular orientation.

    For each bond, da1 and da2 are the spatial displacement components from
    `start` to `end` in the Fourier convention exp(+i k.d).  Gauge variables
    stored in the opposite orientation are converted only by mu_ji=-mu_ij.

    This list is the independent geometry source used by the generic RAW
    Fourier reconstruction and by the finite real-space chain.
    """
    validate_gauge(gauge)
    return (
        RealSpaceBond(1, 2, "x", mu_from_stored(gauge, "mu12x", 1, 2), +1, 0),
        RealSpaceBond(1, 2, "w", mu_from_stored(gauge, "mu12w", 1, 2), -1, 0),
        RealSpaceBond(1, 3, "t", mu_from_stored(gauge, "mu13t", 1, 3), +1, +1),
        RealSpaceBond(1, 3, "z", mu_from_stored(gauge, "mu31z", 1, 3), -1, +1),
        RealSpaceBond(1, 4, "y", mu_from_stored(gauge, "mu41y", 1, 4), 0, +1),
        RealSpaceBond(2, 3, "y", mu_from_stored(gauge, "mu23y", 2, 3), 0, +1),
        RealSpaceBond(2, 4, "t", mu_from_stored(gauge, "mu24t", 2, 4), +1, +1),
        RealSpaceBond(2, 4, "z", mu_from_stored(gauge, "mu42z", 2, 4), -1, +1),
        RealSpaceBond(3, 4, "x", mu_from_stored(gauge, "mu34x", 3, 4), -1, 0),
        RealSpaceBond(3, 4, "w", mu_from_stored(gauge, "mu43w", 3, 4), +1, 0),
    )


def _numeric_coupling_map(
    Kx_: float,
    Ky_: float,
    Kz_: float,
    Kt_: float,
    Kw_: float,
) -> dict[str, float]:
    return {
        "x": float(Kx_),
        "y": float(Ky_),
        "z": float(Kz_),
        "t": float(Kt_),
        "w": float(Kw_),
    }


def raw_hamiltonian_from_bonds_symbolic(
    gauge: Mapping[str, int] = RIGHT_GAUGE,
) -> sp.Matrix:
    """
    Derive H_raw(q1,q2) generically from the explicit real-space bond list.

    Each oriented bond contributes
        i K_gamma mu_ij exp[i(da1 q1 + da2 q2)]
    to the corresponding upper-triangular matrix element; Hermiticity supplies
    the reverse entry.  No closed f_mn formula is reused here.
    """
    H = sp.zeros(4)
    for bond in real_space_bonds(gauge):
        phase = sp.exp(sp.I * (bond.da1 * q1 + bond.da2 * q2))
        amplitude = sp.I * COUPLING_SYMBOLS[bond.gamma] * bond.mu * phase
        i = bond.start - 1
        j = bond.end - 1
        H[i, j] += amplitude
        H[j, i] += sp.conjugate(amplitude)
    return H.applyfunc(sp.simplify)


def raw_hamiltonian_from_bonds_numeric(
    q1_: float,
    q2_: float,
    Kx_: float,
    Ky_: float,
    Kz_: float,
    Kt_: float,
    Kw_: float,
    gauge: Mapping[str, int] = RIGHT_GAUGE,
) -> np.ndarray:
    """Numerical generic bond-list reconstruction of H_raw(q1,q2)."""
    couplings = _numeric_coupling_map(Kx_, Ky_, Kz_, Kt_, Kw_)
    H = np.zeros((4, 4), dtype=np.complex128)
    for bond in real_space_bonds(gauge):
        phase = np.exp(1j * (bond.da1 * q1_ + bond.da2 * q2_))
        amplitude = 1j * couplings[bond.gamma] * bond.mu * phase
        i = bond.start - 1
        j = bond.end - 1
        H[i, j] += amplitude
        H[j, i] += np.conj(amplitude)
    return H


# =============================================================================
# 6. Closed RAW Fourier formulas — BEFORE Eq.(27)
# =============================================================================


def raw_fourier_coefficients_symbolic(
    gauge: Mapping[str, int] = RIGHT_GAUGE,
) -> dict[str, sp.Expr]:
    """
    Six independent coefficients of H_raw(q1,q2).

    The correct Fourier convention has a common exp(+i q2) factor in all four
    lower/upper cross-chain coefficients f13, f14, f23 and f24.
    """
    validate_gauge(gauge)
    phase2 = sp.exp(+sp.I * q2)
    return {
        "f12": sp.I
        * (
            gauge["mu12x"] * Kx * sp.exp(+sp.I * q1)
            + gauge["mu12w"] * Kw * sp.exp(-sp.I * q1)
        ),
        "f13": phase2
        * (
            sp.I * gauge["mu13t"] * Kt * sp.exp(+sp.I * q1)
            - sp.I * gauge["mu31z"] * Kz * sp.exp(-sp.I * q1)
        ),
        "f14": phase2 * (-sp.I * gauge["mu41y"] * Ky),
        "f23": phase2 * (+sp.I * gauge["mu23y"] * Ky),
        "f24": phase2
        * (
            sp.I * gauge["mu24t"] * Kt * sp.exp(+sp.I * q1)
            - sp.I * gauge["mu42z"] * Kz * sp.exp(-sp.I * q1)
        ),
        "f34": sp.I
        * (
            gauge["mu34x"] * Kx * sp.exp(-sp.I * q1)
            - gauge["mu43w"] * Kw * sp.exp(+sp.I * q1)
        ),
    }


def raw_fourier_coefficients_numeric(
    q1_: float,
    q2_: float,
    Kx_: float,
    Ky_: float,
    Kz_: float,
    Kt_: float,
    Kw_: float,
    gauge: Mapping[str, int] = RIGHT_GAUGE,
) -> dict[str, complex]:
    """Numerical closed RAW coefficients with the required exp(+i q2) phase."""
    validate_gauge(gauge)
    phase2 = np.exp(+1j * q2_)
    return {
        "f12": 1j
        * (
            gauge["mu12x"] * Kx_ * np.exp(+1j * q1_)
            + gauge["mu12w"] * Kw_ * np.exp(-1j * q1_)
        ),
        "f13": phase2
        * (
            1j * gauge["mu13t"] * Kt_ * np.exp(+1j * q1_)
            - 1j * gauge["mu31z"] * Kz_ * np.exp(-1j * q1_)
        ),
        "f14": phase2 * (-1j * gauge["mu41y"] * Ky_),
        "f23": phase2 * (+1j * gauge["mu23y"] * Ky_),
        "f24": phase2
        * (
            1j * gauge["mu24t"] * Kt_ * np.exp(+1j * q1_)
            - 1j * gauge["mu42z"] * Kz_ * np.exp(-1j * q1_)
        ),
        "f34": 1j
        * (
            gauge["mu34x"] * Kx_ * np.exp(-1j * q1_)
            - gauge["mu43w"] * Kw_ * np.exp(+1j * q1_)
        ),
    }


def raw_fourier_hamiltonian_symbolic(
    gauge: Mapping[str, int] = RIGHT_GAUGE,
) -> sp.Matrix:
    """Closed-form H_raw immediately after Fourier transform and before Eq.(27)."""
    c = raw_fourier_coefficients_symbolic(gauge)
    return ham_qsl_symbolic(c["f12"], c["f13"], c["f14"], c["f23"], c["f24"], c["f34"])


def h_raw_numeric(
    q1_: float,
    q2_: float,
    Kx_: float,
    Ky_: float,
    Kz_: float,
    Kt_: float,
    Kw_: float,
    gauge: Mapping[str, int] = RIGHT_GAUGE,
) -> np.ndarray:
    """Closed-form numerical H_raw for an arbitrary valid gauge."""
    c = raw_fourier_coefficients_numeric(
        q1_, q2_, Kx_, Ky_, Kz_, Kt_, Kw_, gauge
    )
    return ham_qsl_numeric(c["f12"], c["f13"], c["f14"], c["f23"], c["f24"], c["f34"])


def h_right_raw(
    q1_: float,
    q2_: float,
    Kx_: float,
    Ky_: float,
    Kz_: float,
    Kt_: float,
    Kw_: float,
) -> np.ndarray:
    """H_raw for the physical right-oriented gauge."""
    return h_raw_numeric(q1_, q2_, Kx_, Ky_, Kz_, Kt_, Kw_, RIGHT_GAUGE)


RAW_COEFFS_SYMBOLIC = {
    name: sp.simplify(expr)
    for name, expr in raw_fourier_coefficients_symbolic(RIGHT_GAUGE).items()
}
H_RAW_SYMBOLIC = raw_fourier_hamiltonian_symbolic(RIGHT_GAUGE)


# =============================================================================
# 7. Eq.(27): active Bloch-basis transformation Theta' = U Theta
# =============================================================================


def bloch_basis_matrix_symbolic(q2_: sp.Expr = q2) -> sp.Matrix:
    """Eq.(27): U = diag(1,1,exp(-i q2),exp(-i q2))."""
    return sp.diag(1, 1, sp.exp(-sp.I * q2_), sp.exp(-sp.I * q2_))


def bloch_basis_matrix_numeric(q2_: float) -> np.ndarray:
    """Numerical Eq.(27): U = diag(1,1,e^{-iq2},e^{-iq2})."""
    phase = np.exp(-1j * q2_)
    return np.diag([1.0, 1.0, phase, phase]).astype(np.complex128)


def apply_bloch_basis_transform_symbolic(
    H_raw: sp.Matrix,
    q2_: sp.Expr = q2,
) -> sp.Matrix:
    """
    Apply Eq.(27) in the active convention Theta' = U Theta.

    Starting from Theta^dagger H_raw Theta and inserting Theta' = U Theta gives

        H_Bloch = U^dagger H_raw U.
    """
    U = bloch_basis_matrix_symbolic(q2_)
    return (U.conjugate().T * H_raw * U).applyfunc(sp.simplify)


def apply_bloch_basis_transform_numeric(H_raw: np.ndarray, q2_: float) -> np.ndarray:
    """Numerical active Eq.(27): H_Bloch = U.conj().T @ H_raw @ U."""
    U = bloch_basis_matrix_numeric(q2_)
    return U.conj().T @ H_raw @ U


def h_bloch_for_gauge(
    q1_: float,
    Kx_: float,
    Ky_: float,
    Kz_: float,
    Kt_: float,
    Kw_: float,
    q2_: float = 0.0,
    gauge: Mapping[str, int] = RIGHT_GAUGE,
) -> np.ndarray:
    """Final Bloch Hamiltonian obtained strictly as U^dagger H_raw U."""
    H_raw = h_raw_numeric(q1_, q2_, Kx_, Ky_, Kz_, Kt_, Kw_, gauge)
    return apply_bloch_basis_transform_numeric(H_raw, q2_)


def h_right(
    q1_: float,
    Kx_: float,
    Ky_: float,
    Kz_: float,
    Kt_: float,
    Kw_: float,
    q2_: float = 0.0,
) -> np.ndarray:
    """Final q2-independent Bloch Hamiltonian for the physical RIGHT_GAUGE."""
    return h_bloch_for_gauge(
        q1_, Kx_, Ky_, Kz_, Kt_, Kw_, q2_=q2_, gauge=RIGHT_GAUGE
    )


def h_old(
    q1_: float,
    Kx_: float,
    Ky_: float,
    Kz_: float,
    Kt_: float,
    Kw_: float,
    q2_: float = 0.0,
) -> np.ndarray:
    """Old zero-flux representative used only for differential audit tests."""
    return h_bloch_for_gauge(
        q1_, Kx_, Ky_, Kz_, Kt_, Kw_, q2_=q2_, gauge=OLD_GAUGE
    )


H_RIGHT_SYMBOLIC = apply_bloch_basis_transform_symbolic(H_RAW_SYMBOLIC, q2)


def extract_independent_coefficients_symbolic(H: sp.Matrix) -> dict[str, sp.Expr]:
    """Extract the six independent upper-triangular entries of a 4x4 matrix."""
    return {
        "f12": sp.simplify(H[0, 1]),
        "f13": sp.simplify(H[0, 2]),
        "f14": sp.simplify(H[0, 3]),
        "f23": sp.simplify(H[1, 2]),
        "f24": sp.simplify(H[1, 3]),
        "f34": sp.simplify(H[2, 3]),
    }


def extract_independent_coefficients_numeric(H: np.ndarray) -> dict[str, complex]:
    """Extract the six independent upper-triangular entries of numerical H."""
    return {
        "f12": complex(H[0, 1]),
        "f13": complex(H[0, 2]),
        "f14": complex(H[0, 3]),
        "f23": complex(H[1, 2]),
        "f24": complex(H[1, 3]),
        "f34": complex(H[2, 3]),
    }


RIGHT_COEFFS_SYMBOLIC = extract_independent_coefficients_symbolic(H_RIGHT_SYMBOLIC)


# =============================================================================
# 8. Final closed Bloch coefficients
# =============================================================================


def expected_bloch_coefficients_symbolic(
    gauge: Mapping[str, int] = RIGHT_GAUGE,
) -> dict[str, sp.Expr]:
    """Closed post-Eq.(27) formulas used as an audit target, not as RAW proof."""
    validate_gauge(gauge)
    return {
        "f12": sp.I
        * (
            gauge["mu12x"] * Kx * sp.exp(+sp.I * q1)
            + gauge["mu12w"] * Kw * sp.exp(-sp.I * q1)
        ),
        "f13": (
            sp.I * gauge["mu13t"] * Kt * sp.exp(+sp.I * q1)
            - sp.I * gauge["mu31z"] * Kz * sp.exp(-sp.I * q1)
        ),
        "f14": -sp.I * gauge["mu41y"] * Ky,
        "f23": +sp.I * gauge["mu23y"] * Ky,
        "f24": (
            sp.I * gauge["mu24t"] * Kt * sp.exp(+sp.I * q1)
            - sp.I * gauge["mu42z"] * Kz * sp.exp(-sp.I * q1)
        ),
        "f34": sp.I
        * (
            gauge["mu34x"] * Kx * sp.exp(-sp.I * q1)
            - gauge["mu43w"] * Kw * sp.exp(+sp.I * q1)
        ),
    }


def expected_bloch_coefficients_numeric(
    q1_: float,
    Kx_: float,
    Ky_: float,
    Kz_: float,
    Kt_: float,
    Kw_: float,
    gauge: Mapping[str, int] = RIGHT_GAUGE,
) -> dict[str, complex]:
    """Numerical closed post-Eq.(27) coefficient formulas."""
    validate_gauge(gauge)
    return {
        "f12": 1j
        * (
            gauge["mu12x"] * Kx_ * np.exp(+1j * q1_)
            + gauge["mu12w"] * Kw_ * np.exp(-1j * q1_)
        ),
        "f13": (
            1j * gauge["mu13t"] * Kt_ * np.exp(+1j * q1_)
            - 1j * gauge["mu31z"] * Kz_ * np.exp(-1j * q1_)
        ),
        "f14": -1j * gauge["mu41y"] * Ky_,
        "f23": +1j * gauge["mu23y"] * Ky_,
        "f24": (
            1j * gauge["mu24t"] * Kt_ * np.exp(+1j * q1_)
            - 1j * gauge["mu42z"] * Kz_ * np.exp(-1j * q1_)
        ),
        "f34": 1j
        * (
            gauge["mu34x"] * Kx_ * np.exp(-1j * q1_)
            - gauge["mu43w"] * Kw_ * np.exp(+1j * q1_)
        ),
    }


def bloch_coefficients_numeric(
    q1_: float,
    Kx_: float,
    Ky_: float,
    Kz_: float,
    Kt_: float,
    Kw_: float,
    q2_: float = 0.0,
) -> dict[str, complex]:
    """Final Bloch coefficients extracted from the transformed Hamiltonian."""
    return extract_independent_coefficients_numeric(
        h_right(q1_, Kx_, Ky_, Kz_, Kt_, Kw_, q2_=q2_)
    )


EXPECTED_RIGHT_COEFFS_SYMBOLIC = {
    name: sp.simplify(expr)
    for name, expr in expected_bloch_coefficients_symbolic(RIGHT_GAUGE).items()
}


# =============================================================================
# 9. Reduced matrix and characteristic polynomial
# =============================================================================


def right_reduced_hamiltonian() -> sp.Matrix:
    """Reduced structure with f23=f14, f24=f13, f34=conjugate(f12)."""
    return ham_qsl_symbolic(f12, f13, f14, f14, f13, sp.conjugate(f12))


H_RIGHT_REDUCED = right_reduced_hamiltonian()


def right_reduced_characteristic_polynomial() -> sp.Expr:
    """Exact characteristic polynomial of the reduced symbolic matrix."""
    return sp.factor((lambda_ * sp.eye(4) - H_RIGHT_REDUCED).det())


def right_reduced_characteristic_polynomial_compact() -> sp.Expr:
    """
    Compact P(lambda) for a=f12, b=f13, c=f14 and f34=conjugate(a).

    P = lambda^4 - 2 A lambda^2 - 8 Re(a) Re(b c*) lambda
        + A^2 - 4[|a|^2|c|^2 + Re(a)^2|b|^2 + Re(b c*)^2],
    A = |a|^2 + |b|^2 + |c|^2.
    """
    a, b, c = f12, f13, f14
    A = sp.Abs(a) ** 2 + sp.Abs(b) ** 2 + sp.Abs(c) ** 2
    re_a = sp.re(a)
    re_bc = sp.re(b * sp.conjugate(c))
    return sp.expand(
        lambda_**4
        - 2 * A * lambda_**2
        - 8 * re_a * re_bc * lambda_
        + A**2
        - 4
        * (
            sp.Abs(a) ** 2 * sp.Abs(c) ** 2
            + re_a**2 * sp.Abs(b) ** 2
            + re_bc**2
        )
    )


def characteristic_polynomial_real_component_test() -> sp.Expr:
    """
    Return the fully simplified symbolic difference P_exact-P_compact.

    The six component symbols are explicitly real, avoiding assumptions about
    generic complex SymPy symbols.  A correct implementation returns exactly 0.
    """
    ar, ai, br, bi, cr, ci = sp.symbols(
        "ar ai br bi cr ci", real=True
    )
    a = ar + sp.I * ai
    b = br + sp.I * bi
    c = cr + sp.I * ci
    H = ham_qsl_symbolic(a, b, c, c, b, sp.conjugate(a))
    P_exact = sp.expand((lambda_ * sp.eye(4) - H).det())
    A = sp.Abs(a) ** 2 + sp.Abs(b) ** 2 + sp.Abs(c) ** 2
    P_compact = (
        lambda_**4
        - 2 * A * lambda_**2
        - 8 * sp.re(a) * sp.re(b * sp.conjugate(c)) * lambda_
        + A**2
        - 4
        * (
            sp.Abs(a) ** 2 * sp.Abs(c) ** 2
            + sp.re(a) ** 2 * sp.Abs(b) ** 2
            + sp.re(b * sp.conjugate(c)) ** 2
        )
    )
    return sp.simplify(sp.expand_complex(sp.expand(P_exact - P_compact)))


# =============================================================================
# 10. Diagonalization and bands
# =============================================================================


def eigensystem_numeric(H: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Hermitian eigenvalues/eigenvectors, with eigenvalues in ascending order."""
    H = np.asarray(H, dtype=np.complex128)
    if H.shape != (4, 4):
        raise ValueError(f"Expected a 4x4 Hamiltonian, got shape {H.shape}.")
    if not np.allclose(H, H.conj().T, atol=1e-12, rtol=1e-12):
        raise ValueError("Hamiltonian is not Hermitian within numerical tolerance.")
    return np.linalg.eigh(H)


def energies_right(
    q1_: float,
    Kx_: float,
    Ky_: float,
    Kz_: float,
    Kt_: float,
    Kw_: float,
) -> np.ndarray:
    """Four sorted real bands E_1(q) <= ... <= E_4(q)."""
    return np.linalg.eigvalsh(h_right(q1_, Kx_, Ky_, Kz_, Kt_, Kw_))


def band_structure_right(
    q_values: np.ndarray,
    Kx_: float,
    Ky_: float,
    Kz_: float,
    Kt_: float,
    Kw_: float,
) -> np.ndarray:
    """Return a band array with shape (len(q_values), 4)."""
    q_values = np.asarray(q_values, dtype=float)
    if q_values.ndim != 1 or q_values.size == 0:
        raise ValueError("q_values must be a non-empty one-dimensional array.")
    return np.vstack(
        [energies_right(q, Kx_, Ky_, Kz_, Kt_, Kw_) for q in q_values]
    )


def occupied_energy_sum(
    q_: float,
    Kx_: float,
    Ky_: float,
    Kz_: float,
    Kt_: float,
    Kw_: float,
) -> float:
    """F(q) = sum of all strictly negative eigenvalues of h(q)."""
    eigs = energies_right(q_, Kx_, Ky_, Kz_, Kt_, Kw_)
    return float(eigs[eigs < 0.0].sum())


def _validated_integration_grid(
    q_values: np.ndarray | None,
    lower: float,
    upper: float,
    num_points: int,
    name: str,
) -> np.ndarray:
    if q_values is None:
        if num_points < 3:
            raise ValueError("num_points must be at least 3.")
        return np.linspace(lower, upper, num_points, dtype=float)

    q = np.asarray(q_values, dtype=float)
    if q.ndim != 1 or q.size < 3:
        raise ValueError(f"{name}: q_values must be a 1D array with at least 3 points.")
    if not np.all(np.isfinite(q)):
        raise ValueError(f"{name}: q_values contains non-finite entries.")
    if not np.all(np.diff(q) > 0.0):
        raise ValueError(f"{name}: q_values must be strictly increasing.")
    if q[0] < lower - 1e-12 or q[-1] > upper + 1e-12:
        raise ValueError(f"{name}: q_values must stay inside [{lower}, {upper}].")
    if not np.isclose(q[0], lower, atol=1e-12, rtol=0.0) or not np.isclose(
        q[-1], upper, atol=1e-12, rtol=0.0
    ):
        raise ValueError(
            f"{name}: integration grid must span the full closed numerical "
            f"representative [{lower}, {upper}]."
        )
    return q


def ground_state_energy_per_cell_physical_bz(
    Kx_: float,
    Ky_: float,
    Kz_: float,
    Kt_: float,
    Kw_: float,
    q_values: np.ndarray | None = None,
    num_points: int = 4001,
) -> float:
    """
    Ground-state energy per four-site cell over the physical BZ.

        E_cell = (1/pi) integral_{-pi/2}^{pi/2} F(q) dq,
        F(q) = sum_{E_n(q)<0} E_n(q).

    The endpoint pi/2 is included only as a measure-zero numerical quadrature
    representative of the half-open physical BZ [-pi/2, pi/2).
    """
    q = _validated_integration_grid(
        q_values, -np.pi / 2.0, np.pi / 2.0, num_points,
        "ground_state_energy_per_cell_physical_bz"
    )
    F = np.array(
        [occupied_energy_sum(x, Kx_, Ky_, Kz_, Kt_, Kw_) for x in q],
        dtype=float,
    )
    return float(np.trapezoid(F, q) / np.pi)


def ground_state_energy_per_cell_extended_bz(
    Kx_: float,
    Ky_: float,
    Kz_: float,
    Kt_: float,
    Kw_: float,
    q_values: np.ndarray | None = None,
    num_points: int = 8001,
) -> float:
    """
    Equivalent extended-zone energy per cell.

        E_cell = (1/(2 pi)) integral_{-pi}^{pi} F(q) dq.
    """
    q = _validated_integration_grid(
        q_values, -np.pi, np.pi, num_points,
        "ground_state_energy_per_cell_extended_bz"
    )
    F = np.array(
        [occupied_energy_sum(x, Kx_, Ky_, Kz_, Kt_, Kw_) for x in q],
        dtype=float,
    )
    return float(np.trapezoid(F, q) / (2.0 * np.pi))


def ground_state_energy_per_site(
    Kx_: float,
    Ky_: float,
    Kz_: float,
    Kt_: float,
    Kw_: float,
    q_values: np.ndarray | None = None,
    num_points: int = 4001,
) -> float:
    """Ground-state energy per site: E_site = E_cell / 4."""
    return ground_state_energy_per_cell_physical_bz(
        Kx_, Ky_, Kz_, Kt_, Kw_, q_values=q_values, num_points=num_points
    ) / 4.0


# =============================================================================
# 11. Independent finite periodic Majorana chain
# =============================================================================


def _cell_shift_from_a1_displacement(bond: RealSpaceBond) -> int:
    """
    Convert physical da1 into a cell-index shift for a cell of length 2 a1.

    Internal coordinates are r=(0,1,1,0), so
        da1 = 2 Delta_cell + r_end - r_start.
    """
    r_i = int(INTERNAL_A1_POSITIONS[bond.start - 1])
    r_j = int(INTERNAL_A1_POSITIONS[bond.end - 1])
    numerator = bond.da1 - (r_j - r_i)
    if numerator % 2 != 0:
        raise ValueError(f"Bond {bond} is incompatible with the 2 a1 cell geometry.")
    return numerator // 2


def finite_majorana_A(
    ncell: int,
    Kx_: float,
    Ky_: float,
    Kz_: float,
    Kt_: float,
    Kw_: float,
    gauge: Mapping[str, int] = RIGHT_GAUGE,
) -> np.ndarray:
    """
    Construct the real antisymmetric Majorana matrix A for a periodic chain.

    Convention:
        H = (i/4) gamma^T A gamma,
        h = i A / 2.

    Consequently each oriented real-space bond contributes
        A_ij = 2 K_gamma mu_ij,
        A_ji = -A_ij.
    """
    if not isinstance(ncell, int) or ncell < 3:
        raise ValueError("ncell must be an integer >= 3.")
    couplings = _numeric_coupling_map(Kx_, Ky_, Kz_, Kt_, Kw_)
    A = np.zeros((4 * ncell, 4 * ncell), dtype=float)

    for cell in range(ncell):
        for bond in real_space_bonds(gauge):
            shift = _cell_shift_from_a1_displacement(bond)
            target_cell = (cell + shift) % ncell
            i = 4 * cell + (bond.start - 1)
            j = 4 * target_cell + (bond.end - 1)
            value = 2.0 * couplings[bond.gamma] * bond.mu
            A[i, j] += value
            A[j, i] -= value

    if not np.allclose(A.T, -A, atol=1e-14, rtol=0.0):
        raise AssertionError("Constructed finite Majorana matrix is not antisymmetric.")
    return A


def finite_chain_bloch_block(
    A: np.ndarray,
    q_: float,
    ncell: int,
) -> np.ndarray:
    """
    Fourier-project h=iA/2 onto one allowed momentum of the finite chain.

    The basis uses the actual internal coordinates x=2R+r_m, so this test does
    not insert any closed f_mn expression or any G_pi relation by hand.
    """
    A = np.asarray(A, dtype=float)
    if A.shape != (4 * ncell, 4 * ncell):
        raise ValueError(
            f"A must have shape {(4*ncell, 4*ncell)}, got {A.shape}."
        )
    h_real = 0.5j * A.astype(np.complex128)
    F = np.zeros((4 * ncell, 4), dtype=np.complex128)
    norm = np.sqrt(float(ncell))
    for cell in range(ncell):
        for m in range(4):
            x = 2 * cell + int(INTERNAL_A1_POSITIONS[m])
            F[4 * cell + m, m] = np.exp(1j * q_ * x) / norm
    return F.conj().T @ h_real @ F


def finite_chain_bloch_audit(
    ncell: int = 17,
    couplings: tuple[float, float, float, float, float] = (1.13, 0.79, 1.31, 0.92, 0.67),
) -> float:
    """Compare finite-chain Fourier blocks with H_old(q) and H_new(q)."""
    A_new = finite_majorana_A(ncell, *couplings, gauge=RIGHT_GAUGE)
    A_old = finite_majorana_A(ncell, *couplings, gauge=OLD_GAUGE)
    max_error = 0.0

    # For odd ncell, these are the ncell momenta in the physical BZ.
    half = ncell // 2
    momentum_indices = range(-half, half + 1)
    if len(tuple(momentum_indices)) != ncell:
        raise AssertionError("Finite-chain momentum construction requires odd ncell.")

    for ell in momentum_indices:
        q_ = np.pi * ell / ncell
        H_new_from_A = finite_chain_bloch_block(A_new, q_, ncell)
        H_old_from_A = finite_chain_bloch_block(A_old, q_, ncell)
        err_new = np.max(np.abs(H_new_from_A - h_right(q_, *couplings)))
        err_old = np.max(np.abs(H_old_from_A - h_old(q_, *couplings)))
        max_error = max(max_error, float(err_new), float(err_old))

    return max_error


# =============================================================================
# 12. Gauge-equivalence and differential old/new audits
# =============================================================================


def diagonal_z2_gauge_equivalence(
    q_values: tuple[float, ...] = (-1.21, -0.47, 0.19, 0.73, 1.38),
    coupling_sets: tuple[tuple[float, float, float, float, float], ...] = (
        (1.11, 0.83, 1.27, 0.91, 0.64),
        (0.72, 1.34, 0.88, 1.19, 1.51),
    ),
    atol: float = 1e-11,
) -> tuple[bool, tuple[int, int, int, int] | None]:
    """Try all 16 diagonal Z2 matrices G and test H_new=G H_old G^dagger."""
    for signs in product((-1, +1), repeat=4):
        G = np.diag(signs).astype(np.complex128)
        works = True
        for couplings in coupling_sets:
            for q_ in q_values:
                lhs = h_right(q_, *couplings)
                rhs = G @ h_old(q_, *couplings) @ G.conj().T
                if not np.allclose(lhs, rhs, atol=atol, rtol=atol):
                    works = False
                    break
            if not works:
                break
        if works:
            return True, tuple(int(s) for s in signs)
    return False, None


def old_new_differential_error(
    q_: float,
    couplings: tuple[float, float, float, float, float],
) -> float:
    """Return the largest deviation from the required old-vs-new coefficient change."""
    c_new = extract_independent_coefficients_numeric(h_right(q_, *couplings))
    c_old = extract_independent_coefficients_numeric(h_old(q_, *couplings))
    errors = [
        abs(c_new[name] - c_old[name])
        for name in ("f12", "f13", "f14", "f23", "f24")
    ]
    errors.append(abs(c_new["f34"] + c_old["f34"]))
    return float(max(errors))


# =============================================================================
# 13. Scientific self-test
# =============================================================================


def _matrix_error(A: np.ndarray, B: np.ndarray) -> float:
    return float(np.max(np.abs(np.asarray(A) - np.asarray(B))))


def self_test() -> dict[str, float]:
    """
    Execute the full scientific audit suite and return maximum numerical errors.

    The randomized section uses 20 asymmetric positive coupling sets and 50
    independent (q1,q2) pairs per set: exactly 1000 random cases.
    """
    atol = 1e-11
    rng = np.random.default_rng(2026)

    maxima = {
        "raw_fourier": 0.0,
        "hermiticity": 0.0,
        "unitary_transform": 0.0,
        "finite_chain_bloch": 0.0,
        "majorana_symmetry": 0.0,
        "energy_normalization": 0.0,
        "q2_cancellation": 0.0,
        "final_coefficients": 0.0,
        "periodicity_2pi": 0.0,
        "periodicity_pi_spectrum": 0.0,
    }

    # [1] Gauge validity.
    validate_gauge(RIGHT_GAUGE)
    validate_gauge(OLD_GAUGE)
    assert RIGHT_GAUGE["mu34x"] == -1 and RIGHT_GAUGE["mu43w"] == -1

    # [2] x-block fluxes.
    x_flux = x_block_fluxes(RIGHT_GAUGE)
    assert x_flux == {"W124": +1, "W234": -1, "W123": +1, "W134": -1}

    # [3] w-block fluxes.
    w_flux = w_block_fluxes(RIGHT_GAUGE)
    assert w_flux == {"W213": +1, "W143": -1, "W214": +1, "W243": -1}

    # [4] Symbolic unitarity of U.
    U_sym = bloch_basis_matrix_symbolic(q2)
    assert (U_sym.conjugate().T * U_sym - sp.eye(4)).applyfunc(sp.simplify) == sp.zeros(4)

    # [5] Independent symbolic RAW: formulas vs geometry-derived bonds.
    H_raw_from_bonds_sym = raw_hamiltonian_from_bonds_symbolic(RIGHT_GAUGE)
    raw_symbolic_difference = (H_RAW_SYMBOLIC - H_raw_from_bonds_sym).applyfunc(sp.simplify)
    assert raw_symbolic_difference == sp.zeros(4)

    # [6] Eq.(27) direction U^dagger H_raw U is explicit and exact.
    H_eq27_sym = (U_sym.conjugate().T * H_RAW_SYMBOLIC * U_sym).applyfunc(sp.simplify)
    assert (H_RIGHT_SYMBOLIC - H_eq27_sym).applyfunc(sp.simplify) == sp.zeros(4)

    # [7] Symbolic elimination of q2 from every final matrix element.
    for m in range(4):
        for n in range(4):
            assert sp.simplify(sp.diff(H_RIGHT_SYMBOLIC[m, n], q2)) == 0
            assert q2 not in H_RIGHT_SYMBOLIC[m, n].free_symbols

    # [8] Final H vs supplied physical formulas.
    for name, expected_expr in EXPECTED_RIGHT_COEFFS_SYMBOLIC.items():
        assert sp.simplify(RIGHT_COEFFS_SYMBOLIC[name] - expected_expr) == 0

    # [9] Reduced coefficient identities.
    assert sp.simplify(RIGHT_COEFFS_SYMBOLIC["f23"] - RIGHT_COEFFS_SYMBOLIC["f14"]) == 0
    assert sp.simplify(RIGHT_COEFFS_SYMBOLIC["f24"] - RIGHT_COEFFS_SYMBOLIC["f13"]) == 0
    assert sp.simplify(
        RIGHT_COEFFS_SYMBOLIC["f34"] - sp.conjugate(RIGHT_COEFFS_SYMBOLIC["f12"])
    ) == 0

    # [10]-[13] deterministic symmetry sanity check before random stress test.
    couplings0 = (1.10, 0.80, 1.20, 0.90, 0.70)
    for q_ in np.linspace(-np.pi, np.pi, 33):
        H = h_right(q_, *couplings0)
        maxima["hermiticity"] = max(maxima["hermiticity"], _matrix_error(H, H.conj().T))
        assert np.allclose(H, H.conj().T, atol=atol, rtol=atol)

        H_minus = h_right(-q_, *couplings0)
        maj_err = _matrix_error(H_minus, -H.conj())
        maxima["majorana_symmetry"] = max(maxima["majorana_symmetry"], maj_err)
        assert maj_err < atol
        assert np.allclose(
            np.linalg.eigvalsh(H_minus),
            -np.linalg.eigvalsh(H)[::-1],
            atol=atol,
            rtol=atol,
        )

        H_2pi = h_right(q_ + 2.0 * np.pi, *couplings0)
        err_2pi = _matrix_error(H_2pi, H)
        maxima["periodicity_2pi"] = max(maxima["periodicity_2pi"], err_2pi)
        assert err_2pi < atol

        H_pi = h_right(q_ + np.pi, *couplings0)
        H_pi_gauge = G_PI_NUMERIC @ H @ G_PI_NUMERIC.conj().T
        assert _matrix_error(H_pi, H_pi_gauge) < atol
        spec_pi_err = float(
            np.max(np.abs(np.linalg.eigvalsh(H_pi) - np.linalg.eigvalsh(H)))
        )
        maxima["periodicity_pi_spectrum"] = max(
            maxima["periodicity_pi_spectrum"], spec_pi_err
        )
        assert spec_pi_err < atol

    # [14] Truly symbolic characteristic-polynomial identity.
    charpoly_difference = characteristic_polynomial_real_component_test()
    assert charpoly_difference == 0, charpoly_difference

    # [15] Finite periodic real-space Majorana chain.
    finite_err = finite_chain_bloch_audit(ncell=17)
    maxima["finite_chain_bloch"] = max(maxima["finite_chain_bloch"], finite_err)
    assert finite_err < atol

    # [16] Differential old vs new: only f34 changes sign.
    for q_ in np.linspace(-1.43, 1.43, 21):
        diff_err = old_new_differential_error(q_, couplings0)
        assert diff_err < atol

    # Different plaquette fluxes already establish gauge inequivalence.
    assert x_block_fluxes(OLD_GAUGE) != x_block_fluxes(RIGHT_GAUGE)
    assert w_block_fluxes(OLD_GAUGE) != w_block_fluxes(RIGHT_GAUGE)

    # [17] Exhaustive diagonal Z2 gauge search: no solution among 16 matrices.
    equivalent, witness = diagonal_z2_gauge_equivalence()
    assert not equivalent and witness is None

    # [18] Energy normalization: physical BZ and extended BZ agree.
    q_phys = np.linspace(-np.pi / 2.0, np.pi / 2.0, 2001)
    q_ext = np.linspace(-np.pi, np.pi, 4001)
    E_phys = ground_state_energy_per_cell_physical_bz(
        *couplings0, q_values=q_phys
    )
    E_ext = ground_state_energy_per_cell_extended_bz(
        *couplings0, q_values=q_ext
    )
    energy_err = abs(E_phys - E_ext)
    maxima["energy_normalization"] = max(maxima["energy_normalization"], energy_err)
    assert energy_err < 5e-10
    assert np.isclose(
        ground_state_energy_per_site(*couplings0, q_values=q_phys),
        E_phys / 4.0,
        atol=1e-13,
        rtol=1e-13,
    )

    # [19] Isotropic benchmark at the requested momenta.
    isotropic_expected = {
        0.0: np.array([-3.0, -1.0, 1.0, 3.0]),
        np.pi / 4.0: np.array([-1.8284271247461903, -1.0, -1.0, 3.8284271247461903]),
        np.pi / 2.0: np.array([-3.0, -1.0, 1.0, 3.0]),
        np.pi: np.array([-3.0, -1.0, 1.0, 3.0]),
        -np.pi / 4.0: np.array([-3.8284271247461903, 1.0, 1.0, 1.8284271247461903]),
    }
    for q_, expected_eigs in isotropic_expected.items():
        eigs = energies_right(q_, 1.0, 1.0, 1.0, 1.0, 1.0)
        assert np.allclose(eigs, expected_eigs, atol=2e-12, rtol=2e-12)

    # [20] 1000 deterministic random tests: 20 couplings x 50 (q1,q2).
    random_case_count = 0
    for _ in range(20):
        couplings = tuple(float(x) for x in rng.uniform(0.17, 2.41, size=5))
        # Generic continuous RNG makes exact accidental equality negligibly likely;
        # this assertion protects the intended asymmetric stress test explicitly.
        assert len({round(x, 12) for x in couplings}) == 5

        for _ in range(50):
            q1_ = float(rng.uniform(-np.pi, np.pi))
            q2_ = float(rng.uniform(-np.pi, np.pi))
            random_case_count += 1

            # (1) RAW formula vs independent bond-list reconstruction.
            H_raw_formula = h_raw_numeric(q1_, q2_, *couplings, RIGHT_GAUGE)
            H_raw_bonds = raw_hamiltonian_from_bonds_numeric(
                q1_, q2_, *couplings, RIGHT_GAUGE
            )
            raw_err = _matrix_error(H_raw_formula, H_raw_bonds)
            maxima["raw_fourier"] = max(maxima["raw_fourier"], raw_err)
            assert raw_err < atol

            # (2) Hermiticity.
            herm_raw = _matrix_error(H_raw_formula, H_raw_formula.conj().T)
            maxima["hermiticity"] = max(maxima["hermiticity"], herm_raw)
            assert herm_raw < atol

            # (3) U U^dagger = I.
            U = bloch_basis_matrix_numeric(q2_)
            unitary_err = _matrix_error(U @ U.conj().T, np.eye(4))
            maxima["unitary_transform"] = max(
                maxima["unitary_transform"], unitary_err
            )
            assert unitary_err < atol

            # (4) H_Bloch = U^dagger H_raw U.
            H_bloch = h_right(q1_, *couplings, q2_=q2_)
            herm_bloch = _matrix_error(H_bloch, H_bloch.conj().T)
            maxima["hermiticity"] = max(maxima["hermiticity"], herm_bloch)
            assert herm_bloch < atol

            explicit_transform = U.conj().T @ H_raw_formula @ U
            transform_err = _matrix_error(H_bloch, explicit_transform)
            maxima["unitary_transform"] = max(
                maxima["unitary_transform"], transform_err
            )
            assert transform_err < atol

            # (5) Unitariy preserves the spectrum.
            assert np.allclose(
                np.linalg.eigvalsh(H_raw_formula),
                np.linalg.eigvalsh(H_bloch),
                atol=atol,
                rtol=atol,
            )

            # (6) q2 cancels in the final Bloch matrix.
            H_q2_zero = h_right(q1_, *couplings, q2_=0.0)
            q2_err = _matrix_error(H_bloch, H_q2_zero)
            maxima["q2_cancellation"] = max(maxima["q2_cancellation"], q2_err)
            assert q2_err < atol

            # (7) Final coefficient formulas.
            coeffs = extract_independent_coefficients_numeric(H_bloch)
            expected = expected_bloch_coefficients_numeric(q1_, *couplings)
            coeff_err = max(abs(coeffs[name] - expected[name]) for name in expected)
            maxima["final_coefficients"] = max(
                maxima["final_coefficients"], float(coeff_err)
            )
            assert coeff_err < atol

            # (8)-(10) Reduced identities.
            assert abs(coeffs["f23"] - coeffs["f14"]) < atol
            assert abs(coeffs["f24"] - coeffs["f13"]) < atol
            assert abs(coeffs["f34"] - np.conj(coeffs["f12"])) < atol

            # (11) Majorana symmetry.
            H_minus = h_right(-q1_, *couplings)
            maj_err = _matrix_error(H_minus, -H_q2_zero.conj())
            maxima["majorana_symmetry"] = max(
                maxima["majorana_symmetry"], maj_err
            )
            assert maj_err < atol
            assert np.allclose(
                np.linalg.eigvalsh(H_minus),
                -np.linalg.eigvalsh(H_q2_zero)[::-1],
                atol=atol,
                rtol=atol,
            )

            # (12) 2pi matrix periodicity.
            H_2pi = h_right(q1_ + 2.0 * np.pi, *couplings)
            err_2pi = _matrix_error(H_2pi, H_q2_zero)
            maxima["periodicity_2pi"] = max(maxima["periodicity_2pi"], err_2pi)
            assert err_2pi < atol

            # (13) pi gauge-periodicity and spectral periodicity.
            H_pi = h_right(q1_ + np.pi, *couplings)
            H_pi_expected = G_PI_NUMERIC @ H_q2_zero @ G_PI_NUMERIC.conj().T
            assert _matrix_error(H_pi, H_pi_expected) < atol
            spec_pi_err = float(
                np.max(
                    np.abs(
                        np.linalg.eigvalsh(H_pi)
                        - np.linalg.eigvalsh(H_q2_zero)
                    )
                )
            )
            maxima["periodicity_pi_spectrum"] = max(
                maxima["periodicity_pi_spectrum"], spec_pi_err
            )
            assert spec_pi_err < atol

    assert random_case_count >= 1000

    print("Scientific self-test PASSED")
    print(f"random cases tested: {random_case_count}")
    print(f"x-block fluxes: {x_flux}")
    print(f"w-block fluxes: {w_flux}")
    print(f"maximum raw-Fourier error: {maxima['raw_fourier']:.3e}")
    print(f"maximum Hermiticity error: {maxima['hermiticity']:.3e}")
    print(f"maximum unitary-transform error: {maxima['unitary_transform']:.3e}")
    print(f"maximum finite-chain/Bloch error: {maxima['finite_chain_bloch']:.3e}")
    print(f"maximum Majorana-symmetry error: {maxima['majorana_symmetry']:.3e}")
    print(
        "maximum energy-normalization discrepancy: "
        f"{maxima['energy_normalization']:.3e}"
    )
    return maxima


# =============================================================================
# 14. Command-line demonstration
# =============================================================================


if __name__ == "__main__":
    print("Corrected physical gauge:")
    print(RIGHT_GAUGE)

    print("\nx-block fluxes:")
    print(x_block_fluxes())
    print("w-block fluxes:")
    print(w_block_fluxes())

    print("\nEq.(27) basis matrix U(q2):")
    sp.pprint(bloch_basis_matrix_symbolic(q2))
    print("Transformation convention: H_Bloch = U^dagger H_raw U")

    print("\nRAW Fourier coefficients (correct exp(+i q2) cross-chain phase):")
    for name, expr in RAW_COEFFS_SYMBOLIC.items():
        print(f"{name}_raw(q1,q2) =")
        sp.pprint(expr)

    print("\nFinal q2-independent Bloch coefficients:")
    for name, expr in RIGHT_COEFFS_SYMBOLIC.items():
        print(f"{name}(q1) =")
        sp.pprint(expr)

    print("\nRunning scientific consistency suite...")
    self_test()

    print("\nIsotropic spectrum at q=pi/4:")
    print(energies_right(np.pi / 4.0, 1.0, 1.0, 1.0, 1.0, 1.0))
