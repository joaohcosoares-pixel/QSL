"""
Quadrupolar spin liquid in 1D — right-oriented upper-chain gauge,
with the Fourier -> Eq.(27) Bloch-basis transformation implemented explicitly.

Physical gauge used here
------------------------
Arrow i -> j  <=>  mu_ij = +1, with mu_ji = -mu_ij.
The upper horizontal x and w bonds are reversed relative to the original
zero-flux representative so that both point to the right:

    mu_12,x = +1      mu_34,x = -1
    mu_23,y = +1      mu_41,y = -1
    mu_42,z = +1      mu_31,z = +1
    mu_43,w = -1      mu_12,w = -1
    mu_13,t = -1      mu_24,t = -1

Derivation pipeline implemented in this file
--------------------------------------------
    bond variables mu_ij
        -> Fourier-space raw matrix H_raw(k)
        -> Eq.(27) basis matrix U(k)
        -> H_Bloch(k) = U(k) H_raw(k) U(k)^dagger
        -> eigenvalues / bands

The Eq.(27) matrix is

    U(k) = diag(1, 1, exp(-i k.a2), exp(-i k.a2)).

After this transformation all k.a2 phases cancel from the final Bloch
Hamiltonian, which depends only on q = k.a1 (identified with ka in the 1D
band plots).  For the new gauge one obtains

    f23(q) = f14(q)
    f24(q) = f13(q)
    f34(q) = +conjugate(f12(q)).

Therefore the old compact pointwise +/- spectrum derived for
f34 = -conjugate(f12) must not be reused for this gauge.
"""

from __future__ import annotations

from typing import Mapping

import numpy as np
import sympy as sp


# =============================================================================
# 1. Symbolic variables
# =============================================================================

f12, f13, f14, f23, f24, f34 = sp.symbols(
    "f12 f13 f14 f23 f24 f34", complex=True
)
lambda_ = sp.symbols("lambda", real=True)

# q1 = k.a1 is the dispersive 1D crystal momentum; q2 = k.a2 is retained
# only to reproduce explicitly the pre-Eq.(27) Fourier expressions.
q1, q2 = sp.symbols("q1 q2", real=True)
Kx, Ky, Kz, Kt, Kw = sp.symbols("Kx Ky Kz Kt Kw", real=True)


# =============================================================================
# 2. Gauge configuration: upper horizontal bonds point to the right
# =============================================================================

RIGHT_GAUGE: dict[str, int] = {
    "mu12x": +1,
    "mu34x": -1,   # upper x bond reversed
    "mu23y": +1,
    "mu41y": -1,
    "mu42z": +1,
    "mu31z": +1,
    "mu43w": -1,   # upper w bond reversed
    "mu12w": -1,
    "mu13t": -1,
    "mu24t": -1,
}


def validate_gauge(gauge: Mapping[str, int]) -> None:
    """Require every stored Z2 bond variable to be exactly +1 or -1."""
    required = set(RIGHT_GAUGE)
    missing = required.difference(gauge)
    if missing:
        raise ValueError(f"Gauge is missing bond variables: {sorted(missing)}")

    invalid = {key: value for key, value in gauge.items() if value not in (-1, +1)}
    if invalid:
        raise ValueError(f"Gauge variables must be +/-1. Invalid entries: {invalid}")


# =============================================================================
# 3. Generic 4x4 Hermitian Hamiltonian
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
            [0,                    f12_,                    f13_,                    f14_],
            [sp.conjugate(f12_),   0,                       f23_,                    f24_],
            [sp.conjugate(f13_),   sp.conjugate(f23_),      0,                       f34_],
            [sp.conjugate(f14_),   sp.conjugate(f24_),      sp.conjugate(f34_),      0],
        ]
    )


H_GENERIC = ham_qsl_symbolic(f12, f13, f14, f23, f24, f34)


def generic_characteristic_polynomial() -> sp.Expr:
    """det(lambda I - H) for the fully generic Hermitian matrix."""
    return sp.factor((lambda_ * sp.eye(4) - H_GENERIC).det())


# =============================================================================
# 4. Fluxes of the x block for the new gauge
# =============================================================================


def x_block_fluxes(gauge: Mapping[str, int] = RIGHT_GAUGE) -> dict[str, int]:
    """
    Plaquette circulations for the four triangles of the x block.

    Paths:
        1 -> 2 -> 4 -> 1
        2 -> 3 -> 4 -> 2
        1 -> 2 -> 3 -> 1
        1 -> 3 -> 4 -> 1

    Whenever a path traverses a stored link backwards, mu_ji = -mu_ij is used.
    """
    validate_gauge(gauge)

    mu24z = -gauge["mu42z"]
    mu31t = -gauge["mu13t"]

    return {
        "W124": gauge["mu12x"] * mu24z * gauge["mu41y"],
        "W234": gauge["mu23y"] * gauge["mu34x"] * gauge["mu42z"],
        "W123": gauge["mu12x"] * gauge["mu23y"] * mu31t,
        "W134": gauge["mu13t"] * gauge["mu34x"] * gauge["mu41y"],
    }


# =============================================================================
# 5. Fourier-space RAW coefficients, BEFORE Eq.(27)
# =============================================================================


def raw_fourier_coefficients_symbolic(
    gauge: Mapping[str, int] = RIGHT_GAUGE,
) -> dict[str, sp.Expr]:
    """
    Six independent coefficients of H_raw(k), i.e. after Fourier transform but
    BEFORE the Bloch-basis change of Eq.(27).

    q1 = k.a1
    q2 = k.a2

    Links connecting lower sublattices (1,2) to upper sublattices (3,4)
    retain a common exp(-i q2) factor at this stage.  The 1-2 and 3-4
    horizontal links do not.
    """
    validate_gauge(gauge)

    return {
        "f12": sp.I * (
            gauge["mu12x"] * Kx * sp.exp(sp.I * q1)
            + gauge["mu12w"] * Kw * sp.exp(-sp.I * q1)
        ),
        "f13": sp.exp(-sp.I * q2) * (
            sp.I * gauge["mu13t"] * Kt * sp.exp(sp.I * q1)
            - sp.I * gauge["mu31z"] * Kz * sp.exp(-sp.I * q1)
        ),
        "f14": sp.exp(-sp.I * q2) * (-sp.I * gauge["mu41y"] * Ky),
        "f23": sp.exp(-sp.I * q2) * (sp.I * gauge["mu23y"] * Ky),
        "f24": sp.exp(-sp.I * q2) * (
            sp.I * gauge["mu24t"] * Kt * sp.exp(sp.I * q1)
            - sp.I * gauge["mu42z"] * Kz * sp.exp(-sp.I * q1)
        ),
        "f34": sp.I * (
            gauge["mu34x"] * Kx * sp.exp(-sp.I * q1)
            - gauge["mu43w"] * Kw * sp.exp(sp.I * q1)
        ),
    }


RAW_COEFFS_SYMBOLIC = {
    name: sp.simplify(expr)
    for name, expr in raw_fourier_coefficients_symbolic(RIGHT_GAUGE).items()
}


def raw_fourier_hamiltonian_symbolic() -> sp.Matrix:
    """H_raw(k) immediately after Fourier transform and before Eq.(27)."""
    c = RAW_COEFFS_SYMBOLIC
    return ham_qsl_symbolic(c["f12"], c["f13"], c["f14"], c["f23"], c["f24"], c["f34"])


H_RAW_SYMBOLIC = raw_fourier_hamiltonian_symbolic()


# =============================================================================
# 6. Eq.(27): explicit Bloch-basis transformation
# =============================================================================


def bloch_basis_matrix_symbolic(q2_: sp.Expr = q2) -> sp.Matrix:
    """
    Eq.(27): U(k) = diag(1, 1, exp(-i k.a2), exp(-i k.a2)).
    """
    return sp.diag(1, 1, sp.exp(-sp.I * q2_), sp.exp(-sp.I * q2_))


def apply_bloch_basis_transform_symbolic(
    H_raw: sp.Matrix,
    q2_: sp.Expr = q2,
) -> sp.Matrix:
    """
    Apply Eq.(27) with the convention Theta' = U Theta:

        H_Bloch = U H_raw U^dagger.
    """
    U = bloch_basis_matrix_symbolic(q2_)
    return sp.simplify(U * H_raw * U.conjugate().T)


H_RIGHT_SYMBOLIC = apply_bloch_basis_transform_symbolic(H_RAW_SYMBOLIC, q2)


def extract_independent_coefficients_symbolic(H: sp.Matrix) -> dict[str, sp.Expr]:
    """Extract the six independent upper-triangular entries from a 4x4 matrix."""
    return {
        "f12": sp.simplify(H[0, 1]),
        "f13": sp.simplify(H[0, 2]),
        "f14": sp.simplify(H[0, 3]),
        "f23": sp.simplify(H[1, 2]),
        "f24": sp.simplify(H[1, 3]),
        "f34": sp.simplify(H[2, 3]),
    }


RIGHT_COEFFS_SYMBOLIC = extract_independent_coefficients_symbolic(H_RIGHT_SYMBOLIC)


# =============================================================================
# 7. Expected post-Eq.(27) formulas and reduced symbolic structure
# =============================================================================


def expected_bloch_coefficients_symbolic(
    gauge: Mapping[str, int] = RIGHT_GAUGE,
) -> dict[str, sp.Expr]:
    """Closed expressions expected after Eq.(27); used only as an audit target."""
    validate_gauge(gauge)
    return {
        "f12": sp.I * (
            gauge["mu12x"] * Kx * sp.exp(sp.I * q1)
            + gauge["mu12w"] * Kw * sp.exp(-sp.I * q1)
        ),
        "f13": sp.I * gauge["mu13t"] * Kt * sp.exp(sp.I * q1)
        - sp.I * gauge["mu31z"] * Kz * sp.exp(-sp.I * q1),
        "f14": -sp.I * gauge["mu41y"] * Ky,
        "f23": sp.I * gauge["mu23y"] * Ky,
        "f24": sp.I * gauge["mu24t"] * Kt * sp.exp(sp.I * q1)
        - sp.I * gauge["mu42z"] * Kz * sp.exp(-sp.I * q1),
        "f34": sp.I * (
            gauge["mu34x"] * Kx * sp.exp(-sp.I * q1)
            - gauge["mu43w"] * Kw * sp.exp(sp.I * q1)
        ),
    }


EXPECTED_RIGHT_COEFFS_SYMBOLIC = {
    name: sp.simplify(expr)
    for name, expr in expected_bloch_coefficients_symbolic(RIGHT_GAUGE).items()
}


def right_reduced_hamiltonian() -> sp.Matrix:
    """
    Abstract reduced matrix for the new gauge:
        f23 = f14, f24 = f13, f34 = +conjugate(f12).
    """
    return ham_qsl_symbolic(f12, f13, f14, f14, f13, sp.conjugate(f12))


H_RIGHT_REDUCED = right_reduced_hamiltonian()


def right_reduced_characteristic_polynomial() -> sp.Expr:
    """Exact characteristic polynomial of the new reduced symbolic matrix."""
    return sp.factor((lambda_ * sp.eye(4) - H_RIGHT_REDUCED).det())


def right_reduced_characteristic_polynomial_compact() -> sp.Expr:
    """
    Compact det(lambda I - H) for a=f12, b=f13, c=f14,
    with f34=conjugate(a).
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
        - 4 * (
            sp.Abs(a) ** 2 * sp.Abs(c) ** 2
            + re_a**2 * sp.Abs(b) ** 2
            + re_bc**2
        )
    )


# =============================================================================
# 8. Numerical matrices: raw Fourier -> Eq.(27) -> Bloch
# =============================================================================


def ham_qsl_numeric(
    f12_: complex,
    f13_: complex,
    f14_: complex,
    f23_: complex,
    f24_: complex,
    f34_: complex,
) -> np.ndarray:
    """Numerical Hermitian 4x4 matrix in the same convention as the symbolic one."""
    return np.array(
        [
            [0.0,             f12_,             f13_,             f14_],
            [np.conj(f12_),   0.0,              f23_,             f24_],
            [np.conj(f13_),   np.conj(f23_),    0.0,              f34_],
            [np.conj(f14_),   np.conj(f24_),    np.conj(f34_),    0.0],
        ],
        dtype=np.complex128,
    )


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
    """Numerical pre-Eq.(27) Fourier coefficients."""
    validate_gauge(gauge)

    phase2 = np.exp(-1j * q2_)
    return {
        "f12": 1j * (
            gauge["mu12x"] * Kx_ * np.exp(1j * q1_)
            + gauge["mu12w"] * Kw_ * np.exp(-1j * q1_)
        ),
        "f13": phase2 * (
            1j * gauge["mu13t"] * Kt_ * np.exp(1j * q1_)
            - 1j * gauge["mu31z"] * Kz_ * np.exp(-1j * q1_)
        ),
        "f14": phase2 * (-1j * gauge["mu41y"] * Ky_),
        "f23": phase2 * (1j * gauge["mu23y"] * Ky_),
        "f24": phase2 * (
            1j * gauge["mu24t"] * Kt_ * np.exp(1j * q1_)
            - 1j * gauge["mu42z"] * Kz_ * np.exp(-1j * q1_)
        ),
        "f34": 1j * (
            gauge["mu34x"] * Kx_ * np.exp(-1j * q1_)
            - gauge["mu43w"] * Kw_ * np.exp(1j * q1_)
        ),
    }


def h_right_raw(
    q1_: float,
    q2_: float,
    Kx_: float,
    Ky_: float,
    Kz_: float,
    Kt_: float,
    Kw_: float,
) -> np.ndarray:
    """H_raw(k): after Fourier transform and before Eq.(27)."""
    c = raw_fourier_coefficients_numeric(
        q1_, q2_, Kx_, Ky_, Kz_, Kt_, Kw_, RIGHT_GAUGE
    )
    return ham_qsl_numeric(c["f12"], c["f13"], c["f14"], c["f23"], c["f24"], c["f34"])


def bloch_basis_matrix_numeric(q2_: float) -> np.ndarray:
    """Numerical Eq.(27): diag(1,1,e^{-iq2},e^{-iq2})."""
    phase = np.exp(-1j * q2_)
    return np.diag([1.0, 1.0, phase, phase]).astype(np.complex128)


def apply_bloch_basis_transform_numeric(H_raw: np.ndarray, q2_: float) -> np.ndarray:
    """Apply H_Bloch = U H_raw U^dagger."""
    U = bloch_basis_matrix_numeric(q2_)
    return U @ H_raw @ U.conj().T


def h_right(
    q1_: float,
    Kx_: float,
    Ky_: float,
    Kz_: float,
    Kt_: float,
    Kw_: float,
    q2_: float = 0.0,
) -> np.ndarray:
    """
    Final Bloch Hamiltonian for the right-oriented upper-chain gauge.

    q1 = k.a1 (the 1D dispersive variable, later q=ka)
    q2 = k.a2, kept only to verify explicitly that Eq.(27) removes it.

    The default q2=0 keeps the same convenient public API used previously.
    """
    H_raw = h_right_raw(q1_, q2_, Kx_, Ky_, Kz_, Kt_, Kw_)
    return apply_bloch_basis_transform_numeric(H_raw, q2_)


def extract_independent_coefficients_numeric(H: np.ndarray) -> dict[str, complex]:
    """Extract the six independent upper-triangular entries of a numerical H."""
    return {
        "f12": H[0, 1],
        "f13": H[0, 2],
        "f14": H[0, 3],
        "f23": H[1, 2],
        "f24": H[1, 3],
        "f34": H[2, 3],
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
    """Final post-Eq.(27) Bloch coefficients, extracted from H_Bloch."""
    return extract_independent_coefficients_numeric(
        h_right(q1_, Kx_, Ky_, Kz_, Kt_, Kw_, q2_)
    )


def expected_bloch_coefficients_numeric(
    q1_: float,
    Kx_: float,
    Ky_: float,
    Kz_: float,
    Kt_: float,
    Kw_: float,
    gauge: Mapping[str, int] = RIGHT_GAUGE,
) -> dict[str, complex]:
    """Closed post-Eq.(27) formulas used as an independent numerical audit."""
    validate_gauge(gauge)
    return {
        "f12": 1j * (
            gauge["mu12x"] * Kx_ * np.exp(1j * q1_)
            + gauge["mu12w"] * Kw_ * np.exp(-1j * q1_)
        ),
        "f13": 1j * gauge["mu13t"] * Kt_ * np.exp(1j * q1_)
        - 1j * gauge["mu31z"] * Kz_ * np.exp(-1j * q1_),
        "f14": -1j * gauge["mu41y"] * Ky_,
        "f23": 1j * gauge["mu23y"] * Ky_,
        "f24": 1j * gauge["mu24t"] * Kt_ * np.exp(1j * q1_)
        - 1j * gauge["mu42z"] * Kz_ * np.exp(-1j * q1_),
        "f34": 1j * (
            gauge["mu34x"] * Kx_ * np.exp(-1j * q1_)
            - gauge["mu43w"] * Kw_ * np.exp(1j * q1_)
        ),
    }


# =============================================================================
# 9. Diagonalization, bands and occupied-band integral
# =============================================================================


def eigensystem_numeric(H: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Hermitian eigenvalues/eigenvectors, eigenvalues in ascending order."""
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
    """Return an array with shape (len(q_values), 4)."""
    q_values = np.asarray(q_values, dtype=float)
    return np.vstack(
        [energies_right(x, Kx_, Ky_, Kz_, Kt_, Kw_) for x in q_values]
    )


def occupied_band_integral_right(
    q_values: np.ndarray,
    Kx_: float,
    Ky_: float,
    Kz_: float,
    Kt_: float,
    Kw_: float,
) -> float:
    """
    Compute
        integral_{-pi}^{pi} dq/(2*pi) sum_{E_n(q)<0} E_n(q)
    by trapezoidal integration.

    This is still only the occupied-band integral in the full-BZ convention.
    The overall Majorana BZ/BZ/2 normalization must be fixed separately before
    interpreting it as the final physical ground-state energy per site/cell.
    """
    q_values = np.asarray(q_values, dtype=float)
    bands = band_structure_right(q_values, Kx_, Ky_, Kz_, Kt_, Kw_)
    occupied_sum = np.where(bands < 0.0, bands, 0.0).sum(axis=1)
    return float(np.trapezoid(occupied_sum, q_values) / (2.0 * np.pi))


# =============================================================================
# 10. Consistency tests
# =============================================================================


def self_test() -> None:
    """Checks the gauge, Eq.(27), final coefficients and Majorana structure."""
    validate_gauge(RIGHT_GAUGE)

    # Gauge-flux pattern for the x block.
    fluxes = x_block_fluxes(RIGHT_GAUGE)
    expected_fluxes = {"W124": +1, "W234": -1, "W123": +1, "W134": -1}
    assert fluxes == expected_fluxes, (fluxes, expected_fluxes)

    couplings = (1.10, 0.80, 1.20, 0.90, 0.70)  # Kx, Ky, Kz, Kt, Kw

    # Test several independent q2 values: the raw matrix depends on q2, while
    # the post-Eq.(27) Bloch matrix must not.
    q2_values = (-1.17, -0.41, 0.0, 0.63, 1.31)

    for q1_ in np.linspace(-np.pi, np.pi, 31):
        H_reference = h_right(q1_, *couplings, q2_=0.0)
        expected = expected_bloch_coefficients_numeric(q1_, *couplings)

        for q2_ in q2_values:
            U = bloch_basis_matrix_numeric(q2_)
            H_raw = h_right_raw(q1_, q2_, *couplings)
            H_bloch = h_right(q1_, *couplings, q2_=q2_)
            coeffs = extract_independent_coefficients_numeric(H_bloch)

            # Eq.(27) matrix is unitary.
            assert np.allclose(U @ U.conj().T, np.eye(4), atol=1e-12, rtol=1e-12)

            # Explicit implementation of H_Bloch = U H_raw U^dagger.
            assert np.allclose(
                H_bloch,
                U @ H_raw @ U.conj().T,
                atol=1e-12,
                rtol=1e-12,
            )

            # A unitary basis change cannot alter the spectrum.
            assert np.allclose(
                np.linalg.eigvalsh(H_raw),
                np.linalg.eigvalsh(H_bloch),
                atol=1e-12,
                rtol=1e-12,
            )

            # Eq.(27) must remove all q2 dependence from the final H.
            assert np.allclose(H_bloch, H_reference, atol=1e-12, rtol=1e-12)

            # Final coefficients must agree with the closed Eqs.(21)-(26)-type
            # formulas for the NEW gauge.
            for name in expected:
                assert np.allclose(
                    coeffs[name], expected[name], atol=1e-12, rtol=1e-12
                ), (name, coeffs[name], expected[name])

            # New reduced relations.
            assert np.allclose(coeffs["f23"], coeffs["f14"], atol=1e-12, rtol=1e-12)
            assert np.allclose(coeffs["f24"], coeffs["f13"], atol=1e-12, rtol=1e-12)
            assert np.allclose(
                coeffs["f34"], np.conj(coeffs["f12"]), atol=1e-12, rtol=1e-12
            )

            # Hermiticity and trace.
            assert np.allclose(H_bloch, H_bloch.conj().T, atol=1e-12, rtol=1e-12)
            assert np.isclose(np.trace(H_bloch), 0.0, atol=1e-12)

        # 2*pi periodicity in q1.
        assert np.allclose(
            h_right(q1_ + 2.0 * np.pi, *couplings),
            H_reference,
            atol=1e-12,
            rtol=1e-12,
        )

        # Majorana relation for the final Bloch Hamiltonian.
        H_minus = h_right(-q1_, *couplings)
        assert np.allclose(H_minus, -H_reference.conj(), atol=1e-12, rtol=1e-12)

        e_plus = np.linalg.eigvalsh(H_reference)
        e_minus = np.linalg.eigvalsh(H_minus)
        assert np.allclose(e_minus, -e_plus[::-1], atol=1e-12, rtol=1e-12)

    # Symbolic Eq.(27) cancellation and coefficient checks.
    for name in EXPECTED_RIGHT_COEFFS_SYMBOLIC:
        assert sp.simplify(
            RIGHT_COEFFS_SYMBOLIC[name] - EXPECTED_RIGHT_COEFFS_SYMBOLIC[name]
        ) == 0

    assert sp.simplify(RIGHT_COEFFS_SYMBOLIC["f23"] - RIGHT_COEFFS_SYMBOLIC["f14"]) == 0
    assert sp.simplify(RIGHT_COEFFS_SYMBOLIC["f24"] - RIGHT_COEFFS_SYMBOLIC["f13"]) == 0
    assert sp.simplify(
        RIGHT_COEFFS_SYMBOLIC["f34"] - sp.conjugate(RIGHT_COEFFS_SYMBOLIC["f12"])
    ) == 0

    # Isotropic band sanity check.
    q_grid = np.linspace(-np.pi, np.pi, 101)
    bands = band_structure_right(q_grid, 1.0, 1.0, 1.0, 1.0, 1.0)
    assert bands.shape == (101, 4)
    assert np.all(np.isfinite(bands))
    assert np.all(np.isreal(bands))

    print("Self-test passed for the NEW gauge with Eq.(27) explicit.")
    print("Gauge:", RIGHT_GAUGE)
    print("x-block fluxes:", fluxes)


# =============================================================================
# 11. Command-line demonstration
# =============================================================================


if __name__ == "__main__":
    print("NEW gauge: upper horizontal x and w bonds point to the right")
    print(RIGHT_GAUGE)

    print("\nFluxes in the x block:")
    print(x_block_fluxes())

    print("\nEq.(27) Bloch-basis matrix U(k):")
    sp.pprint(bloch_basis_matrix_symbolic(q2))

    print("\nRaw Fourier coefficients (before Eq.(27)):")
    for name, expr in RAW_COEFFS_SYMBOLIC.items():
        print(f"{name}_raw(q1,q2) =")
        sp.pprint(expr)

    print("\nFinal Bloch coefficients (after Eq.(27)):")
    for name, expr in RIGHT_COEFFS_SYMBOLIC.items():
        print(f"{name}(q1) =")
        sp.pprint(expr)

    print("\nKey symbolic checks:")
    print(
        "f23 - f14 =",
        sp.simplify(RIGHT_COEFFS_SYMBOLIC["f23"] - RIGHT_COEFFS_SYMBOLIC["f14"]),
    )
    print(
        "f24 - f13 =",
        sp.simplify(RIGHT_COEFFS_SYMBOLIC["f24"] - RIGHT_COEFFS_SYMBOLIC["f13"]),
    )
    print(
        "f34 - conjugate(f12) =",
        sp.simplify(
            RIGHT_COEFFS_SYMBOLIC["f34"]
            - sp.conjugate(RIGHT_COEFFS_SYMBOLIC["f12"])
        ),
    )

    print("\nReduced characteristic polynomial for the NEW gauge:")
    sp.pprint(right_reduced_characteristic_polynomial_compact())

    print("\nRunning consistency tests...")
    self_test()

    print("\nSample isotropic spectrum at q1 = pi/4:")
    print(energies_right(np.pi / 4.0, 1.0, 1.0, 1.0, 1.0, 1.0))
