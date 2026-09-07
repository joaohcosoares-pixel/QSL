"""
Python/SymPy translation of the Mathematica notebook

The original notebook:
1) defines a generic 4x4 Hermitian Bloch Hamiltonian with zero diagonal;
2) asks Mathematica for symbolic eigenvalues;
3) imposes the zero-flux relations
      f23 = f14,  f24 = f13,  f34 = -conjugate(f12)
   and obtains the compact particle-hole-symmetric spectrum.

This Python version keeps the symbolic analysis in SymPy and adds a robust
NumPy path for numerical diagonalization at each momentum later in the project.
"""

from __future__ import annotations

import numpy as np
import sympy as sp


# -----------------------------------------------------------------------------
# 1. Symbolic variables
# -----------------------------------------------------------------------------

# The Mathematica notebook assumes these quantities are complex.
f12, f13, f14, f23, f24, f34 = sp.symbols(
    "f12 f13 f14 f23 f24 f34", complex=True
)
lambda_ = sp.symbols("lambda", real=True)


# -----------------------------------------------------------------------------
# 2. Generic 4x4 Hermitian Hamiltonian
# -----------------------------------------------------------------------------

def ham_qsl_symbolic(
    f12_: sp.Expr,
    f13_: sp.Expr,
    f14_: sp.Expr,
    f23_: sp.Expr,
    f24_: sp.Expr,
    f34_: sp.Expr,
) -> sp.Matrix:
    """Return the generic Hermitian 4x4 QSL matrix used in the notebook."""
    return sp.Matrix(
        [
            [0,                 f12_,                 f13_,                 f14_],
            [sp.conjugate(f12_), 0,                  f23_,                 f24_],
            [sp.conjugate(f13_), sp.conjugate(f23_), 0,                    f34_],
            [sp.conjugate(f14_), sp.conjugate(f24_), sp.conjugate(f34_),   0],
        ]
    )


H_generic = ham_qsl_symbolic(f12, f13, f14, f23, f24, f34)


def generic_characteristic_polynomial():
    return sp.factor((lambda_ * sp.eye(4) - H_generic).det())


# -----------------------------------------------------------------------------
# 3. Zero-flux specialization used in the Mathematica notebook
# -----------------------------------------------------------------------------

def zero_flux_hamiltonian() -> sp.Matrix:
    """
    Impose exactly the relations used in the notebook:

        f23 = f14
        f24 = f13
        f34 = -conjugate(f12)
    """
    return ham_qsl_symbolic(
        f12,
        f13,
        f14,
        f14,
        f13,
        -sp.conjugate(f12),
    )


H_zero_flux = zero_flux_hamiltonian()


def zero_flux_closed_form() -> tuple[sp.Expr, sp.Expr, tuple[sp.Expr, ...]]:
    """
    Compact closed-form spectrum corresponding to Mathematica Out[1415].

    Eigenvalues:
        -sqrt(A + sqrt(B)),
        -sqrt(A - sqrt(B)),
        +sqrt(A - sqrt(B)),
        +sqrt(A + sqrt(B)).
    """
    A = sp.Abs(f12) ** 2 + sp.Abs(f13) ** 2 + sp.Abs(f14) ** 2

    B = (
        f14 * sp.conjugate(f13) + f13 * sp.conjugate(f14)
    ) ** 2 + 4 * f13 * sp.conjugate(f13) * sp.im(f12) ** 2

    e_inner = sp.sqrt(A - sp.sqrt(B))
    e_outer = sp.sqrt(A + sp.sqrt(B))

    spectrum = (-e_outer, -e_inner, e_inner, e_outer)
    return sp.simplify(A), sp.simplify(B), spectrum


# -----------------------------------------------------------------------------
# 4. Numerical version (recommended for dispersions E_n(q))
# -----------------------------------------------------------------------------

def ham_qsl_numeric(
    f12_: complex,
    f13_: complex,
    f14_: complex,
    f23_: complex,
    f24_: complex,
    f34_: complex,
) -> np.ndarray:
    """Numerical Hermitian matrix, same convention as ham_qsl_symbolic."""
    return np.array(
        [
            [0.0,              f12_,              f13_,              f14_],
            [np.conj(f12_),    0.0,               f23_,              f24_],
            [np.conj(f13_),    np.conj(f23_),     0.0,               f34_],
            [np.conj(f14_),    np.conj(f24_),     np.conj(f34_),     0.0],
        ],
        dtype=np.complex128,
    )


def eigensystem_numeric(H: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    Eigenvalues/eigenvectors for a Hermitian matrix.

    np.linalg.eigh is preferable to np.linalg.eig because it exploits
    Hermiticity and returns real eigenvalues in ascending order.
    """
    if not np.allclose(H, H.conj().T, atol=1e-12, rtol=1e-12):
        raise ValueError("Hamiltonian is not Hermitian within numerical tolerance.")

    eigenvalues, eigenvectors = np.linalg.eigh(H)
    return eigenvalues, eigenvectors


def zero_flux_numeric(f12_: complex, f13_: complex, f14_: complex) -> np.ndarray:
    """Numerical zero-flux matrix with the same reductions as the notebook."""
    return ham_qsl_numeric(
        f12_,
        f13_,
        f14_,
        f14_,
        f13_,
        -np.conj(f12_),
    )


def zero_flux_closed_form_numeric(
    f12_: complex, f13_: complex, f14_: complex
) -> np.ndarray:
    """Numerical evaluation of the analytic zero-flux eigenvalue formula."""
    A = abs(f12_) ** 2 + abs(f13_) ** 2 + abs(f14_) ** 2

    # The first quantity below is real analytically; np.real_if_close removes
    # tiny floating-point imaginary residues.
    cross = f14_ * np.conj(f13_) + f13_ * np.conj(f14_)
    B = cross**2 + 4.0 * abs(f13_) ** 2 * (np.imag(f12_) ** 2)
    B = float(np.real_if_close(B))

    # Clamp tiny negative roundoff before square roots.
    sqrt_B = np.sqrt(max(B, 0.0))
    inner_minus = max(A - sqrt_B, 0.0)
    inner_plus = max(A + sqrt_B, 0.0)

    e_inner = np.sqrt(inner_minus)
    e_outer = np.sqrt(inner_plus)

    return np.array([-e_outer, -e_inner, e_inner, e_outer], dtype=float)


# -----------------------------------------------------------------------------
# 5. Basic consistency tests
# -----------------------------------------------------------------------------

def self_test() -> None:
    """Compare numerical diagonalization against the closed-form zero-flux spectrum."""
    a = 0.7 + 0.4j
    b = -0.2 + 0.8j
    c = 0.5 - 0.1j

    H = zero_flux_numeric(a, b, c)
    evals, _ = eigensystem_numeric(H)
    evals_formula = zero_flux_closed_form_numeric(a, b, c)

    assert np.allclose(evals, evals_formula, atol=1e-12, rtol=1e-12)
    assert np.allclose(H, H.conj().T, atol=1e-12, rtol=1e-12)
    assert np.isclose(np.trace(H), 0.0)

    print("Self-test passed.")
    print("Numerical eigenvalues:", evals)
    print("Closed-form eigenvalues:", evals_formula)


if __name__ == "__main__":
    print("Generic Hamiltonian H =")
    sp.pprint(H_generic)

    print("\nCharacteristic polynomial det(lambda I - H) =")
    sp.pprint(generic_characteristic_polynomial())

    A, B, spectrum = zero_flux_closed_form()
    print("\nZero-flux reduction:")
    print("A =")
    sp.pprint(A)
    print("B =")
    sp.pprint(B)
    print("Spectrum =")
    sp.pprint(spectrum)

    print("\nRunning numerical consistency test...")
    self_test()
