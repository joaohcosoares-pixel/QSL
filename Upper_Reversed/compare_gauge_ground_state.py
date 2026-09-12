#!/usr/bin/env python3
"""
Compare OLD and NEW Z2-gauge sectors of the audited 1D quadrupolar spin liquid.

This is a compact results script.  The Bloch coefficients are copied from the
already-audited implementation; no Fourier derivation, gauge convention, or
Hamiltonian physics is changed here.

Outputs
-------
- dispersion_old_vs_new.png
- Fq_old_vs_new.png
- ground_state_energy_comparison.csv
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


# =============================================================================
# User-editable parameters
# =============================================================================

Kx = 1.0
Ky = 1.0
Kz = 1.0
Kt = 1.0
Kw = 1.0

NQ_PHYSICAL = 4001
NQ_EXTENDED = 8001

HERMITICITY_TOL = 1.0e-12
GAP_TOL = 1.0e-8
ENERGY_TOL = 1.0e-10
NORMALIZATION_TOL = 5.0e-10

DISPERSION_FILE = "dispersion_old_vs_new.png"
FQ_FILE = "Fq_old_vs_new.png"
CSV_FILE = "ground_state_energy_comparison.csv"


# =============================================================================
# Audited Bloch coefficients: q = k a1
# =============================================================================


def bloch_coefficients_old(
    q: float,
    Kx_: float,
    Ky_: float,
    Kz_: float,
    Kt_: float,
    Kw_: float,
) -> dict[str, complex]:
    """OLD gauge: original upper horizontal x and w orientations."""
    ep = np.exp(+1j * q)
    em = np.exp(-1j * q)

    f12 = 1j * (Kx_ * ep - Kw_ * em)
    f13 = -1j * (Kt_ * ep + Kz_ * em)
    f14 = 1j * Ky_
    f23 = 1j * Ky_
    f24 = -1j * (Kt_ * ep + Kz_ * em)
    f34 = 1j * (Kx_ * em - Kw_ * ep)

    return {
        "f12": f12,
        "f13": f13,
        "f14": f14,
        "f23": f23,
        "f24": f24,
        "f34": f34,
    }


def bloch_coefficients_new(
    q: float,
    Kx_: float,
    Ky_: float,
    Kz_: float,
    Kt_: float,
    Kw_: float,
) -> dict[str, complex]:
    """NEW gauge: upper horizontal x and w links inverted."""
    ep = np.exp(+1j * q)
    em = np.exp(-1j * q)

    f12 = 1j * (Kx_ * ep - Kw_ * em)
    f13 = -1j * (Kt_ * ep + Kz_ * em)
    f14 = 1j * Ky_
    f23 = 1j * Ky_
    f24 = -1j * (Kt_ * ep + Kz_ * em)
    f34 = -1j * (Kx_ * em - Kw_ * ep)

    return {
        "f12": f12,
        "f13": f13,
        "f14": f14,
        "f23": f23,
        "f24": f24,
        "f34": f34,
    }


def _hamiltonian_from_coefficients(c: dict[str, complex]) -> np.ndarray:
    """4x4 Hermitian Bloch matrix in the audited (1,2,3,4) basis."""
    f12 = c["f12"]
    f13 = c["f13"]
    f14 = c["f14"]
    f23 = c["f23"]
    f24 = c["f24"]
    f34 = c["f34"]

    return np.array(
        [
            [0.0, f12, f13, f14],
            [np.conj(f12), 0.0, f23, f24],
            [np.conj(f13), np.conj(f23), 0.0, f34],
            [np.conj(f14), np.conj(f24), np.conj(f34), 0.0],
        ],
        dtype=np.complex128,
    )


def hamiltonian_old(
    q: float,
    Kx_: float,
    Ky_: float,
    Kz_: float,
    Kt_: float,
    Kw_: float,
) -> np.ndarray:
    return _hamiltonian_from_coefficients(
        bloch_coefficients_old(q, Kx_, Ky_, Kz_, Kt_, Kw_)
    )


def hamiltonian_new(
    q: float,
    Kx_: float,
    Ky_: float,
    Kz_: float,
    Kt_: float,
    Kw_: float,
) -> np.ndarray:
    return _hamiltonian_from_coefficients(
        bloch_coefficients_new(q, Kx_, Ky_, Kz_, Kt_, Kw_)
    )


# =============================================================================
# Eigenvalues, bands, F(q), gap, and ground-state energies
# =============================================================================


def _sorted_eigenvalues(H: np.ndarray) -> np.ndarray:
    """Real eigenvalues of a Hermitian matrix, sorted in ascending order."""
    if not np.allclose(H, H.conj().T, atol=HERMITICITY_TOL, rtol=HERMITICITY_TOL):
        raise ValueError("Bloch Hamiltonian is not Hermitian within tolerance.")
    eigs = np.linalg.eigvalsh(H)
    if not np.isrealobj(eigs) or not np.all(np.isfinite(eigs)):
        raise FloatingPointError("Non-real or non-finite eigenvalues detected.")
    return eigs


def bands_sorted(
    hamiltonian,
    q_values: np.ndarray,
    couplings: tuple[float, float, float, float, float],
) -> np.ndarray:
    """Sorted eigenvalues for plotting and diagnostics; shape=(Nq,4)."""
    return np.vstack(
        [_sorted_eigenvalues(hamiltonian(float(q), *couplings)) for q in q_values]
    )


def F_old(
    q: float,
    Kx_: float,
    Ky_: float,
    Kz_: float,
    Kt_: float,
    Kw_: float,
) -> float:
    """F_OLD(q) = sum of the strictly negative eigvalsh eigenvalues."""
    eigs = _sorted_eigenvalues(hamiltonian_old(q, Kx_, Ky_, Kz_, Kt_, Kw_))
    return float(np.sum(eigs[eigs < 0.0]))


def F_new(
    q: float,
    Kx_: float,
    Ky_: float,
    Kz_: float,
    Kt_: float,
    Kw_: float,
) -> float:
    """F_NEW(q) = sum of the strictly negative eigvalsh eigenvalues."""
    eigs = _sorted_eigenvalues(hamiltonian_new(q, Kx_, Ky_, Kz_, Kt_, Kw_))
    return float(np.sum(eigs[eigs < 0.0]))


def _trapezoid(y: np.ndarray, x: np.ndarray) -> float:
    """Use np.trapezoid when available, with a compatibility fallback."""
    if hasattr(np, "trapezoid"):
        return float(np.trapezoid(y, x))
    return float(np.trapz(y, x))


def _F_array(F_function, q_values: np.ndarray, couplings: tuple[float, ...]) -> np.ndarray:
    values = np.array([F_function(float(q), *couplings) for q in q_values], dtype=float)
    if not np.all(np.isfinite(values)):
        raise FloatingPointError("F(q) contains non-finite values.")
    return values


def ground_state_energy_from_F(
    F_values: np.ndarray,
    q_values: np.ndarray,
    extended_bz: bool = False,
) -> tuple[float, float]:
    """Return (energy per four-site cell, energy per site)."""
    integral = _trapezoid(F_values, q_values)
    prefactor_cell = 1.0 / (2.0 * np.pi) if extended_bz else 1.0 / np.pi
    energy_per_cell = prefactor_cell * integral
    energy_per_site = energy_per_cell / 4.0
    return float(energy_per_cell), float(energy_per_site)


def spectral_gap_from_bands(bands: np.ndarray, gap_tol: float = GAP_TOL) -> float:
    """
    Zero-energy spectral gap on the continuous physical BZ.

    A finite mesh may miss an exact zero between neighboring q points.  If the
    number of negative eigenvalues changes between adjacent points, continuity
    of a Hermitian spectrum guarantees a zero-energy crossing in that interval,
    so the gap is exactly zero.  Otherwise the global positive/negative band
    edges determine the Fermi-level gap.
    """
    if not np.all(np.isfinite(bands)):
        raise FloatingPointError("Band array contains non-finite values.")

    if np.any(np.abs(bands) <= gap_tol):
        return 0.0

    negative_count = np.sum(bands < 0.0, axis=1)
    if np.any(np.diff(negative_count) != 0):
        return 0.0

    positive = bands[bands > 0.0]
    negative = bands[bands < 0.0]
    if positive.size == 0 or negative.size == 0:
        return 0.0

    gap = float(np.min(positive) - np.max(negative))
    return max(0.0, gap)


def spectral_character(gap: float) -> str:
    return "gapped" if gap > GAP_TOL else "gapless"


def favored_sector(delta_epsilon: float) -> str:
    if delta_epsilon > ENERGY_TOL:
        return "OLD"
    if delta_epsilon < -ENERGY_TOL:
        return "NEW"
    return "DEGENERATE"


# =============================================================================
# Minimal numerical checks
# =============================================================================


def run_minimal_tests(
    couplings: tuple[float, float, float, float, float],
    q_physical: np.ndarray,
    q_extended: np.ndarray,
) -> dict[str, float]:
    """Run the requested compact consistency checks before final results."""
    assert len(couplings) == 5
    assert np.all(np.isfinite(couplings))

    # OLD and NEW receive the exact same coupling tuple by construction.
    test_q = np.linspace(-np.pi / 2.0, np.pi / 2.0, 21)
    max_hermiticity_error = 0.0

    for q in test_q:
        for H in (
            hamiltonian_old(float(q), *couplings),
            hamiltonian_new(float(q), *couplings),
        ):
            herm_err = float(np.max(np.abs(H - H.conj().T)))
            max_hermiticity_error = max(max_hermiticity_error, herm_err)
            assert herm_err < HERMITICITY_TOL

            eigs = np.linalg.eigvalsh(H)
            assert np.isrealobj(eigs)
            assert np.all(np.isfinite(eigs))

    F_old_phys = _F_array(F_old, q_physical, couplings)
    F_new_phys = _F_array(F_new, q_physical, couplings)
    assert np.all(np.isfinite(F_old_phys))
    assert np.all(np.isfinite(F_new_phys))

    Ecell_old, Esite_old = ground_state_energy_from_F(F_old_phys, q_physical)
    Ecell_new, Esite_new = ground_state_energy_from_F(F_new_phys, q_physical)
    assert np.isfinite(Ecell_old) and np.isfinite(Esite_old)
    assert np.isfinite(Ecell_new) and np.isfinite(Esite_new)
    assert np.isclose(Esite_old, Ecell_old / 4.0, atol=1e-13, rtol=1e-13)
    assert np.isclose(Esite_new, Ecell_new / 4.0, atol=1e-13, rtol=1e-13)

    # Extended-BZ normalization check for both sectors.
    F_old_ext = _F_array(F_old, q_extended, couplings)
    F_new_ext = _F_array(F_new, q_extended, couplings)
    _, Esite_old_ext = ground_state_energy_from_F(F_old_ext, q_extended, extended_bz=True)
    _, Esite_new_ext = ground_state_energy_from_F(F_new_ext, q_extended, extended_bz=True)

    old_norm_err = abs(Esite_old - Esite_old_ext)
    new_norm_err = abs(Esite_new - Esite_new_ext)
    assert old_norm_err < NORMALIZATION_TOL
    assert new_norm_err < NORMALIZATION_TOL

    return {
        "max_hermiticity_error": max_hermiticity_error,
        "old_normalization_error": old_norm_err,
        "new_normalization_error": new_norm_err,
    }


# =============================================================================
# Output helpers
# =============================================================================


def _format_q_axis(ax) -> None:
    ticks = [-np.pi / 2.0, -np.pi / 4.0, 0.0, np.pi / 4.0, np.pi / 2.0]
    labels = [r"$-\pi/2$", r"$-\pi/4$", "0", r"$\pi/4$", r"$\pi/2$"]
    ax.set_xlim(-np.pi / 2.0, np.pi / 2.0)
    ax.set_xticks(ticks, labels)
    ax.grid(True, alpha=0.25)


def plot_dispersions(
    q_values: np.ndarray,
    bands_old: np.ndarray,
    bands_new: np.ndarray,
    output_path: str | Path = DISPERSION_FILE,
) -> None:
    """Save side-by-side OLD/NEW dispersions over the physical BZ."""
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.6), sharey=True)

    for band in range(4):
        axes[0].plot(q_values, bands_old[:, band], linewidth=1.5)
        axes[1].plot(q_values, bands_new[:, band], linewidth=1.5)

    for ax, title in zip(axes, ("OLD gauge", "NEW gauge")):
        ax.axhline(0.0, linewidth=0.9, linestyle="--")
        ax.set_xlabel(r"$q = k a_1$")
        ax.set_title(title)
        _format_q_axis(ax)

    axes[0].set_ylabel(r"$E_n(q)$")
    fig.suptitle("Band dispersions: OLD vs NEW gauge")
    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_Fq(
    q_values: np.ndarray,
    F_old_values: np.ndarray,
    F_new_values: np.ndarray,
    output_path: str | Path = FQ_FILE,
) -> None:
    """Save F_OLD(q) and F_NEW(q) on the same physical-BZ panel."""
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    ax.plot(q_values, F_old_values, label="OLD", linewidth=1.8)
    ax.plot(q_values, F_new_values, label="NEW", linewidth=1.8)
    ax.set_xlabel(r"$q = k a_1$")
    ax.set_ylabel(r"$F(q)$")
    ax.set_title(r"Occupied-band sum $F(q)$: OLD vs NEW")
    ax.legend()
    _format_q_axis(ax)
    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def save_results_csv(
    couplings: tuple[float, float, float, float, float],
    epsilon_old: float,
    epsilon_new: float,
    delta_epsilon: float,
    favored: str,
    output_path: str | Path = CSV_FILE,
) -> None:
    fieldnames = [
        "Kx",
        "Ky",
        "Kz",
        "Kt",
        "Kw",
        "epsilon_old",
        "epsilon_new",
        "Delta_epsilon",
        "favored_sector",
    ]
    row = {
        "Kx": couplings[0],
        "Ky": couplings[1],
        "Kz": couplings[2],
        "Kt": couplings[3],
        "Kw": couplings[4],
        "epsilon_old": epsilon_old,
        "epsilon_new": epsilon_new,
        "Delta_epsilon": delta_epsilon,
        "favored_sector": favored,
    }

    with open(output_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerow(row)


# =============================================================================
# Main
# =============================================================================


def main() -> None:
    couplings = (float(Kx), float(Ky), float(Kz), float(Kt), float(Kw))
    q_values = np.linspace(-np.pi / 2.0, np.pi / 2.0, NQ_PHYSICAL)
    q_extended = np.linspace(-np.pi, np.pi, NQ_EXTENDED)

    checks = run_minimal_tests(couplings, q_values, q_extended)

    # Sorted eigvalsh bands are used directly.  No band tracking enters F(q) or energy.
    bands_old_values = bands_sorted(hamiltonian_old, q_values, couplings)
    bands_new_values = bands_sorted(hamiltonian_new, q_values, couplings)

    F_old_values = _F_array(F_old, q_values, couplings)
    F_new_values = _F_array(F_new, q_values, couplings)

    Ecell_old, epsilon_old = ground_state_energy_from_F(F_old_values, q_values)
    Ecell_new, epsilon_new = ground_state_energy_from_F(F_new_values, q_values)

    delta_epsilon = epsilon_new - epsilon_old
    favored = favored_sector(delta_epsilon)

    gap_old = spectral_gap_from_bands(bands_old_values)
    gap_new = spectral_gap_from_bands(bands_new_values)
    character_old = spectral_character(gap_old)
    character_new = spectral_character(gap_new)

    plot_dispersions(q_values, bands_old_values, bands_new_values)
    plot_Fq(q_values, F_old_values, F_new_values)
    save_results_csv(couplings, epsilon_old, epsilon_new, delta_epsilon, favored)

    print("Minimal numerical checks: PASSED")
    print(f"max Hermiticity error = {checks['max_hermiticity_error']:.3e}")
    print(f"OLD physical/extended BZ energy discrepancy = {checks['old_normalization_error']:.3e}")
    print(f"NEW physical/extended BZ energy discrepancy = {checks['new_normalization_error']:.3e}")

    print("\nCouplings:")
    print(f"Kx = {Kx:.12g}")
    print(f"Ky = {Ky:.12g}")
    print(f"Kz = {Kz:.12g}")
    print(f"Kt = {Kt:.12g}")
    print(f"Kw = {Kw:.12g}")

    print("\nSpectral character:")
    print(f"OLD = {character_old}  (gap = {gap_old:.12g})")
    print(f"NEW = {character_new}  (gap = {gap_new:.12g})")

    print("\nGround-state energy per site:")
    print(f"epsilon_old = {epsilon_old:.12f}")
    print(f"epsilon_new = {epsilon_new:.12f}")
    print(f"Delta_epsilon = {delta_epsilon:.12f}")

    print("\n-------------------------------------------------------")
    print(f"{'Quantity':<34} {'Value':>18}")
    print("-------------------------------------------------------")
    print(f"{'epsilon_old':<34} {epsilon_old:>18.12f}")
    print(f"{'epsilon_new':<34} {epsilon_new:>18.12f}")
    print(f"{'Delta_epsilon':<34} {delta_epsilon:>18.12f}")
    print("-------------------------------------------------------")

    print(f"\nOLD spectral character: {character_old}")
    print(f"NEW spectral character: {character_new}")
    print(f"\nEnergetically favored sector: {favored}")

    print("\nGenerated files:")
    print(f"- {DISPERSION_FILE}")
    print(f"- {FQ_FILE}")
    print(f"- {CSV_FILE}")

    # Explicitly retain the per-cell quantities in the calculation and verify /4.
    assert np.isclose(epsilon_old, Ecell_old / 4.0, atol=1e-13, rtol=1e-13)
    assert np.isclose(epsilon_new, Ecell_new / 4.0, atol=1e-13, rtol=1e-13)


if __name__ == "__main__":
    main()
