#!/usr/bin/env python3
"""Compare OLD and NEW gauge sectors using the validated final Bloch Hamiltonians.

The Hamiltonians and ground-state energy definition are frozen.  Band tracking is
used only for visualization so that plotted lines follow eigenvector continuity
through crossings instead of reconnecting independently sorted eigenvalues.
"""

from __future__ import annotations

import csv
from itertools import permutations
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


# ============================================================================
# PARAMETERS TO EDIT
# ============================================================================
Kx = 1.0
Ky = 1.0
Kz = 1.0
Kt = 1.0
Kw = 1.0

NUM_Q_POINTS = 4001
FIGURE_FILE = Path("dispersion_old_vs_new.png")
CSV_FILE = Path("energy_comparison.csv")

# Numerical tolerances used only by validation / visualization tracking.
TRACKING_ATOL = 1e-11
DEGENERACY_ATOL = 1e-10

# All 4! assignments are cheap enough to test explicitly.
FOUR_BAND_PERMUTATIONS = tuple(permutations(range(4)))


def bloch_coefficients_old(q, Kx_, Ky_, Kz_, Kt_, Kw_):
    """Validated OLD-sector coefficients."""
    f12 = 1j * (Kx_ * np.exp(+1j * q) - Kw_ * np.exp(-1j * q))
    f13 = -1j * (Kt_ * np.exp(+1j * q) + Kz_ * np.exp(-1j * q))
    f14 = 1j * Ky_
    f23 = 1j * Ky_
    f24 = -1j * (Kt_ * np.exp(+1j * q) + Kz_ * np.exp(-1j * q))
    f34 = +1j * (Kx_ * np.exp(-1j * q) - Kw_ * np.exp(+1j * q))
    return f12, f13, f14, f23, f24, f34


def bloch_coefficients_new(q, Kx_, Ky_, Kz_, Kt_, Kw_):
    """Validated NEW-sector coefficients; only f34 changes relative to OLD."""
    f12, f13, f14, f23, f24, f34_old = bloch_coefficients_old(
        q, Kx_, Ky_, Kz_, Kt_, Kw_
    )
    f34 = -f34_old
    return f12, f13, f14, f23, f24, f34


def _hamiltonian_from_coefficients(coeffs):
    f12, f13, f14, f23, f24, f34 = coeffs
    return np.array(
        [
            [0.0, f12, f13, f14],
            [np.conj(f12), 0.0, f23, f24],
            [np.conj(f13), np.conj(f23), 0.0, f34],
            [np.conj(f14), np.conj(f24), np.conj(f34), 0.0],
        ],
        dtype=np.complex128,
    )


def hamiltonian_old(q, Kx_, Ky_, Kz_, Kt_, Kw_):
    return _hamiltonian_from_coefficients(
        bloch_coefficients_old(q, Kx_, Ky_, Kz_, Kt_, Kw_)
    )


def hamiltonian_new(q, Kx_, Ky_, Kz_, Kt_, Kw_):
    return _hamiltonian_from_coefficients(
        bloch_coefficients_new(q, Kx_, Ky_, Kz_, Kt_, Kw_)
    )


# ============================================================================
# A) SORTED SPECTRUM FOR PHYSICS / ENERGY -- DO NOT BAND-TRACK THIS SECTION
# ============================================================================
def bands_old(q_values, couplings):
    """Four sorted real bands E1 <= E2 <= E3 <= E4 for OLD."""
    return np.vstack(
        [np.linalg.eigvalsh(hamiltonian_old(q, *couplings)) for q in q_values]
    )


def bands_new(q_values, couplings):
    """Four sorted real bands E1 <= E2 <= E3 <= E4 for NEW."""
    return np.vstack(
        [np.linalg.eigvalsh(hamiltonian_new(q, *couplings)) for q in q_values]
    )


def _energy_per_site_from_bands(q_values, band_values):
    """epsilon = (1/(4*pi)) integral_BZ dq sum_{E_n(q)<0} E_n(q)."""
    occupied_sum = np.where(band_values < 0.0, band_values, 0.0).sum(axis=1)
    return float(np.trapezoid(occupied_sum, q_values) / (4.0 * np.pi))


def energy_per_site_old(q_values, couplings):
    return _energy_per_site_from_bands(q_values, bands_old(q_values, couplings))


def energy_per_site_new(q_values, couplings):
    return _energy_per_site_from_bands(q_values, bands_new(q_values, couplings))


# ============================================================================
# B) TRACKED SPECTRUM FOR VISUALIZATION ONLY
# ============================================================================
def _best_overlap_permutation(previous_vectors, current_vectors):
    """Return the permutation maximizing total squared eigenvector overlap."""
    overlap = np.abs(previous_vectors.conj().T @ current_vectors) ** 2

    best_perm = None
    best_score = -np.inf
    for perm in FOUR_BAND_PERMUTATIONS:
        score = sum(overlap[m, perm[m]] for m in range(4))
        if score > best_score:
            best_score = score
            best_perm = perm

    return np.asarray(best_perm, dtype=int)


def _degenerate_groups(values, atol=DEGENERACY_ATOL):
    """Return label groups whose eigenvalues are numerically degenerate."""
    values = np.asarray(values, dtype=float)
    scale = max(1.0, float(np.max(np.abs(values))))
    threshold = atol * scale

    unused = set(range(values.size))
    groups = []
    while unused:
        seed = unused.pop()
        group = {seed}
        changed = True
        while changed:
            changed = False
            for i in list(unused):
                if any(abs(values[i] - values[j]) <= threshold for j in group):
                    group.add(i)
                    unused.remove(i)
                    changed = True
        groups.append(tuple(sorted(group)))
    return groups


def _parallel_transport_degenerate_subspaces(previous_vectors, current_vectors, current_values):
    """Stabilize eigenvector labels inside exactly/nearly degenerate subspaces.

    At a degeneracy, individual eigenvectors are not unique.  Rotating vectors
    within the degenerate subspace does not alter any eigenvalue.  We choose the
    unitary rotation that best aligns the current subspace with the previously
    tracked states (orthogonal/Procrustes parallel transport).  This prevents an
    arbitrary LAPACK basis choice at the degeneracy from contaminating the next
    tracking step.
    """
    transported = current_vectors.copy()

    for group in _degenerate_groups(current_values):
        if len(group) <= 1:
            continue

        idx = np.asarray(group, dtype=int)
        previous_block = previous_vectors[:, idx]
        current_block = transported[:, idx]

        # Solve current_block @ Q ~= previous_block with unitary Q.
        cross = current_block.conj().T @ previous_block
        U, _, Vh = np.linalg.svd(cross)
        Q = U @ Vh
        transported[:, idx] = current_block @ Q

    return transported


def tracked_band_structure(hamiltonian_function, q_values, couplings):
    """Track four plotted bands by eigenvector continuity.

    The first q point uses the energy ordering returned by np.linalg.eigh.
    At every subsequent point, the 4! permutations are searched and the one
    maximizing total squared overlap with the previous tracked eigenvectors is
    selected.  Degenerate subspaces are parallel-transported after assignment.

    IMPORTANT: this function is for plotting only.  Ground-state energies remain
    computed from independently sorted np.linalg.eigvalsh spectra.
    """
    q_values = np.asarray(q_values, dtype=float)
    if q_values.ndim != 1 or q_values.size == 0:
        raise ValueError("q_values must be a non-empty one-dimensional array.")

    tracked_values = np.empty((q_values.size, 4), dtype=float)

    first_values, first_vectors = np.linalg.eigh(
        hamiltonian_function(q_values[0], *couplings)
    )
    tracked_values[0] = first_values
    previous_vectors = first_vectors

    for i, q in enumerate(q_values[1:], start=1):
        current_values, current_vectors = np.linalg.eigh(
            hamiltonian_function(q, *couplings)
        )

        perm = _best_overlap_permutation(previous_vectors, current_vectors)
        current_values = current_values[perm]
        current_vectors = current_vectors[:, perm]

        # At a true degeneracy, only the subspace is physical; choose a basis
        # that continues smoothly from the previous step without moving levels.
        current_vectors = _parallel_transport_degenerate_subspaces(
            previous_vectors, current_vectors, current_values
        )

        tracked_values[i] = current_values
        previous_vectors = current_vectors

    return tracked_values


# ============================================================================
# VALIDATION
# ============================================================================
def run_consistency_checks(couplings):
    """Original physics checks; any failure aborts before scientific output."""
    atol = 1e-12
    q_checks = (-np.pi / 2, -0.73, 0.0, 0.41, np.pi / 2)

    if len(couplings) != 5 or not all(np.isfinite(couplings)):
        raise ValueError("Kx, Ky, Kz, Kt, Kw must be five finite real numbers.")

    for q in q_checks:
        old = bloch_coefficients_old(q, *couplings)
        new = bloch_coefficients_new(q, *couplings)
        H_old = hamiltonian_old(q, *couplings)
        H_new = hamiltonian_new(q, *couplings)

        # 1-2. Hermiticity.
        if not np.allclose(H_old, H_old.conj().T, atol=atol, rtol=atol):
            raise AssertionError(f"H_old is not Hermitian at q={q}.")
        if not np.allclose(H_new, H_new.conj().T, atol=atol, rtol=atol):
            raise AssertionError(f"H_new is not Hermitian at q={q}.")

        # 3. Four real eigenvalues from eigvalsh.
        for sector, H in (("OLD", H_old), ("NEW", H_new)):
            eigs = np.linalg.eigvalsh(H)
            if eigs.shape != (4,) or not np.isrealobj(eigs) or not np.all(np.isfinite(eigs)):
                raise AssertionError(f"Invalid {sector} eigenvalues at q={q}: {eigs}")

        # 4-5. Same K_gamma; only f34 changes.
        if not np.allclose(old[:5], new[:5], atol=atol, rtol=atol):
            raise AssertionError(f"A coefficient other than f34 changed at q={q}.")
        if not np.allclose(new[5], -old[5], atol=atol, rtol=atol):
            raise AssertionError(f"f34 did not change sign at q={q}.")

        diff = H_new - H_old
        mask = np.ones((4, 4), dtype=bool)
        mask[2, 3] = mask[3, 2] = False
        if not np.allclose(diff[mask], 0.0, atol=atol, rtol=atol):
            raise AssertionError(f"Matrix entries other than f34 changed at q={q}.")

        # 6-7. Required OLD/NEW identities.
        f12_old, f34_old = old[0], old[5]
        f12_new, f34_new = new[0], new[5]
        if not np.allclose(f34_old, -np.conj(f12_old), atol=atol, rtol=atol):
            raise AssertionError(f"OLD: f34 != -conjugate(f12) at q={q}.")
        if not np.allclose(f34_new, +np.conj(f12_new), atol=atol, rtol=atol):
            raise AssertionError(f"NEW: f34 != +conjugate(f12) at q={q}.")

    # Isotropic NEW sanity checks at +/- pi/4.
    isotropic_checks = {
        +np.pi / 4: np.array([-1.8284271247461903, -1.0, -1.0, 3.8284271247461903]),
        -np.pi / 4: np.array([-3.8284271247461903, 1.0, 1.0, 1.8284271247461903]),
    }
    for q, expected in isotropic_checks.items():
        actual = np.linalg.eigvalsh(hamiltonian_new(q, 1, 1, 1, 1, 1))
        if not np.allclose(actual, expected, atol=2e-12, rtol=2e-12):
            raise AssertionError(
                f"Isotropic NEW sanity check failed at q={q}: "
                f"got {actual}, expected {expected}."
            )


def _analytic_new_isotropic_branches(q_values):
    """Analytic NEW isotropic branches, used only as a numerical sanity check."""
    q = np.asarray(q_values, dtype=float)
    return np.column_stack(
        [
            -2.0 * np.sin(q) + 2.0 * np.cos(q) - 1.0,
            +2.0 * np.sin(q) + 2.0 * np.cos(q) + 1.0,
            +2.0 * np.sin(q) - 2.0 * np.cos(q) - 1.0,
            -2.0 * np.sin(q) - 2.0 * np.cos(q) + 1.0,
        ]
    )


def run_band_tracking_checks(
    q_values,
    couplings,
    sorted_old,
    sorted_new,
    tracked_old,
    tracked_new,
    epsilon_old_before,
    epsilon_new_before,
    delta_before,
):
    """Prove that tracking only reorders the spectrum and cannot change energy."""
    atol = TRACKING_ATOL

    # Fundamental test at every q: tracking may only permute the same eigenvalues.
    if not np.allclose(
        np.sort(tracked_old, axis=1), sorted_old, atol=atol, rtol=atol
    ):
        raise AssertionError("OLD tracked bands do not reproduce the sorted spectrum.")
    if not np.allclose(
        np.sort(tracked_new, axis=1), sorted_new, atol=atol, rtol=atol
    ):
        raise AssertionError("NEW tracked bands do not reproduce the sorted spectrum.")

    # Explicitly repeat the eigvalsh comparison at every q as a hard guard.
    for i, q in enumerate(q_values):
        if not np.allclose(
            np.sort(tracked_old[i]),
            np.linalg.eigvalsh(hamiltonian_old(q, *couplings)),
            atol=atol,
            rtol=atol,
        ):
            raise AssertionError(f"OLD tracking changed the spectrum at q={q}.")
        if not np.allclose(
            np.sort(tracked_new[i]),
            np.linalg.eigvalsh(hamiltonian_new(q, *couplings)),
            atol=atol,
            rtol=atol,
        ):
            raise AssertionError(f"NEW tracking changed the spectrum at q={q}.")

    # Recompute the physical energy after tracking from fresh SORTED spectra only.
    # Tracked bands never enter this calculation.
    sorted_old_after = bands_old(q_values, couplings)
    sorted_new_after = bands_new(q_values, couplings)
    epsilon_old_after = _energy_per_site_from_bands(q_values, sorted_old_after)
    epsilon_new_after = _energy_per_site_from_bands(q_values, sorted_new_after)
    delta_after = epsilon_new_after - epsilon_old_after

    if not np.isclose(epsilon_old_before, epsilon_old_after, atol=1e-14, rtol=1e-14):
        raise AssertionError("epsilon_old changed after band tracking was introduced.")
    if not np.isclose(epsilon_new_before, epsilon_new_after, atol=1e-14, rtol=1e-14):
        raise AssertionError("epsilon_new changed after band tracking was introduced.")
    if not np.isclose(delta_before, delta_after, atol=1e-14, rtol=1e-14):
        raise AssertionError("Delta_epsilon changed after band tracking was introduced.")

    # Special isotropic NEW check against the known analytic branch set.
    if np.allclose(couplings, (1.0, 1.0, 1.0, 1.0, 1.0), atol=0.0, rtol=0.0):
        analytic = _analytic_new_isotropic_branches(q_values)
        if not np.allclose(
            np.sort(tracked_new, axis=1),
            np.sort(analytic, axis=1),
            atol=2e-11,
            rtol=2e-11,
        ):
            raise AssertionError(
                "Isotropic tracked NEW spectrum disagrees with analytic branch set."
            )

    return epsilon_old_after, epsilon_new_after, delta_after


# ============================================================================
# OUTPUT
# ============================================================================
def save_band_figure(q_values, old_tracked_bands, new_tracked_bands):
    """Plot eigenvector-tracked bands; labels denote continuity, not energy rank."""
    fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.2), sharex=True, sharey=True)

    for n in range(4):
        label = f"Tracked band {n + 1}"
        axes[0].plot(q_values, old_tracked_bands[:, n], lw=1.2, label=label)
        axes[1].plot(q_values, new_tracked_bands[:, n], lw=1.2, label=label)

    y_min = float(min(old_tracked_bands.min(), new_tracked_bands.min()))
    y_max = float(max(old_tracked_bands.max(), new_tracked_bands.max()))
    pad = 0.05 * max(y_max - y_min, 1.0)

    for ax, title in zip(
        axes, ("Original gauge (OLD)", "Upper bonds reversed (NEW)")
    ):
        ax.axhline(0.0, color="black", lw=0.8, alpha=0.65)
        ax.set_title(title)
        ax.set_xlabel(r"$q = k a_1$")
        ax.set_xlim(-np.pi / 2, np.pi / 2)
        ax.set_ylim(y_min - pad, y_max + pad)
        ax.grid(alpha=0.20)
        ax.legend(frameon=False, fontsize=8, ncol=2)

    axes[0].set_ylabel(r"$E(q)$")
    fig.tight_layout()
    fig.savefig(FIGURE_FILE, dpi=300, bbox_inches="tight")


def save_energy_csv(couplings, epsilon_old, epsilon_new, delta_epsilon):
    header = [
        "Kx", "Ky", "Kz", "Kt", "Kw",
        "epsilon_old", "epsilon_new", "Delta_epsilon",
    ]
    row = [*couplings, epsilon_old, epsilon_new, delta_epsilon]
    with CSV_FILE.open("w", newline="", encoding="utf-8") as file:
        writer = csv.writer(file)
        writer.writerow(header)
        writer.writerow(row)


def favored_sector(delta_epsilon, couplings):
    tolerance = 1e-10 * max(1.0, *(abs(K) for K in couplings))
    if abs(delta_epsilon) <= tolerance:
        return "DEGENERATE"
    return "OLD" if delta_epsilon > 0.0 else "NEW"


def print_results(couplings, epsilon_old, epsilon_new, delta_epsilon):
    print("Couplings used:")
    for name, value in zip(("Kx", "Ky", "Kz", "Kt", "Kw"), couplings):
        print(f"{name} = {value:.12g}")

    print("\n----------------------------------------------------")
    print(f"{'Quantity':<24}Value")
    print("----------------------------------------------------")
    print(f"{'epsilon_old':<24}{epsilon_old:.12f}")
    print(f"{'epsilon_new':<24}{epsilon_new:.12f}")
    print(f"{'Delta_epsilon':<24}{delta_epsilon:.12f}")
    print("----------------------------------------------------\n")
    print(f"Energetically favored sector: {favored_sector(delta_epsilon, couplings)}")


def main():
    couplings = (float(Kx), float(Ky), float(Kz), float(Kt), float(Kw))
    run_consistency_checks(couplings)

    q_values = np.linspace(-np.pi / 2, np.pi / 2, NUM_Q_POINTS)

    # Physics / energy: keep the original independently sorted eigvalsh spectra.
    sorted_old = bands_old(q_values, couplings)
    sorted_new = bands_new(q_values, couplings)
    epsilon_old = _energy_per_site_from_bands(q_values, sorted_old)
    epsilon_new = _energy_per_site_from_bands(q_values, sorted_new)
    delta_epsilon = epsilon_new - epsilon_old

    # Visualization only: track eigenstates continuously with np.linalg.eigh.
    tracked_old = tracked_band_structure(hamiltonian_old, q_values, couplings)
    tracked_new = tracked_band_structure(hamiltonian_new, q_values, couplings)

    # Hard guards: same spectrum at every q and exactly unchanged energies.
    epsilon_old_after, epsilon_new_after, delta_after = run_band_tracking_checks(
        q_values,
        couplings,
        sorted_old,
        sorted_new,
        tracked_old,
        tracked_new,
        epsilon_old,
        epsilon_new,
        delta_epsilon,
    )

    save_band_figure(q_values, tracked_old, tracked_new)
    save_energy_csv(couplings, epsilon_old_after, epsilon_new_after, delta_after)
    print_results(couplings, epsilon_old_after, epsilon_new_after, delta_after)
    plt.show()


if __name__ == "__main__":
    main()
