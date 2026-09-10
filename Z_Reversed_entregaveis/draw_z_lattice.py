"""Original-PDF geometry, eight oriented triangle audits, and Z-reversed drawing.

The drawing is reconstructed from the original arrows. The previous plotting
script's upper horizontal arrows are not a physical source. All positions and
the original a1=(2,0), a2=(1,2) vectors are preserved. Only z arrows differ
between the two edge sets constructed here. Run with ``python3 -W error``.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from itertools import product
import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/qsl_z_lattice_mpl")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch


ORIGINAL_GAUGE = {
    "u12x": 1, "u34x": 1, "u23y": 1, "u41y": -1,
    "u42z": 1, "u31z": 1, "u43w": 1, "u12w": -1,
    "u13t": -1, "u24t": -1,
}
Z_REVERSED_GAUGE = {
    key: (-value if key.endswith("z") else value)
    for key, value in ORIGINAL_GAUGE.items()
}


@dataclass(frozen=True)
class CellBond:
    start: int
    end: int
    end_cell_shift: int
    kind: str
    key: str


# Both ends are literal sites of the original Eq. (11). The first end lies
# in cell n and the second in cell n+end_cell_shift, with translation 2 a1.
CELL_BONDS = (
    CellBond(1, 2, 0, "x", "u12x"),
    CellBond(3, 4, 0, "x", "u34x"),
    CellBond(2, 3, 0, "y", "u23y"),
    CellBond(4, 1, 0, "y", "u41y"),
    CellBond(4, 2, 0, "z", "u42z"),
    CellBond(3, 1, 1, "z", "u31z"),
    CellBond(4, 3, -1, "w", "u43w"),
    CellBond(1, 2, -1, "w", "u12w"),
    CellBond(1, 3, 0, "t", "u13t"),
    CellBond(2, 4, 1, "t", "u24t"),
)

# Site tuple = (sublattice, cell). Every path is closed and anticlockwise.
PLAQUETTES = {
    "W124": ((1, 0), (2, 0), (4, 0), (1, 0)),
    "W234": ((2, 0), (3, 0), (4, 0), (2, 0)),
    "W123": ((1, 0), (2, 0), (3, 0), (1, 0)),
    "W134": ((1, 0), (3, 0), (4, 0), (1, 0)),
    "W213": ((2, -1), (1, 0), (3, -1), (2, -1)),
    "W143": ((1, 0), (4, 0), (3, -1), (1, 0)),
    "W214": ((2, -1), (1, 0), (4, 0), (2, -1)),
    "W243": ((2, -1), (4, 0), (3, -1), (2, -1)),
}


def position(site: tuple[int, int]) -> tuple[float, float]:
    sublattice, cell = site
    internal = {1: (0., 0.), 2: (2., 0.), 3: (3., 2.), 4: (1., 2.)}
    x, y = internal[sublattice]
    return x + 4. * cell, y


def audit_fluxes(gauge: dict[str, int]) -> dict:
    """Find every directed path edge in the translated real bond graph."""
    directed = {}
    for cell in range(-2, 3):
        for bond in CELL_BONDS:
            a = (bond.start, cell)
            b = (bond.end, cell + bond.end_cell_shift)
            assert (a, b) not in directed and (b, a) not in directed
            directed[a, b] = (gauge[bond.key], bond.key)
            directed[b, a] = (-gauge[bond.key], "-" + bond.key)
    result = {}
    for name, path in PLAQUETTES.items():
        assert path[0] == path[-1]
        coordinates = [position(site) for site in path]
        signed_area = 0.5 * sum(
            a[0] * b[1] - b[0] * a[1]
            for a, b in zip(coordinates[:-1], coordinates[1:])
        )
        assert signed_area > 0, f"Clockwise path: {name}"
        values, expressions = zip(*(directed[a, b] for a, b in zip(path[:-1], path[1:])))
        flux = 1
        for value in values:
            assert value in (-1, 1)
            flux *= value
        result[name] = {
            "path": path, "factors": values, "product": " ".join(f"({e})" for e in expressions),
            "flux": flux, "signed_area": signed_area,
        }
    return result


def geometry_edges(reverse_z: bool) -> list[tuple]:
    b = {i: (float(2 * i), 0.) for i in range(6)}
    t = {i: (float(2 * i + 1), 2.) for i in range(6)}
    edges = []
    for i in range(5):
        kind = "x" if i % 2 == 0 else "w"
        edges.extend(((b[i], b[i + 1], kind), (t[i + 1], t[i], kind)))
    edges.extend((b[i], t[i], "y") for i in range(6))
    edges.extend(
        (b[i + 1], t[i], "z") if reverse_z else (t[i], b[i + 1], "z")
        for i in range(5)
    )
    edges.extend((t[i + 1], b[i], "t") for i in range(5))
    return edges


def draw_lattice(output_path: Path, dpi: int = 300) -> None:
    plt.rcParams.update({"font.family": "DejaVu Sans", "mathtext.fontset": "dejavusans"})
    colors = {"x": "#2474B5", "w": "#D88915", "y": "#CC3333", "z": "#669F2A", "t": "#7A4299"}
    fig = plt.figure(figsize=(11.5, 4.25), dpi=dpi, facecolor="white")
    ax = fig.add_axes((.04, .25, .92, .61))
    for start, end, kind in geometry_edges(True):
        ax.add_patch(FancyArrowPatch(
            start, end, arrowstyle="-|>", mutation_scale=17,
            linewidth=1.65, color=colors[kind], shrinkA=4.7, shrinkB=4.7,
            zorder=2 if kind == "t" else 3,
        ))
    for i in range(6):
        ax.plot(2*i, 0, "o", color="#1B2533", markersize=6.5, zorder=6)
        ax.plot(2*i+1, 2, "o", color="#1B2533", markersize=6.5, zorder=6)
    for label, point, offset, va in (
        ("1", (4, 0), -10, "top"), ("2", (6, 0), -10, "top"),
        ("4", (5, 2), 10, "bottom"), ("3", (7, 2), 10, "bottom"),
    ):
        ax.annotate(label, point, xytext=(0, offset), textcoords="offset points",
                    ha="center", va=va, fontsize=14, fontweight="bold")
    ax.text(-.75, 1., "···", fontsize=21, ha="center", va="center")
    ax.text(11.7, 1., "···", fontsize=21, ha="center", va="center")
    ax.set(xlim=(-1.0, 12.), ylim=(-.65, 2.75), aspect="equal")
    ax.axis("off")
    fig.suptitle("z bonds reversed (Z-REVERSED)", x=.5, y=.975, fontsize=16, weight="bold")
    fig.text(.50, .885, "Original positions, sublattices and bond families", ha="center", fontsize=10.5, color="#566070")
    fig.legend(
        handles=[Line2D([0], [0], color=colors[k], linewidth=2.3, label=k) for k in ("x", "y", "z", "t", "w")],
        loc="lower left", bbox_to_anchor=(.07, .10), ncol=5, frameon=False,
        fontsize=12, columnspacing=1.8, handlelength=2.4,
    )
    fig.text(.08, .071, "uᵢⱼ,γ = +1  ⇔  i → j", fontsize=11.5)
    fig.text(.43, .071, "u₄₂,z = u₃₁,z = −1", fontsize=11.5, color=colors["z"])
    vec = fig.add_axes((.81, .025, .14, .20))
    for end, label, label_xy in (((2, 0), "a₁", (2.05, -.02)), ((1, 2), "a₂", (1.1, 1.94))):
        vec.add_patch(FancyArrowPatch((0, 0), end, arrowstyle="-|>", mutation_scale=11, linewidth=1.25, color="#566070"))
        vec.text(*label_xy, label, fontsize=10, color="#566070")
    vec.set(xlim=(-.3, 2.8), ylim=(-.4, 2.3), aspect="equal")
    vec.axis("off")
    fig.savefig(output_path, dpi=dpi, facecolor="white")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for old, new in zip(geometry_edges(False), geometry_edges(True)):
        if old[2] == "z":
            assert new == (old[1], old[0], "z")
        else:
            assert old == new, "A non-z arrow changed"
    # Independent correspondence of visual arrows to the translated gauge
    # bond list: omit only endpoints outside the six visible sites per chain.
    visible = {position((m, n)) for m in range(1, 5) for n in range(3)}
    for reversed_z, gauge in ((False, ORIGINAL_GAUGE), (True, Z_REVERSED_GAUGE)):
        reconstructed = set()
        for n in range(-1, 4):
            for bond in CELL_BONDS:
                a, b = position((bond.start, n)), position((bond.end, n + bond.end_cell_shift))
                if a in visible and b in visible:
                    reconstructed.add((a, b, bond.kind) if gauge[bond.key] == 1 else (b, a, bond.kind))
        assert reconstructed == set(geometry_edges(reversed_z)), "Gauge/drawing mismatch"
    gauge_candidates = list(product((-1, 1), repeat=4))
    solutions = [
        signs for signs in gauge_candidates
        if all(signs[b.start-1] * ORIGINAL_GAUGE[b.key] * signs[b.end-1] == Z_REVERSED_GAUGE[b.key]
               for b in CELL_BONDS)
    ]
    result = {
        "ORIGINAL": audit_fluxes(ORIGINAL_GAUGE), "Z_REVERSED": audit_fluxes(Z_REVERSED_GAUGE),
        "local_gauges_tested": len(gauge_candidates), "local_gauge_solutions": solutions,
        "non_z_geometry_unchanged": True, "drawing_matches_real_bond_graph": True,
    }
    assert len(gauge_candidates) == 16
    assert not solutions
    output = args.output_dir / "qsl_lattice_z_reversed.png"
    draw_lattice(output)
    json_path = args.output_dir / "geometry_flux_audit.json"
    json_path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    print(f"Saved {output} at 300 dpi")


if __name__ == "__main__":
    main()
