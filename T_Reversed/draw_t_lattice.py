"""Draw and audit the PDF geometry before and after reversing only t arrows.

The production graph imports the scientific core's bond list. A separate,
literal transcription of PDF Eq.(11) and an explicit drawing-edge rule audit
that graph. Neither verification obtains its expected result from the drawing.
"""

from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import asdict, dataclass
from itertools import product
import json
import math
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/qsl_t_lattice_mpl")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch

from quadrupolar_spin_liquid_t_reversed import (
    GEOMETRY,
    INTERNAL_A1,
    INTERNAL_A2,
    ORIGINAL_GAUGE,
    PLAQUETTES,
    T_REVERSED_GAUGE,
    cell_shift,
    metadata,
    real_space_bonds,
    require,
)


# Independent PDF pp.3-4, Eq.(11) transcription: i,j,gamma,key,da1,da2,shift.
# Keep this redundant evidence explicit: changing the core's geometry must not
# silently change the expected source transcription used by this auditor.
PAPER_BONDS = (
    (1, 2, "x", "u12x", +1, 0, 0),
    (3, 4, "x", "u34x", -1, 0, 0),
    (2, 3, "y", "u23y", 0, +1, 0),
    (4, 1, "y", "u41y", 0, -1, 0),
    (4, 2, "z", "u42z", +1, -1, 0),
    (3, 1, "z", "u31z", +1, -1, +1),
    (4, 3, "w", "u43w", -1, 0, -1),
    (1, 2, "w", "u12w", -1, 0, -1),
    (1, 3, "t", "u13t", +1, +1, 0),
    (2, 4, "t", "u24t", +1, +1, +1),
)

PAPER_GAUGE = {
    "u12x": +1, "u34x": +1, "u23y": +1, "u41y": -1,
    "u42z": +1, "u31z": +1, "u43w": +1, "u12w": -1,
    "u13t": -1, "u24t": -1,
}

# Nodes are (cell, sublattice). This transcription records every CCW path,
# separately from the core's implementation; numerical fluxes are not stored.
PAPER_PATHS = {
    "W124": (((0, 1), (0, 2), (0, 4), (0, 1)), "xzy"),
    "W234": (((0, 2), (0, 3), (0, 4), (0, 2)), "yxz"),
    "W123": (((0, 1), (0, 2), (0, 3), (0, 1)), "xyt"),
    "W134": (((0, 1), (0, 3), (0, 4), (0, 1)), "txy"),
    "W213": (((-1, 2), (0, 1), (-1, 3), (-1, 2)), "wzy"),
    "W143": (((0, 1), (0, 4), (-1, 3), (0, 1)), "ywz"),
    "W214": (((-1, 2), (0, 1), (0, 4), (-1, 2)), "wyt"),
    "W243": (((-1, 2), (0, 4), (-1, 3), (-1, 2)), "twy"),
}

COLORS = {"x": "#2474B5", "y": "#CB3335", "z": "#689F2C", "t": "#7C429B", "w": "#DE901C"}


@dataclass(frozen=True)
class Edge:
    start: tuple[float, float]
    end: tuple[float, float]
    kind: str
    key: str
    stored_start_cell: int


def position(site: tuple[int, int]) -> tuple[float, float]:
    cell, sublattice = site
    a1, a2 = INTERNAL_A1[sublattice - 1], INTERNAL_A2[sublattice - 1]
    return float(4 * cell + 2 * a1 + a2), float(2 * a2)


def geometry_edges(gauge: dict[str, int]) -> tuple[Edge, ...]:
    """Visible arrows from the production real bond graph, in stable order."""
    visible = {position((n, m)) for n in range(3) for m in range(1, 5)}
    result = []
    for n in range(-1, 4):
        for spec, bond in zip(GEOMETRY, real_space_bonds(gauge), strict=True):
            a = position((n, bond.start))
            b = position((n + cell_shift(bond), bond.end))
            if a in visible and b in visible:
                start, end = (a, b) if bond.u == +1 else (b, a)
                result.append(Edge(start, end, bond.gamma, spec[3], n))
    return tuple(result)


def paper_figure_edges(reverse_t: bool) -> set[tuple]:
    """Independent arrow rules read from the figure on PDF page 4."""
    lower = {i: (float(2 * i), 0.) for i in range(6)}
    upper = {i: (float(2 * i + 1), 2.) for i in range(6)}
    result = []
    for i in range(5):
        kind = "x" if i % 2 == 0 else "w"
        pair = ((lower[i], lower[i + 1], kind), (upper[i + 1], upper[i], kind))
        result.extend(pair)
    result.extend((lower[i], upper[i], "y") for i in range(6))
    result.extend((upper[i], lower[i + 1], "z") for i in range(5))
    result.extend((lower[i], upper[i + 1], "t") if reverse_t else
                  (upper[i + 1], lower[i], "t") for i in range(5))
    return set(result)


def audit_fluxes(gauge: dict[str, int]) -> dict:
    """Evaluate each step by lookup in an independent translated edge graph."""
    directed = {}
    for n in range(-2, 3):
        for i, j, kind, key, d1, d2, shift in PAPER_BONDS:
            a, b = (n, i), (n + shift, j)
            require((a, b, kind) not in directed and (b, a, kind) not in directed,
                    f"Duplicate PDF graph edge: {a}, {b}, {kind}")
            directed[a, b, kind] = (gauge[key], key)
            directed[b, a, kind] = (-gauge[key], "-" + key)
    result = {}
    for name, (path, kinds) in PAPER_PATHS.items():
        require(path[0] == path[-1], f"Open plaquette: {name}")
        xy = [position(site) for site in path]
        area = .5 * sum(a[0] * b[1] - b[0] * a[1] for a, b in zip(xy[:-1], xy[1:]))
        require(area > 0, f"Clockwise path: {name}")
        factors, expressions = zip(*(directed[a, b, kind] for a, b, kind in zip(path[:-1], path[1:], kinds)))
        require(all(value in (-1, +1) for value in factors), f"Nonbinary flux factor: {name}")
        result[name] = {
            "path": path,
            "kinds": kinds,
            "factors": factors,
            "product": " ".join(f"({value})" for value in expressions),
            "flux": math.prod(factors),
            "signed_area": area,
        }
    return result


def audit_geometry() -> dict:
    initial_provenance = metadata()
    require(tuple(INTERNAL_A1) == (0, 1, 1, 0), "Core internal a1 positions differ from PDF")
    require(tuple(INTERNAL_A2) == (0, 0, 1, 1), "Core internal a2 positions differ from PDF")
    require(ORIGINAL_GAUGE == PAPER_GAUGE, "Core OLD differs from PDF Eq.(15)")
    require(len(GEOMETRY) == len(PAPER_BONDS) and set(GEOMETRY) == {row[:6] for row in PAPER_BONDS},
            "Core bond list differs from PDF Eq.(11)")
    require(PLAQUETTES == PAPER_PATHS, "Core paths differ from audited CCW paths")
    expected_shifts = {row[3]: row[6] for row in PAPER_BONDS}
    for spec, bond in zip(GEOMETRY, real_space_bonds(ORIGINAL_GAUGE), strict=True):
        require(cell_shift(bond) == expected_shifts[spec[3]], f"Incorrect cell shift: {spec[3]}")
    changed_keys = {key for key in ORIGINAL_GAUGE if ORIGINAL_GAUGE[key] != T_REVERSED_GAUGE[key]}
    require(changed_keys == {"u13t", "u24t"}, f"Gauge contamination: {changed_keys}")
    for key in ("u34x", "u43w", "u42z", "u31z"):
        require(T_REVERSED_GAUGE[key] == +1, f"Forbidden gauge change: {key}")
    old_edges, new_edges = geometry_edges(ORIGINAL_GAUGE), geometry_edges(T_REVERSED_GAUGE)
    changed_edges = []
    for old, new in zip(old_edges, new_edges, strict=True):
        require((old.kind, old.key, old.stored_start_cell) == (new.kind, new.key, new.stored_start_cell),
                "Drawing changed bond family, identity or cell")
        if old.kind == "t":
            require((new.start, new.end) == (old.end, old.start), "A t arrow was not reversed exactly")
            changed_edges.append({"old": asdict(old), "t_reversed": asdict(new)})
        else:
            require(old == new, "A non-t arrow changed")
    for reverse, edges in ((False, old_edges), (True, new_edges)):
        require({(edge.start, edge.end, edge.kind) for edge in edges} == paper_figure_edges(reverse),
                "Core-derived drawing differs from independent PDF figure rules")
    candidates = list(product((-1, +1), repeat=4))
    solutions = [
        signs for signs in candidates
        if all(signs[i - 1] * ORIGINAL_GAUGE[key] * signs[j - 1] == T_REVERSED_GAUGE[key]
               for i, j, kind, key, d1, d2, shift in PAPER_BONDS)
    ]
    require(len(candidates) == 16, "Diagonal gauge search is incomplete")
    require(not solutions, "Unexpected diagonal gauge equivalence")
    old_flux, new_flux = audit_fluxes(ORIGINAL_GAUGE), audit_fluxes(T_REVERSED_GAUGE)
    for name, (_, kinds) in PAPER_PATHS.items():
        require(new_flux[name]["flux"] == (-1 if "t" in kinds else +1) * old_flux[name]["flux"],
                f"Derived flux violates t-edge parity: {name}")
    final_provenance = metadata()
    for kind in ("script_sha256", "source_sha256", "local_source_sha256"):
        require(initial_provenance[kind] == final_provenance[kind], f"Sources changed during geometry audit: {kind}")
    return {
        "metadata": final_provenance,
        "source_geometry": "PDF pages 3-4, Eq.(11), original arrows page 4; gauge Eq.(15), page 5",
        "changed_keys": sorted(changed_keys),
        "visible_edge_count": len(old_edges),
        "visible_edge_counts_by_kind": dict(Counter(edge.kind for edge in old_edges)),
        "reversed_t_edge_count": len(changed_edges),
        "non_t_edges_unchanged": True,
        "positions_and_vectors_unchanged": True,
        "drawing_matches_independent_PDF_figure_rules": True,
        "core_bonds_match_independent_Eq11_transcription": True,
        "changed_edges": changed_edges,
        "OLD_EDGES": [asdict(edge) for edge in old_edges],
        "T_REVERSED_EDGES": [asdict(edge) for edge in new_edges],
        "OLD": old_flux,
        "T_REVERSED": new_flux,
        "local_gauges_tested": len(candidates),
        "local_gauge_solutions": solutions,
    }


def draw_lattice(output_path: Path, dpi: int = 300) -> None:
    plt.rcParams.update({"font.family": "DejaVu Sans", "mathtext.fontset": "dejavusans"})
    fig, axes = plt.subplots(2, 1, figsize=(11.6, 7.1), dpi=dpi, facecolor="white")
    fig.subplots_adjust(left=.055, right=.97, top=.85, bottom=.18, hspace=.52)
    for ax, gauge, title in (
        (axes[0], ORIGINAL_GAUGE, "OLD / ZERO-FLUX"),
        (axes[1], T_REVERSED_GAUGE, "T-REVERSED"),
    ):
        for edge in geometry_edges(gauge):
            ax.add_patch(FancyArrowPatch(
                edge.start, edge.end, arrowstyle="-|>", mutation_scale=16,
                linewidth=1.75 if edge.kind != "t" else 2.25,
                color=COLORS[edge.kind], shrinkA=5, shrinkB=5,
                zorder=2 if edge.kind == "t" else 3,
            ))
        for n in range(3):
            for m in range(1, 5):
                x, y = position((n, m))
                ax.plot(x, y, "o", color="#1D2837", markersize=6.1, zorder=6)
        for m, offset, va in ((1, -10, "top"), (2, -10, "top"), (3, 10, "bottom"), (4, 10, "bottom")):
            ax.annotate(str(m), position((1, m)), xytext=(0, offset), textcoords="offset points",
                        ha="center", va=va, fontsize=13, fontweight="bold")
        ax.text(-.62, 1., "···", fontsize=20, ha="center", va="center")
        ax.text(11.6, 1., "···", fontsize=20, ha="center", va="center")
        ax.set(xlim=(-1., 12.), ylim=(-.62, 2.62), aspect="equal")
        ax.axis("off")
        ax.set_title(title, loc="left", fontsize=12.5, pad=13, color="#29394D", weight="bold")
    fig.suptitle("Inversão exclusiva das ligações t", fontsize=17, y=.972, weight="bold", color="#1D2837")
    fig.text(.5, .913, "Mesmas posições, sub-redes, vetores e demais ligações nos dois painéis",
             ha="center", fontsize=10.8, color="#566070")
    fig.legend(handles=[Line2D([0], [0], color=COLORS[kind], linewidth=2.6, label=kind)
                        for kind in ("x", "y", "z", "t", "w")],
               loc="lower left", bbox_to_anchor=(.075, .081), ncol=5, frameon=False,
               fontsize=12, columnspacing=1.8, handlelength=2.4)
    fig.text(.08, .058, "i → j  ⇔  uᵢⱼ,γ = +1", fontsize=11.5)
    fig.text(.39, .058, "u₁₃,t: -1 → +1,    u₂₄,t: -1 → +1", fontsize=11.5, color=COLORS["t"])
    vector_ax = fig.add_axes((.825, .035, .13, .125))
    for end, label, label_xy in (((2, 0), "a₁", (2.12, -.10)), ((1, 2), "a₂", (1.1, 1.94))):
        vector_ax.add_patch(FancyArrowPatch((0, 0), end, arrowstyle="-|>", mutation_scale=11,
                                          linewidth=1.25, color="#566070"))
        vector_ax.text(*label_xy, label, fontsize=10, color="#566070")
    vector_ax.set(xlim=(-.25, 2.8), ylim=(-.35, 2.35), aspect="equal")
    vector_ax.axis("off")
    fig.savefig(output_path, dpi=dpi, facecolor="white")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    result = audit_geometry()
    output = args.output_dir / "qsl_lattice_t_reversed.png"
    draw_lattice(output)
    current = metadata()
    for kind in ("script_sha256", "source_sha256", "local_source_sha256"):
        require(result["metadata"][kind] == current[kind], f"Sources changed while rendering geometry: {kind}")
    audit_file = args.output_dir / "geometry_t_audit.json"
    audit_file.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"visible_edges": result["visible_edge_count"],
                      "reversed_t_edges": result["reversed_t_edge_count"],
                      "gauge_solutions": result["local_gauge_solutions"],
                      "output": str(output), "audit": str(audit_file)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
