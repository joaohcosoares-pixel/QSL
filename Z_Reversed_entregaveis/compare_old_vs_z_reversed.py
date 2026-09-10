#!/usr/bin/env python3
"""Calculate physical bands, gaps and energies from the audited bond-list core.

Run: python3 compare_old_vs_z_reversed.py
Optional: --couplings Kx Ky Kz Kt Kw --output-dir PATH
Band tracking is used only for the dispersion figure; energy uses eigvalsh.
gap_old / gap_z are FULL zero-energy band gaps, not distances from zero.
"""
from __future__ import annotations

import argparse
import csv
from itertools import permutations
import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/qsl_z_comparison_mpl")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize_scalar
from scipy.special import ellipe

import quadrupolar_spin_liquid_z_reversed as model
from draw_z_lattice import draw_lattice

PERMUTATIONS = tuple(permutations(range(4)))
CONVERGENCE_NQ = (1001, 2001, 4001, 8001, 16001)
GAP_TOL = 1e-8
ENERGY_TOL = 1e-10


def tracked_bands(q_values, couplings, gauge):
    """Maximum squared-overlap assignment, with degenerate subspace transport.

    This function returns a visualization only. It is never used for occupation,
    gaps, energy integration or favored-sector decisions.
    """
    result = np.empty((len(q_values), 4))
    values, previous = np.linalg.eigh(model.h_bloch(q_values[0],couplings,gauge))
    result[0] = values
    for k, q in enumerate(q_values[1:], start=1):
        values, current = np.linalg.eigh(model.h_bloch(q,couplings,gauge))
        overlap = np.abs(previous.conj().T@current)**2
        best = max(PERMUTATIONS, key=lambda p: sum(overlap[m,p[m]] for m in range(4)))
        indices = np.array(best)
        values, current = values[indices], current[:,indices]
        remaining = set(range(4))
        threshold = 1e-11*max(1.,float(np.max(abs(values))))
        while remaining:
            first = min(remaining)
            group = sorted(i for i in remaining if abs(values[i]-values[first]) < threshold)
            remaining.difference_update(group)
            if len(group)>1:
                cross = current[:,group].conj().T@previous[:,group]
                left, singular_values, right = np.linalg.svd(cross)
                current[:,group] = current[:,group]@(left@right)
        result[k] = values
        previous = current
    return result


def zero_energy_band_gap(couplings,gauge):
    """Continuous-BZ numerical search. Full gap=2 min|E| for proven +/- pairing.

    Coarse mesh brackets all visible minima; bounded refinement improves their
    location. Isotropic globality is additionally established analytically in
    the report. The generic routine does not claim a certified global bound.
    """
    grid = np.linspace(-np.pi/2,np.pi/2,4001)
    values = model.bands_sorted(grid,couplings,gauge)
    distance = np.min(abs(values),axis=1)
    def objective(q):
        return float(np.min(abs(np.linalg.eigvalsh(model.h_bloch(q,couplings,gauge)))))
    candidates = [(float(distance[0]),float(grid[0])),(float(distance[-1]),float(grid[-1]))]
    for i in range(1,len(grid)-1):
        if distance[i] <= distance[i-1] and distance[i] <= distance[i+1]:
            candidates.append((float(distance[i]),float(grid[i])))
            if distance[i-1] == distance[i] == distance[i+1]:
                continue
            optimum = minimize_scalar(objective,bounds=(grid[i-1],grid[i+1]),
                                      method="bounded",options={"xatol":1e-13})
            if not optimum.success:
                raise RuntimeError("Gap minimization failed")
            candidates.append((float(optimum.fun),float(optimum.x)))
    d, qmin = min(candidates)
    gap = 0. if d < GAP_TOL else 2*d
    return {"zero_energy_band_gap": gap,"min_abs_E":d,"q_at_minimum":qmin,
            "classification":"gapless" if gap == 0 else "gapped",
            "definition":"min(E>0)-max(E<0)=2*min(abs(E)) for this paired spectrum"}


def configure_axis(ax):
    ax.set_xlim(-np.pi/2,np.pi/2)
    ax.set_xticks([-np.pi/2,-np.pi/4,0,np.pi/4,np.pi/2])
    ax.set_xticklabels(["−π/2","−π/4","0","π/4","π/2"])
    ax.set_xlabel("q = k·a₁")
    ax.grid(alpha=.18,linewidth=.6)
    ax.spines[["top","right"]].set_visible(False)


def write_csv(path,rows,fieldnames):
    with path.open("w",newline="",encoding="utf-8") as stream:
        writer = csv.DictWriter(stream,fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--couplings",type=float,nargs=5,default=[1.]*5)
    parser.add_argument("--output-dir",type=Path,default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    K = tuple(args.couplings)
    model.coupling_map(K)
    model.assert_only_z_reversed()
    out = args.output_dir
    out.mkdir(parents=True,exist_ok=True)
    sectors = model.GAUGES
    convergence = []
    for nq in CONVERGENCE_NQ:
        q = np.linspace(-np.pi/2,np.pi/2,nq)
        energies = [model.energy_per_site(q,model.bands_sorted(q,K,gauge)) for gauge in sectors.values()]
        convergence.append({"Nq":nq,"epsilon_old":energies[0],"epsilon_z":energies[1],
                            "Delta_epsilon":energies[1]-energies[0]})
    last = convergence[-1]
    stable_decimals = {}
    for name in ("epsilon_old","epsilon_z","Delta_epsilon"):
        stable_decimals[name] = max(d for d in range(13)
                                  if len({round(row[name],d) for row in convergence[-3:]}) == 1)
    gaps = {name:zero_energy_band_gap(K,gauge) for name,gauge in sectors.items()}
    delta = last["Delta_epsilon"]
    favored = "OLD" if delta > ENERGY_TOL else "Z-REVERSED" if delta < -ENERGY_TOL else "TIE"
    row = dict(zip(("Kx","Ky","Kz","Kt","Kw"),K))
    row.update({name:last[name] for name in ("epsilon_old","epsilon_z","Delta_epsilon")})
    row.update({"gap_old":gaps["OLD"]["zero_energy_band_gap"],
                "gap_z":gaps["Z-REVERSED"]["zero_energy_band_gap"],"favored_sector":favored})

    # Independent quadrature samples freshly diagonalized physical matrices.
    adaptive = {}
    bz_errors = {}
    for name,gauge in sectors.items():
        def integrand(q):
            return float(model.occupied_sum(np.linalg.eigvalsh(model.h_bloch(q,K,gauge))))
        integral,error = quad(integrand,-np.pi/2,np.pi/2,epsabs=1e-12,epsrel=1e-12,limit=200)
        adaptive[name] = {"epsilon":integral/(4*np.pi),"estimated_error_per_site":error/(4*np.pi)}
        centered = np.linspace(-np.pi/2,np.pi/2,8001)
        shifted = centered+.371
        extended = np.linspace(-np.pi,np.pi,16001)
        eps = model.energy_per_site(centered,model.bands_sorted(centered,K,gauge))
        shift_eps = model.energy_per_site(shifted,model.bands_sorted(shifted,K,gauge))
        ext_eps = model.energy_per_site(extended,model.bands_sorted(extended,K,gauge),extended=True)
        bz_errors[name] = {"shifted_BZ":abs(shift_eps-eps),"extended_BZ":abs(ext_eps-eps),
                           "adaptive_vs_mesh":abs(adaptive[name]["epsilon"]-eps)}
        assert max(bz_errors[name].values()) < 2e-8

    energy_uncertainty = max(ENERGY_TOL,
        sum(abs(last[key]-adaptive[name]["epsilon"])+adaptive[name]["estimated_error_per_site"]
            for name,key in (("OLD","epsilon_old"),("Z-REVERSED","epsilon_z"))),
        sum(abs(convergence[-1][key]-convergence[-2][key]) for key in ("epsilon_old","epsilon_z")))
    favored = ("OLD" if delta > energy_uncertainty else
               "Z-REVERSED" if delta < -energy_uncertainty else
               "TIE" if K == (1.,)*5 else "UNRESOLVED")
    row["favored_sector"] = favored

    qplot = np.linspace(-np.pi/2,np.pi/2,2001)
    sorted_bands = {name:model.bands_sorted(qplot,K,gauge) for name,gauge in sectors.items()}
    tracked = {name:tracked_bands(qplot,K,gauge) for name,gauge in sectors.items()}
    tracking_errors = {name:float(np.max(abs(np.sort(tracked[name],axis=1)-sorted_bands[name]))) for name in sectors}
    assert max(tracking_errors.values()) < 1e-11
    F = {name:model.occupied_sum(values) for name,values in sorted_bands.items()}
    isotropic = None
    if K == (1.,)*5:
        positive_large = np.sqrt(5+4*np.cos(qplot))
        positive_small = np.sqrt(5-4*np.cos(qplot))
        closed_old = np.stack([-positive_large,-positive_small,positive_small,positive_large],axis=1)
        positive_z = np.sqrt(1+8*np.sin(qplot)**2)
        closed_z = np.stack([-positive_z,-positive_z,positive_z,positive_z],axis=1)
        closed_errors = {"OLD":float(np.max(abs(sorted_bands["OLD"]-closed_old))),
                         "Z-REVERSED":float(np.max(abs(sorted_bands["Z-REVERSED"]-closed_z)))}
        assert max(closed_errors.values()) < 1e-11
        exact_energy = float(-3/np.pi*ellipe(8/9))
        assert abs(last["epsilon_old"]-exact_energy) < 1e-12
        assert abs(last["epsilon_z"]-exact_energy) < 1e-12
        assert favored == "TIE"
        isotropic = {"elliptic_energy":exact_energy,"closed_band_max_errors":closed_errors,
                     "exact_equality_proven":True,"zero_crossings":False,
                     "negative_band_counts":{name:sorted(set(np.sum(E<0,axis=1).tolist())) for name,E in sorted_bands.items()}}
    special_points = {name:{label:np.linalg.eigvalsh(model.h_bloch(q,K,gauge)).tolist()
                           for label,q in (("-pi/2",-np.pi/2),("-pi/4",-np.pi/4),("0",0.),("pi/4",np.pi/4),("pi/2",np.pi/2))}
                      for name,gauge in sectors.items()}

    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":11,"axes.titlesize":13})
    fig,axes = plt.subplots(1,2,figsize=(11.5,4.3),sharey=True,layout="constrained")
    colors = ["#24478f","#dd7132","#2d8b69","#8c459a"]
    scale = max(float(np.max(abs(E))) for E in sorted_bands.values())
    for ax,(name,values),title in zip(axes,tracked.items(),("Original gauge (OLD)","z bonds reversed (Z-REVERSED)")):
        for band in range(4):
            ax.plot(qplot,values[:,band],color=colors[band],linewidth=1.65,
                    linestyle="--" if band in (1,3) else "-")
        ax.axhline(0,color="#555555",linewidth=.8)
        ax.set_ylim(-1.12*scale,1.12*scale)
        configure_axis(ax)
        ax.set_title(title)
        if isotropic is not None and name == "Z-REVERSED":
            ax.text(.5,.51,"Each branch has degeneracy 2",transform=ax.transAxes,ha="center",fontsize=9,color="#555555")
    axes[0].set_ylabel("Band energy E")
    fig.suptitle("Kx = Ky = Kz = Kt = Kw = 1" if isotropic else "Couplings (x,y,z,t,w): "+str(K),fontsize=11)
    fig.savefig(out/"dispersion_old_vs_z_reversed.png",dpi=300)
    plt.close(fig)
    fig,ax = plt.subplots(figsize=(8.2,4.5),layout="constrained")
    for name,color in (("OLD","#24478f"),("Z-REVERSED","#be4b31")):
        ax.plot(qplot,F[name],label=name,color=color,linewidth=2.)
    configure_axis(ax)
    ax.set_ylabel("F(q) = sum of negative band energies")
    ax.set_title("Occupied-band sum · OLD versus Z-REVERSED")
    ax.legend(frameon=False)
    fig.savefig(out/"Fq_old_vs_z_reversed.png",dpi=300)
    plt.close(fig)
    draw_lattice(out/"qsl_lattice_z_reversed.png")
    write_csv(out/"ground_state_energy_old_vs_z.csv",[row],list(row))
    write_csv(out/"energy_convergence_old_vs_z.csv",convergence,list(convergence[0]))
    results = {"couplings":dict(zip(("Kx","Ky","Kz","Kt","Kw"),K)),
               "convergence":convergence,"stable_decimals_last_three_meshes_capped_at_12":stable_decimals,
               "gaps":gaps,"summary":row,"adaptive_quadrature":adaptive,"BZ_invariance_errors":bz_errors,
               "tracking_spectrum_errors":tracking_errors,"special_points":special_points,"isotropic":isotropic,
               "energy_tie_tolerance":ENERGY_TOL,"gap_zero_tolerance":GAP_TOL,
               "effective_energy_comparison_tolerance":energy_uncertainty,
               "tracking_used_for_energy":False,"figures_dpi":300}
    (out/"comparison_z_reversed.json").write_text(json.dumps(results,indent=2)+"\n")
    for filename in ("qsl_lattice_z_reversed.png","dispersion_old_vs_z_reversed.png","Fq_old_vs_z_reversed.png","ground_state_energy_old_vs_z.csv"):
        assert (out/filename).stat().st_size > 0
    print(json.dumps(results,indent=2))


if __name__ == "__main__":
    main()
