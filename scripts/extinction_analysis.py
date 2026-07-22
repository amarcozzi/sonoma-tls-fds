#!/usr/bin/env python
"""Throwaway diagnostic: why don't c4_p25 / c8_p25 sustain spread?

The three runs are IDENTICAL except for the tree bulk-density file (verified
by diff), so the answer must lie in canopy structure. This script quantifies,
per run:
  1. Initial canopy structure from PARTICLE BULK DENSITY (prt5): mass per
     metre of fetch and vertical profile, within the |y| <= 10 m band.
  2. Sub-canopy wind from the U-VELOCITY PBY=0 slice, averaged over the
     pre-fire spin-up window: how much does each canopy slow the wind?
  3. Canopy heating state from PARTICLE TEMPERATURE: mass fraction actively
     drying (~100 C plateau) vs pyrolyzing (>250 C) over time.
  4. Canopy consumption: initial vs final canopy mass along x.
  5. Total HRR and water-vapor MLR (drying heat sink) from CHID_hrr.csv.

Particle files are read with a minimal seek-based prt5 parser (fdsreader
would load all 41 GB of c4_p25); only the 16 centre-row meshes are read and
frames are kept every DT_KEEP seconds.
"""

import os
import struct

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

import fdsreader

SIM_ROOT = "/Volumes/T7 Shield/sonoma-tls-fds/simulations"
OUT_DIR = "/Volumes/T7 Shield/sonoma-tls-fds/postprocess/extinction_analysis"
RUNS = ["c4_p10", "c4_p25", "c8_p25"]
COLORS = {"c4_p10": "tab:blue", "c4_p25": "tab:orange", "c8_p25": "tab:green"}

Y_BAND = 10.0          # analysis band, matches ROS band / ROI half-width
DT_KEEP = 2.0          # particle frame subsampling (s)
DX = DY = 0.25         # horizontal grid spacing
X_BINS = np.arange(-80.0, 40.0 + 1e-9, 2.0)   # 2 m bins along spread
Z_BINS = np.arange(0.0, 30.0 + 1e-9, 1.0)     # vertical profile bins
# centre-row meshes (j=3,4 of the 8x8 MULT layout; mesh number = j*8+i+1)
CENTER_MESH_NUMS = [j * 8 + i + 1 for j in (3, 4) for i in range(8)]
PREFIRE_WINDOW = (15.0, 29.0)  # ignition ramp starts at t=30 s


# ---------------------------------------------------------------------------
# Minimal prt5 reader (FDS particle format: Fortran sequential records)
# ---------------------------------------------------------------------------

def _rec(f):
    """Read one Fortran record, return raw payload bytes."""
    head = f.read(4)
    if len(head) < 4:
        return None
    (n,) = struct.unpack("<i", head)
    payload = f.read(n)
    f.seek(4, 1)  # trailing length
    return payload


def _skip_rec(f):
    head = f.read(4)
    if len(head) < 4:
        return False
    (n,) = struct.unpack("<i", head)
    f.seek(n + 4, 1)
    return True


def read_prt5(path, dt_keep=0.0):
    """Yield (time, positions (n,3), quantities dict name->(n,) f32) frames.

    Frames closer than dt_keep to the previously yielded frame are seeked
    past without decoding (headers only), which is what makes scanning tens
    of GB feasible.
    """
    with open(path, "rb") as f:
        _rec(f)  # ONE (endianness check)
        _rec(f)  # FDS version
        n_classes = struct.unpack("<i", _rec(f)[:4])[0]
        class_quants = []
        for _ in range(n_classes):
            n_q = struct.unpack("<ii", _rec(f))[0]
            names = []
            for _ in range(n_q):
                names.append(_rec(f).decode().strip())
                _rec(f)  # unit
            class_quants.append(names)

        t_last = -1e30
        while True:
            payload = _rec(f)
            if payload is None:
                return
            t = struct.unpack("<f", payload[:4])[0]
            keep = t - t_last >= dt_keep
            for names in class_quants:
                n = struct.unpack("<i", _rec(f)[:4])[0]
                # FDS writes the position/tag/data records even when n == 0,
                # so every record must be consumed on every frame
                if keep and n > 0:
                    xyz = np.frombuffer(_rec(f), dtype="<f4").reshape(
                        (n, 3), order="F")
                    _skip_rec(f)  # tags
                    if names:
                        q = np.frombuffer(_rec(f), dtype="<f4").reshape(
                            (n, len(names)), order="F")
                        quants = {nm: q[:, i] for i, nm in enumerate(names)}
                    else:
                        quants = {}
                else:
                    _skip_rec(f)  # positions
                    _skip_rec(f)  # tags
                    if names:
                        _skip_rec(f)
                    xyz, quants = np.empty((0, 3)), {}
            if keep:
                t_last = t
                yield t, xyz, quants


# ---------------------------------------------------------------------------
# Per-run analysis
# ---------------------------------------------------------------------------

def analyze_particles(run, z_nodes):
    """Aggregate canopy particle mass/temperature stats over time.

    Particle mass is approximated as bulk_density * cell volume, with the
    cell height dz taken from the (stretched) grid at the particle height.
    """
    dz_cells = np.diff(z_nodes)

    def cell_dz(z):
        idx = np.clip(np.searchsorted(z_nodes, z) - 1, 0, len(dz_cells) - 1)
        return dz_cells[idx]

    frames = {}  # t -> accumulators
    sim_dir = os.path.join(SIM_ROOT, run)
    for mesh_num in CENTER_MESH_NUMS:
        path = os.path.join(sim_dir, f"out_{run}_cat_{mesh_num}.prt5")
        for t, xyz, quants in read_prt5(path, dt_keep=DT_KEEP):
            key = round(t, 1)
            acc = frames.setdefault(key, {
                "mass": 0.0, "drying": 0.0, "pyro": 0.0,
                "mass_x": np.zeros(X_BINS.size - 1),
                "mass_z": np.zeros(Z_BINS.size - 1),
            })
            if xyz.shape[0] == 0:
                continue
            sel = np.abs(xyz[:, 1]) <= Y_BAND
            if not sel.any():
                continue
            xyz = xyz[sel]
            bd = quants["PARTICLE BULK DENSITY"][sel].astype(np.float64)
            temp = quants["PARTICLE TEMPERATURE"][sel]
            mass = bd * DX * DY * cell_dz(xyz[:, 2])
            acc["mass"] += mass.sum()
            acc["drying"] += mass[(temp >= 90) & (temp <= 110)].sum()
            acc["pyro"] += mass[temp > 250].sum()
            acc["mass_x"] += np.histogram(xyz[:, 0], X_BINS, weights=mass)[0]
            acc["mass_z"] += np.histogram(xyz[:, 2], Z_BINS, weights=mass)[0]

    times = np.array(sorted(frames))
    out = {
        "t": times,
        "mass": np.array([frames[t]["mass"] for t in times]),
        "drying": np.array([frames[t]["drying"] for t in times]),
        "pyro": np.array([frames[t]["pyro"] for t in times]),
        "mass_x0": frames[times[0]]["mass_x"],
        "mass_x_end": frames[times[-1]]["mass_x"],
        "mass_z0": frames[times[0]]["mass_z"],
    }
    return out


def subcanopy_wind(sim):
    """Pre-fire mean U at z<=3 m vs x from the PBY=0 U-VELOCITY slice."""
    slc = next(s for s in sim.slices
               if s.quantity.name == "U-VELOCITY" and s.orientation == 2
               and abs(s.extent.y_start) < 0.5)
    prof_x, prof_u = [], []
    for sub in slc._subslices.values():
        mesh = sub.mesh
        xn, zn = mesh.coordinates["x"], mesh.coordinates["z"]
        xc = 0.5 * (xn[:-1] + xn[1:])
        zc = 0.5 * (zn[:-1] + zn[1:])
        data = np.squeeze(sub.data)  # (n_t, n_x, n_z)
        tsel = (slc.times >= PREFIRE_WINDOW[0]) & (slc.times <= PREFIRE_WINDOW[1])
        zsel = zc <= 3.0
        dz = np.diff(zn)[zsel]
        u = np.average(data[tsel][:, :, zsel].mean(axis=0), axis=1, weights=dz)
        prof_x.append(xc)
        prof_u.append(u)
    slc.clear_cache()
    order = np.argsort([x[0] for x in prof_x])
    return (np.concatenate([prof_x[i] for i in order]),
            np.concatenate([prof_u[i] for i in order]))


def read_hrr(run):
    d = np.genfromtxt(os.path.join(SIM_ROOT, run, f"out_{run}_cat_hrr.csv"),
                      delimiter=",", skip_header=2)
    return d[:, 0], d[:, 1], d[:, 14]  # t, HRR, MLR_WATER VAPOR


# ---------------------------------------------------------------------------

def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    xc_bins = 0.5 * (X_BINS[:-1] + X_BINS[1:])
    zc_bins = 0.5 * (Z_BINS[:-1] + Z_BINS[1:])
    bin_w = np.diff(X_BINS)[0]

    results = {}
    for run in RUNS:
        print(f"--- {run}")
        sim = fdsreader.Simulation(os.path.join(SIM_ROOT, run))
        z_nodes = sim.meshes[0].coordinates["z"]
        res = analyze_particles(run, z_nodes)
        res["wind_x"], res["wind_u"] = subcanopy_wind(sim)
        res["hrr_t"], res["hrr"], res["mlr_w"] = read_hrr(run)
        results[run] = res
        print(f"    initial canopy mass in band: {res['mass'][0]/1e3:.1f} t, "
              f"final: {res['mass'][-1]/1e3:.1f} t")

    # 1. canopy structure
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    for run in RUNS:
        r = results[run]
        axes[0].plot(xc_bins, r["mass_x0"] / bin_w / (2 * Y_BAND),
                     color=COLORS[run], label=run)
        axes[1].plot(r["mass_z0"] / np.diff(Z_BINS)[0], zc_bins,
                     color=COLORS[run], label=run)
    axes[0].axvspan(-10, 10, color="red", alpha=0.06)
    for xv in (-10, 10):
        axes[0].axvline(xv, color="red", ls="--", lw=1)
    axes[0].set_xlabel("x (m)")
    axes[0].set_ylabel("canopy load (kg/m²)")
    axes[0].set_title("Initial canopy load along spread direction (|y| ≤ 10 m)")
    axes[1].set_xlabel("canopy mass per m height (kg/m)")
    axes[1].set_ylabel("z (m)")
    axes[1].set_title("Initial vertical canopy profile")
    for ax in axes:
        ax.legend()
        ax.grid(alpha=0.3)
    fig.savefig(os.path.join(OUT_DIR, "canopy_structure.png"), dpi=200,
                bbox_inches="tight")
    plt.close(fig)

    # 2. sub-canopy wind
    fig, ax = plt.subplots(figsize=(8, 4.5))
    for run in RUNS:
        r = results[run]
        ax.plot(r["wind_x"], r["wind_u"], color=COLORS[run], label=run)
    ax.axhline(0, color="k", lw=0.5)
    for xv in (-10, 10):
        ax.axvline(xv, color="red", ls="--", lw=1)
    ax.axvline(-60, color="gray", ls=":", lw=1)
    ax.annotate("ignition line", (-60, 0.95), xycoords=("data", "axes fraction"),
                rotation=90, va="top", fontsize=8, color="gray")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("mean U below 3 m (m/s)")
    ax.set_title(f"Pre-fire sub-canopy wind (t = {PREFIRE_WINDOW[0]:.0f}–"
                 f"{PREFIRE_WINDOW[1]:.0f} s, y = 0 slice)")
    ax.legend()
    ax.grid(alpha=0.3)
    fig.savefig(os.path.join(OUT_DIR, "subcanopy_wind.png"), dpi=200,
                bbox_inches="tight")
    plt.close(fig)

    # 3. canopy heating state over time
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    for run in RUNS:
        r = results[run]
        axes[0].plot(r["t"], 100 * r["drying"] / r["mass"][0],
                     color=COLORS[run], label=run)
        axes[1].plot(r["t"], 100 * r["pyro"] / r["mass"][0],
                     color=COLORS[run], label=run)
    axes[0].set_ylabel("% of initial canopy mass at 90–110 °C (drying)")
    axes[1].set_ylabel("% of initial canopy mass > 250 °C (pyrolyzing)")
    for ax in axes:
        ax.set_xlabel("time (s)")
        ax.legend()
        ax.grid(alpha=0.3)
    fig.suptitle("Canopy heating state (|y| ≤ 10 m band)")
    fig.savefig(os.path.join(OUT_DIR, "canopy_heating.png"), dpi=200,
                bbox_inches="tight")
    plt.close(fig)

    # 4. canopy mass: time series + consumption along x
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    for run in RUNS:
        r = results[run]
        axes[0].plot(r["t"], r["mass"] / 1e3, color=COLORS[run], label=run)
        consumed = np.where(r["mass_x0"] > 0,
                            100 * (1 - r["mass_x_end"] / np.maximum(r["mass_x0"], 1e-9)),
                            np.nan)
        axes[1].plot(xc_bins, consumed, color=COLORS[run], label=run)
    axes[0].set_xlabel("time (s)")
    axes[0].set_ylabel("canopy mass in band (t)")
    axes[0].set_title("Canopy mass remaining")
    for xv in (-10, 10):
        axes[1].axvline(xv, color="red", ls="--", lw=1)
    axes[1].set_xlabel("x (m)")
    axes[1].set_ylabel("canopy mass consumed (%)")
    axes[1].set_title("Canopy consumption at end of run")
    for ax in axes:
        ax.legend()
        ax.grid(alpha=0.3)
    fig.savefig(os.path.join(OUT_DIR, "canopy_consumption.png"), dpi=200,
                bbox_inches="tight")
    plt.close(fig)

    # 5. HRR + water MLR
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    for run in RUNS:
        r = results[run]
        axes[0].plot(r["hrr_t"], r["hrr"] / 1e6, color=COLORS[run], label=run)
        axes[1].plot(r["hrr_t"], r["mlr_w"], color=COLORS[run], label=run)
    axes[0].set_ylabel("total HRR (GW)")
    axes[1].set_ylabel("water vapor MLR (kg/s)")
    for ax in axes:
        ax.set_xlabel("time (s)")
        ax.set_xlim(0, 130)
        ax.legend()
        ax.grid(alpha=0.3)
    fig.suptitle("Fire intensity vs moisture evaporation (drying heat sink)")
    fig.savefig(os.path.join(OUT_DIR, "hrr_water.png"), dpi=200,
                bbox_inches="tight")
    plt.close(fig)

    # summary numbers
    with open(os.path.join(OUT_DIR, "summary.txt"), "w") as f:
        for run in RUNS:
            r = results[run]
            roi = (xc_bins >= -10) & (xc_bins <= 10)
            fetch = (xc_bins >= -60) & (xc_bins <= 10)
            f.write(
                f"{run}: canopy mass in band {r['mass'][0]/1e3:.1f} t "
                f"(ROI {r['mass_x0'][roi].sum()/1e3:.1f} t); "
                f"mean canopy load over fetch "
                f"{r['mass_x0'][fetch].sum()/(fetch.sum()*bin_w*2*Y_BAND):.2f} kg/m2; "
                f"mean sub-canopy U over fetch "
                f"{np.interp(xc_bins[fetch], r['wind_x'], r['wind_u']).mean():.2f} m/s; "
                f"peak drying {100*r['drying'].max()/r['mass'][0]:.1f} %; "
                f"canopy consumed {100*(1-r['mass'][-1]/r['mass'][0]):.1f} %\n")
    print(open(os.path.join(OUT_DIR, "summary.txt")).read())


if __name__ == "__main__":
    main()
