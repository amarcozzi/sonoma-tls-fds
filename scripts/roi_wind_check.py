#!/usr/bin/env python
"""Is the near-zero ROI sub-canopy wind an equilibrium state or a spin-up
artifact? Track mean u over the ROI (x in [-10,10], z < 3 m) for the FULL
record of each run. In c4_p25 (t_end=320 s, fire dead ~110 s) and c8_p25
(t_end=175 s, fire dead ~105 s), the late fire-free period shows the true
developed wind. Also split the ROI into upwind/downwind halves.
"""

import os
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

import fdsreader

SIM_ROOT = "/Volumes/T7 Shield/sonoma-tls-fds/simulations"
OUT_DIR = "/Volumes/T7 Shield/sonoma-tls-fds/postprocess/wind_checks"
RUNS = ["c4_p10", "c4_p25", "c8_p25"]
COLORS = {"c4_p10": "tab:blue", "c4_p25": "tab:orange", "c8_p25": "tab:green"}
FIRE_DEAD = {"c4_p10": None, "c4_p25": 115.0, "c8_p25": 110.0}

fig, axes = plt.subplots(1, 3, figsize=(14, 4.2), sharey=True)

for k, run in enumerate(RUNS):
    sim = fdsreader.Simulation(os.path.join(SIM_ROOT, run))
    slices = [s for s in sim.slices
              if s.quantity.name == "U-VELOCITY" and s.orientation == 2]
    # average over all 5 PBY slices to reduce single-plane bias
    acc = {}
    times = None
    for slc in slices:
        for sub in slc._subslices.values():
            mesh = sub.mesh
            xn, zn = mesh.coordinates["x"], mesh.coordinates["z"]
            xc = 0.5 * (xn[:-1] + xn[1:])
            zc = 0.5 * (zn[:-1] + zn[1:])
            xsel = (xc >= -10) & (xc <= 10)
            if not xsel.any():
                continue
            zsel = zc <= 3
            data = np.squeeze(sub.data)[:, xsel][:, :, zsel]
            key = "u"
            if key not in acc:
                acc[key] = np.zeros(data.shape[0])
                acc["n"] = 0.0
            acc[key][:data.shape[0]] += data.sum(axis=(1, 2))
            acc["n"] += data.shape[1] * data.shape[2]
        times = slc.times
        slc.clear_cache()
    u_roi = acc["u"] / acc["n"]

    ax = axes[k]
    ax.plot(times[:u_roi.size], u_roi, color=COLORS[run], lw=1)
    ax.axvline(30, color="gray", ls=":", lw=1)
    ax.annotate("ignition", (30, 0.97), xycoords=("data", "axes fraction"),
                rotation=90, va="top", ha="right", fontsize=7, color="gray")
    pre = (times >= 15) & (times <= 29)
    ax.hlines(u_roi[pre[:u_roi.size]].mean(), 15, 29, color="k", lw=2.5)
    ax.annotate(f"pre-fire {u_roi[pre[:u_roi.size]].mean():.2f}",
                (15, u_roi[pre[:u_roi.size]].mean()), fontsize=8,
                va="bottom", ha="left")
    td = FIRE_DEAD[run]
    if td is not None:
        late = times[:u_roi.size] >= td + 30
        m = u_roi[late].mean()
        ax.hlines(m, td + 30, times[:u_roi.size][-1], color="red", lw=2.5)
        ax.annotate(f"fire-free developed {m:.2f}",
                    (td + 30, m), fontsize=8, va="bottom", color="red")
        ax.axvspan(30, td, color="orange", alpha=0.08)
    else:
        ax.axvspan(30, times[:u_roi.size][-1], color="orange", alpha=0.08)
    ax.set_title(f"{run} (shaded: fire active)")
    ax.set_xlabel("t (s)")
    ax.grid(alpha=0.3)
    if k == 0:
        ax.set_ylabel("mean u in ROI, z<3 m (m/s)")
    print(run, "pre-fire ROI u =", round(float(u_roi[pre[:u_roi.size]].mean()), 3),
          "| late fire-free u =",
          None if td is None else round(float(u_roi[times[:u_roi.size] >= td + 30].mean()), 3))

fig.suptitle("ROI sub-canopy wind: pre-fire (15–29 s) vs fully-developed fire-free state")
fig.savefig(os.path.join(OUT_DIR, "roi_wind_equilibrium.png"), dpi=170,
            bbox_inches="tight")
