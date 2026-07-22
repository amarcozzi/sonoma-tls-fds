#!/usr/bin/env python
"""Test two claims about the wind field using existing slice data:

Claim A (entrance effect): the inflow/stand-edge adjustment produces a
sub-canopy jet near XMIN. Question: does it decay before the ROI, or does it
contaminate the flow where ROS/flame height are measured?
  -> time-averaged pre-fire u(x,z) cross-section (PBY=0), per run, plus
     vertical profiles of u at x = -70, -40, -20, 0 to test self-similarity.

Claim B (spin-up): 30 s is not enough for the flow to equilibrate.
  -> stationarity: u(t) averaged over sub-canopy region (from the PBY=0
     slice) in windows 5-15 s vs 15-29 s, plus full time series of
     spatial-mean u at z=2, 10, 20 m (PBZ slices) for t < 30 s.
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


def assemble_xz(slc, tsel):
    """Time-averaged u(x,z) from a vertical slice, stitched across meshes."""
    pieces = []
    for sub in slc._subslices.values():
        mesh = sub.mesh
        xn, zn = mesh.coordinates["x"], mesh.coordinates["z"]
        xc = 0.5 * (xn[:-1] + xn[1:])
        zc = 0.5 * (zn[:-1] + zn[1:])
        data = np.squeeze(sub.data)  # (nt, nx, nz)
        pieces.append((xc, zc, data[tsel].mean(axis=0)))
    pieces.sort(key=lambda p: p[0][0])
    # meshes overlapping y=0 produce two copies per x -- average duplicates
    xs, out = [], []
    for xc, zc, u in pieces:
        if xs and abs(xs[-1][0] - xc[0]) < 1e-6:
            out[-1] = 0.5 * (out[-1] + u)
        else:
            xs.append(xc)
            out.append(u)
    return np.concatenate(xs), zc, np.vstack(out)


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    fig_cs, axes_cs = plt.subplots(3, 1, figsize=(11, 10), sharex=True)
    fig_pr, axes_pr = plt.subplots(1, 3, figsize=(13, 4.5), sharey=True)
    fig_st, axes_st = plt.subplots(1, 3, figsize=(13, 4), sharey=True)

    for k, run in enumerate(RUNS):
        sim = fdsreader.Simulation(os.path.join(SIM_ROOT, run))
        slc = next(s for s in sim.slices
                   if s.quantity.name == "U-VELOCITY" and s.orientation == 2
                   and abs(s.extent.y_start) < 0.5)

        # --- Claim A: cross-section, pre-fire average
        tsel = (slc.times >= 15) & (slc.times <= 29)
        x, zc, u_xz = assemble_xz(slc, tsel)
        ax = axes_cs[k]
        zmask = zc <= 40
        pc = ax.pcolormesh(x, zc[zmask], u_xz[:, zmask].T, cmap="RdBu_r",
                           vmin=-8, vmax=8, shading="nearest")
        ax.contour(x, zc[zmask], u_xz[:, zmask].T, levels=[1, 2, 4],
                   colors="k", linewidths=0.5)
        for xv in (-10, 10):
            ax.axvline(xv, color="red", ls="--", lw=1)
        ax.axvline(-60, color="gray", ls=":", lw=1)
        ax.set_ylabel("z (m)")
        ax.set_title(f"{run}: mean u(x,z), t=15–29 s "
                     f"(contours at 1/2/4 m/s; red: ROI, gray: ignition)")
        if k == 2:
            ax.set_xlabel("x (m)")
        fig_cs.colorbar(pc, ax=ax, label="u (m/s)")

        # --- Claim A: vertical profile self-similarity
        axp = axes_pr[k]
        for x0, c in zip((-70, -40, -20, 0),
                         ("tab:red", "tab:purple", "tab:brown", "k")):
            i = np.argmin(np.abs(x - x0))
            band = u_xz[max(0, i - 4):i + 5].mean(axis=0)  # +-1 m in x
            axp.plot(band[zmask], zc[zmask], color=c, label=f"x = {x0} m")
        axp.set_title(run)
        axp.set_xlabel("u (m/s)")
        axp.grid(alpha=0.3)
        if k == 0:
            axp.set_ylabel("z (m)")
            axp.legend(fontsize=8)

        # --- Claim B: stationarity of the sub-canopy flow
        axs = axes_st[k]
        # sub-canopy mean over the whole fetch, every timestep to t=35
        zsel = zc <= 3
        fetch = None
        ts, us = [], []
        for it, t in enumerate(slc.times):
            if t > 35:
                break
            ts.append(t)
        ts = np.array(ts)
        # recompute cheaply: stitched instantaneous mean over x in fetch, z<3
        u_t = np.zeros(ts.size)
        w_t = np.zeros(ts.size)
        for sub in slc._subslices.values():
            mesh = sub.mesh
            xn, zn = mesh.coordinates["x"], mesh.coordinates["z"]
            xc = 0.5 * (xn[:-1] + xn[1:])
            zcs = 0.5 * (zn[:-1] + zn[1:])
            data = np.squeeze(sub.data)
            xsel = (xc >= -55) & (xc <= 10)
            zs = zcs <= 3
            if not xsel.any():
                continue
            u_t += data[:ts.size][:, xsel][:, :, zs].sum(axis=(1, 2))
            w_t += xsel.sum() * zs.sum()
        axs.plot(ts, u_t / w_t, color=COLORS[run])
        for w0, w1, c in ((5, 15, "0.6"), (15, 29, "0.3")):
            m = u_t[(ts >= w0) & (ts <= w1)].sum() / w_t[(ts >= w0) & (ts <= w1)].sum()
            axs.hlines(m, w0, w1, color=c, lw=2)
            axs.annotate(f"{m:.2f}", (w1, m), fontsize=8, va="bottom")
        axs.axvline(30, color="red", ls="--", lw=1)
        axs.set_title(f"{run}: mean u (x∈[-55,10], z<3 m)")
        axs.set_xlabel("t (s)")
        axs.grid(alpha=0.3)
        if k == 0:
            axs.set_ylabel("u (m/s)")
        slc.clear_cache()
        print(run, "done")

    fig_cs.savefig(os.path.join(OUT_DIR, "u_cross_sections.png"), dpi=170,
                   bbox_inches="tight")
    fig_pr.suptitle("Vertical u profiles at several fetch positions (t=15–29 s mean)")
    fig_pr.savefig(os.path.join(OUT_DIR, "u_profiles_vs_x.png"), dpi=170,
                   bbox_inches="tight")
    fig_st.suptitle("Pre-ignition stationarity (ignition at t=30 s, red line)")
    fig_st.savefig(os.path.join(OUT_DIR, "stationarity.png"), dpi=170,
                   bbox_inches="tight")


if __name__ == "__main__":
    main()
