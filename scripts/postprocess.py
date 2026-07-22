#!/usr/bin/env python
"""Postprocess sonoma-tls-fds simulations with fdsreader.

Computes, per simulation:
  1. Fire front position and rate of spread (ROS) over time, from the
     cumulative 'FIRE ARRIVAL TIME' boundary file on the ground surface
     (FDS UG Sec. 22.10.34: arrival = first time near-surface HRRPUV exceeds
     the Smokeview fire-rendering threshold, UG Eq. 22.18).
  2. Flame height h_F over time from the vertical HRRPUV slices aligned with
     the spread direction (PBY = -2..2). Flame height follows Fig. 3 of
     Alexander (1982) / IJWF: the VERTICAL distance from the fuel bed to the
     highest point of flame, i.e. max z of cells with HRRPUV >= threshold
     (not the along-axis flame length L).

The default HRRPUV threshold is the same one FDS uses internally for
FIRE ARRIVAL TIME: qdot_min = min(200, 20/dx) kW/m3 (UG Eq. 22.18), which is
80 kW/m3 for the 0.25 m near-ground grid, so ROS and flame height are based
on a consistent flame definition.

ROI markers: the region of interest spans x in (-10, 10) (and y in (-10, 10)).
The front-position plot gets dashed horizontal red lines at x = -10/+10; the
time-series plots (ROS, flame height) get dashed red vertical lines at the
times the front entered/exited the ROI.

Usage:
    python scripts/postprocess.py "/Volumes/T7 Shield/sonoma-tls-fds/simulations"/*
    python scripts/postprocess.py --root "/Volumes/T7 Shield/sonoma-tls-fds/simulations"
"""

import argparse
import glob
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import fdsreader

ROI_X = (-10.0, 10.0)
ROI_Y = (-10.0, 10.0)
UNBURNED = 1.0e5  # arrival times >= this are "never burned" (FDS writes 1e6)


def ignition_time(sim_dir, chid, margin=2.0):
    """Time the fire actually starts, from CHID_hrr.csv.

    Spurious pre-ignition combustion (numerical startup transients of tens of
    kW) trips the FIRE ARRIVAL TIME threshold at scattered nodes long before
    the ignition burner fires, so arrivals before this time are masked.
    Ignition is taken as the first time total HRR exceeds 1 % of its peak
    (the burner is O(10 MW), the transients are O(10 kW)), minus a small
    margin. Returns 0.0 if the file is missing.
    """
    path = os.path.join(sim_dir, f"{chid}_hrr.csv")
    if not os.path.isfile(path):
        return 0.0
    data = np.genfromtxt(path, delimiter=",", skip_header=2)
    t, hrr = data[:, 0], data[:, 1]
    burning = hrr > 0.01 * hrr.max()
    return max(0.0, float(t[burning][0]) - margin) if burning.any() else 0.0


# ----------------------------------------------------------------------------
# Fire arrival time / ROS
# ----------------------------------------------------------------------------

def assemble_arrival_map(sim):
    """Stitch the per-mesh ground-face FIRE ARRIVAL TIME patches into one map.

    Returns (x_nodes, y_nodes, arrival) where arrival[i, j] is the arrival
    time (s) at (x_nodes[i], y_nodes[j]); unburned nodes hold ~1e6.
    The BNDF is cumulative, so only the final frame is needed.
    """
    patches = []  # (x_nodes, y_nodes, final_frame)
    for mesh in sim.meshes:
        try:
            bndf = mesh.get_boundary_data("FIRE ARRIVAL TIME")
        except StopIteration:
            continue
        for orientation, patch in bndf.data.items():
            if orientation != 3:  # ground face (+z normal)
                continue
            # a patch may cover only part of the mesh footprint (e.g. the
            # ignition VENT splits the ground), so use the patch extent
            nx, ny = patch.shape
            xm = np.round(np.linspace(patch.extent.x_start, patch.extent.x_end, nx), 6)
            ym = np.round(np.linspace(patch.extent.y_start, patch.extent.y_end, ny), 6)
            patches.append((xm, ym, patch.data[-1]))
        bndf.clear_cache()

    if not patches:
        raise RuntimeError("no ground FIRE ARRIVAL TIME patches found")

    x_nodes = np.unique(np.concatenate([p[0] for p in patches]))
    y_nodes = np.unique(np.concatenate([p[1] for p in patches]))
    arrival = np.full((x_nodes.size, y_nodes.size), 1.0e6)
    for xm, ym, data in patches:
        i0 = np.searchsorted(x_nodes, xm[0])
        j0 = np.searchsorted(y_nodes, ym[0])
        sl = arrival[i0:i0 + xm.size, j0:j0 + ym.size]
        # shared nodes on mesh interfaces: keep the earliest arrival
        np.minimum(sl, data, out=sl)
    return x_nodes, y_nodes, arrival


def front_arrival_curve(x_nodes, y_nodes, arrival, y_band, stat="median",
                        t_min=0.0):
    """Reduce the 2-D arrival map to tau(x), the front arrival time at each x.

    stat='median': time at which half the y-band at x has ignited (robust,
    represents the mean fireline). stat='min': leading edge (first finger).
    Returns (x, tau, tau_lead, tau_trail) with NaN where the front never
    arrived; tau is forced monotonic in x (spread is +x).
    """
    sel = np.abs(y_nodes) <= y_band
    band = arrival[:, sel]
    burned = (band < UNBURNED) & (band >= t_min)

    with np.errstate(invalid="ignore"):
        tau_lead = np.where(burned.any(axis=1),
                            np.min(np.where(burned, band, np.inf), axis=1),
                            np.nan)
        tau_trail = np.where(burned.all(axis=1),
                             band.max(axis=1), np.nan)
        # median with unburned treated as +inf (censored)
        band_inf = np.where(burned, band, np.inf)
        med = np.median(band_inf, axis=1)
        tau_med = np.where(np.isfinite(med), med, np.nan)

    tau = tau_med if stat == "median" else tau_lead
    # enforce monotonic arrival along +x (removes noise from fingering)
    valid = ~np.isnan(tau)
    tau_mono = tau.copy()
    tau_mono[valid] = np.maximum.accumulate(tau[valid])
    return x_nodes, tau_mono, tau_lead, tau_trail


def crossing_time(x, tau, x_target):
    """Time at which the front (tau vs x) reaches x_target, or NaN."""
    valid = np.isfinite(tau)
    if not valid.any() or x[valid].max() < x_target or x[valid].min() > x_target:
        return np.nan
    return float(np.interp(x_target, x[valid], tau[valid]))


def ros_timeseries(x, tau, dt=0.5, smooth_s=5.0):
    """Invert tau(x) to X_front(t) on a uniform grid and differentiate.

    Returns (t, x_front, ros). x_front is smoothed with a moving average of
    width smooth_s before differentiating, since the 0.25 m arrival grid makes
    the raw derivative staircase-noisy.
    """
    valid = np.isfinite(tau)
    xv, tv = x[valid], tau[valid]
    t = np.arange(tv.min(), tv.max(), dt)
    x_front = np.interp(t, tv, xv)

    n = max(3, int(round(smooth_s / dt)) | 1)  # odd window
    kernel = np.ones(n) / n
    pad = np.pad(x_front, n // 2, mode="edge")
    x_smooth = np.convolve(pad, kernel, mode="valid")
    ros = np.gradient(x_smooth, t)
    return t, x_front, ros


# ----------------------------------------------------------------------------
# Flame height from HRRPUV slices
# ----------------------------------------------------------------------------

def flame_height_from_slices(sim, threshold, x_limits=None):
    """Flame height h_F(t) from each vertical (orientation-2) HRRPUV slice.

    h_F = max cell-center height z where HRRPUV >= threshold anywhere in the
    slice (vertical distance above the fuel bed, per Fig. 3), assembled
    manually per mesh because Slice.to_global() cannot handle the stretched
    TRNZ grid. Returns (times, heights[n_slice, n_t], y_positions).
    """
    slices = sorted(
        (s for s in sim.slices
         if s.quantity.name == "HRRPUV" and s.orientation == 2),
        key=lambda s: s.extent.y_start,
    )
    if not slices:
        raise RuntimeError("no vertical HRRPUV slices found")

    ref_times = slices[0].times
    heights = np.full((len(slices), ref_times.size), np.nan)
    y_pos = []

    for k, slc in enumerate(slices):
        y_pos.append(slc.extent.y_start)
        h = np.full(slc.times.size, 0.0)
        for sub in slc._subslices.values():
            mesh = sub.mesh
            xn = mesh.coordinates["x"]
            zn = mesh.coordinates["z"]
            xc = 0.5 * (xn[:-1] + xn[1:])
            zc = 0.5 * (zn[:-1] + zn[1:])
            data = np.squeeze(sub.data)  # (n_t, n_x, n_z), cell centered
            if x_limits is not None:
                keep = (xc >= x_limits[0]) & (xc <= x_limits[1])
                if not keep.any():
                    continue
                data = data[:, keep, :]
            mask = data >= threshold  # flame per UG Eq. 22.18 threshold
            h_sub = np.max(np.where(mask.any(axis=1), zc[None, :], 0.0), axis=1)
            np.maximum(h, h_sub, out=h)
        slc.clear_cache()
        # slices share DT but interpolate defensively onto the reference grid
        heights[k] = np.interp(ref_times, slc.times, h)

    return np.asarray(ref_times), heights, y_pos


# ----------------------------------------------------------------------------
# Output
# ----------------------------------------------------------------------------

def write_csv(path, header, columns):
    data = np.column_stack(columns)
    np.savetxt(path, data, delimiter=",", header=",".join(header), comments="")


def process_simulation(sim_dir, out_root, args):
    name = os.path.basename(os.path.normpath(sim_dir))
    out_dir = os.path.join(out_root, name)
    os.makedirs(out_dir, exist_ok=True)
    print(f"=== {name} ===")

    sim = fdsreader.Simulation(sim_dir)

    # grid size near the ground sets the FDS fire threshold (UG Eq. 22.18)
    dx = float(np.diff(sim.meshes[0].coordinates["x"]).min())
    threshold = args.hrrpuv_threshold
    if threshold is None:
        threshold = min(200.0, 20.0 / dx)
    print(f"  HRRPUV flame threshold: {threshold:.0f} kW/m3 (dx={dx:.2f} m)")

    # --- ROS from arrival time ---
    t_ign = ignition_time(sim_dir, sim.chid)
    print(f"  ignition detected at t={t_ign:.1f} s (arrivals before this masked)")
    x_nodes, y_nodes, arrival = assemble_arrival_map(sim)
    x, tau, tau_lead, tau_trail = front_arrival_curve(
        x_nodes, y_nodes, arrival, args.y_band, stat=args.front_stat,
        t_min=t_ign)

    t_enter = crossing_time(x, tau, ROI_X[0])
    t_exit = crossing_time(x, tau, ROI_X[1])
    t_enter_lead = crossing_time(x, np.where(np.isfinite(tau_lead), tau_lead, np.nan), ROI_X[0])
    roi_ros = (ROI_X[1] - ROI_X[0]) / (t_exit - t_enter) \
        if np.isfinite(t_enter) and np.isfinite(t_exit) else np.nan

    t, x_front, ros = ros_timeseries(x, tau, smooth_s=args.smooth)

    # --- flame height ---
    times_f, heights, y_pos = flame_height_from_slices(
        sim, threshold, x_limits=(args.xmin, args.xmax))
    h_mean = heights.mean(axis=0)
    h_max = heights.max(axis=0)

    # --- CSVs ---
    write_csv(os.path.join(out_dir, "front_trajectory.csv"),
              ["x_m", f"t_{args.front_stat}_s", "t_leading_s", "t_trailing_s"],
              [x, tau, tau_lead, tau_trail])
    write_csv(os.path.join(out_dir, "ros_time.csv"),
              ["t_s", "x_front_m", "ros_m_per_s"], [t, x_front, ros])
    write_csv(os.path.join(out_dir, "flame_height.csv"),
              ["t_s", "hF_mean_m", "hF_max_m"]
              + [f"hF_y{y:+.2f}_m" for y in y_pos],
              [times_f, h_mean, h_max] + [heights[k] for k in range(len(y_pos))])

    in_roi = (t >= t_enter) & (t <= t_exit) if np.isfinite(t_enter) and np.isfinite(t_exit) else np.zeros_like(t, bool)
    in_roi_f = (times_f >= t_enter) & (times_f <= t_exit) if in_roi.any() else np.zeros_like(times_f, bool)
    with open(os.path.join(out_dir, "summary.csv"), "w") as f:
        f.write("quantity,value,units\n")
        f.write(f"hrrpuv_threshold,{threshold:.1f},kW/m3\n")
        f.write(f"t_enter_roi_{args.front_stat},{t_enter:.2f},s\n")
        f.write(f"t_exit_roi_{args.front_stat},{t_exit:.2f},s\n")
        f.write(f"t_enter_roi_leading_edge,{t_enter_lead:.2f},s\n")
        f.write(f"mean_ros_in_roi,{roi_ros:.4f},m/s\n")
        if in_roi.any():
            f.write(f"mean_ros_in_roi_timeseries,{ros[in_roi].mean():.4f},m/s\n")
        if in_roi_f.any():
            f.write(f"mean_flame_height_in_roi,{h_mean[in_roi_f].mean():.2f},m\n")
            f.write(f"max_flame_height_in_roi,{h_max[in_roi_f].max():.2f},m\n")

    # --- plots ---
    label_stat = "median" if args.front_stat == "median" else "leading-edge"

    fig, ax = plt.subplots(figsize=(8, 5))
    valid = np.isfinite(tau_lead) & np.isfinite(tau_trail)
    ax.fill_betweenx(x[valid], tau_lead[valid], tau_trail[valid],
                     alpha=0.2, color="tab:orange",
                     label="leading/trailing edge envelope")
    ax.plot(tau, x, color="tab:red", lw=2, label=f"front position ({label_stat})")
    for xr, lbl in zip(ROI_X, ("ROI start (x=-10 m)", "ROI end (x=+10 m)")):
        ax.axhline(xr, color="red", ls="--", lw=1.2)
        ax.annotate(lbl, (0.01, xr), xycoords=("axes fraction", "data"),
                    va="bottom", fontsize=8, color="red")
    if np.isfinite(t_enter):
        ax.plot([t_enter], [ROI_X[0]], "ro", ms=5)
        ax.annotate(f"enter {t_enter:.1f} s", (t_enter, ROI_X[0]),
                    textcoords="offset points", xytext=(6, -12), fontsize=8)
    if np.isfinite(t_exit):
        ax.plot([t_exit], [ROI_X[1]], "ro", ms=5)
        ax.annotate(f"exit {t_exit:.1f} s", (t_exit, ROI_X[1]),
                    textcoords="offset points", xytext=(6, 6), fontsize=8)
    ax.set_xlabel("time (s)")
    ax.set_ylabel("fire front position x (m)")
    ax.set_title(f"{name}: fire front from FIRE ARRIVAL TIME (|y| ≤ {args.y_band:g} m)")
    ax.legend(loc="upper left", fontsize=8)
    ax.grid(alpha=0.3)
    fig.savefig(os.path.join(out_dir, "front_position.png"), dpi=200,
                bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(t, ros, color="tab:blue", lw=1.5)
    for tc, lbl in ((t_enter, "fire enters ROI"), (t_exit, "fire exits ROI")):
        if np.isfinite(tc):
            ax.axvline(tc, color="red", ls="--", lw=1.2)
            ax.annotate(lbl, (tc, 0.98), xycoords=("data", "axes fraction"),
                        rotation=90, va="top", ha="right", fontsize=8, color="red")
    if np.isfinite(roi_ros):
        ax.axhline(roi_ros, color="k", ls=":", lw=1,
                   label=f"ROI mean ROS = {roi_ros:.3f} m/s")
        ax.legend(fontsize=8)
    ax.set_xlabel("time (s)")
    ax.set_ylabel("rate of spread (m/s)")
    ax.set_title(f"{name}: ROS ({label_stat} front, {args.smooth:g} s smoothing)")
    ax.grid(alpha=0.3)
    fig.savefig(os.path.join(out_dir, "ros.png"), dpi=200, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 5))
    for k, y in enumerate(y_pos):
        ax.plot(times_f, heights[k], lw=0.7, alpha=0.4,
                label=f"y = {y:+.2f} m")
    ax.plot(times_f, h_mean, color="k", lw=2, label="mean of slices")
    for tc, lbl in ((t_enter, "fire enters ROI"), (t_exit, "fire exits ROI")):
        if np.isfinite(tc):
            ax.axvline(tc, color="red", ls="--", lw=1.2)
            ax.annotate(lbl, (tc, 0.98), xycoords=("data", "axes fraction"),
                        rotation=90, va="top", ha="right", fontsize=8, color="red")
    ax.set_xlabel("time (s)")
    ax.set_ylabel("flame height $h_F$ (m)")
    ax.set_title(f"{name}: flame height (HRRPUV ≥ {threshold:.0f} kW/m³)")
    ax.legend(fontsize=7, ncol=2)
    ax.grid(alpha=0.3)
    fig.savefig(os.path.join(out_dir, "flame_height.png"), dpi=200,
                bbox_inches="tight")
    plt.close(fig)

    print(f"  t_enter={t_enter:.1f} s  t_exit={t_exit:.1f} s  "
          f"ROI ROS={roi_ros:.3f} m/s")
    return {
        "name": name, "t": t, "x_front": x_front, "ros": ros,
        "tau": tau, "x": x, "times_f": times_f, "h_mean": h_mean,
        "t_enter": t_enter, "t_exit": t_exit, "roi_ros": roi_ros,
    }


def comparison_plots(results, out_root):
    fig, ax = plt.subplots(figsize=(8, 5))
    for r in results:
        ax.plot(r["tau"], r["x"], lw=1.8, label=r["name"])
    for xr in ROI_X:
        ax.axhline(xr, color="red", ls="--", lw=1.2)
    ax.set_xlabel("time (s)")
    ax.set_ylabel("fire front position x (m)")
    ax.set_title("Fire front position (dashed red: ROI bounds)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.savefig(os.path.join(out_root, "comparison_front_position.png"),
                dpi=200, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 5))
    for r in results:
        ax.plot(r["times_f"], r["h_mean"], lw=1.8, label=r["name"])
    ax.set_xlabel("time (s)")
    ax.set_ylabel("flame height $h_F$ (m)")
    ax.set_title("Mean flame height")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.savefig(os.path.join(out_root, "comparison_flame_height.png"),
                dpi=200, bbox_inches="tight")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("sim_dirs", nargs="*", help="simulation directories")
    ap.add_argument("--root", help="process every subdirectory of this dir")
    ap.add_argument("--out", default="postprocess",
                    help="output directory (default: ./postprocess)")
    ap.add_argument("--hrrpuv-threshold", type=float, default=None,
                    help="flame HRRPUV threshold in kW/m3 "
                         "(default: min(200, 20/dx) per FDS UG Eq. 22.18)")
    ap.add_argument("--y-band", type=float, default=10.0,
                    help="half-width of y band for the front (default 10 m = ROI)")
    ap.add_argument("--front-stat", choices=("median", "min"), default="median",
                    help="front position statistic across the y band "
                         "(median = mean fireline, min = leading edge)")
    ap.add_argument("--smooth", type=float, default=5.0,
                    help="smoothing window (s) for the ROS derivative")
    ap.add_argument("--xmin", type=float, default=None,
                    help="restrict flame-height search to x >= xmin")
    ap.add_argument("--xmax", type=float, default=None,
                    help="restrict flame-height search to x <= xmax")
    args = ap.parse_args()

    sim_dirs = list(args.sim_dirs)
    if args.root:
        sim_dirs += sorted(d for d in glob.glob(os.path.join(args.root, "*"))
                           if os.path.isdir(d))
    if not sim_dirs:
        ap.error("no simulation directories given (positional args or --root)")
    if args.xmin is None:
        args.xmin = -np.inf
    if args.xmax is None:
        args.xmax = np.inf

    os.makedirs(args.out, exist_ok=True)
    results = []
    for d in sim_dirs:
        try:
            results.append(process_simulation(d, args.out, args))
        except Exception as e:
            print(f"  ERROR processing {d}: {e}", file=sys.stderr)
    if len(results) > 1:
        comparison_plots(results, args.out)
    print(f"Done. Outputs in {os.path.abspath(args.out)}")


if __name__ == "__main__":
    main()
