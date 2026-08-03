#!/usr/bin/env python
"""Per-simulation postprocessing for sonoma-tls-fds (runs on the HPC).

Designed to run one Slurm array task per simulation, next to the raw FDS
output, so only a few MB of CSV/PNG per run ever leaves Ceres.
`aggregate_metrics.py` then merges the per-simulation `metrics.json` files
into the cross-simulation tables and comparison figures.

What is computed
----------------
1. **Fire front and rate of spread (ROS)** from the cumulative
   'FIRE ARRIVAL TIME' boundary file on the ground surface (FDS UG
   Sec. 22.10.34: arrival = first time near-surface HRRPUV exceeds the
   Smokeview fire-rendering threshold, UG Eq. 22.18). Two ROS estimators:

   * `ros`      - d/dt of the smoothed front position X_front(t). Dense in
                  time, but successive samples are correlated by the
                  smoothing window.
   * `ros_seg`  - chord ROS over fixed `--segment` metre intervals of x,
                  i.e. seg / (tau(x+seg) - tau(x)). Sparser but the samples
                  are close to independent, so its spread is the more
                  defensible measure of ROS variability.

2. **Flame geometry** from the vertical HRRPUV slices spanning the spread
   direction (PBY = -2..2), per Fig. 3 of Alexander (1982) / IJWF:

   * `flame_height` h_F - vertical distance from the fuel bed to the highest
                          point of the flame.
   * `flame_length` L   - straight-line distance from the *middle of the
                          flame base* to the flame tip. For a wind-tilted
                          flame L > h_F; this is the quantity that feeds
                          Byram-style fireline-intensity relations.
   * `flame_tilt`       - angle of that base->tip vector from vertical (deg).
   * `flame_depth`      - along-wind extent of the flaming zone at the base.

   Each slice is assembled across meshes by hand (fdsreader's
   `Slice.to_global()` cannot handle the stretched TRNZ grid), the flaming
   cells are restricted to a +/-`--flame-window` metre band around the tracked
   front so residual burning near the igniter cannot be mistaken for the head
   flame, and the flame is taken to be the ground-connected component of that
   masked field nearest the front. Detached flamelets are therefore excluded.

3. **Distributional statistics** of every time series above, over three
   analysis windows:

   * `roi`    - front inside the region of interest, x in (-10, 10).
   * `steady` - front inside `--steady-x` (default -30..20 m): past the
                stand-edge entrance jet, short of the kill devices at x=+30.
   * `full`   - the whole tracked propagation record.

   Statistics: n, mean, std, min, p02, p05, p25, p50, p75, p95, p98, max.
   Flame statistics pool all slices and all in-window times (so n is ~5x the
   number of frames); ROS statistics use the single front track.

Usage
-----
    python scripts/postprocess.py simulations/c4_p10 --out postprocess
    python scripts/postprocess.py --root simulations --out postprocess
    python scripts/postprocess.py --root simulations --no-figures   # tables only
"""

import argparse
import glob
import json
import math
import os
import sys
import traceback
import warnings

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import ndimage

import fdsreader

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import plotstyle

ROI_X = (-10.0, 10.0)
UNBURNED = 1.0e5   # arrival times >= this are "never burned" (FDS writes 1e6)

STAT_NAMES = ("n", "mean", "std", "min", "p02", "p05", "p25", "p50", "p75",
              "p95", "p98", "max")
PERCENTILES = (2, 5, 25, 50, 75, 95, 98)

# metric key -> (units, human description) ; drives the data dictionary
METRICS = {
    "ros": ("m/s", "ROS from d/dt of the smoothed front position"),
    "ros_seg": ("m/s", "chord ROS over fixed-length x segments"),
    "flame_length": ("m", "base-midpoint to tip distance (tilted flame length)"),
    "flame_height": ("m", "vertical fuel bed to flame tip distance"),
    "flame_tilt": ("deg", "flame axis angle from vertical"),
    "flame_depth": ("m", "along-wind extent of the flaming zone at the base"),
}
WINDOWS = ("roi", "steady", "full")


# ---------------------------------------------------------------------------
# statistics helpers
# ---------------------------------------------------------------------------

def describe(values):
    """n / mean / std / min / percentiles / max of the finite entries."""
    v = np.asarray(values, dtype=float).ravel()
    v = v[np.isfinite(v)]
    if v.size == 0:
        return {k: (0 if k == "n" else float("nan")) for k in STAT_NAMES}
    q = np.percentile(v, PERCENTILES)
    out = {
        "n": int(v.size),
        "mean": float(v.mean()),
        "std": float(v.std(ddof=1)) if v.size > 1 else float("nan"),
        "min": float(v.min()),
        "max": float(v.max()),
    }
    for p, val in zip(PERCENTILES, q):
        out[f"p{p:02d}"] = float(val)
    return out


def mask_window(times, bounds):
    """Boolean mask of `times` inside the (t0, t1) window; empty if degenerate."""
    t0, t1 = bounds
    if not (np.isfinite(t0) and np.isfinite(t1)) or t1 <= t0:
        return np.zeros(np.shape(times), dtype=bool)
    return (np.asarray(times) >= t0) & (np.asarray(times) <= t1)


# ---------------------------------------------------------------------------
# ignition time
# ---------------------------------------------------------------------------

def ignition_time(sim_dir, chid, margin=2.0):
    """Time the fire actually starts, from CHID_hrr.csv.

    Spurious pre-ignition combustion (numerical startup transients of tens of
    kW) trips the FIRE ARRIVAL TIME threshold at scattered nodes long before
    the ignition burner fires, so arrivals before this time are masked.
    Ignition is the first time total HRR exceeds 1 % of its peak (the burner
    is O(10 MW), the transients are O(10 kW)), minus a small margin.
    """
    path = os.path.join(sim_dir, f"{chid}_hrr.csv")
    if not os.path.isfile(path):
        return 0.0
    try:
        data = np.genfromtxt(path, delimiter=",", skip_header=2)
    except Exception:
        return 0.0
    if data.ndim != 2 or data.shape[0] == 0:
        return 0.0
    t, hrr = data[:, 0], data[:, 1]
    burning = hrr > 0.01 * np.nanmax(hrr)
    return max(0.0, float(t[burning][0]) - margin) if burning.any() else 0.0


# ---------------------------------------------------------------------------
# fire arrival time -> front -> ROS
# ---------------------------------------------------------------------------

def assemble_arrival_map(sim):
    """Stitch the per-mesh ground-face FIRE ARRIVAL TIME patches into one map.

    Returns (x_nodes, y_nodes, arrival) where arrival[i, j] is the arrival time
    (s) at (x_nodes[i], y_nodes[j]); unburned nodes hold ~1e6. The BNDF is
    cumulative, so only the final frame is needed.
    """
    patches = []
    for mesh in sim.meshes:
        try:
            bndf = mesh.get_boundary_data("FIRE ARRIVAL TIME")
        except (StopIteration, KeyError):
            continue
        for orientation, patch in bndf.data.items():
            if orientation != 3:      # ground face (+z normal)
                continue
            # a patch may cover only part of the mesh footprint (the ignition
            # VENT splits the ground), so use the patch extent, not the mesh
            nx, ny = patch.shape
            xm = np.round(np.linspace(patch.extent.x_start, patch.extent.x_end, nx), 6)
            ym = np.round(np.linspace(patch.extent.y_start, patch.extent.y_end, ny), 6)
            patches.append((xm, ym, np.asarray(patch.data[-1], dtype=float)))
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


def ignition_line_x(x_nodes, tau_lead, tol=1.5):
    """x of the ignition strip: the upwind edge of the earliest-arriving band.

    Everything upwind of this is the BACKING fire, which spreads the other way
    and therefore arrives late at low x. Sorted by x those points come first,
    so leaving them in front_arrival_curve's monotonic pass lets a 55 s backing
    arrival clamp the entire forward front to 55 s (a flat tau(x), i.e. an
    apparently infinite ROS). They must be dropped before the pass, not after.
    """
    if not np.isfinite(tau_lead).any():
        return -np.inf
    t0 = np.nanmin(tau_lead)
    return float(x_nodes[np.isfinite(tau_lead) & (tau_lead <= t0 + tol)].min())


def front_arrival_curve(x_nodes, y_nodes, arrival, y_band, stat="median",
                        t_min=0.0, x_ignition=None):
    """Reduce the 2-D arrival map to tau(x), the front arrival time at each x.

    stat='median': time at which half the y-band at x has ignited (robust,
    represents the mean fireline, and immune to the scattered spot-fire nodes
    that make the leading edge non-monotonic). stat='min': leading edge.

    Only x >= the ignition line is kept (see ignition_line_x); over that range
    tau is forced monotonic in x, since spread is +x and the raw curve is
    staircase-noisy at the 0.25 m arrival resolution.

    Returns (x, tau, tau_raw, tau_lead, tau_trail, x_ignition) where tau_raw is
    the same statistic BEFORE the monotonic pass, so the effect of that pass
    stays auditable in the output CSV.
    """
    sel = np.abs(y_nodes) <= y_band
    band = arrival[:, sel]
    burned = (band < UNBURNED) & (band >= t_min)

    with np.errstate(invalid="ignore"):
        tau_lead = np.where(burned.any(axis=1),
                            np.min(np.where(burned, band, np.inf), axis=1),
                            np.nan)
        tau_trail = np.where(burned.all(axis=1), band.max(axis=1), np.nan)
        band_inf = np.where(burned, band, np.inf)      # unburned = censored
        med = np.median(band_inf, axis=1)
        tau_med = np.where(np.isfinite(med), med, np.nan)

    if x_ignition is None:
        x_ignition = ignition_line_x(x_nodes, tau_lead)

    tau_raw = np.where(x_nodes >= x_ignition,
                       tau_med if stat == "median" else tau_lead, np.nan)
    valid = ~np.isnan(tau_raw)
    tau_mono = tau_raw.copy()
    tau_mono[valid] = np.maximum.accumulate(tau_raw[valid])   # spread is +x
    return x_nodes, tau_mono, tau_raw, tau_lead, tau_trail, x_ignition


def flattened_fraction(tau_raw, tau_mono, tol=0.5):
    """Fraction of the tracked front the monotonic pass had to lift.

    A healthy front is a few percent. A large value means tau(x) is dominated
    by one late upwind arrival being smeared downwind, and every ROS derived
    from it is meaningless -- so it is reported, not silently accepted.
    """
    valid = np.isfinite(tau_raw) & np.isfinite(tau_mono)
    if not valid.any():
        return float("nan")
    return float(np.mean(tau_mono[valid] > tau_raw[valid] + tol))


def crossing_time(x, tau, x_target):
    """Time at which the front (tau vs x) reaches x_target, or NaN."""
    valid = np.isfinite(tau)
    if not valid.any():
        return float("nan")
    xv, tv = x[valid], tau[valid]
    if xv.max() < x_target or xv.min() > x_target:
        return float("nan")
    return float(np.interp(x_target, xv, tv))


def ros_timeseries(x, tau, dt=0.5, smooth_s=5.0):
    """Invert tau(x) to X_front(t) on a uniform grid and differentiate.

    x_front is smoothed with a moving average of width smooth_s before
    differentiating, since the 0.25 m arrival grid makes the raw derivative
    staircase-noisy.
    """
    valid = np.isfinite(tau)
    if valid.sum() < 3:
        empty = np.array([])
        return empty, empty, empty
    xv, tv = x[valid], tau[valid]
    if tv.max() - tv.min() < dt:
        empty = np.array([])
        return empty, empty, empty
    t = np.arange(tv.min(), tv.max(), dt)
    x_front = np.interp(t, tv, xv)

    n = max(3, int(round(smooth_s / dt)) | 1)      # odd window
    kernel = np.ones(n) / n
    pad = np.pad(x_front, n // 2, mode="edge")
    x_smooth = np.convolve(pad, kernel, mode="valid")
    return t, x_front, np.gradient(x_smooth, t)


def segment_ros(x, tau, seg_len=5.0):
    """Chord ROS over non-overlapping seg_len intervals of x.

    Returns (x_mid, t_mid, ros). Nearly independent samples, so this is the
    estimator to quote spread from.
    """
    valid = np.isfinite(tau)
    if valid.sum() < 2:
        empty = np.array([])
        return empty, empty, empty
    xv, tv = x[valid], tau[valid]
    lo = math.ceil(xv.min() / seg_len) * seg_len
    hi = math.floor(xv.max() / seg_len) * seg_len
    if hi - lo < seg_len:
        empty = np.array([])
        return empty, empty, empty
    edges = np.arange(lo, hi + 0.5 * seg_len, seg_len)
    t_edge = np.interp(edges, xv, tv)
    dt = np.diff(t_edge)
    with np.errstate(divide="ignore", invalid="ignore"):
        ros = np.where(dt > 0, seg_len / dt, np.nan)
    return 0.5 * (edges[:-1] + edges[1:]), 0.5 * (t_edge[:-1] + t_edge[1:]), ros


def front_at_times(x, tau, times):
    """X_front interpolated onto arbitrary times; NaN outside the record."""
    valid = np.isfinite(tau)
    if valid.sum() < 2:
        return np.full(np.shape(times), np.nan)
    return np.interp(times, tau[valid], x[valid],
                     left=np.nan, right=np.nan)


# ---------------------------------------------------------------------------
# flame geometry from vertical HRRPUV slices
# ---------------------------------------------------------------------------

def _sub_axes(sub, n_x, n_z):
    """Cell-centre x and z axes for a subslice, tolerant of centring."""
    coords = sub.mesh.coordinates
    xn, zn = np.asarray(coords["x"]), np.asarray(coords["z"])
    xc, zc = 0.5 * (xn[:-1] + xn[1:]), 0.5 * (zn[:-1] + zn[1:])
    if xc.size == n_x and zc.size == n_z:
        return xc, zc
    if xn.size == n_x and zn.size == n_z:
        return xn, zn
    ext = sub.extent      # last resort: even spacing over the subslice extent
    return (np.linspace(ext.x_start, ext.x_end, n_x),
            np.linspace(ext.z_start, ext.z_end, n_z))


def assemble_slice(slc):
    """Assemble a vertical (PBY) slice across meshes.

    Returns (times, xc, zc, data[n_t, n_x, n_z]) as float32. Written by hand
    because fdsreader 1.11.x `Slice.to_global()` fails on the stretched TRNZ
    grid.
    """
    pieces = []
    for sub in slc._subslices.values():
        raw = np.squeeze(np.asarray(sub.data))
        if raw.ndim != 3:
            continue
        xc, zc = _sub_axes(sub, raw.shape[1], raw.shape[2])
        pieces.append((np.round(xc, 6), np.round(zc, 6), raw.astype(np.float32)))
    if not pieces:
        raise RuntimeError(f"slice at y={slc.extent.y_start} has no usable data")

    zc = pieces[0][1]
    for _, z, _ in pieces[1:]:
        if z.size != zc.size or not np.allclose(z, zc):
            raise RuntimeError("meshes disagree on the z grid; cannot assemble "
                               "the slice (expected an identical TRNZ on every "
                               "mesh)")

    x_all = np.unique(np.concatenate([p[0] for p in pieces]))
    n_t = min(p[2].shape[0] for p in pieces)
    data = np.zeros((n_t, x_all.size, zc.size), dtype=np.float32)
    for xc_p, _, d in pieces:
        i0 = int(np.searchsorted(x_all, xc_p[0]))
        block = data[:, i0:i0 + xc_p.size, :]
        np.maximum(block, d[:n_t], out=block)
    return np.asarray(slc.times[:n_t], dtype=float), x_all, zc, data


def flame_geometry_frames(times, xc, zc, data, threshold, x_front,
                          z_base=0.5, window=15.0):
    """Head-fire flame geometry per frame from one assembled vertical slice.

    At each time the flaming field (HRRPUV >= threshold) is restricted to
    |x - X_front(t)| <= window, connected components are labelled, and the
    ground-touching component whose base midpoint is nearest the front is
    taken to be the head flame. Frames with no tracked front, no flame, or no
    ground-connected flame yield NaN.
    """
    n_t = times.size
    keys = ("flame_length", "flame_height", "flame_tilt", "flame_depth",
            "x_tip", "z_tip", "x_base")
    out = {k: np.full(n_t, np.nan) for k in keys}

    base_rows = zc <= z_base
    if not base_rows.any():
        base_rows = np.zeros(zc.size, dtype=bool)
        base_rows[0] = True
    struct = np.ones((3, 3), dtype=bool)     # 8-connectivity in (x, z)

    for it in range(n_t):
        xf = x_front[it]
        if not np.isfinite(xf):
            continue
        keep = np.abs(xc - xf) <= window
        if not keep.any():
            continue
        xk = xc[keep]
        flaming = data[it][keep] >= threshold
        if not flaming.any():
            continue

        labels, n_lab = ndimage.label(flaming, structure=struct)
        if n_lab == 0:
            continue
        grounded = np.unique(labels[:, base_rows])
        grounded = grounded[grounded > 0]
        if grounded.size == 0:
            continue

        best, best_dist, best_base = None, np.inf, None
        for lab in grounded:
            cells = labels == lab
            base_x = xk[cells[:, base_rows].any(axis=1)]
            if base_x.size == 0:
                continue
            dist = abs(0.5 * (base_x.min() + base_x.max()) - xf)
            if dist < best_dist:
                best, best_dist, best_base = lab, dist, base_x
        if best is None:
            continue

        ix, iz = np.nonzero(labels == best)
        k_tip = int(iz.max())
        z_tip = float(zc[k_tip])
        x_tip = float(np.median(xk[ix[iz == k_tip]]))
        x_lo, x_hi = float(best_base.min()), float(best_base.max())
        x_base = 0.5 * (x_lo + x_hi)
        run = abs(x_tip - x_base)

        out["flame_height"][it] = z_tip
        out["flame_length"][it] = float(math.hypot(run, z_tip))
        out["flame_tilt"][it] = float(math.degrees(math.atan2(run, z_tip)))
        out["flame_depth"][it] = x_hi - x_lo
        out["x_tip"][it] = x_tip
        out["z_tip"][it] = z_tip
        out["x_base"][it] = x_base
    return out


def flame_geometry(sim, threshold, x, tau, z_base, window, x_limits=None):
    """Flame geometry from every vertical HRRPUV slice.

    Returns (times, per_slice, y_positions) where per_slice[key] has shape
    (n_slice, n_t).
    """
    slices = sorted(
        (s for s in sim.slices
         if s.quantity.name == "HRRPUV" and s.orientation == 2),
        key=lambda s: s.extent.y_start,
    )
    if not slices:
        raise RuntimeError("no vertical HRRPUV slices found")

    ref_times = None
    collected, y_pos = [], []
    for slc in slices:
        try:
            times, xc, zc, data = assemble_slice(slc)
        except RuntimeError as exc:
            print(f"    slice y={slc.extent.y_start:+.1f}: {exc}", file=sys.stderr)
            slc.clear_cache()
            continue
        if x_limits is not None:
            keep = (xc >= x_limits[0]) & (xc <= x_limits[1])
            if keep.any():
                xc, data = xc[keep], data[:, keep, :]
        res = flame_geometry_frames(times, xc, zc, data, threshold,
                                    front_at_times(x, tau, times),
                                    z_base=z_base, window=window)
        slc.clear_cache()
        del data
        if ref_times is None:
            ref_times = times
        collected.append((times, res))
        y_pos.append(float(slc.extent.y_start))

    if ref_times is None:
        raise RuntimeError("no vertical HRRPUV slice could be assembled")

    keys = ("flame_length", "flame_height", "flame_tilt", "flame_depth")
    per_slice = {k: np.full((len(collected), ref_times.size), np.nan)
                 for k in keys}
    for i, (times, res) in enumerate(collected):
        for k in keys:
            if times.size == ref_times.size and np.allclose(times, ref_times):
                per_slice[k][i] = res[k]
            else:
                # NaN-safe resample onto the reference clock
                good = np.isfinite(res[k])
                if good.sum() >= 2:
                    per_slice[k][i] = np.interp(ref_times, times[good],
                                                res[k][good],
                                                left=np.nan, right=np.nan)
    return ref_times, per_slice, y_pos


# ---------------------------------------------------------------------------
# output
# ---------------------------------------------------------------------------

def write_csv(path, header, columns):
    if not columns or np.size(columns[0]) == 0:
        with open(path, "w") as fh:
            fh.write(",".join(header) + "\n")
        return
    np.savetxt(path, np.column_stack(columns), delimiter=",",
               header=",".join(header), comments="", fmt="%.6g")


def parse_sim_id(name):
    """'c4_p10' -> ('c4', 'p10'). Neutral labels: these are TLS canopy fields."""
    parts = name.split("_")
    group = parts[0] if parts else ""
    plot = parts[1] if len(parts) > 1 else ""
    return group, plot


def find_chid(sim_dir):
    smv = glob.glob(os.path.join(sim_dir, "*.smv"))
    return os.path.splitext(os.path.basename(smv[0]))[0] if smv else None


def fds_completed(sim_dir, chid):
    """Did FDS shut down cleanly?

    The .end file is the usual marker, but a run stopped by the KILL control
    does not always leave one, so the .out banner is checked as well. Note
    that CATF renames the run: the outputs carry CHID_cat, which is what
    fdsreader reports, so both spellings are tried.
    """
    stems = [chid]
    if chid.endswith("_cat"):
        stems.append(chid[:-4])
    for stem in stems:
        if os.path.isfile(os.path.join(sim_dir, f"{stem}.end")):
            return True
    for stem in stems:
        out = os.path.join(sim_dir, f"{stem}.out")
        if not os.path.isfile(out):
            continue
        try:
            with open(out, "rb") as fh:
                fh.seek(max(0, os.path.getsize(out) - 4096))
                tail = fh.read().decode("utf-8", "replace")
        except OSError:
            continue
        if "completed successfully" in tail or "STOP: FDS completed" in tail:
            return True
    return False


def process_simulation(sim_dir, out_root, args):
    name = os.path.basename(os.path.normpath(sim_dir))
    out_dir = os.path.join(out_root, name)
    os.makedirs(out_dir, exist_ok=True)
    print(f"=== {name} ===", flush=True)

    sim = fdsreader.Simulation(sim_dir)
    chid = sim.chid
    completed = fds_completed(sim_dir, chid)
    print(f"  chid={chid}  FDS completed={completed}", flush=True)

    dx = float(np.diff(sim.meshes[0].coordinates["x"]).min())
    threshold = args.hrrpuv_threshold
    if threshold is None:
        threshold = min(200.0, 20.0 / dx)      # FDS UG Eq. 22.18
    print(f"  HRRPUV flame threshold {threshold:.0f} kW/m3 (dx={dx:.2f} m)",
          flush=True)

    # --- front and ROS -----------------------------------------------------
    t_ign = ignition_time(sim_dir, chid)
    print(f"  ignition at t={t_ign:.1f} s (earlier arrivals masked)", flush=True)
    x_nodes, y_nodes, arrival = assemble_arrival_map(sim)
    x, tau, tau_raw, tau_lead, tau_trail, x_ign_line = front_arrival_curve(
        x_nodes, y_nodes, arrival, args.y_band, stat=args.front_stat,
        t_min=t_ign, x_ignition=args.x_ignition)
    flat_frac = flattened_fraction(tau_raw, tau)
    print(f"  ignition line at x={x_ign_line:.2f} m "
          f"(backing fire upwind of it excluded)", flush=True)
    if np.isfinite(flat_frac) and flat_frac > 0.2:
        print(f"  WARNING monotonic pass lifted {flat_frac:.0%} of the front — "
              f"tau(x) is largely flat, so the ROS below is NOT trustworthy; "
              f"inspect front_position.png and arrival_map.png",
              file=sys.stderr)

    finite_tau = tau[np.isfinite(tau)]
    x_front_max = float(x[np.isfinite(tau)].max()) if finite_tau.size else float("nan")
    t_track_end = float(finite_tau.max()) if finite_tau.size else float("nan")

    t_enter = crossing_time(x, tau, ROI_X[0])
    t_exit = crossing_time(x, tau, ROI_X[1])
    t_enter_lead = crossing_time(x, tau_lead, ROI_X[0])
    ros_bulk_roi = ((ROI_X[1] - ROI_X[0]) / (t_exit - t_enter)
                    if np.isfinite(t_enter) and np.isfinite(t_exit)
                    and t_exit > t_enter else float("nan"))

    t_ros, x_front, ros = ros_timeseries(x, tau, dt=args.dt, smooth_s=args.smooth)
    x_seg, t_seg, ros_seg = segment_ros(x, tau, seg_len=args.segment)

    if np.isfinite(t_exit):
        status = "propagated"
    elif np.isfinite(t_enter):
        status = "stalled_in_roi"
    elif finite_tau.size:
        status = "no_roi_arrival"
    else:
        status = "no_spread"
    print(f"  status={status}  x_front_max={x_front_max:.1f} m", flush=True)

    # --- windows -----------------------------------------------------------
    steady_x0, steady_x1 = args.steady_x
    bounds = {
        "roi": (t_enter, t_exit),
        "steady": (crossing_time(x, tau, steady_x0),
                   crossing_time(x, tau, steady_x1)),
        "full": (t_ign, t_track_end),
    }

    # --- flame geometry ----------------------------------------------------
    flame_times = np.array([])
    per_slice, y_pos = {}, []
    try:
        flame_times, per_slice, y_pos = flame_geometry(
            sim, threshold, x, tau, z_base=args.z_base,
            window=args.flame_window, x_limits=(args.xmin, args.xmax))
        print(f"  flame geometry from {len(y_pos)} slice(s) at y={y_pos}",
              flush=True)
    except RuntimeError as exc:
        print(f"  WARNING flame geometry unavailable: {exc}", file=sys.stderr)

    with warnings.catch_warnings():
        # all-NaN columns are expected before ignition and after burnout
        warnings.simplefilter("ignore", RuntimeWarning)
        flame_mean = {k: np.nanmean(v, axis=0) if v.size else np.array([])
                      for k, v in per_slice.items()}

    # --- statistics --------------------------------------------------------
    stats = {m: {} for m in METRICS}
    with np.errstate(invalid="ignore"):
        for win in WINDOWS:
            b = bounds[win]
            stats["ros"][win] = describe(ros[mask_window(t_ros, b)])
            stats["ros_seg"][win] = describe(ros_seg[mask_window(t_seg, b)])
            fmask = mask_window(flame_times, b)
            for key in ("flame_length", "flame_height", "flame_tilt",
                        "flame_depth"):
                arr = per_slice.get(key)
                # pool every slice at every in-window time
                stats[key][win] = describe(arr[:, fmask] if arr is not None
                                           and arr.size else [])

    # --- per-simulation CSVs ----------------------------------------------
    write_csv(os.path.join(out_dir, "front_trajectory.csv"),
              ["x_m", f"t_{args.front_stat}_s", f"t_{args.front_stat}_raw_s",
               "t_leading_s", "t_trailing_s"],
              [x, tau, tau_raw, tau_lead, tau_trail])
    write_csv(os.path.join(out_dir, "ros_time.csv"),
              ["t_s", "x_front_m", "ros_m_per_s", "in_roi", "in_steady"],
              [t_ros, x_front, ros,
               mask_window(t_ros, bounds["roi"]).astype(int),
               mask_window(t_ros, bounds["steady"]).astype(int)])
    write_csv(os.path.join(out_dir, "ros_segments.csv"),
              ["x_mid_m", "t_mid_s", "ros_m_per_s", "in_roi", "in_steady"],
              [x_seg, t_seg, ros_seg,
               mask_window(t_seg, bounds["roi"]).astype(int),
               mask_window(t_seg, bounds["steady"]).astype(int)])
    if flame_times.size:
        cols = [flame_times,
                mask_window(flame_times, bounds["roi"]).astype(int),
                mask_window(flame_times, bounds["steady"]).astype(int)]
        head = ["t_s", "in_roi", "in_steady"]
        for key in ("flame_length", "flame_height", "flame_tilt", "flame_depth"):
            cols.append(flame_mean[key])
            head.append(f"{key}_mean")
        for i, y in enumerate(y_pos):
            cols.append(per_slice["flame_length"][i])
            head.append(f"flame_length_y{y:+.0f}")
            cols.append(per_slice["flame_height"][i])
            head.append(f"flame_height_y{y:+.0f}")
        write_csv(os.path.join(out_dir, "flame_geometry.csv"), head, cols)

    # --- metrics.json (the unit aggregate_metrics.py consumes) -------------
    group, plot = parse_sim_id(name)
    metrics = {
        "sim_id": name,
        "group": group,
        "plot": plot,
        "status": status,
        "fds_completed": bool(completed),
        "meta": {
            "chid": chid,
            "dx_m": dx,
            "hrrpuv_threshold_kw_m3": float(threshold),
            "front_stat": args.front_stat,
            "y_band_m": args.y_band,
            "ros_dt_s": args.dt,
            "ros_smooth_s": args.smooth,
            "ros_segment_m": args.segment,
            "flame_window_m": args.flame_window,
            "flame_z_base_m": args.z_base,
            "slice_y_m": y_pos,
            "roi_x_m": list(ROI_X),
            "steady_x_m": list(args.steady_x),
            "fdsreader_version": getattr(fdsreader, "__version__", "unknown"),
        },
        "scalars": {
            "t_ignition_s": float(t_ign),
            "t_enter_roi_s": float(t_enter),
            "t_exit_roi_s": float(t_exit),
            "t_enter_roi_leading_edge_s": float(t_enter_lead),
            "t_track_end_s": float(t_track_end),
            "x_front_max_m": x_front_max,
            "ros_bulk_roi_m_per_s": float(ros_bulk_roi),
            "x_ignition_line_m": float(x_ign_line),
            "front_flattened_fraction": float(flat_frac),
        },
        "windows": {w: {"t0": float(bounds[w][0]), "t1": float(bounds[w][1])}
                    for w in WINDOWS},
        "stats": stats,
    }
    with open(os.path.join(out_dir, "metrics.json"), "w") as fh:
        json.dump(metrics, fh, indent=2, allow_nan=True)

    # --- figures -----------------------------------------------------------
    if not args.no_figures:
        make_figures(out_dir, name, args, x, tau, tau_lead, tau_trail,
                     tau_raw, x_ign_line, x_nodes, y_nodes, arrival,
                     t_ros, ros, x_seg, ros_seg, flame_times, flame_mean,
                     per_slice, y_pos, bounds, ros_bulk_roi, threshold)

    print(f"  ROI: enter {t_enter:.1f} s  exit {t_exit:.1f} s  "
          f"bulk ROS {ros_bulk_roi:.3f} m/s", flush=True)
    return metrics


# ---------------------------------------------------------------------------
# per-simulation figures
# ---------------------------------------------------------------------------

def make_figures(out_dir, name, args, x, tau, tau_lead, tau_trail,
                 tau_raw, x_ign_line, x_nodes, y_nodes, arrival,
                 t_ros, ros, x_seg, ros_seg, flame_times, flame_mean,
                 per_slice, y_pos, bounds, ros_bulk_roi, threshold):
    label_stat = "median" if args.front_stat == "median" else "leading-edge"
    S = plotstyle

    # 1. front position ----------------------------------------------------
    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    ok = np.isfinite(tau_lead) & np.isfinite(tau_trail)
    if ok.any():
        ax.fill_betweenx(x[ok], tau_lead[ok], tau_trail[ok], alpha=0.18,
                         color=S.ACCENT, lw=0,
                         label="leading-trailing edge envelope")
    # the raw curve alongside the monotonic one: where they separate, the
    # monotonic pass is inventing the front rather than de-noising it
    ax.plot(tau_raw, x, color=S.MUTED, lw=1.0, ls=":",
            label=f"{label_stat}, before monotonic pass")
    ax.plot(tau, x, color=S.PRIMARY, lw=2,
            label=f"front position ({label_stat})")
    if np.isfinite(x_ign_line) and abs(x_ign_line) < 1e5:
        ax.axhline(x_ign_line, color=S.MUTED, ls="-.", lw=0.9, zorder=0)
        ax.annotate(f"ignition line x={x_ign_line:.0f} m (backing fire below)",
                    (0.005, x_ign_line), xycoords=("axes fraction", "data"),
                    va="top", fontsize=7, color=S.MUTED)
    S.roi_marker(ax, ROI_X, axis="y", label_fmt="ROI x={:+.0f} m")
    for tc, xc_, txt in ((bounds["roi"][0], ROI_X[0], "enter"),
                         (bounds["roi"][1], ROI_X[1], "exit")):
        if np.isfinite(tc):
            ax.plot([tc], [xc_], "o", ms=5, color=S.CRITICAL,
                    mec=S.SURFACE, mew=1.5, zorder=5)
            ax.annotate(f"{txt} {tc:.0f} s", (tc, xc_), fontsize=8,
                        color=S.INK_2, textcoords="offset points",
                        xytext=(7, -3))
    ax.set_xlabel("time (s)")
    ax.set_ylabel("front position x (m)")
    ax.set_title(f"{name} — fire front from FIRE ARRIVAL TIME (|y| ≤ {args.y_band:g} m)")
    ax.legend(loc="upper left")
    fig.savefig(os.path.join(out_dir, "front_position.png"))
    plt.close(fig)

    # 2. arrival map -------------------------------------------------------
    fig, ax = plt.subplots(figsize=(7.0, 4.0))
    masked = np.ma.masked_where(arrival >= UNBURNED, arrival)
    cmap = plt.get_cmap("magma_r").copy()
    cmap.set_bad(S.GRID)
    im = ax.pcolormesh(x_nodes, y_nodes, masked.T, cmap=cmap, shading="auto")
    cb = fig.colorbar(im, ax=ax, pad=0.02)
    cb.set_label("fire arrival time (s)", color=S.INK_2)
    cb.outline.set_visible(False)
    for xr in ROI_X:
        ax.axvline(xr, color=S.CRITICAL, ls="--", lw=1.0)
    ax.axhline(-args.y_band, color=S.CRITICAL, ls=":", lw=0.9)
    ax.axhline(args.y_band, color=S.CRITICAL, ls=":", lw=0.9)
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_title(f"{name} — arrival map (grey = never burned; "
                 f"dashed = ROI, dotted = front y-band)")
    ax.grid(False)
    fig.savefig(os.path.join(out_dir, "arrival_map.png"))
    plt.close(fig)

    # 3. ROS ---------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    if t_ros.size:
        ax.plot(t_ros, ros, color=S.PRIMARY, lw=1.4,
                label=f"instantaneous ({args.smooth:g} s smoothing)")
    ax.set_xlabel("time (s)")
    ax.set_ylabel("rate of spread (m/s)")
    for tc, lbl in zip(bounds["roi"], ("enters ROI", "exits ROI")):
        if np.isfinite(tc):
            ax.axvline(tc, color=S.CRITICAL, ls="--", lw=1.0, zorder=0)
            ax.annotate(lbl, (tc, 0.98), xycoords=("data", "axes fraction"),
                        rotation=90, va="top", ha="right", fontsize=7,
                        color=S.CRITICAL)
    if np.isfinite(ros_bulk_roi):
        ax.axhline(ros_bulk_roi, color=S.INK_2, ls=":", lw=1.1,
                   label=f"ROI bulk ROS = {ros_bulk_roi:.3f} m/s")
    ax.set_title(f"{name} — rate of spread ({label_stat} front)")
    ax.legend(loc="upper left")
    fig.savefig(os.path.join(out_dir, "ros.png"))
    plt.close(fig)

    # 4. segment ROS vs x --------------------------------------------------
    if x_seg.size:
        fig, ax = plt.subplots(figsize=(6.5, 3.6))
        ax.step(x_seg, ros_seg, where="mid", color=S.PRIMARY, lw=1.6)
        ax.plot(x_seg, ros_seg, "o", ms=4, color=S.PRIMARY,
                mec=S.SURFACE, mew=1.2)
        S.roi_marker(ax, ROI_X, axis="x", label_fmt="ROI x={:+.0f} m")
        ax.set_xlabel("x (m)")
        ax.set_ylabel("chord ROS (m/s)")
        ax.set_title(f"{name} — ROS over {args.segment:g} m segments "
                     f"(near-independent samples)")
        fig.savefig(os.path.join(out_dir, "ros_segments.png"))
        plt.close(fig)

    # 5. flame geometry ----------------------------------------------------
    if flame_times.size and per_slice:
        fig, axes = plt.subplots(2, 1, figsize=(6.5, 6.0), sharex=True)
        ax = axes[0]
        for i, y in enumerate(y_pos):
            ax.plot(flame_times, per_slice["flame_length"][i], lw=0.7,
                    alpha=0.35, color=S.MUTED,
                    label="individual slices" if i == 0 else None)
        ax.plot(flame_times, flame_mean["flame_length"], color=S.PRIMARY, lw=2,
                label="flame length $L$ (slice mean)")
        ax.plot(flame_times, flame_mean["flame_height"], color=S.ACCENT, lw=2,
                label="flame height $h_F$ (slice mean)")
        ax.set_ylabel("length / height (m)")
        ax.set_title(f"{name} — head-fire flame geometry "
                     f"(HRRPUV ≥ {threshold:.0f} kW/m³)")
        ax.legend(loc="upper left")

        ax = axes[1]
        ax.plot(flame_times, flame_mean["flame_tilt"], color=S.SERIES[2], lw=1.8,
                label="tilt from vertical")
        ax.set_ylabel("tilt (deg)")
        ax.set_xlabel("time (s)")
        ax.legend(loc="upper left")

        for a in axes:
            for tc, lbl in zip(bounds["roi"], ("enters ROI", "exits ROI")):
                if np.isfinite(tc):
                    a.axvline(tc, color=S.CRITICAL, ls="--", lw=1.0, zorder=0)
        fig.savefig(os.path.join(out_dir, "flame_geometry.png"))
        plt.close(fig)


# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Per-simulation ROS and flame-geometry postprocessing.")
    ap.add_argument("sim_dirs", nargs="*", help="simulation directories")
    ap.add_argument("--root", help="process every subdirectory of this dir")
    ap.add_argument("--out", default="postprocess",
                    help="output directory (default: ./postprocess)")
    ap.add_argument("--hrrpuv-threshold", type=float, default=None,
                    help="flame HRRPUV threshold in kW/m3 "
                         "(default: min(200, 20/dx) per FDS UG Eq. 22.18)")
    ap.add_argument("--y-band", type=float, default=10.0,
                    help="half-width of the y band for the front (default 10 m)")
    ap.add_argument("--x-ignition", type=float, default=None,
                    help="x (m) of the ignition line; the backing fire upwind "
                         "of it is excluded from the front. Default: detected "
                         "from the earliest arrival")
    ap.add_argument("--front-stat", choices=("median", "min"), default="median",
                    help="front statistic across the y band "
                         "(median = mean fireline, min = leading edge)")
    ap.add_argument("--dt", type=float, default=0.5,
                    help="resampling interval (s) for X_front(t)")
    ap.add_argument("--smooth", type=float, default=5.0,
                    help="smoothing window (s) for the ROS derivative")
    ap.add_argument("--segment", type=float, default=2.5,
                    help="segment length (m) for the chord ROS estimator. "
                         "2.5 m gives 8 samples across the 20 m ROI and 20 "
                         "across the default steady window; raising it trades "
                         "sample count for independence")
    ap.add_argument("--steady-x", type=float, nargs=2, default=(-30.0, 20.0),
                    metavar=("X0", "X1"),
                    help="x bounds of the 'steady' analysis window "
                         "(default -30 20: past the entrance jet, short of "
                         "the kill devices at x=+30)")
    ap.add_argument("--flame-window", type=float, default=15.0,
                    help="half-width (m) of the band around the front searched "
                         "for the head flame (default 15)")
    ap.add_argument("--z-base", type=float, default=0.5,
                    help="height (m) below which flaming cells count as the "
                         "flame base (default 0.5)")
    ap.add_argument("--xmin", type=float, default=None,
                    help="restrict the flame search to x >= xmin")
    ap.add_argument("--xmax", type=float, default=None,
                    help="restrict the flame search to x <= xmax")
    ap.add_argument("--no-figures", action="store_true",
                    help="skip figures; write tables only")
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

    plotstyle.apply()
    os.makedirs(args.out, exist_ok=True)

    failures = 0
    for d in sim_dirs:
        try:
            process_simulation(d, args.out, args)
        except Exception as exc:
            failures += 1
            name = os.path.basename(os.path.normpath(d))
            print(f"  ERROR processing {name}: {exc}", file=sys.stderr)
            traceback.print_exc()
            fail_dir = os.path.join(args.out, name)
            os.makedirs(fail_dir, exist_ok=True)
            with open(os.path.join(fail_dir, "metrics.json"), "w") as fh:
                json.dump({"sim_id": name,
                           "group": parse_sim_id(name)[0],
                           "plot": parse_sim_id(name)[1],
                           "status": "error",
                           "error": f"{type(exc).__name__}: {exc}"},
                          fh, indent=2)

    print(f"Done. Outputs in {os.path.abspath(args.out)}"
          f"{f' ({failures} failed)' if failures else ''}")
    return 1 if failures and len(sim_dirs) == failures else 0


if __name__ == "__main__":
    sys.exit(main())
