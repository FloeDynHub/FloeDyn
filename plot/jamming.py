"""Jamming analysis of FloeDyn channel/funnel ensembles.

For each output .h5, computes the average floe speed inside a y-window over time, then:
  - jamming time  : cumulated SIMULATED time where that average is below --vmin while the system still
                    holds real ice to discharge (the "blocked group" observable), over a FIXED horizon
                    (the run, or --horizon) — a restricted-mean measure, comparable across runs even when
                    some never unblock;
  - drain time    : last time the system still had ice to discharge or was discharging;
  - outcome       : drained (emptied before the horizon) / arrested (no outflow at all over the last
                    --arrest-frac of the horizon, ice still held back: permanent clog) / unresolved
                    (neither: still intermittent at the horizon — the only truly censored case);
  - frac_evacuated: fraction of the pack past the outlet at the horizon (bounded, no infinity issue).

Outputs (all optional, combinable — designed to run headless on the cluster and be driven by a UI):
  --csv PATH        one row per simulation (machine-readable)
  --configs-csv P   one row per configuration (width, shape): P(arrest), outcome counts, T_jam and
                    evacuated-fraction statistics over the draws
  --plots DIR       one PNG per simulation (avg speed curve + jammed intervals highlighted)
  --grid PATH       small multiples rows = widths, cols = shapes: ensemble beam of the avg speed over the
                    draws (median + interquartile band, each draw as a thin grey trace). Once a draw has
                    drained, its (empty-window) speed is dropped, not counted as 0.
  --states PATH     same layout: fraction of the draws flowing / jammed / drained at each time
  --discharge PATH  one panel per width: evacuated fraction vs time, median over the draws per shape
  --summary PATH    vs width, one series per shape: T_jam, P(arrest), evacuated fraction; each draw's
                    marker shows its outcome
  --table PATH      dashboard: width x shape tables of P(arrest), median evacuated fraction, median T_jam
  --scaling PATH    P(arrest) and evacuated fraction vs outlet width / floe size MEASURED on the actual floes
                    (equivalent diameter and min width): does the size ratio explain clogging whatever the shape?
  Floe speeds come from the displacement between outputs (--speed disp), not the recorded velocity.
  Runs whose output does not start at t=0 (overwritten by a restart) or ends early are flagged
  ('complete' = 0 in jam.csv) and left out of every ensemble statistic.
  Names must expose wXX, the shape (pN|real) and fileN (e.g. in_chenalv4_w17_p4_file2) for the
  ensemble figures.

Usage examples:
  python3 plot/jamming.py io/outputs/out_chenalv4_*.h5 --out-dir jam_out
  python3 plot/jamming.py --manifest runs.csv --out-dir jam_out     # UI-driven: file,label per line
  python3 plot/jamming.py --manifest runs.csv --out-dir jam_14d --horizon 1209600   # horizon check
"""
import argparse
import csv
import math
import os
import re
import sys

import h5py
import numpy as np

DAY = 86400.0


# ---------------------------------------------------------------- data extraction

def load_series(path, ymin, ymax, outlet, skip_first, horizon=None, speed_mode="disp"):
    """Return (times, avg_speed, count, n_floes, exited) for floes inside the y-window.
    Obstacles are excluded: the first `skip_first` shapes if given, else auto-detected as the
    floes whose speed is exactly 0 over the whole run (walls never move; even frozen floes
    carry a small phantom drag speed, so this is unambiguous).
    speed_mode 'disp' (default): speed = displacement of the floe centre to the next output / elapsed
    time. The recorded velocity of a floe at rest is not 0 but one drag step (a*dt), which grows with the
    time step once a pack freezes and was seen to hover right at vmin in clogged runs, making the jammed
    state flicker; a blocked floe does not move, whatever dt. 'recorded': the velocity stored in the file.
    `horizon` (s) truncates the series, to test how conclusions depend on the analysis horizon."""
    with h5py.File(path, "r") as f:
        times = f["time"][:]
        nt = len(times) if horizon is None else int(np.searchsorted(times, horizon, side="right"))
        times = times[:nt]
        states = f["floe_states"][:nt]  # (T, N, 11): 0,1=pos 3,4=speed
    y = states[:, :, 1]
    recorded = np.hypot(states[:, :, 3], states[:, :, 4])
    if skip_first is None:
        moving = (recorded.max(axis=0) > 0)       # auto: walls have exactly 0 speed forever
    else:
        moving = np.ones(recorded.shape[1], bool)
        moving[:skip_first] = False
    if speed_mode == "disp" and len(times) > 1:
        # forward difference: speed[k] = motion between outputs k and k+1 (the interval the jam test uses)
        d = np.hypot(np.diff(states[:, :, 0], axis=0), np.diff(y, axis=0)) / np.diff(times)[:, None]
        speed = np.vstack([d, d[-1:]])
    else:
        speed = recorded
    inside = (y > ymin) & (y < ymax) & moving[None, :]
    count = inside.sum(axis=1)
    total = np.where(inside, speed, 0.0).sum(axis=1)
    avg = np.divide(total, count, out=np.zeros_like(total), where=count > 0)
    # Cumulative discharge: floes that have passed the outlet (y <= outlet). Monotone (below the outlet
    # they keep drifting away), so a rise between two snapshots = floes actually leaving = live outflow.
    exited = ((y <= outlet) & moving[None, :]).sum(axis=1)
    return times, avg, count, int(moving.sum()), exited


def _convex_hull(p):
    """Monotone-chain convex hull (counter-clockwise), for the min width of non-convex outlines."""
    p = sorted(map(tuple, p))
    if len(p) < 3:
        return np.array(p)
    def half(pts):
        h = []
        for q in pts:
            while len(h) >= 2 and ((h[-1][0] - h[-2][0]) * (q[1] - h[-2][1])
                                   - (h[-1][1] - h[-2][1]) * (q[0] - h[-2][0])) <= 0:
                h.pop()
            h.append(q)
        return h
    lower, upper = half(p), half(p[::-1])
    return np.array(lower[:-1] + upper[:-1])


def floe_sizes(path):
    """Measured size of the pack's floes (obstacles excluded), from the outlines stored in the output:
    returns (area-weighted mean equivalent diameter, area-weighted mean min width), in m.
      equivalent diameter = diameter of the disc of same area;
      min width           = smallest caliper width (narrowest way through, what an arch has to span).
    Area-weighted, so the big floes — the ones that build arches — dominate. Measuring the actual
    floes makes the comparison independent of how each shape was sized at generation."""
    deq, wmin, area = [], [], []
    with h5py.File(path, "r") as f:
        if "floe_shapes" not in f:
            return float("nan"), float("nan")
        for ds in f["floe_shapes"].values():
            if int(ds.attrs.get("obstacle", 0)):
                continue
            p = ds[:]
            x, y = p[:, 0], p[:, 1]
            a = 0.5 * abs(np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1)))
            h = _convex_hull(p)
            e = np.diff(np.vstack([h, h[:1]]), axis=0)
            nrm = np.c_[-e[:, 1], e[:, 0]] / np.hypot(e[:, 0], e[:, 1])[:, None]
            area.append(a)
            deq.append(2 * math.sqrt(a / math.pi))
            wmin.append(np.ptp(h @ nrm.T, axis=0).min())
    if not area:
        return float("nan"), float("nan")
    return float(np.average(deq, weights=area)), float(np.average(wmin, weights=area))


def jamming_metrics(times, avg, count, exited, n_floes, vmin, min_floes, drain_frac, arrest_frac):
    """Jamming time = cumulated time the channel is BLOCKED: the group in the window is slow
    (avg < vmin) WHILE there is still substantial ice left to discharge, over the working phase only.
      remaining(t) = n_floes - exited(t) : floes still in the system (above the outlet, not evacuated).
                     This — NOT the window population — is what tells a flow-blocking clog (lots of ice
                     held back, often backed up ABOVE the window) apart from a few residual floes lodged
                     in the funnel corners after the reservoir emptied.
      has_ice      : remaining >= max(min_floes, drain_frac * n_floes) — real load still to evacuate.
      discharging  : floes still cross the outlet this interval (exited rises) = live outflow.
      t_start      : first time the group reached flow speed (avg >= vmin) — excludes the initial
                     at-rest / acceleration phase, which is not a jam.
      drain_t      : LAST time the system still had ice to discharge OR was discharging.
    Regimes: live outflow -> not drained; no outflow + little ice -> drained (residual tail dropped);
    no outflow + much ice held back -> stays jammed (the clog we must count, wherever it sits).
    Outcome at the horizon (last sample):
      drained    : the system emptied before the horizon (T_jam is exact and final);
      arrested   : not drained, and not a single floe left over the last `arrest_frac` of the horizon
                   (permanent clog: T_jam ~ horizon - t_arrest, it would grow with a longer run);
      unresolved : still intermittent at the horizon (the only truly censored case)."""
    dt = np.diff(times)
    thresh = max(float(min_floes), drain_frac * n_floes)   # ice (floes) that must remain to count as loaded
    remaining = n_floes - exited                           # floes not yet discharged (monotone non-increasing)
    has_ice = remaining >= thresh
    discharging = np.zeros(count.shape, bool)
    discharging[:-1] = np.diff(exited) > 0                  # a floe left through the outlet during interval k
    active = has_ice | discharging                         # channel not yet drained
    last_active = np.nonzero(active)[0]
    drain_t = float(times[last_active[-1]]) if len(last_active) else float(times[-1])
    moved = np.nonzero(avg >= vmin)[0]
    t_start = float(times[moved[0]]) if len(moved) else 0.0
    jammed = ((avg[:-1] < vmin) & has_ice[:-1]
              & (times[:-1] >= t_start) & (times[:-1] <= drain_t))
    t_jam = float(dt[jammed].sum())
    horizon = min(drain_t, float(times[-1])) - t_start
    n_jams = int(np.sum(np.diff(np.r_[0, jammed.astype(int)]) == 1))   # distinct jammed episodes

    t_end = float(times[-1])
    outflows = np.nonzero(np.diff(exited) > 0)[0]
    t_last_out = float(times[outflows[-1] + 1]) if len(outflows) else float(times[0])
    if not active[-1]:
        outcome = "drained"
    elif t_last_out <= t_end - arrest_frac * (t_end - float(times[0])):
        outcome = "arrested"
    else:
        outcome = "unresolved"
    return {
        "t_jam": t_jam,
        "t_end": t_end,
        "t_start": t_start,
        "t_drain": drain_t,
        "jam_frac": t_jam / horizon if horizon > 0 else 0.0,
        "n_floes": n_floes,
        "drain_thresh": thresh,
        "outcome": outcome,
        "t_arrest": t_last_out if outcome == "arrested" else float("nan"),
        "frac_evac": float(exited[-1]) / n_floes if n_floes else float("nan"),
        "n_jams": n_jams,
        "jam_mask": np.append(jammed, False),  # per-sample, for plotting
        "remaining": remaining, "exited": exited, "count": count,  # for --trace
    }


# ---------------------------------------------------------------- naming / ensemble structure

def parse_name(name):
    """Extract (width, shape, draw) from a file/label name, each independently and best-effort:
    'in_chenalv4_w17_p4_file2' -> (17, 'p4', 2); 'out_chenal_1_4_real_w25' -> (25, 'real', None)."""
    w = re.search(r"w(\d+)", name)
    s = re.search(r"(p\d+|real)", name)
    d = re.search(r"file(\d+)", name)
    return (int(w.group(1)) if w else None,
            s.group(1) if s else None,
            int(d.group(1)) if d else None)


def shape_key(s):
    """Physical order: polygons by number of sides (p3, p4, p6, p25), then real outlines."""
    m = re.fullmatch(r"p(\d+)", s or "")
    return (0, int(m.group(1))) if m else (1, s or "")


# ---------------------------------------------------------------- styling

# Polygons are ordinal (more sides = rounder), so they share one blue ramp light -> dark; the realistic
# outlines take a distinct hue. Every shape also carries its own marker (secondary encoding).
SHAPE_COLORS = {"p3": "#86b6ef", "p4": "#3987e5", "p6": "#1c5cab", "p25": "#0d366b", "real": "#eb6834"}
EXTRA_COLORS = ["#1baf7a", "#4a3aa7", "#e87ba4"]   # unexpected shapes, in fixed order
STATE_COLORS = {"flowing": "#86b6ef", "jammed": "#d03b3b", "drained": "#c9c7c1"}
INK, INK_2, SURFACE = "#0b0b0b", "#52514e", "#fcfcfb"


def shape_color(shape, shapes):
    if shape in SHAPE_COLORS:
        return SHAPE_COLORS[shape]
    extras = [s for s in shapes if s not in SHAPE_COLORS]
    return EXTRA_COLORS[extras.index(shape) % len(EXTRA_COLORS)]


def shape_marker(shape):
    """Matplotlib marker echoing the floe shape: p3=triangle, p4=square, p6=hexagon, p25=circle,
    real=irregular blob."""
    from matplotlib.path import Path
    if shape == "p3":
        return "^"
    if shape == "p4":
        return "s"
    if shape == "p6":
        return "h"
    if shape == "p25":
        return "o"
    # "real" (or anything else): a fixed irregular polygon, like a realistic floe outline
    rng = np.random.default_rng(4)
    ang = np.linspace(0, 2 * np.pi, 9, endpoint=False)
    r = 0.65 + 0.5 * rng.random(9)
    pts = np.c_[r * np.cos(ang), r * np.sin(ang)]
    return Path(np.vstack([pts, pts[0]]))


def setup_mpl():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
        "axes.edgecolor": "#b5b3ad", "axes.labelcolor": INK_2, "xtick.color": INK_2, "ytick.color": INK_2,
        "text.color": INK, "grid.color": "#e4e2dd", "grid.linewidth": 0.6,
        "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False,
    })
    return plt


# ---------------------------------------------------------------- ensemble helpers

def rolling_nanmean(a, w):
    """Centred rolling mean over `w` samples ignoring NaN (edge-aware: the window shrinks at the ends).
    Samples that are NaN in the input stay NaN, so a drained draw's masked tail is not smeared back in.
    Pure numpy (cumsum) so no pandas/scipy dependency on the cluster venv."""
    a = np.asarray(a, float)
    n = a.size
    if w <= 1 or n == 0:
        return a.copy()
    valid = ~np.isnan(a)
    c = np.cumsum(np.insert(np.where(valid, a, 0.0), 0, 0.0))
    cv = np.cumsum(np.insert(valid.astype(float), 0, 0.0))
    half = w // 2
    idx = np.arange(n)
    lo = np.clip(idx - half, 0, n)
    hi = np.clip(idx + half + 1, 0, n)
    cnt = cv[hi] - cv[lo]
    out = np.divide(c[hi] - c[lo], cnt, out=np.full(n, np.nan), where=cnt > 0)
    out[~valid] = np.nan
    return out


def ensemble_band(stack, min_draws):
    """Per-time median and interquartile range over the draws (rows), ignoring NaN; NaN where fewer
    than `min_draws` draws are still defined (a beam resting on one or two runs would mislead)."""
    n_valid = (~np.isnan(stack)).sum(axis=0)
    ok = n_valid >= max(1, min_draws)
    med, q25, q75 = (np.full(stack.shape[1], np.nan) for _ in range(3))
    if ok.any():
        med[ok] = np.nanmedian(stack[:, ok], axis=0)
        q25[ok], q75[ok] = np.nanpercentile(stack[:, ok], [25, 75], axis=0)
    return med, q25, q75


def group_by_config(rows):
    """{(width, shape): [rows]} for complete rows whose name exposes both; widths and shapes in plotting
    order. Incomplete runs (output not starting at t=0) are left out of every ensemble statistic."""
    groups = {}
    for r in rows:
        if r["w"] is not None and r["shape"] is not None and r["complete"]:
            groups.setdefault((r["w"], r["shape"]), []).append(r)
    widths = sorted({k[0] for k in groups})
    shapes = sorted({k[1] for k in groups}, key=shape_key)
    return groups, widths, shapes


def time_grid(rows):
    """Common time axis for the ensemble figures: regular, from 0 to the longest complete run, with the
    typical output step. Runs are put on it by time, not by sample index, so a run with missing or
    extra outputs cannot shift or cut the others."""
    ok = [r for r in rows if r["complete"] and len(r["times"]) > 1]
    if not ok:
        return None
    step = float(np.median([np.median(np.diff(r["times"])) for r in ok]))
    return np.arange(0.0, max(r["times"][-1] for r in ok) + step / 2, step)


def on_grid(r, key, grid):
    """Series `key` of run r sampled on `grid`: nearest output, NaN where the run has no output within
    3/4 of a grid step (a hole in the file, or after its end) — no data is invented."""
    t = r["times"]
    v = np.asarray(r[key], float)
    j = np.clip(np.searchsorted(t, grid), 1, len(t) - 1)
    pick = np.where(np.abs(t[j] - grid) < np.abs(t[j - 1] - grid), j, j - 1)
    out = v[pick]
    out[np.abs(t[pick] - grid) > 0.75 * (grid[1] - grid[0])] = np.nan
    return out


def stacked(rs, key, grid):
    return np.vstack([on_grid(r, key, grid) for r in rs])


def outcome_tally(rs):
    n = len(rs)
    k = {o: sum(r["outcome"] == o for r in rs) for o in ("arrested", "drained", "unresolved")}
    txt = f"arrested {k['arrested']}/{n} · drained {k['drained']}/{n}"
    if k["unresolved"]:
        txt += f" · unres. {k['unresolved']}"
    return txt


def wilson(k, n, z=1.96):
    """Wilson 95% interval for a binomial proportion (well-behaved at k = 0 or n, unlike +-sqrt(p(1-p)/n))."""
    if n == 0:
        return float("nan"), float("nan")
    p = k / n
    den = 1 + z * z / n
    mid = (p + z * z / (2 * n)) / den
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return max(0.0, mid - half), min(1.0, mid + half)


# ---------------------------------------------------------------- plots

def plot_single(times, avg, m, vmin, label, path):
    plt = setup_mpl()
    td = times / DAY
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(td, avg, lw=0.9, color=SHAPE_COLORS["p4"])
    ax.axhline(vmin, color=STATE_COLORS["jammed"], ls="--", lw=0.8, label=f"vmin = {vmin} m/s")
    ax.fill_between(td, 0, max(avg.max(), vmin) * 1.05, where=m["jam_mask"], color=STATE_COLORS["jammed"],
                    alpha=0.12, lw=0, label=f"jammed  (T_jam = {m['t_jam'] / DAY:.1f} d)")
    if m["t_start"] > times[0]:
        ax.axvline(m["t_start"] / DAY, color="0.5", ls=":", lw=1, label="flow start")
    if m["t_drain"] < times[-1]:
        ax.axvline(m["t_drain"] / DAY, color="0.3", ls=":", lw=1, label="drained")
    ax.set_xlabel("time (days)"); ax.set_ylabel("avg speed in window (m/s)")
    ax.set_title(f"{label}  —  {m['outcome']}, {100 * m['frac_evac']:.0f}% evacuated", fontsize=10)
    ax.legend(fontsize=8); ax.grid(alpha=0.6)
    fig.savefig(path, dpi=110, bbox_inches="tight"); plt.close(fig)


def plot_grid(rows, grid, path, vmin, smooth, min_draws):
    """Small multiples rows = widths, cols = shapes. Each panel is the ensemble of draws: every draw as a
    thin grey trace (keeps the intermittency, which averaging would smear since jams do not happen at the
    same time in each draw), the median as the line and the interquartile range as the band. A drained
    draw's empty window is NaN, so it leaves the statistics instead of pulling them to 0."""
    groups, widths, shapes = group_by_config(rows)
    if not groups:
        print("grid: names do not expose wXX/shape -> skipped"); return
    plt = setup_mpl()
    fig, axes = plt.subplots(len(widths), len(shapes), figsize=(3.3 * len(shapes), 2.5 * len(widths)),
                             sharex=True, sharey=True, squeeze=False, constrained_layout=True)
    for i, w in enumerate(widths):
        for j, s in enumerate(shapes):
            ax = axes[i][j]
            ax.axhline(vmin, color=INK_2, ls="--", lw=0.7)
            ax.grid(alpha=0.6); ax.tick_params(labelsize=8); ax.margins(x=0)
            if i == 0: ax.set_title(s, fontsize=11)
            if j == 0: ax.set_ylabel(f"w{w}\navg speed (m/s)", fontsize=9)
            if i == len(widths) - 1: ax.set_xlabel("time (days)", fontsize=9)
            rs = groups.get((w, s))
            if not rs:
                ax.text(0.5, 0.5, "no run", transform=ax.transAxes, ha="center", color=INK_2, fontsize=8)
                continue
            stack = stacked(rs, "avg_plot", grid)
            td = grid / DAY
            stack = np.vstack([rolling_nanmean(a, smooth) for a in stack])
            for a in stack:
                ax.plot(td, a, lw=0.5, color="0.55", alpha=0.35)
            med, q25, q75 = ensemble_band(stack, min_draws)
            col = shape_color(s, shapes)
            ax.fill_between(td, q25, q75, color=col, alpha=0.30, lw=0)
            ax.plot(td, med, lw=1.6, color=col)
            ax.text(0.98, 0.95, outcome_tally(rs), transform=ax.transAxes, ha="right", va="top",
                    fontsize=7, color=INK_2)
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    fig.legend([Line2D([], [], color=INK_2, lw=1.6), Patch(color=INK_2, alpha=0.3),
                Line2D([], [], color="0.55", lw=0.6), Line2D([], [], color=INK_2, ls="--", lw=0.7)],
               ["median of the draws", "interquartile range", "one draw", f"vmin = {vmin} m/s"],
               loc="outside upper center", ncol=4, fontsize=9,
               title=f"smoothed over {smooth} samples; drained draws leave the statistics "
                     f"(beam shown while >= {min_draws} draws remain)")
    fig.savefig(path, dpi=160, bbox_inches="tight"); plt.close(fig)
    print(f"Saved {path}")


def plot_states(rows, grid, path):
    """Same layout as the grid: fraction of the draws flowing / jammed / drained at each time.
    It is the 'probability of being blocked' along the run, and shows where the draws disagree."""
    groups, widths, shapes = group_by_config(rows)
    if not groups:
        print("states: names do not expose wXX/shape -> skipped"); return
    plt = setup_mpl()
    fig, axes = plt.subplots(len(widths), len(shapes), figsize=(3.3 * len(shapes), 2.3 * len(widths)),
                             sharex=True, sharey=True, squeeze=False, constrained_layout=True)
    names = ["flowing", "jammed", "drained"]
    for i, w in enumerate(widths):
        for j, s in enumerate(shapes):
            ax = axes[i][j]
            ax.tick_params(labelsize=8); ax.margins(x=0); ax.set_ylim(0, 1)
            if i == 0: ax.set_title(s, fontsize=11)
            if j == 0: ax.set_ylabel(f"w{w}\nfraction of draws", fontsize=9)
            if i == len(widths) - 1: ax.set_xlabel("time (days)", fontsize=9)
            rs = groups.get((w, s))
            if not rs:
                continue
            st = stacked(rs, "state", grid)
            n_def = (~np.isnan(st)).sum(axis=0)       # draws with an output at that time
            fr = [np.divide((st == k).sum(axis=0), n_def, out=np.full(len(grid), np.nan), where=n_def > 0)
                  for k in range(3)]
            # no edge lines: a draw flipping state at every output drew one white edge per flip, which at
            # this density painted whole bands white
            ax.stackplot(grid / DAY, *fr, colors=[STATE_COLORS[n] for n in names], lw=0)
    from matplotlib.patches import Patch
    fig.legend([Patch(color=STATE_COLORS[n]) for n in names], names, loc="outside upper center",
               ncol=3, fontsize=9, title="state of each draw at time t")
    fig.savefig(path, dpi=160, bbox_inches="tight"); plt.close(fig)
    print(f"Saved {path}")


def plot_discharge(rows, grid, path, min_draws):
    """One panel per width: evacuated fraction (floes past the outlet / pack) vs time, the median over the
    draws for each shape. Monotone curves, so the medians overlay legibly (five overlapping bands would
    not: the spread over the draws is in the summary, panel c); a plateau below 1 is a permanent clog,
    its height the ice that got through."""
    groups, widths, shapes = group_by_config(rows)
    if not groups:
        print("discharge: names do not expose wXX/shape -> skipped"); return
    plt = setup_mpl()
    fig, axes = plt.subplots(1, len(widths), figsize=(4.2 * len(widths), 3.6), sharey=True,
                             squeeze=False, constrained_layout=True)
    for i, w in enumerate(widths):
        ax = axes[0][i]
        ax.set_title(f"w{w}", fontsize=11); ax.grid(alpha=0.6); ax.margins(x=0)
        ax.set_xlabel("time (days)", fontsize=9); ax.set_ylim(0, 1.02)
        if i == 0: ax.set_ylabel("evacuated fraction of the pack", fontsize=9)
        for s in shapes:
            rs = groups.get((w, s))
            if not rs:
                continue
            stack = stacked(rs, "evac", grid)
            med, _, _ = ensemble_band(stack, min_draws)
            col = shape_color(s, shapes)
            td = grid / DAY
            ax.plot(td, med, lw=1.8, color=col, label=s)
            k = np.nonzero(~np.isnan(med))[0]
            if len(k):   # shape marker at the end of each median: identity is not colour-alone
                ax.scatter([td[k[-1]]], [med[k[-1]]], marker=shape_marker(s), s=60, color=col,
                           edgecolors=SURFACE, linewidths=1.5, zorder=3, clip_on=False)
    from matplotlib.lines import Line2D
    fig.legend([Line2D([], [], color=shape_color(s, shapes), marker=shape_marker(s), lw=1.6, ms=8)
                for s in shapes], shapes, loc="outside upper center", ncol=len(shapes), fontsize=9,
               title="shape — median of the draws (spread over the draws: summary, panel c)")
    fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)
    print(f"Saved {path}")


def plot_summary(rows, path):
    """Three views vs the biggest-floe width, one series per shape (dodged around each width):
      (a) jamming time over the horizon, (b) probability of permanent arrest k/n with its Wilson 95%
      interval, (c) evacuated fraction at the horizon. In (a) and (c) every draw is a small marker whose
      style tells its outcome (filled = drained, x = arrested, hollow = unresolved/censored); the large
      marker + bar are the median and the interquartile range."""
    groups, widths, shapes = group_by_config(rows)
    if not groups:
        print("summary: names do not expose wXX/shape -> skipped"); return
    plt = setup_mpl()
    from matplotlib.lines import Line2D
    fig, axes = plt.subplots(1, 3, figsize=(16, 5), constrained_layout=True)
    offs = np.linspace(-0.36, 0.36, len(shapes)) if len(shapes) > 1 else [0.0]
    panels = [("t_jam_d", "jamming time over the horizon (days)"),
              ("p_arrest", "P(permanent arrest)"),
              ("frac_evac", "evacuated fraction at the horizon")]
    for ax, (key, ylabel) in zip(axes, panels):
        for s, off in zip(shapes, offs):
            col, mk = shape_color(s, shapes), shape_marker(s)
            xs_line, ys_line = [], []
            for iw, w in enumerate(widths):
                rs = groups.get((w, s))
                if not rs:
                    continue
                x = iw + off
                if key == "p_arrest":
                    k = sum(r["outcome"] == "arrested" for r in rs)
                    lo, hi = wilson(k, len(rs))
                    ax.plot([x, x], [lo, hi], color=col, lw=2, solid_capstyle="round")
                    ax.scatter([x], [k / len(rs)], marker=mk, s=70, color=col, edgecolors=SURFACE,
                               linewidths=1.5, zorder=3)
                    xs_line.append(x); ys_line.append(k / len(rs))
                    continue
                vals = np.array([r["t_jam"] / DAY if key == "t_jam_d" else r["frac_evac"] for r in rs])
                jit = np.linspace(-0.05, 0.05, len(rs))
                for r, v, dj in zip(rs, vals, jit):
                    if r["outcome"] == "arrested":
                        ax.scatter([x + dj], [v], marker="x", s=16, color=col, linewidths=0.9, alpha=0.5)
                    else:
                        filled = r["outcome"] == "drained"
                        ax.scatter([x + dj], [v], marker=mk, s=14, alpha=0.5, linewidths=0.8,
                                   facecolors=col if filled else "none", edgecolors=col)
                q25, med, q75 = np.percentile(vals, [25, 50, 75])
                ax.plot([x, x], [q25, q75], color=col, lw=2.5, solid_capstyle="round", alpha=0.9)
                ax.scatter([x], [med], marker=mk, s=80, color=col, edgecolors=SURFACE, linewidths=1.5,
                           zorder=3)
                xs_line.append(x); ys_line.append(med)
            ax.plot(xs_line, ys_line, lw=1, alpha=0.5, color=col)
        ax.set_xticks(range(len(widths)), [f"w{w}" for w in widths])
        for iw in range(len(widths) - 1):     # separate the width groups, so each dodged cluster reads as one
            ax.axvline(iw + 0.5, color="#e4e2dd", lw=0.8, zorder=0)
        ax.set_xlim(-0.5, len(widths) - 0.5)
        ax.set_xlabel("biggest floe width (% of channel width)"); ax.set_ylabel(ylabel)
        ax.grid(axis="y", alpha=0.6)
        if key != "t_jam_d":
            ax.set_ylim(-0.03, 1.03)
    axes[0].set_title("(a) jamming time", fontsize=10, loc="left")
    axes[1].set_title("(b) probability of permanent arrest (95% interval)", fontsize=10, loc="left")
    axes[2].set_title("(c) evacuated fraction", fontsize=10, loc="left")
    shape_h = [Line2D([], [], color=shape_color(s, shapes), marker=shape_marker(s), lw=1, ms=8) for s in shapes]
    out_h = [Line2D([], [], color=INK_2, marker="o", ls="", ms=5),
             Line2D([], [], color=INK_2, marker="x", ls="", ms=5),
             Line2D([], [], color=INK_2, marker="o", ls="", ms=5, markerfacecolor="none")]
    fig.legend(shape_h + out_h, shapes + ["drained", "arrested", "unresolved (censored)"],
               loc="outside upper center", ncol=len(shapes) + 3, fontsize=9,
               title="shape (median, bar = interquartile range)  ·  one small marker per draw, style = outcome")
    fig.savefig(path, dpi=130, bbox_inches="tight"); plt.close(fig)
    print(f"Saved {path}")


def plot_table(rows, path):
    """Dashboard: three width x shape tables, each cell coloured by its value and labelled with it —
    P(permanent arrest) (with k/n), median evacuated fraction, median jamming time. Reads the trends at a
    glance: along a row the effect of the shape, down a column the effect of the floe size."""
    groups, widths, shapes = group_by_config(rows)
    if not groups:
        print("table: names do not expose wXX/shape -> skipped"); return
    plt = setup_mpl()
    from matplotlib.colors import LinearSegmentedColormap
    cmap = LinearSegmentedColormap.from_list("ramp", ["#e8f1fc", "#86b6ef", "#2a78d6", "#104281"])
    cmap.set_bad(SURFACE)
    horizon_d = max(r["t_end"] for rs in groups.values() for r in rs) / DAY

    def cell(rs, key):
        if key == "p_arrest":
            k = sum(r["outcome"] == "arrested" for r in rs)
            return k / len(rs), f"{k / len(rs):.0%}\n{k}/{len(rs)}"
        if key == "frac_evac":
            v = float(np.median([r["frac_evac"] for r in rs]))
            return v, f"{v:.0%}"
        v = float(np.median([r["t_jam"] for r in rs])) / DAY
        return v, f"{v:.1f} d"

    # colour = how blocked, in every panel (darker = more jamming): the evacuated fraction is coloured by
    # the ice held back (1 - fraction), its label still shows the fraction itself
    panels = [("p_arrest", "P(permanent arrest)", 1.0),
              ("frac_evac", "evacuated fraction at the horizon (median)", 1.0),
              ("t_jam", f"jamming time (median, days, over {horizon_d:.1f} d)", horizon_d)]
    fig, axes = plt.subplots(1, 3, figsize=(4.0 * 3 + 1.5, 0.75 * len(widths) + 1.6), constrained_layout=True)
    for ax, (key, title, vmax) in zip(axes, panels):
        val = np.full((len(widths), len(shapes)), np.nan)
        txt = {}
        for (w, s), rs in groups.items():
            i, j = widths.index(w), shapes.index(s)
            val[i, j], txt[i, j] = cell(rs, key)
        if key == "frac_evac":
            val = 1.0 - val
        ax.imshow(np.ma.masked_invalid(val), cmap=cmap, vmin=0, vmax=vmax, aspect="auto")
        for (i, j), t in txt.items():
            dark = val[i, j] / vmax > 0.55
            ax.text(j, i, t, ha="center", va="center", fontsize=8.5, color="white" if dark else INK)
        ax.set_xticks(range(len(shapes)), shapes)
        ax.set_yticks(range(len(widths)), [f"w{w}" for w in widths])
        ax.tick_params(length=0)
        for sp in ax.spines.values():
            sp.set_visible(False)
        ax.set_xticks(np.arange(len(shapes) + 1) - 0.5, minor=True)   # 2px surface gap between cells
        ax.set_yticks(np.arange(len(widths) + 1) - 0.5, minor=True)
        ax.grid(which="minor", color=SURFACE, linewidth=2)
        ax.tick_params(which="minor", length=0)
        ax.set_title(title, fontsize=10, loc="left")
    axes[0].set_ylabel("biggest floe width\n(% of channel width)", fontsize=9)
    fig.suptitle("darker = more jamming", fontsize=9, color=INK_2)
    fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)
    print(f"Saved {path}")


def plot_scaling(rows, path, channel):
    """Scaling law: probability of permanent arrest and evacuated fraction vs the outlet-to-floe size
    ratio, one point per (width, shape) configuration, sizes MEASURED on the actual floes (left: equivalent
    diameter, right: min width). If the shapes collapse on one curve, the floe size relative to the outlet
    explains the clogging whatever the shape — the classical hopper/silo result; the residual spread
    between shapes at equal ratio is the genuine shape effect (here still mixed with the initial ice
    concentration, which differs between shapes)."""
    groups, widths, shapes = group_by_config(rows)
    if not groups:
        print("scaling: names do not expose wXX/shape -> skipped"); return
    plt = setup_mpl()
    from matplotlib.lines import Line2D
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), sharex="col", constrained_layout=True)
    sizes = [("d_eq", "equivalent diameter"), ("w_min", "min width")]
    for j, (skey, sname) in enumerate(sizes):
        for s in shapes:
            col, mk = shape_color(s, shapes), shape_marker(s)
            pts = []
            for w in widths:
                rs = groups.get((w, s))
                if not rs:
                    continue
                x = channel / np.median([r[skey] for r in rs])
                if not np.isfinite(x):
                    continue
                k = sum(r["outcome"] == "arrested" for r in rs)
                lo, hi = wilson(k, len(rs))
                fe = np.array([r["frac_evac"] for r in rs])
                pts.append((x, k / len(rs), lo, hi, *np.percentile(fe, [25, 50, 75])))
            pts.sort()
            if not pts:
                continue
            x, p, lo, hi, q25, med, q75 = map(np.array, zip(*pts))
            for ax, yv, ylo, yhi in ((axes[0][j], p, lo, hi), (axes[1][j], med, q25, q75)):
                ax.vlines(x, ylo, yhi, color=col, lw=1.6, alpha=0.6)
                ax.plot(x, yv, lw=1, color=col, alpha=0.5)
                ax.scatter(x, yv, marker=mk, s=70, color=col, edgecolors=SURFACE, linewidths=1.5, zorder=3)
        axes[0][j].set_title(f"floe size = {sname} (area-weighted, measured)", fontsize=10, loc="left")
        axes[1][j].set_xlabel(f"outlet width / floe {sname}")
        for ax in (axes[0][j], axes[1][j]):
            ax.grid(alpha=0.6); ax.set_ylim(-0.03, 1.03)
    axes[0][0].set_ylabel("P(permanent arrest)\n(bar = 95% interval)")
    axes[1][0].set_ylabel("evacuated fraction at the horizon\n(median, bar = interquartile range)")
    fig.legend([Line2D([], [], color=shape_color(s, shapes), marker=shape_marker(s), lw=1, ms=8) for s in shapes],
               shapes, loc="outside upper center", ncol=len(shapes), fontsize=9,
               title=f"one point per (width, shape); outlet width = {channel / 1e3:g} km")
    fig.savefig(path, dpi=140, bbox_inches="tight"); plt.close(fig)
    print(f"Saved {path}")


def write_configs_csv(rows, path):
    """One row per (width, shape): outcome counts, P(arrest) with its Wilson interval, and T_jam /
    evacuated-fraction statistics over the draws (cv = std/mean: how much the draw matters)."""
    groups, widths, shapes = group_by_config(rows)
    with open(path, "w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["width", "shape", "n_draws", "n_arrested", "n_drained", "n_unresolved",
                     "p_arrest", "p_arrest_lo95", "p_arrest_hi95",
                     "t_jam_median_s", "t_jam_mean_s", "t_jam_q25_s", "t_jam_q75_s", "t_jam_cv",
                     "frac_evac_median", "frac_evac_q25", "frac_evac_q75",
                     "t_arrest_median_s", "n_jams_median", "d_eq_m", "w_min_m"])
        for w in widths:
            for s in shapes:
                rs = groups.get((w, s))
                if not rs:
                    continue
                n = len(rs)
                k = {o: sum(r["outcome"] == o for r in rs) for o in ("arrested", "drained", "unresolved")}
                lo, hi = wilson(k["arrested"], n)
                tj = np.array([r["t_jam"] for r in rs])
                fe = np.array([r["frac_evac"] for r in rs])
                ta = np.array([r["t_arrest"] for r in rs if r["outcome"] == "arrested"])
                wr.writerow([w, s, n, k["arrested"], k["drained"], k["unresolved"],
                             f"{k['arrested'] / n:.3f}", f"{lo:.3f}", f"{hi:.3f}",
                             f"{np.median(tj):.1f}", f"{tj.mean():.1f}",
                             f"{np.percentile(tj, 25):.1f}", f"{np.percentile(tj, 75):.1f}",
                             f"{tj.std() / tj.mean():.3f}" if tj.mean() > 0 else "",
                             f"{np.median(fe):.4f}", f"{np.percentile(fe, 25):.4f}", f"{np.percentile(fe, 75):.4f}",
                             f"{np.median(ta):.1f}" if len(ta) else "",
                             f"{np.median([r['n_jams'] for r in rs]):.1f}",
                             f"{np.median([r['d_eq'] for r in rs]):.0f}", f"{np.median([r['w_min'] for r in rs]):.0f}"])
    print(f"Saved {path}")


# ---------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("files", nargs="*", help="output .h5 files (globs expanded by the shell)")
    ap.add_argument("--manifest", help="CSV with columns file[,label] (alternative to positional files, for UI driving)")
    ap.add_argument("--ymin", type=float, default=0.0, help="window lower y bound (default 0)")
    ap.add_argument("--ymax", type=float, default=60000.0, help="window upper y bound (default 60000)")
    ap.add_argument("--vmin", type=float, default=0.03, help="avg-speed threshold defining 'jammed' (m/s, default 0.03 = 10%% of the 0.3 m/s free drift)")
    ap.add_argument("--min-floes", type=int, default=10, help="absolute floor for the window to count as 'populated' (default 10)")
    ap.add_argument("--drain-frac", type=float, default=0.1, help="the channel is 'drained' only once fewer than this FRACTION of the pack's floes REMAIN in the system (not yet discharged past the outlet) AND the outflow has stopped (default 0.1). Discards a few residual floes lodged in the funnel corners, while a flow-blocking clog that holds back a lot of ice (even backed up above the window) keeps counting as jammed.")
    ap.add_argument("--outlet", type=float, default=0.0, help="outlet y coordinate: a floe counts as discharged once y <= outlet (default 0). Drives both the live-outflow signal and the 'ice remaining' load.")
    ap.add_argument("--arrest-frac", type=float, default=0.2, help="a non-drained run is 'arrested' (permanent clog) when not a single floe left through the outlet over this last FRACTION of the horizon (default 0.2, i.e. ~5.8 days of a 28.9-day run); otherwise it is 'unresolved' (censored).")
    ap.add_argument("--horizon", type=float, default=None, help="analysis horizon in seconds: truncate every run there (default: the whole run). Re-run with e.g. 1209600 (14 d) to check that the conclusions do not depend on the horizon.")
    ap.add_argument("--speed", choices=["disp", "recorded"], default="disp", help="floe speed used for the window average: 'disp' (default) = displacement between consecutive outputs / elapsed time, immune to the phantom drag speed of floes at rest (a*dt, which grows with dt in frozen packs); 'recorded' = the velocity stored in the file (previous behaviour)")
    ap.add_argument("--min-draws", type=int, default=3, help="ensemble beams are drawn only where at least this many draws are still defined (not drained) (default 3)")
    ap.add_argument("--trace", help="dump the per-time diagnostic series (count, exited, remaining, avg, jammed) of the sim whose label CONTAINS this substring, to <trace>_trace.csv — to validate a single case (e.g. --trace simu_179).")
    ap.add_argument("--skip-first", type=int, default=None, help="treat the first N shapes as obstacles (default: auto-detect = shapes with exactly zero speed all run long)")
    ap.add_argument("--smooth", type=int, default=21, help="grid smoothing window in samples (rolling mean per draw before the ensemble statistics; default 21, 1 = raw)")
    ap.add_argument("--out-dir", help="write EVERY artifact into this one folder (created if needed), with default names (jam.csv, configs.csv, table.png, scaling.png, grid.png, states.png, discharge.png, summary.png, per_sim/, <trace>_trace.csv) unless overridden below. Lets you rsync/cp the whole folder in one go.")
    ap.add_argument("--csv", help="write the per-simulation table to this CSV (default jam.csv under --out-dir)")
    ap.add_argument("--configs-csv", help="write the per-configuration table to this CSV (default configs.csv under --out-dir)")
    ap.add_argument("--plots", help="write one PNG per simulation into this directory (default per_sim/ under --out-dir)")
    ap.add_argument("--grid", help="write the ensemble avg-speed small multiples to this path (default grid.png under --out-dir)")
    ap.add_argument("--states", help="write the flowing/jammed/drained fractions figure to this path (default states.png under --out-dir)")
    ap.add_argument("--discharge", help="write the evacuated-fraction-vs-time figure to this path (default discharge.png under --out-dir)")
    ap.add_argument("--summary", help="write the per-width summary figure to this path (default summary.png under --out-dir)")
    ap.add_argument("--scaling", help="write the scaling-law figure (P(arrest) and evacuated fraction vs outlet width / measured floe size) to this path (default scaling.png under --out-dir)")
    ap.add_argument("--channel-width", type=float, default=40000.0, help="outlet (narrowest channel) width in m, for the scaling figure (default 40000 = chenal v4)")
    ap.add_argument("--table", help="write the width x shape dashboard (P(arrest), evacuated fraction, T_jam) to this path (default table.png under --out-dir)")
    args = ap.parse_args()

    # --out-dir: one folder for everything. Fill in default artifact names, and resolve every relative
    # output path under it (absolute paths are left as given).
    out_dir = args.out_dir
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
        args.csv = args.csv or "jam.csv"
        args.configs_csv = args.configs_csv or "configs.csv"
        args.grid = args.grid or "grid.png"
        args.states = args.states or "states.png"
        args.discharge = args.discharge or "discharge.png"
        args.summary = args.summary or "summary.png"
        args.table = args.table or "table.png"
        args.scaling = args.scaling or "scaling.png"
        args.plots = args.plots or "per_sim"

    def resolve(p):
        return os.path.join(out_dir, p) if (p and out_dir and not os.path.isabs(p)) else p

    csv_path, configs_path, plots_dir = resolve(args.csv), resolve(args.configs_csv), resolve(args.plots)
    grid_path, states_path = resolve(args.grid), resolve(args.states)
    discharge_path, summary_path, table_path = resolve(args.discharge), resolve(args.summary), resolve(args.table)
    scaling_path = resolve(args.scaling)

    entries = []  # (path, label)
    if args.manifest:
        with open(args.manifest) as fh:
            for row in csv.DictReader(fh):
                entries.append((row["file"], row.get("label") or os.path.basename(row["file"])))
    entries += [(p, os.path.splitext(os.path.basename(p))[0]) for p in args.files]
    if not entries:
        ap.error("no input files (positional or --manifest)")

    rows = []
    for path, label in entries:
        try:
            times, avg, count, n_floes, exited = load_series(path, args.ymin, args.ymax, args.outlet,
                                                             args.skip_first, args.horizon, args.speed)
            m = jamming_metrics(times, avg, count, exited, n_floes, args.vmin, args.min_floes,
                                args.drain_frac, args.arrest_frac)
            w, shape, draw = parse_name(label)
            d_eq, w_min = floe_sizes(path)
            drained_after = (times > m["t_drain"]) if m["outcome"] == "drained" else np.zeros(len(times), bool)
            avg_plot = avg.astype(float).copy()
            avg_plot[(count == 0) | drained_after] = np.nan   # an empty/drained window is not a speed of 0
            state = np.where(drained_after, 2, np.where(m["jam_mask"], 1, 0))
            step = float(np.median(np.diff(times))) if len(times) > 1 else 0.0
            rows.append(dict(file=path, label=label, w=w, shape=shape, draw=draw,
                             complete=bool(times[0] <= 1.5 * step), step=step, d_eq=d_eq, w_min=w_min,
                             times=times, avg_plot=avg_plot, state=state,
                             evac=exited / n_floes if n_floes else np.full(len(times), np.nan),
                             t_jam=m["t_jam"], jam_frac=m["jam_frac"], t_start=m["t_start"],
                             t_drain=m["t_drain"], t_end=m["t_end"], n_floes=m["n_floes"],
                             outcome=m["outcome"], t_arrest=m["t_arrest"], frac_evac=m["frac_evac"],
                             n_jams=m["n_jams"]))
            print(f"{label}: {m['outcome']:<10} T_jam={m['t_jam'] / DAY:.2f} d ({100 * m['jam_frac']:.1f}%)  "
                  f"evacuated={100 * m['frac_evac']:.0f}%  drain={m['t_drain'] / DAY:.2f} d  "
                  f"end={m['t_end'] / DAY:.2f} d  ({n_floes} floes, drain<{m['drain_thresh']:.0f})")
            if plots_dir:
                os.makedirs(plots_dir, exist_ok=True)
                plot_single(times, avg, m, args.vmin, label, os.path.join(plots_dir, label + ".png"))
            if args.trace and (args.trace in label or args.trace in path):
                tpath = resolve(f"{args.trace}_trace.csv")
                jm = m["jam_mask"]
                with open(tpath, "w", newline="") as th:
                    tw = csv.writer(th); tw.writerow(["time", "count_in_window", "exited", "remaining", "avg_speed", "jammed"])
                    for k in range(len(times)):
                        tw.writerow([f"{times[k]:.0f}", int(m['count'][k]), int(m['exited'][k]),
                                     int(m['remaining'][k]), f"{avg[k]:.4f}", int(jm[k])])
                print(f"  trace -> {tpath}  (t_start={m['t_start']:.0f}, t_drain={m['t_drain']:.0f}, thresh={m['drain_thresh']:.0f})")
        except Exception as e:
            print(f"ERROR {path}: {e}", file=sys.stderr)

    # Incomplete runs: output starting after t=0 (e.g. a restart that overwrote the file) or ending well
    # before the others. Their outcome and T_jam are not comparable: kept in jam.csv, flagged, and left out
    # of every ensemble figure and of configs.csv.
    if rows:
        t_end_max = max(r["t_end"] for r in rows)
        for r in rows:
            if r["t_end"] < t_end_max - 1.5 * r["step"]:
                r["complete"] = False
        bad = [r for r in rows if not r["complete"]]
        if bad:
            print(f"WARNING {len(bad)} incomplete run(s), excluded from the ensemble statistics:", file=sys.stderr)
            for r in bad:
                print(f"  {r['label']}  ({r['file']}): output from {r['times'][0] / DAY:.2f} d "
                      f"to {r['t_end'] / DAY:.2f} d", file=sys.stderr)
    grid = time_grid(rows)

    if csv_path:
        with open(csv_path, "w", newline="") as fh:
            wr = csv.writer(fh)
            wr.writerow(["file", "label", "width", "shape", "draw", "complete", "outcome",
                         "t_jam_s", "jam_frac", "frac_evacuated", "n_jams", "t_start_s", "t_arrest_s",
                         "t_drain_s", "t_end_s", "n_floes", "d_eq_m", "w_min_m"])
            for r in rows:
                wr.writerow([r["file"], r["label"], r["w"], r["shape"], r["draw"], int(r["complete"]), r["outcome"],
                             f"{r['t_jam']:.1f}", f"{r['jam_frac']:.4f}", f"{r['frac_evac']:.4f}", r["n_jams"],
                             f"{r['t_start']:.1f}", "" if math.isnan(r["t_arrest"]) else f"{r['t_arrest']:.1f}",
                             f"{r['t_drain']:.1f}", f"{r['t_end']:.1f}", r["n_floes"],
                             f"{r['d_eq']:.0f}", f"{r['w_min']:.0f}"])
        print(f"Saved {csv_path}")
    if configs_path:
        write_configs_csv(rows, configs_path)
    if table_path:
        plot_table(rows, table_path)
    if scaling_path:
        plot_scaling(rows, scaling_path, args.channel_width)
    if grid is None:
        print("no complete run -> ensemble figures skipped"); return
    if grid_path:
        plot_grid(rows, grid, grid_path, args.vmin, args.smooth, args.min_draws)
    if states_path:
        plot_states(rows, grid, states_path)
    if discharge_path:
        plot_discharge(rows, grid, discharge_path, args.min_draws)
    if summary_path:
        plot_summary(rows, summary_path)


if __name__ == "__main__":
    main()
