"""Measured mixing fraction across grain boundaries: every pattern near the boundary between two large grains (including the seam of
not-indexed / one-pixel 'grains' the instrument leaves there) is fitted as a*R_1 + b*R_2, where R_1 and R_2 are means of clean patterns
5-8 px inside each grain next to that stretch of boundary (measured references, NOTES 15.8). f = b / (a + b) = fraction of the second grain.
Controls: deeper interior points fitted the same way (threshold = their 99th percentile). Calibration: known fractions of a real interior
pattern of the other grain injected into real interior patterns. Patterns are binned, flat-fielded by the map mean and band-passed
(all linear, so a mixture of patterns stays the same mixture).
usage: gb_mixing.py --src work\hrebsd_Pt6gnd --stem Pt6gnd --out work\gbmix_Pt6 [--dmax 6] [--ref 5 8] [--win 8] [--bin 4] [--min-grain 40]
Coordinates: d_edge = 1 for the last indexed pixel of a large grain, 2 the next one in, ...; s = signed distance (px) from the midline
between the two grains (negative = lower grain id). Writes mixing.npz, mixing.csv, summary.json, mixing.png. Values are fractions, not strain."""
import argparse, csv, json, os, numpy as np, tifffile
from scipy import ndimage as ndi

ap = argparse.ArgumentParser()
ap.add_argument("--src", required=True); ap.add_argument("--stem", required=True); ap.add_argument("--out", required=True)
ap.add_argument("--dmax", type=float, default=6, help="analyse grain pixels up to this d_edge (px)")
ap.add_argument("--ref", nargs=2, type=float, default=[5, 8], help="reference band: distance (px) to the other grain")
ap.add_argument("--win", type=float, default=8); ap.add_argument("--bin", type=int, default=4); ap.add_argument("--min-grain", type=int, default=40)
ap.add_argument("--seed", type=int, default=0)
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True); rng = np.random.default_rng(a.seed)
d = np.load(os.path.join(a.src, "maps.npz")); nr, nc = int(d["nrows"]), int(d["ncols"]); grain = d["grain"]
step = float(d["x_um"][0, 1] - d["x_um"][0, 0])
indexed = ~d["col0"] & ~d["not_indexed"]                       # references: indexed, not scan column 0 (anomalous, NOTES 12.13)
ids = np.unique(grain); big = [int(g) for g in ids if g > 0 and (indexed & (grain == g)).sum() >= a.min_grain]
tif = lambda r, c: os.path.join(a.src, "Corrected_EBSPs", f"{a.stem}_r{r}c{c}.tif")
binp = lambda im: im[: im.shape[0] // a.bin * a.bin, : im.shape[1] // a.bin * a.bin].reshape(im.shape[0] // a.bin, a.bin, -1, a.bin).mean((1, 3))
raw, P = {}, {}
def rawp(r, c):
    if (r, c) not in raw: raw[(r, c)] = binp(tifffile.imread(tif(r, c)).astype(np.float32))
    return raw[(r, c)]
bg = ndi.gaussian_filter(np.mean([rawp(r, c) for r in range(0, nr, 3) for c in range(1, nc, 3)], 0), 2)   # static background
def pre(im):                                                                       # linear: flat-field, then band-pass
    q = im / bg; return (ndi.gaussian_filter(q, 1) - ndi.gaussian_filter(q, 8)).ravel().astype(np.float32)
def pp(r, c):
    if (r, c) not in P: P[(r, c)] = pre(rawp(r, c))
    return P[(r, c)]

rr, cc = np.mgrid[0:nr, 0:nc]
dist = {g: ndi.distance_transform_edt(grain != g) for g in big}                     # px to the nearest pixel of grain g (0 inside g)
D = np.stack([dist[g] for g in big]); order = np.argsort(D, 0)
g1 = np.array(big)[order[0]]; g2 = np.array(big)[order[1]]                          # nearest and second-nearest large grain
d1 = np.take_along_axis(D, order[:1], 0)[0]; d2 = np.take_along_axis(D, order[1:2], 0)[0]
lo, hi = np.minimum(g1, g2), np.maximum(g1, g2)
t = np.where(g1 == lo, d2 - d1, d1 - d2)                                            # >0 on the lower-id side
s = np.zeros((nr, nc))
for gl, gh in {(int(x), int(y)) for x, y in zip(lo.ravel(), hi.ravel())}:
    m = (lo == gl) & (hi == gh); mid = m & (np.abs(t) <= 1)
    if mid.any(): s[m] = (-np.sign(t) * ndi.distance_transform_edt(~mid))[m]
d_edge = np.zeros((nr, nc))
for g in big: d_edge[grain == g] = ndi.distance_transform_edt(grain == g)[grain == g]
cross = np.full((nr, nc), "", object)                                               # direction to the other grain: across rows / columns
for g in big:
    _, ind = ndi.distance_transform_edt(grain != g, return_indices=True)
    m = g2 == g; dy, dx = np.abs(ind[0] - rr), np.abs(ind[1] - cc)
    cross[m] = np.where(dy[m] > dx[m], "across rows", np.where(dx[m] > dy[m], "across columns", "diagonal"))

def refs(r, c, g, h, win):
    """mean preprocessed pattern of indexed grain-g points whose distance to grain h is in the reference band, within win px of (r, c)"""
    m = (grain == g) & indexed & (dist[h] >= a.ref[0]) & (dist[h] <= a.ref[1]) & ((rr - r) ** 2 + (cc - c) ** 2 <= win ** 2)
    pts = list(zip(*np.nonzero(m)))
    return (np.mean([pp(*p) for p in pts], 0), len(pts)) if len(pts) >= 3 else (None, len(pts))
def fit(q, R1, R2):
    X = np.stack([R1, R2], 1); coef, *_ = np.linalg.lstsq(X, q, rcond=None); tot = coef.sum()
    return (coef[1] / tot if tot > 0 else np.nan), float(tot), float(np.linalg.norm(q - X @ coef) / np.linalg.norm(q))

rows, ctrl, ctrl_res, inj = [], [], [], []; levels = [0.02, 0.05, 0.1, 0.2, 0.3, 0.5]; no_ref = 0
fmap = np.full((nr, nc), np.nan)
for r in range(nr):
    for c in range(1, nc):
        ga, gb, da, db = int(g1[r, c]), int(g2[r, c]), d1[r, c], d2[r, c]
        in_grain = grain[r, c] == ga and indexed[r, c]
        seam = grain[r, c] not in big and db <= a.ref[0]
        is_ctrl = in_grain and a.ref[1] + 2 <= db <= a.ref[1] + 6
        if not (seam or is_ctrl or (in_grain and d_edge[r, c] <= a.dmax and db <= a.dmax + 4)): continue
        R1, n1 = refs(r, c, ga, gb, a.win); R2, n2 = refs(r, c, gb, ga, a.win)
        if R1 is None or R2 is None: R1, n1 = refs(r, c, ga, gb, 2 * a.win); R2, n2 = refs(r, c, gb, ga, 2 * a.win)
        if R1 is None or R2 is None: no_ref += 1; continue
        f, tot, res = fit(pp(r, c), R1, R2)
        if is_ctrl:
            ctrl.append(f); ctrl_res.append(res)
            if len(inj) < 60 * len(levels):                                       # inject a real interior pattern of the other grain
                cand = np.argwhere((grain == gb) & indexed & (dist[ga] >= a.ref[1] + 2) & ((rr - r) ** 2 + (cc - c) ** 2 <= (3 * a.win) ** 2))
                if len(cand):
                    pb = rawp(*cand[rng.integers(len(cand))])
                    inj += [(lv, fit(pre((1 - lv) * rawp(r, c) + lv * pb), R1, R2)[0]) for lv in levels]
            continue
        fmap[r, c] = f
        rows.append(dict(row=r, col=c, label=int(grain[r, c]), grain=ga, other=gb, kind="seam" if seam else "grain", d_edge=float(d_edge[r, c]),
                         s_px=float(s[r, c]), f_other=float(f), f_higher_id=float(f if gb > ga else 1 - f), direction=cross[r, c],
                         a_plus_b=tot, resid=res, n_ref=n1, n_ref_other=n2))
print(f"points fitted {len(rows)} (seam {sum(x['kind'] == 'seam' for x in rows)}), controls {len(ctrl)}, without references {no_ref}; large grains {big}")
if not rows: raise SystemExit("nothing fitted: no boundary between two large grains with references on both sides (check --ref / --win / --min-grain)")
ctrl = np.array(ctrl); thr = float(np.nanpercentile(ctrl, 99)) if len(ctrl) else float("nan")
rthr = float(np.nanpercentile(ctrl_res, 99)) if len(ctrl_res) else float("nan")       # unexplained pattern fraction of pure points
for x in rows: x["unreliable"] = int(x["kind"] == "seam" or x["f_other"] > thr or x["resid"] > rthr)
bad = np.zeros((nr, nc), bool)
for x in rows: bad[x["row"], x["col"]] = bool(x["unreliable"])
inj = np.array(inj)
cal = {str(lv): dict(median=float(np.nanmedian(inj[inj[:, 0] == lv, 1])), p10=float(np.nanpercentile(inj[inj[:, 0] == lv, 1], 10)),
                     p90=float(np.nanpercentile(inj[inj[:, 0] == lv, 1], 90)), n=int((inj[:, 0] == lv).sum())) for lv in levels} if len(inj) else {}
G = [x for x in rows if x["kind"] == "grain"]
de = np.array([x["d_edge"] for x in G]); fo = np.array([x["f_other"] for x in G]); di = np.array([x["direction"] for x in G])
ub = np.array([x["unreliable"] for x in G], bool)
gg = np.array([x["grain"] for x in G]); prof = {}
for lab, sel in [("all", np.ones(len(G), bool)), ("across rows", di == "across rows"), ("across columns", di == "across columns")] + [(f"grain {g}", gg == g) for g in big]:
    prof[lab] = []
    for k in range(1, int(a.dmax) + 1):
        m = sel & (np.round(de) == k)
        prof[lab].append(dict(d_edge=k, d_um=k * step, n=int(m.sum()), median=float(np.nanmedian(fo[m])) if m.any() else None,
                              p90=float(np.nanpercentile(fo[m], 90)) if m.any() else None, frac_above_thr=float(np.mean(fo[m] > thr)) if m.any() else None,
                              frac_unreliable=float(np.mean(ub[m])) if m.any() else None))
S_ = np.array([x["s_px"] for x in rows]); FH = np.array([x["f_higher_id"] for x in rows])
across = [dict(s_px=float(k), n=int((np.round(S_) == k).sum()), median_f_higher=float(np.nanmedian(FH[np.round(S_) == k])))
          for k in range(-int(a.dmax) - 2, int(a.dmax) + 3) if (np.round(S_) == k).any()]
seamf = np.array([x["f_other"] for x in rows if x["kind"] == "seam"])
near4 = np.isin(grain, big) & indexed & (d_edge <= 4) & (d2 <= 4 + a.ref[0])         # what a 4 px margin removes from the large grains
summ = dict(src=a.src, step_um=step, large_grains=big, n_grain_points=len(G), n_seam_points=int(len(seamf)), n_controls=int(len(ctrl)),
            control_median=float(np.nanmedian(ctrl)) if len(ctrl) else None, threshold_p99=thr, residual_threshold_p99=rthr,
            grain_points_unreliable=int(ub.sum()), grain_points_within_4px=int(near4.sum()),
            within_4px_but_clean=int((near4 & ~bad).sum()), injection=cal, profile_d_edge=prof,
            profile_across=across, seam_min_side_fraction_median=float(np.nanmedian(np.minimum(seamf, 1 - seamf))) if len(seamf) else None,
            settings=vars(a))
json.dump(summ, open(os.path.join(a.out, "summary.json"), "w"), indent=1)
with open(os.path.join(a.out, "mixing.csv"), "w", newline="") as fh:
    wr = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); wr.writeheader(); wr.writerows(rows)
np.savez_compressed(os.path.join(a.out, "mixing.npz"), f_other=fmap, unreliable=bad, d_edge=d_edge, s_px=s, grain=grain, nearest=g1, second=g2,
                    threshold=thr, residual_threshold=rthr, step_um=step)
print(json.dumps({k: summ[k] for k in ("n_grain_points", "n_seam_points", "n_controls", "control_median", "threshold_p99", "residual_threshold_p99",
                                       "grain_points_unreliable", "grain_points_within_4px", "within_4px_but_clean", "seam_min_side_fraction_median")}, indent=1))
print("calibration (injected -> median recovered [p10, p90]):", {k: (round(v["median"], 3), round(v["p10"], 3), round(v["p90"], 3)) for k, v in cal.items()})
print("per d_edge: (d_edge, n, median fraction, share above threshold, share unreliable)")
for lab, pr in prof.items():
    print(f"{lab:15s}", [(p["d_edge"], p["n"], None if p["median"] is None else round(p["median"], 3), None if p["frac_above_thr"] is None else round(p["frac_above_thr"], 2),
                          None if p["frac_unreliable"] is None else round(p["frac_unreliable"], 2)) for p in pr])
print("across the boundary (s_px, n, median fraction of the higher-id grain):", [(x["s_px"], x["n"], round(x["median_f_higher"], 3)) for x in across])

import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 4, figsize=(21, 4.6))
im = ax[0].imshow(fmap, vmin=0, vmax=0.5, cmap="viridis", interpolation="nearest"); plt.colorbar(im, ax=ax[0], fraction=0.03, label="fraction of the other grain")
ax[0].contour(np.isin(grain, big).astype(float), levels=[0.5], colors="w", linewidths=0.6); ax[0].set_title("measured mixing fraction (white: edge of the large grains)")
ax[1].plot([x["s_px"] * step for x in across], [x["median_f_higher"] for x in across], "o-")
ax[1].set_xlabel("distance from the boundary midline (um)"); ax[1].set_ylabel("median fraction of the higher-id grain"); ax[1].set_title("profile across the boundary")
for lab, mk in (("all", "o-"), ("across rows", "s--"), ("across columns", "^--")):
    pr = [p for p in prof[lab] if p["median"] is not None]
    ax[2].plot([p["d_um"] for p in pr], [p["median"] for p in pr], mk, label=lab)
ax[2].axhline(thr, color="r", lw=1, label=f"control 99th pct ({thr:.3f})"); ax[2].set_xlabel("distance into the grain from its last indexed pixel (um)")
ax[2].set_ylabel("median fraction of the other grain"); ax[2].legend(fontsize=8); ax[2].set_title("grain pixels")
if cal:
    lv = np.array(levels); ax[3].plot(lv, [cal[str(x)]["median"] for x in levels], "o-", label="median recovered")
    ax[3].fill_between(lv, [cal[str(x)]["p10"] for x in levels], [cal[str(x)]["p90"] for x in levels], alpha=0.3, label="10-90%")
    ax[3].plot([0, 0.5], [0, 0.5], "k:", lw=1); ax[3].set_xlabel("injected fraction"); ax[3].set_ylabel("recovered"); ax[3].legend(fontsize=8)
    ax[3].set_title("calibration with real patterns")
fig.tight_layout(); fig.savefig(os.path.join(a.out, "mixing.png"), dpi=105)
