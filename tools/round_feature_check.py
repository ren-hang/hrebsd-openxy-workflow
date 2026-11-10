"""Check a map's pixel geometry with a round feature (hole, indent, particle) when no fiducial grid is available.
A circle that is mapped with row pixels covering a fraction r of the column pixel on the sample appears stretched by 1/r along the rows, and a
shear between the map axes tilts the ellipse. From the second moments C of the feature (pixel units, y down) the metric of the map follows up to
scale: (row pixel / column pixel) = sqrt(Cxx / Cyy), angle between the map axes on the sample = arccos(-Cxy / sqrt(Cxx Cyy)).
This assumes the feature is a circle; check that once in a top-view (0 deg tilt) SEM image, or use many features of different origin (an
individual feature's shape error then averages out, while the scan's distortion does not). NOTES 15.4.
Feature detection: the largest compact region of low IQ (or low CI) that does not touch the map border, at --levels between the 1st percentile
and the median of the map (each level gives one estimate; their spread is reported).
usage: round_feature_check.py MAP.ANG [MAP.ANG ...] --out DIR [--metric iq|ci] [--levels 0.3 0.5 0.7] [--min-area 40]
Output DIR/round_features.csv (one row per map and level), DIR/round_features_summary.json, DIR/round_features.png."""
import argparse, csv, json, os, re, numpy as np
from scipy import ndimage as ndi

ap = argparse.ArgumentParser(); ap.add_argument("ang", nargs="+"); ap.add_argument("--out", required=True); ap.add_argument("--metric", default="iq", choices=["iq", "ci"])
ap.add_argument("--levels", nargs="+", type=float, default=[0.3, 0.5, 0.7]); ap.add_argument("--min-area", type=int, default=40); ap.add_argument("--min-solidity", type=float, default=0.85)
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)

def read_ang(p):
    hdr = {}; rows = []
    for line in open(p, encoding="latin-1"):
        if line.startswith("#"):
            m = re.match(r"#\s*(\w+):\s*(.*)", line)
            if m: hdr[m.group(1).upper()] = m.group(2).strip()
        elif line.strip(): rows.append(line.split())
    d = np.array(rows, float); nc = int(hdr.get("NCOLS_ODD") or hdr.get("NCOLS") or 0); nr = int(hdr.get("NROWS", 0))   # Pathfinder writes NCOLS
    if nc * nr != len(d):                                                   # fall back to the x column
        nc = len(np.unique(np.round(d[:, 3], 4))); nr = len(d) // nc
    return hdr, d, nr, nc

out, figs = [], []
for p in a.ang:
    try: hdr, d, nr, nc = read_ang(p)
    except Exception as e: print("skip", p, e); continue
    if hdr.get("GRID", "").lower().startswith("hex"): print("skip (hex grid)", p); continue
    v = (d[:, 5] if a.metric == "iq" else d[:, 6]).reshape(nr, nc).astype(float); v = ndi.gaussian_filter(v, 0.7)
    lo, med = np.percentile(v, 1), np.median(v); best = None
    for lev in a.levels:
        m = v < lo + lev * (med - lo); lab, n = ndi.label(ndi.binary_opening(m))
        cands = []
        for k in range(1, n + 1):
            r = lab == k; ys, xs = np.nonzero(r)
            if r.sum() < a.min_area or ys.min() == 0 or xs.min() == 0 or ys.max() == nr - 1 or xs.max() == nc - 1: continue
            f = ndi.binary_fill_holes(r); hull_area = f.sum()                 # solidity proxy: area / filled-ellipse area of same moments
            C = np.cov(np.vstack(np.nonzero(f)[::-1])); ell = np.pi * 4 * np.sqrt(np.linalg.det(C))
            sol = min(f.sum() / ell, ell / f.sum()) if ell > 0 else 0
            if sol >= a.min_solidity: cands.append((f.sum(), f, C, sol))
        if not cands: continue
        area, f, C, sol = max(cands, key=lambda c: c[0]); ys, xs = np.nonzero(f)
        ratio = float(np.sqrt(C[0, 0] / C[1, 1])); ang = float(np.degrees(np.arccos(np.clip(-C[0, 1] / np.sqrt(C[0, 0] * C[1, 1]), -1, 1))))
        step = float(hdr.get("XSTEP", "nan")); diam = 2 * np.sqrt(area / np.pi) * step
        out.append(dict(map=os.path.basename(p), level=lev, nrows=nr, ncols=nc, step_um=step, area_px=int(area), eq_diameter_um=round(diam, 2),
                        centre_col=round(float(xs.mean()), 1), centre_row=round(float(ys.mean()), 1), apparent_height_over_width=round(1 / ratio, 4),
                        row_over_col_pixel_if_circle=round(ratio, 4), axis_angle_if_circle_deg=round(ang, 2), ellipse_fill=round(sol, 3)))
        if best is None or lev == 0.5: best = (p, v, f)
    if best: figs.append(best)
with open(os.path.join(a.out, "round_features.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(out[0])); w.writeheader(); w.writerows(out)
r = np.array([o["row_over_col_pixel_if_circle"] for o in out]); g = np.array([o["axis_angle_if_circle_deg"] for o in out])
per_map = {}
for o in out: per_map.setdefault(o["map"], []).append(o["row_over_col_pixel_if_circle"])
pm = np.array([np.median(v) for v in per_map.values()])
summ = dict(n_maps_with_feature=len(per_map), n_estimates=len(out), row_over_col_median_over_maps=round(float(np.median(pm)), 4),
            row_over_col_IQR_over_maps=[round(float(v), 4) for v in np.percentile(pm, [25, 75])], fraction_of_maps_below_0p95=round(float(np.mean(pm < 0.95)), 3),
            axis_angle_median=round(float(np.median(g)), 2), axis_angle_IQR=[round(float(v), 2) for v in np.percentile(g, [25, 75])])
json.dump(summ, open(os.path.join(a.out, "round_features_summary.json"), "w"), indent=1); print(json.dumps(summ, indent=1))
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
nf = len(figs); ncol = min(6, nf); nrow = (nf + ncol - 1) // ncol
fig, axs = plt.subplots(nrow, ncol, figsize=(3 * ncol, 2.4 * nrow), squeeze=False)
for ax_, (p, v, f) in zip(axs.flat, figs):
    ax_.imshow(v, cmap="gray"); ax_.contour(f, [0.5], colors="r", linewidths=0.8); ax_.set_axis_off()
    ax_.set_title("%s\nh/w %.2f" % (os.path.basename(p)[:28], 1 / np.median(per_map[os.path.basename(p)])), fontsize=6)
for ax_ in list(axs.flat)[nf:]: ax_.set_axis_off()
fig.tight_layout(); fig.savefig(os.path.join(a.out, "round_features.png"), dpi=90)
