"""Register a map to a regular grid of fiducials (SECCM deposits, indents, FIB marks, a TEM grid ...) and measure the map's true geometry.
The fiducials must show up as bright dots (or dark, --dark) in an image of the map: the per-pattern mean counts of a .pat (computed and cached
as <out>/meanmap.npy), or any 2-D .npy / .tif / .png. Method (NOTES 15.6):
  1. dots = local maxima of the high-passed image above --sigma robust SDs; the image's own dot lattice is fitted from nearest-neighbour vectors;
  2. the fiducial list (x, y in um; optional present 0/1 column for missing fiducials) is indexed on its own lattice; all 8 rotations/mirrors
     of that index grid and all integer shifts are scored on deposit/no-deposit agreement with the dot lattice (agree - disagree - unexplained
     dots), so missing fiducials, an extra point or an asymmetric marker make the answer unique; the margin to the runner-up is reported;
  3. the affine (um -> px) is refined by least squares on the matched dots.
Reports the orientation (rotation, and whether it is mirrored), RMS residual, and the physical pixel size along the map's columns and rows plus
the angle between the map axes on the sample (a nominal-step check, e.g. for GND, which uses the nominal step in both directions).
usage: register_grid.py --image P.pat|img.npy --points fid.csv [--xcol x --ycol y --present-col present] --step 0.546 --out DIR [--dark]
Output DIR/register.json, DIR/points_in_map.csv, DIR/register.png."""
import argparse, csv, json, os, numpy as np
from scipy import ndimage as ndi
from scipy.spatial import cKDTree

ap = argparse.ArgumentParser(); ap.add_argument("--image", required=True); ap.add_argument("--points", required=True); ap.add_argument("--out", required=True)
ap.add_argument("--xcol", default="x"); ap.add_argument("--ycol", default="y"); ap.add_argument("--present-col", default="present")
ap.add_argument("--step", type=float, required=True, help="nominal map step (um per pixel), only used to seed the scale")
ap.add_argument("--sigma", type=float, default=4.0); ap.add_argument("--min-sep", type=int, default=5); ap.add_argument("--tol", type=float, default=2.5)
ap.add_argument("--dark", action="store_true"); ap.add_argument("--hp", type=float, default=3.0, help="high-pass gaussian sigma (px)")
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)

# ---- image
if a.image.lower().endswith(".pat"):
    fm = os.path.join(a.out, "meanmap.npy")
    if os.path.exists(fm): img = np.load(fm)
    else:
        h = np.fromfile(a.image, dtype="<u4", count=5); NC, NR, PW, PH = (int(v) for v in h[1:5]); rec = PW * PH * 2 + 8
        m = np.memmap(a.image, dtype="u1", mode="r"); img = np.empty(NR * NC)
        for k in range(NR * NC): o = 25 + k * rec; img[k] = np.frombuffer(m[o:o + PW * PH * 2].tobytes(), "<u2").mean()
        img = img.reshape(NR, NC); img[0, 0] = np.median(img); np.save(fm, img)          # record 0 of a Thermo .pat is dark
else:
    if a.image.lower().endswith(".npy"): img = np.load(a.image).astype(float)
    else:
        import tifffile
        from PIL import Image
        img = tifffile.imread(a.image) if a.image.lower().endswith((".tif", ".tiff")) else np.array(Image.open(a.image))
        img = (img[..., 0] if img.ndim == 3 else img).astype(float)
img = -img if a.dark else img; NR, NC = img.shape
hp = img - ndi.gaussian_filter(img, a.hp); rob = 1.4826 * np.median(np.abs(hp - np.median(hp)))
pk = (hp == ndi.maximum_filter(hp, a.min_sep)) & (hp > a.sigma * rob); dots = np.argwhere(pk)[:, ::-1].astype(float)
print("image %d x %d; dots %d" % (NC, NR, len(dots)))

# ---- lattices
def basis(P, lo, hi):
    tr = cKDTree(P); _, nb = tr.query(P, k=min(5, len(P))); V = (P[nb[:, 1:]] - P[:, None]).reshape(-1, 2); r = np.hypot(*V.T); V = V[(r > lo) & (r < hi)]
    hx = V[np.abs(V[:, 0]) > np.abs(V[:, 1])]; hx *= np.sign(hx[:, :1]); vy = V[np.abs(V[:, 0]) <= np.abs(V[:, 1])]; vy *= -np.sign(vy[:, 1:])
    return np.c_[np.median(hx, 0), np.median(vy, 0)]
rd = list(csv.DictReader(open(a.points))); F = np.array([[float(r[a.xcol]), float(r[a.ycol])] for r in rd])
pres = np.array([(r.get(a.present_col, "1") or "1") not in ("0", "False", "false") for r in rd])
dnn = np.median(cKDTree(F).query(F, k=2)[0][:, 1]); Pb = basis(F * [1, -1], 0.5 * dnn, 1.5 * dnn) * [[1], [-1]]   # fiducial lattice (um, y up)
Fi = np.round((F - F[0]) @ np.linalg.inv(Pb).T).astype(int); Fi -= Fi.min(0)
pitch_px = dnn / a.step; L = basis(dots, 0.5 * pitch_px, 1.6 * pitch_px); Li = np.linalg.inv(L)
for _ in range(3):                                                                      # refine the dot lattice on integer-indexed dots
    o = np.median(dots - np.round((dots - dots[0]) @ Li.T) @ L.T, 0); uv = np.round((dots - o) @ Li.T)
    good = np.hypot(*((dots - o) - uv @ L.T).T) < 2; sol, *_ = np.linalg.lstsq(np.c_[uv[good], np.ones(good.sum())], dots[good], rcond=None); L, o = sol[:2].T, sol[2]; Li = np.linalg.inv(L)
print("fiducial pitch %.3f um; dot lattice a1 %s a2 %s px" % (dnn, L[:, 0].round(2), L[:, 1].round(2)))
def node(u, v):
    p = o + L @ [u, v]; c, r = int(round(p[0])), int(round(p[1]))
    return hp[r - 1:r + 2, c - 1:c + 2].max() / rob > a.sigma if 1 <= c < NC - 1 and 1 <= r < NR - 1 else None
span = int(max(NR, NC) / min(np.hypot(*L[:, 0]), np.hypot(*L[:, 1]))) + Fi.max() + 2
E = {(u, v): node(u, v) for u in range(-span, span + 1) for v in range(-span, span + 1)}; E = {k: e for k, e in E.items() if e is not None}
S = {tuple(k): bool(p) for k, p in zip(Fi, pres)}
D8 = {"rot0": ((1, 0), (0, 1)), "rot90": ((0, -1), (1, 0)), "rot180": ((-1, 0), (0, -1)), "rot270": ((0, 1), (-1, 0)),
      "mirror_x": ((-1, 0), (0, 1)), "mirror_y": ((1, 0), (0, -1)), "transpose": ((0, 1), (1, 0)), "antitranspose": ((0, -1), (-1, 0))}
dotset = {q for q, e in E.items() if e}; res = []
for nm, M in D8.items():
    M = np.array(M); Q = {k: tuple(M @ k) for k in S}
    for du in range(-span, span + 1):
        for dv in range(-span, span + 1):
            keys = {k: (q[0] + du, q[1] + dv) for k, q in Q.items()}; n = ag = 0
            for k, q in keys.items():
                if q in E: n += 1; ag += E[q] == S[k]
            if n < 0.3 * len(S): continue
            un = len(dotset - set(keys.values())); res.append((2 * ag - n - un, n, ag, n - ag, un, nm, du, dv))
res.sort(reverse=True); best = res[0]; runner = next((r for r in res[1:] if (r[5], r[6], r[7]) != best[5:]), None)
_, nin, ag, dis, un, nm, du, dv = best; M = np.array(D8[nm])
# fiducial (x, y um) -> index (i, j) = Pb^-1 (F - F0) - min ; -> lattice M (i, j) + (du, dv) ; -> pixel o + L (u, v)
Pbi = np.linalg.inv(Pb); off = (np.round((F - F[0]) @ Pbi.T).astype(int)).min(0)
A = L @ M @ Pbi; t = o + L @ (np.array([du, dv]) - M @ off) - A @ F[0]
for _ in range(8):
    d, j = cKDTree(dots).query(F @ A.T + t); ok = (d < a.tol) & pres
    sol, *_ = np.linalg.lstsq(np.c_[F[ok], np.ones(ok.sum())], dots[j[ok]], rcond=None); A, t = sol[:2].T, sol[2]
FP = F @ A.T + t; d, j = cKDTree(dots).query(FP); inside = (FP[:, 0] >= 0) & (FP[:, 0] <= NC - 1) & (FP[:, 1] >= 0) & (FP[:, 1] <= NR - 1)
matched = (d < a.tol) & inside; Ai = np.linalg.inv(A)                                  # px -> um
px_col, px_row = np.hypot(*Ai[:, 0]), np.hypot(*Ai[:, 1]); axis_angle = np.degrees(np.arccos(np.dot(Ai[:, 0], Ai[:, 1]) / (px_col * px_row)))
rot = np.degrees(np.arctan2(-A[1, 0], A[0, 0]))                                         # direction of fiducial +x in the map, CCW with map y up
out = dict(orientation=nm, mirrored=bool(round(np.linalg.det(M)) < 0), score=int(best[0]), runner_up=dict(score=int(runner[0]), orientation=runner[5]) if runner else None,
           fiducials_inside=int(inside.sum()), present_inside=int((inside & pres).sum()), present_matched=int((matched & pres).sum()),
           absent_with_dot=int((matched & ~pres).sum()), unexplained_dots=int(un), rms_px=round(float(np.sqrt(np.mean(d[matched & pres] ** 2))), 3),
           affine_um_to_px=np.r_[A.ravel(), t].round(6).tolist(), fiducial_x_direction_in_map_deg=round(float(rot), 2),
           pixel_um_along_columns=round(float(px_col), 4), pixel_um_along_rows=round(float(px_row), 4), map_axis_angle_deg=round(float(axis_angle), 2),
           nominal_step_um=a.step, ratio_columns=round(float(px_col / a.step), 4), ratio_rows=round(float(px_row / a.step), 4))
json.dump(out, open(os.path.join(a.out, "register.json"), "w"), indent=1); print(json.dumps(out, indent=1))
with open(os.path.join(a.out, "points_in_map.csv"), "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["row_in_file", "x_um", "y_um", "present", "map_col", "map_row", "dot_dist_px", "matched"])
    for k in range(len(F)): w.writerow([k, *F[k], int(pres[k]), round(FP[k, 0], 2), round(FP[k, 1], 2), round(d[k], 2), int(matched[k])])
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(11, 11 * NR / NC + 0.5)); ax.imshow(hp, cmap="gray", vmin=-3 * rob, vmax=6 * rob); ax.plot(dots[:, 0], dots[:, 1], "r.", ms=3, label="dots")
ax.plot(FP[inside & pres, 0], FP[inside & pres, 1], "o", mfc="none", mec="c", ms=7, label="fiducials (fit)"); ax.plot(FP[inside & ~pres, 0], FP[inside & ~pres, 1], "yx", ms=7, label="absent fiducials")
ax.set_xlim(-0.5, NC - 0.5); ax.set_ylim(NR - 0.5, -0.5); ax.legend(fontsize=8, loc="lower left")
ax.set_title("%s, %d/%d matched, rms %.2f px; pixel %.3f x %.3f um, axes %.1f deg" % (nm, out["present_matched"], out["present_inside"], out["rms_px"], px_col, px_row, axis_angle), fontsize=9)
fig.tight_layout(); fig.savefig(os.path.join(a.out, "register.png"), dpi=90)
