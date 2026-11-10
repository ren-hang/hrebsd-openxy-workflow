"""Affine registration of an EBSD map to a top-view SEM image of the same area by its grain boundaries; measures the map's true pixel geometry.
No fiducials or round features are needed: every grain boundary in the overlap constrains the fit (NOTES 15.4).
EBSD edge map: misorientation to the right and lower neighbour (cubic m-3m by default, orix) scaled to 0-1 at --gb-deg, plus the edge of the
non-indexed region (CI < --min-ci, e.g. a hole). SEM edge map: gradient magnitude of the smoothed image (grain contrast steps, grooves, rims),
block-averaged to --sem-px um (default step/2). The SEM image must be a 0 deg tilt, tilt-correction-off image (square pixels on the sample).
Search: rotation 0..360 by --rot-step (proper rotations only unless --allow-mirror), translation by masked normalised cross-correlation at the
nominal scale (skimage); then all six affine parameters by maximising the NCC (Powell). Uncertainty: the fit repeated on each quadrant of the
EBSD map that lies >= 60% inside the SEM field. "reliable" = NCC beats the best fit >20 deg away by > 0.1 and quadrants agree (ratio SD < 0.06,
angle SD < 3 deg). Use the SEM image taken of this map's area (same mount if possible); a smaller SEM field is fine.
Reports the size of one EBSD pixel along the map's columns and rows on the sample (um), the angle between the map axes, the rotation of the
map in the SEM image, whether the frames are mirrored, and the NCC. --selftest warps the SEM edge map with a known distortion and checks recovery.
usage: register_sem_ebsd.py --ang MAP.ANG --sem SEM.tif --out DIR [--gb-deg 5] [--sem-px 0.63] [--rot-step 2] [--selftest]
Output DIR/register_sem_ebsd.json, DIR/register_sem_ebsd.png."""
import argparse, json, os, re, numpy as np, tifffile
from scipy import ndimage as ndi
from scipy.optimize import minimize
from skimage.registration import phase_cross_correlation

ap = argparse.ArgumentParser(); ap.add_argument("--ang", required=True); ap.add_argument("--sem", required=True); ap.add_argument("--out", required=True)
ap.add_argument("--gb-deg", type=float, default=5.0); ap.add_argument("--min-ci", type=float, default=0.1); ap.add_argument("--sem-px", type=float, default=None)
ap.add_argument("--rot-step", type=float, default=2.0); ap.add_argument("--selftest", action="store_true"); ap.add_argument("--no-quadrants", action="store_true")
ap.add_argument("--allow-mirror", action="store_true", help="also try mirrored frames (both images view the surface from the beam side, so normally not)")
ap.add_argument("--row-scales", nargs="+", type=float, default=[0.78, 0.86, 0.94, 1.02, 1.1, 1.18], help="row/column scales tried in the coarse search")
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)

# ---- EBSD
hdr, rows = {}, []
for line in open(a.ang, encoding="latin-1"):
    if line.startswith("#"):
        m = re.match(r"#\s*(\w+):\s*(.*)", line)
        if m: hdr[m.group(1).upper()] = m.group(2).strip()
    elif line.strip(): rows.append(line.split())
d = np.array(rows, float); NR, step = int(hdr["NROWS"]), float(hdr["XSTEP"])
NC = int(hdr.get("NCOLS_ODD") or hdr.get("NCOLS") or len(d) // NR)                         # Pathfinder writes NCOLS, TSL NCOLS_ODD
from orix.quaternion import Orientation, symmetry
O = Orientation.from_euler(d[:, :3], symmetry.Oh).reshape(NR, NC)
misr = np.zeros((NR, NC)); misd = np.zeros((NR, NC))
misr[:, :-1] = np.degrees(O[:, :-1].angle_with(O[:, 1:])); misd[:-1] = np.degrees(O[:-1].angle_with(O[1:]))
ci = d[:, 6].reshape(NR, NC); bad = ndi.binary_opening(ci < a.min_ci, iterations=1)
edge_e = np.clip(np.maximum(misr, misd) / a.gb_deg, 0, 1); edge_e[bad] = 0
edge_e = np.maximum(edge_e, ndi.binary_dilation(bad) & ~ndi.binary_erosion(bad))
edge_e = ndi.gaussian_filter(edge_e.astype(float), 0.6)

# ---- SEM
t = tifffile.TiffFile(a.sem); p = t.pages[0]; s = str(p.tags[34682].value) if 34682 in p.tags else ""
g = lambda k: (re.search(k + r"'?\s*:\s*'?([-0-9.eE+]+)", s) or [None, None])[1]
pw, ph = float(g("PixelWidth")) * 1e6, float(g("PixelHeight")) * 1e6; tilt = np.degrees(float(g("StageT") or 0))
if abs(pw - ph) > 1e-3 * pw or abs(tilt) > 0.5: print("WARNING: SEM pixels not square or tilted (%.4f/%.4f um, %.2f deg)" % (pw, ph, tilt))
img = p.asarray(); img = (img[..., 0] if img.ndim == 3 else img).astype(np.float32); img = img[:int(img.shape[0] * 0.935)]   # FEI databar
pd = a.sem_px or step / 2; f = max(1, int(round(pd / pw))); pd = f * pw
H, W = img.shape[0] // f, img.shape[1] // f; small = img[:H * f, :W * f].reshape(H, f, W, f).mean((1, 3))
gs = ndi.gaussian_filter(small, 1.0); edge_s = np.hypot(ndi.sobel(gs, 0), ndi.sobel(gs, 1)); edge_s = np.clip(edge_s / np.percentile(edge_s, 99), 0, 1)
print("EBSD %d x %d at %.3f um; SEM %d x %d at %.3f um (binned %dx from %.4f um)" % (NC, NR, step, W, H, pd, f, pw))

def warp(Amat, tvec, src, shape):
    """EBSD image (row, col) -> SEM grid: SEM (x, y) = A (col, row) + t; returns warped image and validity mask."""
    Ai = np.linalg.inv(Amat); M = np.array([[Ai[1, 1], Ai[1, 0]], [Ai[0, 1], Ai[0, 0]]])        # (y, x) -> (row, col)
    off = np.array([-(Ai @ tvec)[1], -(Ai @ tvec)[0]])
    w = ndi.affine_transform(src, M, offset=off, output_shape=shape, order=1, cval=0.0)
    m = ndi.affine_transform(np.ones_like(src), M, offset=off, output_shape=shape, order=0, cval=0.0) > 0.5
    return w, m
def ncc(Amat, tvec, src, ref):
    w, m = warp(Amat, tvec, src, ref.shape)
    if m.sum() < 200: return -1.0
    x, y = w[m] - w[m].mean(), ref[m] - ref[m].mean(); return float((x * y).sum() / np.sqrt((x * x).sum() * (y * y).sum() + 1e-12))
def register(src, ref, s0, quiet=False):
    best = (-2, None); scan = []
    for mir in ((1, -1) if a.allow_mirror else (1,)):
      for fr in a.row_scales:                                                                   # rows can be off by up to ~20% (NOTES 15.13)
        for th in np.arange(0, 360, a.rot_step):
            c, sn = np.cos(np.radians(th)), np.sin(np.radians(th)); A0 = s0 * np.array([[c, -sn], [sn, c]]) @ np.diag([1, mir * fr])
            corners = A0 @ np.array([[0, 0, src.shape[1], src.shape[1]], [0, src.shape[0], 0, src.shape[0]]]); t0 = -corners.min(1)
            w, m = warp(A0, t0, src, ref.shape)
            if m.sum() < 200: continue
            sh = phase_cross_correlation(ref, w, reference_mask=np.ones_like(ref, bool), moving_mask=m, overlap_ratio=0.25)
            sh = np.asarray(sh[0] if isinstance(sh, tuple) else sh, float); t1 = t0 + sh[::-1]
            sc = ncc(A0, t1, src, ref); scan.append((th, mir, sc))
            if sc > best[0]: best = (sc, (A0, t1, th, mir, fr))
    sc0, (A0, t1, th, mir, fr) = best
    other = max([s for t_, m_, s in scan if m_ != mir or abs((t_ - th + 180) % 360 - 180) > 20], default=np.nan)   # best fit elsewhere
    if not quiet: print("coarse: rotation %.0f deg, row scale %.2f, %s, NCC %.3f (best NCC >20 deg away: %.3f)" % (th, fr, "mirrored" if mir < 0 else "proper", sc0, other))
    x0 = np.r_[A0.ravel(), t1]; f_ = lambda x: -ncc(x[:4].reshape(2, 2), x[4:], src, ref)
    r = minimize(f_, x0, method="Powell", options=dict(xtol=1e-4, ftol=1e-6, maxfev=6000))
    register.other = other
    return r.x[:4].reshape(2, 2), r.x[4:], -r.fun, sc0
def geometry(Amat):
    vc, vr = pd * Amat[:, 0], pd * Amat[:, 1]                                                   # um on the sample per EBSD column / row step
    ang = np.degrees(np.arccos(np.dot(vc, vr) / np.linalg.norm(vc) / np.linalg.norm(vr)))
    return dict(pixel_um_along_columns=round(float(np.linalg.norm(vc)), 4), pixel_um_along_rows=round(float(np.linalg.norm(vr)), 4),
                map_axis_angle_deg=round(float(ang), 2), ratio_rows_to_columns=round(float(np.linalg.norm(vr) / np.linalg.norm(vc)), 4),
                columns_direction_in_sem_deg=round(float(np.degrees(np.arctan2(-vc[1], vc[0]))), 2), mirrored=bool(np.linalg.det(Amat) < 0))

s0 = step / pd
if a.selftest:                                                                                 # known distortion: rows 0.83x, axes 86.3 deg, 37 deg rotation
    true = dict(col=step, row=0.83 * step, angle=86.3, rot=37.0)
    vc = true["col"] / pd * np.array([np.cos(np.radians(true["rot"])), -np.sin(np.radians(true["rot"]))])
    ar = np.radians(true["rot"] - true["angle"]); vr = true["row"] / pd * np.array([np.cos(ar), -np.sin(ar)])
    At = np.c_[vc, vr]; tt = np.array([edge_s.shape[1] * 0.2, edge_s.shape[0] * 0.25]); n = (int(edge_s.shape[0] * 0.35 / (true["row"] / pd)), int(edge_s.shape[1] * 0.35 / s0))
    rr, cc = np.mgrid[0:n[0], 0:n[1]]; XY = At @ np.vstack([cc.ravel(), rr.ravel()]) + tt[:, None]
    syn = ndi.map_coordinates(ndi.gaussian_filter(edge_s, 0.5 * s0), [XY[1], XY[0]], order=1).reshape(n)
    Af, tf, sc, _ = register(syn, edge_s, s0); gm = geometry(Af)
    print("selftest truth: col %.3f, row %.3f, angle %.1f, rotation %.1f" % (true["col"], true["row"], true["angle"], true["rot"])); print("selftest fit:", gm, "NCC %.3f" % sc)
    json.dump(dict(selftest_truth=true, selftest_fit=gm, ncc=sc), open(os.path.join(a.out, "selftest.json"), "w"), indent=1); raise SystemExit
Af, tf, sc, sc0 = register(edge_e, edge_s, s0); gm = geometry(Af)
quad = []
if not a.no_quadrants:
    for (r0, r1, c0, c1) in [(0, NR // 2, 0, NC // 2), (0, NR // 2, NC // 2, NC), (NR // 2, NR, 0, NC // 2), (NR // 2, NR, NC // 2, NC)]:
        cr = np.array([[c0, c1, c0, c1], [r0, r0, r1, r1]], float); qc = Af @ cr + tf[:, None]          # quadrant inside the SEM field?
        gq = np.mgrid[r0:r1:8, c0:c1:8].reshape(2, -1)[::-1]; xy = Af @ gq + tf[:, None]
        inside = np.mean((xy[0] >= 0) & (xy[0] < edge_s.shape[1]) & (xy[1] >= 0) & (xy[1] < edge_s.shape[0]))
        if inside < 0.6: continue
        sub = np.zeros_like(edge_e); sub[r0:r1, c0:c1] = edge_e[r0:r1, c0:c1]
        def fq(x, sub=sub, r0=r0, r1=r1, c0=c0, c1=c1):
            w, m = warp(x[:4].reshape(2, 2), x[4:], np.pad(np.ones((r1 - r0, c1 - c0)), ((r0, NR - r1), (c0, NC - c1))), edge_s.shape)
            ww, _ = warp(x[:4].reshape(2, 2), x[4:], sub, edge_s.shape); m = m > 0.5
            if m.sum() < 200: return 1.0
            x_, y_ = ww[m] - ww[m].mean(), edge_s[m] - edge_s[m].mean(); return -float((x_ * y_).sum() / np.sqrt((x_ * x_).sum() * (y_ * y_).sum() + 1e-12))
        r = minimize(fq, np.r_[Af.ravel(), tf], method="Powell", options=dict(xtol=1e-4, ftol=1e-6, maxfev=4000)); quad.append(geometry(r.x[:4].reshape(2, 2)) | dict(ncc=round(-r.fun, 3), inside=round(float(inside), 2)))
out = dict(ebsd=os.path.basename(a.ang), sem=os.path.basename(a.sem), nominal_step_um=step, sem_pixel_um_used=round(pd, 4), ncc=round(sc, 4), coarse_ncc=round(sc0, 4),
           best_ncc_other_rotation=round(float(register.other), 4), **gm,
           ratio_columns_to_nominal=round(gm["pixel_um_along_columns"] / step, 4), ratio_rows_to_nominal=round(gm["pixel_um_along_rows"] / step, 4), quadrants=quad)
if len(quad) >= 2:
    for k in ("ratio_rows_to_columns", "map_axis_angle_deg"):
        v = np.array([q[k] for q in quad]); out[k + "_quadrant_sd"] = round(float(v.std(ddof=1)), 4)
out["reliable"] = bool(sc - register.other > 0.1 and out.get("ratio_rows_to_columns_quadrant_sd", 1) < 0.06 and out.get("map_axis_angle_deg_quadrant_sd", 99) < 3)
json.dump(out, open(os.path.join(a.out, "register_sem_ebsd.json"), "w"), indent=1); print(json.dumps({k: v for k, v in out.items() if k != "quadrants"}, indent=1))
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
w, m = warp(Af, tf, edge_e, edge_s.shape); rgb = np.dstack([np.clip(edge_s * 1.2, 0, 1), np.clip(w * 1.2, 0, 1), np.zeros_like(edge_s)])
fig, ax = plt.subplots(1, 2, figsize=(14, 6)); ax[0].imshow(small, cmap="gray"); ax[0].contour(m, [0.5], colors="y", linewidths=0.8); ax[0].set_title("SEM, EBSD footprint (yellow)", fontsize=9)
ax[1].imshow(rgb); ax[1].set_title("red: SEM edges, green: EBSD boundaries (yellow = overlap); NCC %.3f\nrows/cols %.3f, axes %.1f deg" % (sc, gm["ratio_rows_to_columns"], gm["map_axis_angle_deg"]), fontsize=9)
for x in ax: x.set_axis_off()
fig.tight_layout(); fig.savefig(os.path.join(a.out, "register_sem_ebsd.png"), dpi=90)
