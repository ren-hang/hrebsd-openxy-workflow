"""Measure the true step of an earlier beam scan from the contamination dots it left, in a top-view SEM image (NOTES 15.4).
Long-dwell scans (EBSD/HR-EBSD maps, point grids) leave a dark contamination dot at each beam position. In a 0 deg tilt, tilt-correction-off
SEM image the dot lattice gives the physical step along both scan axes and the angle between them, with no registration and no shape assumption.
Dots: normalised cross-correlation with a dark Gaussian spot of --dot-um diameter on the flattened, block-averaged image; lattice: dominant
nearest-neighbour directions, refined by least squares on integer-indexed dots (residual and on-lattice fraction reported).
EXPERIMENTAL: validated on one image (Deba 15 Jan 2026 hole1_map2_001: 9.91 x 9.71 um, 90.0 deg, 47 of 83 dots, rms 0.27 um, consistent with the
grain-boundary registration of the same session). On images with many other small dark features (Apr 2025 set) it finds no stable lattice:
trust a result only if most dots are on the lattice and the angle is near 90 deg; otherwise crop to one clean dot array first.
usage: scan_dot_lattice.py SEM.tif [SEM.tif ...] [--dot-um 2.2] [--min-sep-um 7.5] [--max-sep-um 13] [--thr 0.35] [--out DIR]
Prints per image: a1, a2 (um, direction in the image, y up), a2/a1, angle between, dots used / found, residual; writes DIR/scan_dot_lattice.json."""
import argparse, json, os, re, numpy as np, tifffile
from scipy import ndimage as ndi
from scipy.spatial import cKDTree
from skimage.feature import match_template, peak_local_max

ap = argparse.ArgumentParser(); ap.add_argument("sem", nargs="+"); ap.add_argument("--dot-um", type=float, default=2.2); ap.add_argument("--min-sep-um", type=float, default=7.5)
ap.add_argument("--max-sep-um", type=float, default=13.0); ap.add_argument("--thr", type=float, default=0.35); ap.add_argument("--out", default=None); a = ap.parse_args()
out = []
for f in a.sem:
    p = tifffile.TiffFile(f).pages[0]; s = str(p.tags[34682].value) if 34682 in p.tags else ""
    m1 = re.search(r"PixelWidth'?\s*:\s*'?([-0-9.eE+]+)", s); tl = re.search(r"StageT'?\s*:\s*'?([-0-9.eE+]+)", s)
    if not m1: print("skip (no pixel size):", f); continue
    px = float(m1.group(1)) * 1e6; tilt = np.degrees(float(tl.group(1))) if tl else 0.0
    img = p.asarray(); img = (img[..., 0] if img.ndim == 3 else img).astype(np.float32); img = img[:int(img.shape[0] * 0.935)]
    k = max(1, int(round(0.08 / px)))                                                       # work at ~80 nm pixels
    H, W = img.shape[0] // k, img.shape[1] // k; b = img[:H * k, :W * k].reshape(H, k, W, k).mean((1, 3)); pk_ = px * k
    sd = a.dot_um / 4 / pk_; b = ndi.gaussian_filter(b, max(0.5, sd / 3)) - ndi.gaussian_filter(b, 6 * sd)
    r_ = int(3 * sd); yy, xx = np.mgrid[-r_:r_ + 1, -r_:r_ + 1]; tpl = -np.exp(-(xx ** 2 + yy ** 2) / (2 * sd ** 2)).astype(np.float32)
    cc = match_template(b, tpl, pad_input=True); P = peak_local_max(cc, min_distance=int(0.6 * a.min_sep_um / pk_), threshold_abs=a.thr)[:, ::-1].astype(float) * pk_
    res = dict(image=os.path.basename(f), pixel_um=round(px, 5), stage_tilt_deg=round(tilt, 2), dots_found=int(len(P)))
    if len(P) < 12: print("%-40s dots %d: too few" % (res["image"], len(P))); out.append(res); continue
    _, j = cKDTree(P).query(P, k=min(7, len(P))); V = (P[j[:, 1:]] - P[:, None]).reshape(-1, 2); r = np.hypot(*V.T); V = V[(r > a.min_sep_um) & (r < a.max_sep_um)]
    ang = np.degrees(np.arctan2(-V[:, 1], V[:, 0])) % 180; h, e = np.histogram(ang, bins=np.arange(0, 181, 3)); order = np.argsort(h)[::-1]
    def vec(a0):
        mm = np.abs(((ang - a0) + 90) % 180 - 90) < 6; v = V[mm]; v = v * np.sign(v @ [np.cos(np.radians(a0)), -np.sin(np.radians(a0))])[:, None]
        rr = np.hypot(*v.T); return np.median(v[np.abs(rr - np.median(rr)) < 1.0], 0)
    a1 = e[order[0]] + 1.5; a2 = next((e[i] + 1.5 for i in order[1:] if abs(((e[i] + 1.5 - a1) + 90) % 180 - 90) > 30), None)
    if a2 is None: print("%-40s no second lattice direction" % res["image"]); out.append(res); continue
    v1, v2 = vec(a1), vec(a2)
    for _ in range(5):
        L = np.c_[v1, v2]; uv = np.linalg.solve(L, (P - P[0]).T).T; kk = np.round(uv); good = np.hypot(*((uv - kk) @ L.T).T) < 0.12 * min(np.linalg.norm(v1), np.linalg.norm(v2))
        if good.sum() < 6: break
        sol, *_ = np.linalg.lstsq(np.c_[kk[good], np.ones(good.sum())], P[good], rcond=None); v1, v2 = sol[0], sol[1]
    l1, l2 = np.linalg.norm(v1), np.linalg.norm(v2); rms = np.sqrt(np.mean(np.hypot(*(P[good] - (kk[good] @ np.vstack([v1, v2]) + sol[2])).T) ** 2))
    res.update(a1_um=round(float(l1), 4), a1_dir_deg=round(float(np.degrees(np.arctan2(-v1[1], v1[0]))), 2), a2_um=round(float(l2), 4), a2_dir_deg=round(float(np.degrees(np.arctan2(-v2[1], v2[0]))), 2),
               a2_over_a1=round(float(l2 / l1), 4), angle_deg=round(float(np.degrees(np.arccos(abs(v1 @ v2) / l1 / l2))), 2), dots_on_lattice=int(good.sum()), residual_um=round(float(rms), 3))
    print("%-40s a1 %.3f um (%.1f deg), a2 %.3f um (%.1f deg), a2/a1 %.4f, angle %.2f deg, %d/%d dots, rms %.3f um, tilt %.1f"
          % (res["image"], l1, res["a1_dir_deg"], l2, res["a2_dir_deg"], l2 / l1, res["angle_deg"], good.sum(), len(P), rms, tilt)); out.append(res)
if a.out: os.makedirs(a.out, exist_ok=True); json.dump(out, open(os.path.join(a.out, "scan_dot_lattice.json"), "w"), indent=1)
