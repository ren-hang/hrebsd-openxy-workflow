"""Whole-map pattern-centre (+ optional detector tilt) fit with the ebsdsim Pt master pattern and kikuchipy projections.
Points: up to --per-grain well-indexed points (instrument CI >= --min-ci, >= --margin px from a grain boundary, column > 0) in each grain, spread
within the grain; orientation fixed to the corrected ang. Prints the fit, the Pathfinder-unit equivalent (inverse Brady) and per-grain scores.
usage: map_pc_fit.py --pat P.pat --ang orig.ANG --corr corrected.ang --grain G.txt --pc0 x y z [--free-elev] [--master npz] [--ncols 64]"""
import argparse, os, sys, numpy as np, kikuchipy as kp
from orix.quaternion import Rotation
from scipy.ndimage import zoom, gaussian_filter, distance_transform_cdt
from scipy.optimize import minimize
sys.path.insert(0, os.path.dirname(__file__))
from pat_az_scan import read_pat

def prep(im, n):
    h, w = im.shape; x0 = (w - h) // 2; im = im[:, x0:x0 + h].astype(float)
    im = zoom(im, n / h, order=1); im = im - gaussian_filter(im, n / 8); return (im - im.mean()) / im.std()

ap = argparse.ArgumentParser()
ap.add_argument("--pat", required=True); ap.add_argument("--ang", required=True); ap.add_argument("--corr", required=True); ap.add_argument("--grain", required=True)
ap.add_argument("--pc0", nargs=3, type=float, required=True); ap.add_argument("--elev", type=float, default=11.04); ap.add_argument("--tilt", type=float, default=70.0)
ap.add_argument("--free-elev", action="store_true"); ap.add_argument("--master", default="work/pt_master/Pt_20kV_hw250.npz"); ap.add_argument("--ncols", type=int, default=64)
ap.add_argument("--per-grain", type=int, default=2); ap.add_argument("--min-ci", type=float, default=0.5); ap.add_argument("--margin", type=int, default=2); ap.add_argument("--size", type=int, default=232)
a = ap.parse_args(); n = a.size; nc = a.ncols
orig = np.loadtxt(a.ang, comments="#"); corr = np.loadtxt(a.corr, comments="#"); gid = np.loadtxt(a.grain, usecols=(8,), dtype=int); nr = len(orig) // nc
lab = gid.reshape(nr, nc); b = np.zeros((nr, nc), bool)
b[:-1] |= lab[:-1] != lab[1:]; b[1:] |= lab[:-1] != lab[1:]; b[:, :-1] |= lab[:, :-1] != lab[:, 1:]; b[:, 1:] |= lab[:, :-1] != lab[:, 1:]
bd = distance_transform_cdt(~b, metric="chessboard"); ci = orig[:, 6].reshape(nr, nc); ok = (bd >= a.margin) & (ci >= a.min_ci); ok[:, 0] = False; ok &= (lab > 0)
pick, pg = [], []
for g in np.unique(lab[ok]):
    ys, xs = np.nonzero(ok & (lab == g)); k0 = np.argmax(ci[ys, xs]); sel = [k0]
    for _ in range(a.per_grain - 1):
        d = np.min([np.hypot(ys - ys[s], xs - xs[s]) for s in sel], axis=0); sel.append(int(np.argmax(d)))
    for s in sel:
        k = int(ys[s] * nc + xs[s])
        if k not in pick: pick.append(k); pg.append(int(g))
print("%d points in %d grains" % (len(pick), len(set(pg))))
pats = [prep(read_pat(a.pat, k)[0], n) for k in pick]; R0 = Rotation.from_euler(corr[pick, :3])
mp = kp.load(a.master); yy, xx = np.mgrid[0:n, 0:n]; mask = ((xx - n / 2) ** 2 + (yy - n / 2) ** 2) < (0.48 * n) ** 2
def scores(p):
    det = kp.detectors.EBSDDetector(shape=(n, n), pc=p[:3], sample_tilt=a.tilt, tilt=(p[3] if a.free_elev else a.elev), convention="tsl")
    sim = mp.get_patterns(rotations=R0, detector=det, energy=20, compute=True, show_progressbar=False).data.reshape(len(pick), n, n).astype(float)
    out = []
    for k in range(len(pick)):
        q = sim[k] - gaussian_filter(sim[k], n / 8); q = (q - q[mask].mean()) / q[mask].std(); out.append((q[mask] * pats[k][mask]).mean())
    return np.array(out)
ref = np.array(list(a.pc0) + ([a.elev] if a.free_elev else [])); s0 = scores(ref)
print("score at given PC: mean %.3f; per grain:" % s0.mean(), {g: round(float(np.mean([s for s, gg in zip(s0, pg) if gg == g])), 3) for g in sorted(set(pg))})
step = np.array([0.02] * 3 + ([1.0] if a.free_elev else [])); sim0 = np.vstack([ref] + [ref + np.diag(step)[i] for i in range(len(ref))])
res = minimize(lambda p: -scores(p).mean(), ref, method="Nelder-Mead", options=dict(initial_simplex=sim0, xatol=1e-4, fatol=1e-6, maxfev=300))
s1 = scores(res.x); x, y, z = res.x[:3]
print("refined: mean %.3f  PC=(%.4f, %.4f, %.4f)%s  diff from given %s" % (s1.mean(), x, y, z, "  elev %.2f" % res.x[3] if a.free_elev else "", np.round(res.x - ref, 4)))
print("Pathfinder units: X %.2f Y %.2f L %.2f mm" % (15.8 + 40 * (x - 0.5), 20 * (1 - y), 20 * z))
print("per grain after:", {g: round(float(np.mean([s for s, gg in zip(s1, pg) if gg == g])), 3) for g in sorted(set(pg))})
