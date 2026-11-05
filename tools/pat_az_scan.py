"""Read a Thermo Pathfinder .pat (reverse-engineered layout), Hough-index (PyEBSDIndex) a few patterns and compare with the .ANG
orientation corrected by Rz(az)*Rcorr (as pei_az_scan.py, but for raw .pat).
layout: 20-byte header = 5 x uint32 (version, ncols, nrows, pw, ph); record k: uint16 LE pixels (pw*ph) at 25 + k*(pw*ph*2+8)
usage: pat_az_scan.py --pat P.pat --ang A.ang --rmat M.mat --pc X Y L [--n 12] [--save-png out.png]"""
import argparse, struct, numpy as np, scipy.io as sio
from pyebsdindex import ebsd_index
from scipy.ndimage import zoom
from scipy.spatial.transform import Rotation as SR


def read_pat(path, k):
    with open(path, "rb") as f:
        _, nc, nr, w, h = struct.unpack("<5I", f.read(20))
        n = w * h * 2
        f.seek(25 + k * (n + 8))
        return np.frombuffer(f.read(n), "<u2").reshape(h, w).astype(np.float32), (nc, nr, w, h)


def brady(X, Y, L):            # Pathfinder X,Y,L (mm) -> TSL x*,y*,z* (Brady E.1-E.3, see NOTES.md)
    return 0.5 + 0.5 * ((X - 15.8) / 20), 1 - Y / 20, L / 20


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--pat", required=True); ap.add_argument("--ang", required=True); ap.add_argument("--rmat", required=True)
    ap.add_argument("--pc", nargs=3, type=float, required=True, help="Pathfinder X Y L")
    ap.add_argument("--elev", type=float, default=11.04); ap.add_argument("--n", type=int, default=12); ap.add_argument("--size", type=int, default=232)
    ap.add_argument("--nbands", type=int, default=9); ap.add_argument("--phase", default="FCC"); ap.add_argument("--seed", type=int, default=0); ap.add_argument("--min-ci", type=float, default=0.0)
    a = ap.parse_args()
    M = np.array(sio.loadmat(a.rmat, squeeze_me=True, struct_as_record=False)["calib"].Rmat, float)
    ang = np.loadtxt(a.ang, comments="#"); sym = SR.create_group("O").as_matrix()
    mis = lambda g1, g2: min(np.degrees(np.arccos(np.clip((np.trace(g1.T @ (S @ g2)) - 1) / 2, -1, 1))) for S in sym)
    _, (nc, nr, w, h) = read_pat(a.pat, 0)
    print("pat grid", nc, "x", nr, " pattern", w, "x", h, " ang rows", len(ang))
    pc = brady(*a.pc); print("PC (TSL) x*,y*,z* =", np.round(pc, 4))
    rng = np.random.default_rng(a.seed)
    ids = [i for i in rng.permutation(nc * nr) if np.all(np.isfinite(ang[i, :3])) and ang[i, 6] >= a.min_ci][: a.n]
    tot = np.zeros(360); rows = []
    for k in ids:
        im, _ = read_pat(a.pat, int(k)); x0 = (w - h) // 2; im = im[:, x0:x0 + h]
        pat = zoom(im, a.size / h, order=1).astype(np.float32)[None]
        d, _ = ebsd_index.index_pats(patsin=pat, PC=list(pc), sampleTilt=70.0, camElev=a.elev, vendor="EDAX", phaselist=[a.phase], nRho=90, nTheta=180, nBands=a.nbands, verbose=0)
        d = d[-1]; q = d["quat"].reshape(-1, 4)[0]
        if not np.any(q): print("idx %5d  Hough indexing failed" % k); continue
        gH = SR.from_quat([q[1], q[2], q[3], q[0]]).as_matrix()
        gA = SR.from_euler("ZXZ", ang[k, :3]).as_matrix().T
        res = np.array([mis(gA @ (SR.from_euler("z", al, degrees=True).as_matrix() @ M).T, gH) for al in range(360)])
        rows.append((int(k), d["fit"][0], d["nmatch"][0], res[0], res[90], res.argmin(), res.min())); tot += res
        print("idx %5d  fit=%.2f nmatch=%d | az=0: %5.1f | az=90: %5.1f | best az=%3d: %5.1f" % rows[-1])
    print("median over patterns: az=0 %.1f  az=90 %.1f ; best common az=%d (median-of-sum %.1f)" % (np.median([r[3] for r in rows]), np.median([r[4] for r in rows]), tot.argmin(), tot.min() / len(rows)))
