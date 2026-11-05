"""Cut an N x N block out of a Pathfinder .pat/.ANG map for OpenXY: writes <stem>_r{r}c{c}.tif (uint16, pattern as stored),
<stem>.ang and an 11-column grain file <stem>.txt (grain id = connected cubic misorientation < tol, phase name given).
The block is placed in the largest grain (or at --r0 --c0).  Orientations are written unmodified (use correct_ang.py afterwards).
usage: make_pat_subset.py --pat P.pat --ang A.ang --n 10 --out DIR --stem ORR10 [--r0 R --c0 C] [--tol 5] [--phase Platinum]"""
import argparse, os, re, sys, numpy as np, tifffile
from scipy.ndimage import label
from scipy.spatial.transform import Rotation as SR
sys.path.insert(0, os.path.dirname(__file__))
from pat_az_scan import read_pat

ap = argparse.ArgumentParser()
ap.add_argument("--pat", required=True); ap.add_argument("--ang", required=True); ap.add_argument("--n", type=int, default=10); ap.add_argument("--full", action="store_true", help="whole map instead of an N x N block")
ap.add_argument("--out", required=True); ap.add_argument("--stem", required=True)
ap.add_argument("--r0", type=int); ap.add_argument("--c0", type=int); ap.add_argument("--tol", type=float, default=5.0)
ap.add_argument("--phase", default="Platinum")
ap.add_argument("--ref-margin", type=int, default=0, help="set IQ=0, CI=0, Fit=99 for points within this many pixels of a grain boundary so reference selection (votes on max IQ, max CI, min Fit) avoids them")
a = ap.parse_args()
hdr, rows = [], []
for L in open(a.ang, errors="replace"):
    (hdr if L.lstrip().startswith("#") else rows).append(L.rstrip("\r\n"))
rows = [r for r in rows if r.strip()]
nc = int(re.search(r"NCOLS:\s*(\d+)", "\n".join(hdr)).group(1)); nr = int(re.search(r"NROWS:\s*(\d+)", "\n".join(hdr)).group(1))
V = np.array([r.split()[:10] for r in rows], float); assert len(V) == nc * nr
E = V[:, :3].reshape(nr, nc, 3); ok = np.isfinite(E).all(2)
g = np.zeros((nr, nc, 3, 3)); g[ok] = SR.from_euler("ZXZ", E[ok]).as_matrix()
sym = SR.create_group("O").as_matrix()
def mis(g1, g2):                      # cubic misorientation (deg) between two orientation matrices
    best = max(np.einsum("ij,jk,ik->", g1, S, g2) for S in sym)
    return np.degrees(np.arccos(np.clip((best - 1) / 2, -1, 1)))
# label grains by flood fill over 4-neighbours with misorientation < tol
lab = -np.ones((nr, nc), int); cur = 0
for r in range(nr):
    for c in range(nc):
        if not ok[r, c] or lab[r, c] >= 0: continue
        stack = [(r, c)]; lab[r, c] = cur
        while stack:
            y, x = stack.pop()
            for dy, dx in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                y2, x2 = y + dy, x + dx
                if 0 <= y2 < nr and 0 <= x2 < nc and ok[y2, x2] and lab[y2, x2] < 0 and mis(g[y, x], g[y2, x2]) < a.tol:
                    lab[y2, x2] = cur; stack.append((y2, x2))
        cur += 1
sizes = np.bincount(lab[lab >= 0]); big = int(sizes.argmax())
from scipy.ndimage import distance_transform_cdt
bnd = np.zeros((nr, nc), bool)
bnd[:-1, :] |= lab[:-1, :] != lab[1:, :]; bnd[1:, :] |= lab[:-1, :] != lab[1:, :]
bnd[:, :-1] |= lab[:, :-1] != lab[:, 1:]; bnd[:, 1:] |= lab[:, :-1] != lab[:, 1:]
bdist = distance_transform_cdt(~bnd, metric="chessboard")
print("grains:", cur, " largest grain id", big, "size", sizes[big])
n = a.n
if a.full: n = None
if a.full: a.r0, a.c0 = 0, 0; nB, mB = nr, nc
else: nB = mB = n
if a.r0 is None:                      # place the block where it fits entirely inside the largest grain, nearest the grain centroid
    ys, xs = np.nonzero(lab == big); cy, cx = ys.mean(), xs.mean(); best = None
    for r in range(nr - n + 1):
        for c in range(nc - n + 1):
            blk = lab[r:r + n, c:c + n]; frac = (blk == big).mean()
            sc = (-frac, np.hypot(r + n / 2 - cy, c + n / 2 - cx))
            if best is None or sc < best[0]: best = (sc, r, c)
    a.r0, a.c0 = best[1], best[2]; print("block fraction in grain %.2f" % -best[0][0])
print("block r0,c0 =", a.r0, a.c0)
os.makedirs(a.out, exist_ok=True)
sel = [(r, c) for r in range(a.r0, a.r0 + nB) for c in range(a.c0, a.c0 + mB)]
x0 = V[a.r0 * nc + a.c0, 3]; y0 = V[a.r0 * nc + a.c0, 4]
new_hdr = [re.sub(r"(NROWS:\s*)\d+", r"\g<1>%d" % nB, re.sub(r"(NCOLS(?:_ODD|_EVEN)?:\s*)\d+", r"\g<1>%d" % mB, h)) for h in hdr]
grain_out, ids, meta = [], {}, []
for i, (r, c) in enumerate(sel):
    k = r * nc + c; v = V[k]
    rr, cc = divmod(i, mB)
    xs_, ys_ = v[3] - x0, v[4] - y0
    gid = ids.setdefault(lab[r, c], len(ids) + 1) if lab[r, c] >= 0 else 0
    gs = int(sizes[lab[r, c]]) if lab[r, c] >= 0 else 0
    near = bdist[r, c] < a.ref_margin      # OpenXY votes on max IQ, max CI and min Fit: make all three unattractive near boundaries
    iq, ci, fit = (0, 0, 99) if near else (v[5], v[6], 1)
    meta.append((k, v[6], c))
    grain_out.append("%.6f %.6f %.6f %.6f %.6f %g %g %g %d 1 %s" % (v[0], v[1], v[2], xs_, ys_, iq, ci, fit, gid, a.phase))
    im, _ = read_pat(a.pat, k)
    tifffile.imwrite(os.path.join(a.out, "%s_r%dc%d.tif" % (a.stem, rr, cc)), im.astype(np.uint16))
# .ang rows: shift x,y so the subset starts at 0 (OpenXY derives step from them)
fixed = []
for i, (r, c) in enumerate(sel):
    t = rows[r * nc + c].split(); rr, cc = divmod(i, mB)
    t[3] = "%.6f" % (V[r * nc + c, 3] - x0); t[4] = "%.6f" % (V[r * nc + c, 4] - y0)
    if bdist[r, c] < a.ref_margin: t[5], t[6], t[-1] = "0", "0", "99"
    fixed.append("  ".join(t))
open(os.path.join(a.out, a.stem + ".ang"), "w").write("\n".join(new_hdr + fixed) + "\n")
open(os.path.join(a.out, a.stem + ".txt"), "w").write("\n".join(grain_out) + "\n")
np.savetxt(os.path.join(a.out, a.stem + "_index.txt"), np.array(meta), fmt="%d %g %d", header="original_index instrument_CI original_column")
print("wrote %d patterns + ang + grain file to %s; grains in block: %d" % (len(sel), a.out, len(ids)))
