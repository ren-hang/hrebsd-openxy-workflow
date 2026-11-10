"""Check that an OpenXY grain file (11-column .txt) lists the points in the same order as the scan .ang, and optionally fix it.
OpenXY reads the grain ID (column 9) and phase (column 11) by ROW NUMBER and ignores the grain file's x,y, so a file written in
another row order puts every grain ID on the wrong pattern (seen 2026-11: the lab MTEX export_GrainFile after
'convertSpatial2EulerReferenceFrame' writes the rows bottom-to-top -> OpenXY's grain boundaries are the real ones flipped top-bottom).
Tests the 8 grid orders (x-first/y-first x flips) by Euler angles (exact, or a constant rotation if one file was corrected) and
checks the grain IDs as OpenXY will read them against the .ang orientations (cubic misorientation, orix).
usage: check_grain_file.py --ang A.ang --grain G.txt [--fix OUT.txt] [--phase Platinum] [--tol 5]
exit 0: same order and IDs consistent; 1: different order / inconsistent (with --fix: OUT written in the .ang order)."""
import argparse, re, sys, numpy as np
from scipy.spatial.transform import Rotation as SR
from orix.quaternion import Orientation, symmetry

ap = argparse.ArgumentParser()
ap.add_argument("--ang", required=True); ap.add_argument("--grain", required=True); ap.add_argument("--fix")
ap.add_argument("--phase", help="phase name for column 11 of the fixed file (default: keep; numeric/missing -> required)")
ap.add_argument("--tol", type=float, default=5.0, help="grain tolerance (deg) for the consistency check")
a = ap.parse_args()

hdr = [L for L in open(a.ang, errors="replace") if L.lstrip().startswith("#")]
h = "".join(hdr)
nc = int((re.search(r"NCOLS_ODD:\s*(\d+)", h) or re.search(r"NCOLS:\s*(\d+)", h)).group(1)); nr = int(re.search(r"NROWS:\s*(\d+)", h).group(1))
A = np.loadtxt(a.ang, comments="#", usecols=range(5)); N = nr * nc
assert len(A) == N, f".ang has {len(A)} rows, header {nc} x {nr}"
glines = [L.rstrip("\r\n") for L in open(a.grain, errors="replace")]
ghdr = [L for L in glines if L.lstrip().startswith("#")]
G = [L.split() for L in glines if L.strip() and not L.lstrip().startswith("#")]
if len(G) != N: sys.exit(f"FAIL grain file has {len(G)} rows, .ang has {N}: OpenXY stops ('different number of points')")
GE = np.array([g[:3] for g in G], float); gid = np.array([int(float(g[8])) for g in G])

r, c = np.divmod(np.arange(N), nc)                      # .ang pixel i = (r, c), x steps first
orders = {"same order (x first)": r * nc + c, "flipped top-bottom": (nr - 1 - r) * nc + c, "flipped left-right": r * nc + (nc - 1 - c),
          "rotated 180": (nr - 1 - r) * nc + (nc - 1 - c), "y first": c * nr + r, "y first, flipped top-bottom": c * nr + (nr - 1 - r),
          "y first, flipped left-right": (nc - 1 - c) * nr + r, "y first, rotated 180": (nc - 1 - c) * nr + (nr - 1 - r)}
ok = np.isfinite(A[:, :3]).all(1)
Ra = SR.from_euler("ZXZ", np.where(ok[:, None], A[:, :3], 0)).as_matrix()
o = Orientation.from_euler(np.where(ok[:, None], A[:, :3], 0), symmetry.Oh).reshape(nr, nc)
mis_r = np.degrees(o[:, :-1].angle_with(o[:, 1:])); mis_d = np.degrees(o[:-1, :].angle_with(o[1:, :]))
okr = ok.reshape(nr, nc); vr = okr[:, :-1] & okr[:, 1:]; vd = okr[:-1, :] & okr[1:, :]

def consistency(ids):                                   # ids: grain ID per .ang pixel (nr x nc)
    same_r, same_d = ids[:, :-1] == ids[:, 1:], ids[:-1, :] == ids[1:, :]
    big = np.concatenate([mis_r[vr] > a.tol, mis_d[vd] > a.tol]); same = np.concatenate([same_r[vr], same_d[vd]])
    return float((same & big).sum() / max(big.sum(), 1)), float((~same & ~big).sum() / max((~big).sum(), 1)), int(big.sum())

print(f".ang {nc} x {nr}; grain file {len(G)} rows, {len(np.unique(gid))} grain IDs")
xy_same = np.allclose(np.array([g[3:5] for g in G], float), A[:, 3:5], atol=0.01)
print("grain-file x,y columns equal the .ang x,y row by row:", xy_same, "(OpenXY ignores them either way)")
res = {}
for name, j in orders.items():
    m = ok & np.isfinite(GE[j]).all(1)
    d = np.abs((GE[j][m] - A[m, :3] + np.pi) % (2 * np.pi) - np.pi).max(1)
    M = SR.from_matrix(SR.from_euler("ZXZ", GE[j][m]).as_matrix() @ np.transpose(Ra[m], (0, 2, 1)))
    Mm = M[np.argmin(np.degrees((M * M.mean().inv()).magnitude()))]          # a typical row's relative rotation
    res[name] = (float((d < 1e-3).mean()), float((np.degrees((M * Mm.inv()).magnitude()) < 0.1).mean()), Mm)
best = max(res, key=lambda k: (res[k][0] > 0.99, res[k][1]))
print(f"{'row order':30s} {'rows with identical Euler':>26s} {'rows with one common rotation (0.1 deg)':>40s}")
for k, (e, s, _) in res.items(): print(f"{k:30s} {100 * e:25.1f}% {100 * s:39.1f}%" + ("   <- match" if k == best else ""))
e, s, Mm = res[best]
if s < 0.99: print("no grid order matches the .ang (Euler angles unrelated): cannot tell which grain ID belongs to which point"); sys.exit(1)
if e < 0.99: print(f"orientations match after a constant rotation of {np.degrees(Mm.magnitude()):.2f} deg (one of the files was corrected)")
merged, split, nb = consistency(gid[orders["same order (x first)"]].reshape(nr, nc))
print(f"grain IDs as OpenXY reads them (row i -> .ang point i): {100 * merged:.1f} % of the {nb} neighbour pairs with misorientation "
      f"> {a.tol:g} deg share an ID (missed boundaries); {100 * split:.1f} % of the low-misorientation pairs are split (false boundaries)")
if best == "same order (x first)":
    print("PASS grain file is in the .ang order"); sys.exit(0 if merged < 0.05 else 1)
merged2, split2, _ = consistency(gid[orders[best]].reshape(nr, nc))
print(f"FAIL grain file rows are '{best}' relative to the .ang; reordered: {100 * merged2:.1f} % missed, {100 * split2:.1f} % false boundaries")
if a.fix:
    j = orders[best]; out = list(ghdr)
    for i in range(N):
        g = list(G[j[i]]) + [""] * max(0, 11 - len(G[j[i]]))
        g[3], g[4] = f"{A[i, 3]:.6f}", f"{A[i, 4]:.6f}"
        if a.phase: g[10] = a.phase
        elif not re.search(r"[A-Za-z]", g[10]): sys.exit("column 11 (phase) is numeric or missing: give --phase (e.g. Platinum)")
        out.append(" ".join(g[:11]))
    open(a.fix, "w").write("\n".join(out) + "\n")
    print("wrote", a.fix, "(rows in the .ang order, x,y from the .ang)")
sys.exit(1)
