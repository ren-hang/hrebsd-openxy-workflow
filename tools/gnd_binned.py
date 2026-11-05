"""GND from the saved rotation field with k x k binning (noise floor falls ~ 1/(k * step)).
Rotation vectors (deg, OpenXY frame, maps.npz rotx/roty/rotz) are averaged over k x k blocks of usable, same-grain points; curvature kappa_ij = d(theta_i)/d(x_j)
(forward difference of neighbouring block means, both blocks in the same grain); Nye alpha_13 = kappa_31, alpha_23 = kappa_32, alpha_33 = -(kappa_11 + kappa_22);
GND = 3 * (|alpha_13| + |alpha_23| + |alpha_33|) / b   (same 3-component form and 30/10 prefactor as OpenXY alpha_total3; b = 0.2775 nm for Pt).
k = 1 is a check against OpenXY's gnd3. Returns (gnd, ok) on the full grid, each block's value written back to its member points.
usage: gnd_binned.py <maps.npz> [--k 1 2 3] [--b 0.2775e-9]  (prints medians per grain)"""
import argparse, numpy as np

def gnd_binned(d, k, b=0.2775e-9):
    nr, nc = int(d["nrows"]), int(d["ncols"]); step = float(d["x_um"][0, 1] - d["x_um"][0, 0]) * 1e-6
    th = np.radians(np.stack([d["rotx"], d["roty"], d["rotz"]], -1)); use = d["usable"] & ~d["gnd_filtered"] if "gnd_filtered" in d else d["usable"]; g = d["grain"]
    nbr, nbc = nr // k, nc // k; M = np.full((nbr, nbc, 3), np.nan); G = np.zeros((nbr, nbc), int)
    for i in range(nbr):
        for j in range(nbc):
            sl = (slice(i * k, (i + 1) * k), slice(j * k, (j + 1) * k)); gg = g[sl]; uu = use[sl]
            if uu.sum() < max(1, (k * k + 1) // 2) or len(np.unique(gg[uu])) != 1: continue
            M[i, j] = th[sl][uu].mean(0); G[i, j] = gg[uu][0]
    h = k * step; out = np.full((nbr, nbc), np.nan)
    for i in range(nbr - 1):
        for j in range(nbc - 1):
            if not np.isfinite(M[i, j, 0]) or not (np.isfinite(M[i, j + 1, 0]) and np.isfinite(M[i + 1, j, 0])) or G[i, j + 1] != G[i, j] or G[i + 1, j] != G[i, j]: continue
            dx = (M[i, j + 1] - M[i, j]) / h; dy = (M[i + 1, j] - M[i, j]) / h          # kappa_i1, kappa_i2
            out[i, j] = 3 * (abs(dx[2]) + abs(dy[2]) + abs(dx[0] + dy[1])) / b          # alpha13 = kappa_31, alpha23 = kappa_32, alpha33 = -(kappa_11 + kappa_22)
    full = np.full((nr, nc), np.nan)
    for i in range(nbr):
        for j in range(nbc): full[i * k:(i + 1) * k, j * k:(j + 1) * k] = out[i, j]
    return full

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("npz"); ap.add_argument("--k", nargs="+", type=int, default=[1, 2, 3]); ap.add_argument("--b", type=float, default=0.2775e-9)
    a = ap.parse_args(); d = np.load(a.npz)
    for k in a.k:
        f = gnd_binned(d, k, a.b); s = "k=%d" % k
        for gi in np.unique(d["grain"][d["usable"]]):
            m = d["usable"] & (d["grain"] == gi) & np.isfinite(f)
            if m.sum() > 50: s += " | grain %d n=%d median %.2e" % (gi, m.sum(), np.median(f[m])) + ("  (OpenXY gnd3 %.2e)" % np.median(d["gnd3"][m]) if k == 1 and "gnd3" in d else "")
        print(s)
