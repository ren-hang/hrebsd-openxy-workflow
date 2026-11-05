"""Synthetic OpenXY validation dataset from a kikuchipy master pattern.

Renders a 10x10 'scan' of one FCC (Ni) orientation with KNOWN deformations, using
OpenXY's own geometry model (CalcF.m / Theoretical_Pixel_Shift.m) for pixel -> ray, and
kikuchipy's master-pattern sampler for intensity. Deformation convention (identical to
Theoretical_Pixel_Shift): a feature at crystal-frame ray d appears at ray F*d, so
I'(ray) = I_ref(F^-1 * ray), F expressed in the crystal frame.

Outputs in --out:  <stem>.ang, <stem>.txt (grain file), Corrected_EBSPs/Ni_r{Y}c{X}.tif,
                   truth.csv (one row per scan point), render.json (geometry used).
Index i = r*10 + c  (x fastest), exactly the order OpenXY expects.
"""
import argparse
import json
from pathlib import Path

import numpy as np
import tifffile
import kikuchipy as kp
from kikuchipy.signals.util._master_pattern import _project_single_pattern_from_master_pattern

N = 10


def euler2gmat(p1, P, p2):
    """TSL/Bunge sample->crystal matrix (same as OpenXY euler2gmat)."""
    c1, s1, c, s, c2, s2 = np.cos(p1), np.sin(p1), np.cos(P), np.sin(P), np.cos(p2), np.sin(p2)
    return np.array([
        [c1 * c2 - s1 * s2 * c, s1 * c2 + c1 * s2 * c, s2 * s],
        [-c1 * s2 - s1 * c2 * c, -s1 * s2 + c1 * c2 * c, c2 * s],
        [s1 * s, -c1 * s, c]])


def rot_axis(axis, ang):
    a = np.asarray(axis, float)
    a /= np.linalg.norm(a)
    K = np.array([[0, -a[2], a[1]], [a[2], 0, -a[0]], [-a[1], a[0], 0]])
    return np.eye(3) + np.sin(ang) * K + (1 - np.cos(ang)) * K @ K


def build_rays_crystal(P, xs, ys, zs, tilt, elev, Qsc, shape=None):
    """Unit rays (crystal frame) for every pixel; pixel (col c, row r) 1-based <-> view coords (xc,yc)=(c,r)."""
    alpha = np.pi / 2 - tilt + elev
    Qps = np.array([[0, -np.cos(alpha), -np.sin(alpha)],
                    [-1, 0, 0],
                    [0, np.sin(alpha), -np.cos(alpha)]])
    rr, cc = np.mgrid[1:P + 1, 1:P + 1].astype(float)
    Xp = xs * P - cc
    Yp = (1 - ys) * P - rr
    pts_p = np.stack([Xp.ravel(), Yp.ravel(), np.full(Xp.size, -zs * P)], axis=1)
    pts_s = pts_p @ Qps.T
    pts_c = pts_s @ Qsc.T
    return pts_c / np.linalg.norm(pts_c, axis=1, keepdims=True)


def make_cases(rng, profile="small"):
    """List of dicts: name, F (crystal frame), shift (dx,dy px), noise (fraction of range)."""
    cases = [dict(group="ref", name="ref", F=np.eye(3), shift=(0, 0), noise=0.0)]
    comps = {"e11": (0, 0), "e22": (1, 1), "e33": (2, 2), "e12": (0, 1), "e13": (0, 2), "e23": (1, 2)}
    for k, (i, j) in comps.items():
        for a in (-5e-3, -1e-3, 1e-3, 5e-3):
            E = np.zeros((3, 3))
            E[i, j] = E[j, i] = a
            cases.append(dict(group="strain", name=f"{k}={a:+.0e}", F=np.eye(3) + E, shift=(0, 0), noise=0.0))
    for ax, v in (("x", (1, 0, 0)), ("y", (0, 1, 0)), ("z", (0, 0, 1))):
        for deg in (-0.5, -0.1, 0.1, 0.5):
            cases.append(dict(group="rot", name=f"rot{ax}={deg:+.1f}deg", F=rot_axis(v, np.radians(deg)), shift=(0, 0), noise=0.0))
    for dx, dy in [(0.5, 0), (-0.5, 0), (2, 0), (-2, 0), (0, 0.5), (0, -0.5), (0, 2), (0, -2),
                   (0.5, 0.5), (-0.5, -0.5), (2, 2), (-2, -2), (2, -2), (-2, 2), (0.1, 0.1), (5, 0)]:
        cases.append(dict(group="shift", name=f"shift=({dx};{dy})px", F=np.eye(3), shift=(dx, dy), noise=0.0))
    for k in range(12):
        E = rng.normal(0, 1.5e-3, (3, 3))
        E = (E + E.T) / 2
        W = rot_axis(rng.normal(size=3), np.radians(rng.uniform(-0.3, 0.3)))
        cases.append(dict(group="mixed", name=f"mixed{k}", F=W @ (np.eye(3) + E), shift=(0, 0), noise=0.0))
    if profile == "bigrot":
        # larger misorientations (typical intra-grain lattice rotation range) to expose pattern-centre sensitivity
        for ax, v in (("x", (1, 0, 0)), ("y", (0, 1, 0)), ("z", (0, 0, 1))):
            for deg in (-4, -2, -1, 1, 2, 4):
                cases.append(dict(group="bigrot", name=f"rot{ax}={deg:+.0f}deg", F=rot_axis(v, np.radians(deg)), shift=(0, 0), noise=0.0))
        for k in range(17):
            E = rng.normal(0, 1e-3, (3, 3))
            E = (E + E.T) / 2
            W = rot_axis(rng.normal(size=3), np.radians(rng.uniform(-3, 3)))
            cases.append(dict(group="bigmix", name=f"bigmix{k}", F=W @ (np.eye(3) + E), shift=(0, 0), noise=0.0))
    else:
        for lev in (0.01, 0.03, 0.05, 0.10, 0.20):
            for k in range(7):
                cases.append(dict(group="noise", name=f"noise{lev:.2f}_{k}", F=np.eye(3), shift=(0, 0), noise=lev))
    assert len(cases) == N * N, len(cases)
    return cases


def fourier_shift(img, dx, dy):
    """Sub-pixel translation by (dx,dy): image content moves right by dx, down by dy."""
    if dx == 0 and dy == 0:
        return img
    ky = np.fft.fftfreq(img.shape[0])[:, None]
    kx = np.fft.fftfreq(img.shape[1])[None, :]
    return np.real(np.fft.ifft2(np.fft.fft2(img) * np.exp(-2j * np.pi * (kx * dx + ky * dy))))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--stem", default="Syn10x10")
    ap.add_argument("--P", type=int, default=1160)
    ap.add_argument("--pc", nargs=3, type=float, default=[0.51275, 0.657, 0.9545], help="x* y* z* (OpenXY convention, fractions of P)")
    ap.add_argument("--tilt", type=float, default=70.0)
    ap.add_argument("--elev", type=float, default=0.0)
    ap.add_argument("--euler", nargs=3, type=float, default=[20.0, 35.0, 10.0], help="Bunge phi1 PHI phi2 (deg)")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--profile", choices=["small", "bigrot"], default="small")
    a = ap.parse_args()

    out = Path(a.out)
    (out / "Corrected_EBSPs").mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(a.seed)

    mp = kp.data.nickel_ebsd_master_pattern_small(projection="lambert", hemisphere="both")
    npx, npy = mp.axes_manager.signal_shape
    scale = (npx - 1) / 2
    mup, mlo = mp._get_master_pattern_arrays_from_energy(None)
    mup, mlo = np.ascontiguousarray(mup, dtype=np.float64), np.ascontiguousarray(mlo, dtype=np.float64)

    eul = np.radians(a.euler)
    Qsc = euler2gmat(*eul)
    xs, ys, zs = a.pc
    tilt, elev = np.radians(a.tilt), np.radians(a.elev)
    rays = build_rays_crystal(a.P, xs, ys, zs, tilt, elev, Qsc)
    ident = np.array([1.0, 0, 0, 0])

    def render(F):
        d = rays @ np.linalg.inv(F).T          # (F^-1 ray) as row vectors
        d /= np.linalg.norm(d, axis=1, keepdims=True)
        v = _project_single_pattern_from_master_pattern(ident, np.ascontiguousarray(d), mup, mlo, int(npx), int(npy), float(scale), False, 1.0, 2.0, np.float64)
        return v.reshape(a.P, a.P)

    cases = make_cases(rng, a.profile)
    ref = render(np.eye(3))
    lo, hi = float(ref.min()), float(ref.max())
    print(f"reference render range {lo:.3g}..{hi:.3g}")

    truth = []
    for idx, cs in enumerate(cases):
        img = ref if np.allclose(cs["F"], np.eye(3)) else render(cs["F"])
        img = fourier_shift(img, *cs["shift"])
        if cs["noise"] > 0:
            img = img + rng.normal(0, cs["noise"] * (hi - lo), img.shape)
        u16 = np.clip((img - lo) / (hi - lo), 0, 1) * 65535
        r, c = divmod(idx, N)
        tifffile.imwrite(out / "Corrected_EBSPs" / f"Ni_r{r}c{c}.tif", np.round(u16).astype(np.uint16))
        truth.append((idx, r, c, cs["group"], cs["name"], *cs["F"].ravel(), *cs["shift"], cs["noise"]))

    hdr = "idx,row,col,group,name," + ",".join(f"F{i+1}{j+1}" for i in range(3) for j in range(3)) + ",dx,dy,noise"
    with open(out / "truth.csv", "w") as f:
        f.write(hdr + "\n")
        for t in truth:
            f.write(",".join(str(x) for x in t) + "\n")

    # ang + grain file (point 0 is the reference: best IQ, CI, Fit)
    step = 1.0
    ang = ["# x-star   0.000000", "# y-star   0.000000", "# z-star   0.000000", "# WorkingDistance 10.0", "#",
           "# Phase 1", "# MaterialName Nickel", "# Formula", "# Info", "# Symmetry 43",
           "# LatticeConstants 3.5236 3.5236 3.5236 90.000 90.000 90.000", "# NumberFamilies 0", "# Categories0 0 0 0 0", "#",
           "# GRID: SqrGrid", f"# XSTEP: {step:.6f}", f"# YSTEP: {step:.6f}", f"# NCOLS_ODD: {N}", f"# NCOLS_EVEN: {N}", f"# NROWS: {N}", "#"]
    grain = []
    for idx in range(N * N):
        r, c = divmod(idx, N)
        iq, ci = (1000.0, 0.99) if idx == 0 else (500.0, 0.90)
        row = f"{eul[0]:.6f} {eul[1]:.6f} {eul[2]:.6f} {c*step:.6f} {r*step:.6f} {iq:.1f} {ci:.3f} 0 0 1.000"
        ang.append(row)
        grain.append(f"{eul[0]:.6f} {eul[1]:.6f} {eul[2]:.6f} {c*step:.6f} {r*step:.6f} {iq:.1f} {ci:.3f} 1 1 0 Nickel")
    (out / f"{a.stem}.ang").write_text("\n".join(ang) + "\n")
    (out / f"{a.stem}.txt").write_text("\n".join(grain) + "\n")
    (out / "render.json").write_text(json.dumps(dict(P=a.P, pc_xyz=a.pc, tilt_deg=a.tilt, camera_elev_deg=a.elev,
                                                     euler_bunge_deg=a.euler, master="kikuchipy nickel_ebsd_master_pattern_small 20 kV",
                                                     pixel_coords="1-based (col,row) == OpenXY view coords (xc,yc)", seed=a.seed), indent=1))
    print("wrote", out)


if __name__ == "__main__":
    main()
